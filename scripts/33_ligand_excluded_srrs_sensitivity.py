from pathlib import Path
import itertools
import math

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
import matplotlib.pyplot as plt


PROJECT = Path.cwd()
DISCOVERY_SCORES = PROJECT / "results" / "srrs_framework" / "GSE165816_SRRS_discovery_sample_scores.tsv"
LOCKED_GENE_SETS = PROJECT / "results" / "srrs_framework" / "SRRS_locked_gene_sets.tsv"
GSE223964_H5AD = PROJECT / "data" / "raw" / "GSE223964" / "GSE223964_Integrated_all_cells.h5ad"
SPATIAL_SPOT_FILE = PROJECT / "results" / "srrs_framework" / "spatial_gse241124" / "GSE241124_spatial_SRRS_spot_scores.tsv"
OUT_DIR = PROJECT / "results" / "srrs_framework" / "ligand_excluded_srrs"
FIG_DIR = PROJECT / "figures" / "srrs_framework" / "ligand_excluded_srrs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

LIGANDS = ["TNC", "IL11", "IL6", "INHBA", "THBS1", "SERPINE1", "CCL20", "WNT5A", "ADAM12", "PTGS2"]

MARKER_SETS = {
    "keratinocyte": ["KRT14", "KRT5", "KRT1", "KRT10", "KRT15", "DSG1", "DSP"],
    "fibroblast_stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A1"],
    "myeloid": ["LST1", "LYZ", "AIF1", "C1QA", "C1QB", "CD68", "FCGR3A"],
    "endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "CLDN5"],
    "pericyte_smc": ["RGS5", "PDGFRB", "MCAM", "CSPG4", "ACTA2", "TAGLN", "MYH11"],
    "t_nk": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"],
    "b_plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1"],
    "melanocyte": ["PMEL", "MLANA", "TYR", "DCT"],
    "schwann": ["MPZ", "PLP1", "SOX10", "S100B"],
}


def z(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.full_like(x, np.nan)
    return (x - np.nanmean(x)) / sd


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    pv = p[ok]
    if pv.size == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * pv.size / (np.arange(pv.size) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty_like(adj)
    tmp[order] = np.clip(adj, 0, 1)
    out[ok] = tmp
    return out


def exact_perm_p(values, labels):
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    ok = np.isfinite(values)
    values = values[ok]
    labels = labels[ok]
    n_h = int(np.sum(labels == "Healer"))
    obs = values[labels == "Healer"].mean() - values[labels == "Non-healer"].mean()
    deltas = []
    for idx in itertools.combinations(range(len(values)), n_h):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(idx)] = True
        deltas.append(values[mask].mean() - values[~mask].mean())
    deltas = np.asarray(deltas)
    return float((np.sum(np.abs(deltas) >= abs(obs)) + 1) / (len(deltas) + 1))


def bootstrap_ci(values, labels, n_boot=10000, seed=1):
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    h = values[(labels == "Healer") & np.isfinite(values)]
    n = values[(labels == "Non-healer") & np.isfinite(values)]
    if len(h) < 2 or len(n) < 2:
        return math.nan, math.nan
    boot = np.array([rng.choice(h, size=len(h), replace=True).mean() - rng.choice(n, size=len(n), replace=True).mean() for _ in range(n_boot)])
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def test_group(df, feature):
    values = df[feature].to_numpy(dtype=float)
    labels = df["healing_status"].to_numpy()
    h = values[labels == "Healer"]
    n = values[labels == "Non-healer"]
    lo, hi = bootstrap_ci(values, labels)
    return {
        "feature": feature,
        "healer_n": int(np.sum(labels == "Healer")),
        "nonhealer_n": int(np.sum(labels == "Non-healer")),
        "healer_mean": float(np.nanmean(h)),
        "nonhealer_mean": float(np.nanmean(n)),
        "delta_healer_minus_nonhealer": float(np.nanmean(h) - np.nanmean(n)),
        "p_mannwhitney": float(stats.mannwhitneyu(h, n, alternative="two-sided").pvalue),
        "exact_permutation_p": exact_perm_p(values, labels),
        "bootstrap_ci_low": lo,
        "bootstrap_ci_high": hi,
    }


def lodo(df, feature):
    full = test_group(df, feature)["delta_healer_minus_nonhealer"]
    rows = []
    for i, row in df.iterrows():
        sub = df.drop(index=i)
        delta = test_group(sub, feature)["delta_healer_minus_nonhealer"]
        rows.append({
            "feature": feature,
            "dropped_sample": row["sample_code"],
            "delta_healer_minus_nonhealer": delta,
            "same_direction_as_full": np.sign(delta) == np.sign(full),
        })
    return pd.DataFrame(rows)


def score_gene_set(x, var_names, genes):
    genes = [g for g in genes if g in var_names]
    if len(genes) < 3:
        return np.full(x.shape[0], np.nan), genes
    idx = [var_names.get_loc(g) for g in genes]
    arr = x[:, idx].toarray().astype(float) if sparse.issparse(x) else np.asarray(x[:, idx], dtype=float)
    means = arr.mean(axis=0)
    stds = arr.std(axis=0)
    stds[stds == 0] = 1.0
    return ((arr - means) / stds).mean(axis=1), genes


def get_matrix(adata):
    return adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(np.asarray(adata.X))


def load_gene_sets():
    tbl = pd.read_csv(LOCKED_GENE_SETS, sep="\t")
    sets = {m: sub["gene"].dropna().astype(str).tolist() for m, sub in tbl.groupby("module", sort=False)}
    repair = [g for g in sets["fibroblast_repair_activation"] if g not in set(LIGANDS)]
    d7 = [g for g in sets["d7_acute_wound_alignment"] if g not in set(LIGANDS)]
    receiver = [g for g in sets["vascular_perivascular_receiver_coupling"] if g not in set(LIGANDS)]
    return sets, {
        "repair_no_ligands": repair,
        "d7_no_ligands": d7,
        "receiver_no_ligands": receiver,
    }


def discovery_ligand_excluded():
    df = pd.read_csv(DISCOVERY_SCORES, sep="\t")
    # Existing sample scores already include z-transformed repair, D7 and receiver modules.
    # Excluding the explicit ligand panel is a first non-circular check.
    df["SRRS_no_ligand_panel_existing"] = df[[
        "z_fibroblast_repair_activation",
        "z_d7_acute_wound_alignment",
        "z_vascular_perivascular_receiver_coupling",
    ]].mean(axis=1)
    tests = [test_group(df, "SRRS_no_ligand_panel_existing")]
    lod = lodo(df, "SRRS_no_ligand_panel_existing")
    tests = pd.DataFrame(tests)
    lod_sum = lod.groupby("feature").agg(
        lodo_min_delta=("delta_healer_minus_nonhealer", "min"),
        lodo_max_delta=("delta_healer_minus_nonhealer", "max"),
        lodo_same_direction=("same_direction_as_full", lambda x: f"{int(x.sum())}/{len(x)}"),
    ).reset_index()
    tests = tests.merge(lod_sum, on="feature", how="left")
    df[["sample_code", "healing_status", "SRRS_no_ligand_panel_existing"]].to_csv(
        OUT_DIR / "GSE165816_ligand_panel_excluded_existing_sample_scores.tsv", sep="\t", index=False
    )
    tests.to_csv(OUT_DIR / "GSE165816_ligand_panel_excluded_existing_tests.tsv", sep="\t", index=False)
    lod.to_csv(OUT_DIR / "GSE165816_ligand_panel_excluded_existing_lodo.tsv", sep="\t", index=False)
    return tests


def annotate_gse223964(obs):
    marker_cols = [f"score_marker_{x}" for x in MARKER_SETS]
    cluster_scores = obs.groupby("leiden", observed=True)[marker_cols].mean()
    cluster_scores["assigned_broad_cell_type"] = [
        col.replace("score_marker_", "") for col in cluster_scores.idxmax(axis=1)
    ]
    obs["broad_cell_type"] = obs["leiden"].map(cluster_scores["assigned_broad_cell_type"].to_dict()).astype(str)
    return obs


def external_gse223964_ligand_excluded():
    sets, le_sets = load_gene_sets()
    adata = ad.read_h5ad(GSE223964_H5AD)
    adata.var_names_make_unique()
    x = get_matrix(adata)
    obs = adata.obs.copy()
    for name, genes in MARKER_SETS.items():
        obs[f"score_marker_{name}"], _ = score_gene_set(x, adata.var_names, genes)
    obs = annotate_gse223964(obs)

    for comp, genes in le_sets.items():
        obs[f"score_{comp}"], _ = score_gene_set(x, adata.var_names, genes)
    fib_mask = obs["broad_cell_type"].eq("fibroblast_stromal").to_numpy()
    fib_obs = obs.loc[fib_mask].copy()
    fib_obs["SRRS_ligand_excluded"] = np.nanmean(
        np.vstack([
            z(fib_obs["score_repair_no_ligands"].to_numpy()),
            z(fib_obs["score_d7_no_ligands"].to_numpy()),
            z(fib_obs["score_receiver_no_ligands"].to_numpy()),
        ]).T,
        axis=1,
    )
    q25 = np.nanquantile(fib_obs["SRRS_ligand_excluded"], 0.25)
    q75 = np.nanquantile(fib_obs["SRRS_ligand_excluded"], 0.75)
    fib_obs["SRRS_ligand_excluded_state"] = np.where(
        fib_obs["SRRS_ligand_excluded"] >= q75,
        "SRRS_high",
        np.where(fib_obs["SRRS_ligand_excluded"] <= q25, "SRRS_low", "intermediate"),
    )

    present_ligands = [g for g in LIGANDS if g in adata.var_names]
    lig_idx = [adata.var_names.get_loc(g) for g in present_ligands]
    lig_expr = pd.DataFrame(
        x[fib_mask, :][:, lig_idx].toarray().astype(float),
        columns=present_ligands,
        index=fib_obs.index,
    )
    lig_expr["candidate_ligand_panel"] = lig_expr.mean(axis=1)
    for col in lig_expr.columns:
        fib_obs[f"expr_{col}"] = lig_expr[col].to_numpy()

    paired_rows = []
    for (condition, sample), sub in fib_obs.groupby(["condition", "sample"], observed=True):
        hi = sub[sub["SRRS_ligand_excluded_state"] == "SRRS_high"]
        lo = sub[sub["SRRS_ligand_excluded_state"] == "SRRS_low"]
        if len(hi) < 5 or len(lo) < 5:
            continue
        for gene in present_ligands + ["candidate_ligand_panel"]:
            col = f"expr_{gene}"
            paired_rows.append({
                "condition": condition,
                "sample": sample,
                "feature": gene,
                "n_high": len(hi),
                "n_low": len(lo),
                "mean_high": hi[col].mean(),
                "mean_low": lo[col].mean(),
                "delta_high_minus_low": hi[col].mean() - lo[col].mean(),
            })
    paired = pd.DataFrame(paired_rows)
    paired.to_csv(OUT_DIR / "GSE223964_ligand_excluded_SRRS_high_low_ligand_deltas.tsv", sep="\t", index=False)

    def wilcoxon_paired_delta(deltas):
        deltas = np.asarray(deltas, dtype=float)
        deltas = deltas[np.isfinite(deltas)]
        if len(deltas) < 3:
            return math.nan
        try:
            return stats.wilcoxon(deltas).pvalue
        except Exception:
            return math.nan

    paired_tests = (
        paired.groupby("feature", observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_delta_high_minus_low=("delta_high_minus_low", "mean"),
            median_delta_high_minus_low=("delta_high_minus_low", "median"),
            positive_samples=("delta_high_minus_low", lambda x: int(np.sum(np.asarray(x) > 0))),
            p_paired_wilcoxon=("delta_high_minus_low", wilcoxon_paired_delta),
        )
        .reset_index()
    )
    paired_tests["padj_bh"] = bh_adjust(paired_tests["p_paired_wilcoxon"])
    paired_tests.to_csv(OUT_DIR / "GSE223964_ligand_excluded_SRRS_high_low_ligand_tests.tsv", sep="\t", index=False)

    # Random panel control using ligand-excluded SRRS states.
    expressed_frac = np.asarray((x[fib_mask, :] > 0).mean(axis=0)).ravel()
    pool = np.asarray(adata.var_names[expressed_frac >= 0.02])
    pool = np.asarray([g for g in pool if g not in set(LIGANDS) and not (g.startswith("MT-") or g.startswith("RPL") or g.startswith("RPS"))])
    observed = paired_tests.loc[paired_tests["feature"] == "candidate_ligand_panel", "mean_delta_high_minus_low"].iloc[0]
    rng = np.random.default_rng(20260522)
    state = fib_obs[["condition", "sample", "SRRS_ligand_excluded_state"]].copy()
    null = []
    for _ in range(1000):
        genes = rng.choice(pool, size=len(present_ligands), replace=False)
        idx = [adata.var_names.get_loc(g) for g in genes]
        arr = x[fib_mask, :][:, idx].toarray().astype(float).mean(axis=1)
        tmp = state.copy()
        tmp["panel"] = arr
        deltas = []
        for _, sub in tmp.groupby(["condition", "sample"], observed=True):
            hi = sub.loc[sub["SRRS_ligand_excluded_state"] == "SRRS_high", "panel"]
            lo = sub.loc[sub["SRRS_ligand_excluded_state"] == "SRRS_low", "panel"]
            if len(hi) >= 5 and len(lo) >= 5:
                deltas.append(hi.mean() - lo.mean())
        null.append(np.mean(deltas))
    null = np.asarray(null)
    pd.DataFrame({
        "observed_candidate_ligand_within_sample_delta": [observed],
        "random_panel_n": [len(null)],
        "random_mean_delta": [np.nanmean(null)],
        "random_sd_delta": [np.nanstd(null, ddof=1)],
        "empirical_p_random_ge_observed": [(np.sum(null >= observed) + 1) / (len(null) + 1)],
    }).to_csv(OUT_DIR / "GSE223964_ligand_excluded_SRRS_random_panel_control.tsv", sep="\t", index=False)
    return paired_tests


def spatial_ligand_excluded():
    spots = pd.read_csv(SPATIAL_SPOT_FILE, sep="\t")
    # Use existing spatial module scores but exclude the explicit ligand panel.
    components = ["fibroblast_repair_activation", "d7_acute_wound_alignment", "vascular_perivascular_receiver_coupling"]
    for c in components:
        spots[f"z_le_{c}"] = z(spots[c].to_numpy())
    spots["spatial_SRRS_ligand_panel_excluded"] = spots[[f"z_le_{c}" for c in components]].mean(axis=1)
    sample = (
        spots.groupby(["Donor", "Condition", "Sample_name"], observed=True)
        .agg(
            spatial_SRRS_ligand_panel_excluded=("spatial_SRRS_ligand_panel_excluded", "mean"),
            spots=("spatial_SRRS_ligand_panel_excluded", "size"),
        )
        .reset_index()
    )
    sample.to_csv(OUT_DIR / "GSE241124_spatial_ligand_panel_excluded_sample_scores.tsv", sep="\t", index=False)
    rows = []
    for cond in ["Wound1", "Wound7", "Wound30"]:
        wide = sample.pivot_table(index="Donor", columns="Condition", values="spatial_SRRS_ligand_panel_excluded", aggfunc="mean")
        dat = wide[["Skin", cond]].dropna()
        diffs = dat[cond] - dat["Skin"]
        rng = np.random.default_rng(1)
        boots = np.array([np.mean(rng.choice(diffs, size=len(diffs), replace=True)) for _ in range(10000)])
        try:
            p = stats.wilcoxon(diffs).pvalue if len(diffs) >= 3 else math.nan
        except Exception:
            p = math.nan
        rows.append({
            "comparison": f"{cond}_vs_Skin",
            "n_donors": len(diffs),
            "mean_delta": float(np.mean(diffs)),
            "bootstrap_ci_low": float(np.quantile(boots, 0.025)),
            "bootstrap_ci_high": float(np.quantile(boots, 0.975)),
            "positive_donors": int(np.sum(diffs > 0)),
            "p_paired_wilcoxon": p,
        })
    tests = pd.DataFrame(rows)
    tests.to_csv(OUT_DIR / "GSE241124_spatial_ligand_panel_excluded_donor_tests.tsv", sep="\t", index=False)
    return tests


def figures(discovery_tests, external_tests, spatial_tests):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), constrained_layout=True)
    ax = axes[0]
    d = discovery_tests.copy()
    ax.errorbar(
        d["delta_healer_minus_nonhealer"],
        [0],
        xerr=[
            d["delta_healer_minus_nonhealer"] - d["bootstrap_ci_low"],
            d["bootstrap_ci_high"] - d["delta_healer_minus_nonhealer"],
        ],
        fmt="o",
        color="#0072B2",
        ecolor="#777777",
        capsize=2,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks([0])
    ax.set_yticklabels(["Discovery\nno ligand panel"])
    ax.set_xlabel("Healer minus non-healer delta")
    ax.set_title("GSE165816")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    e = external_tests.sort_values("mean_delta_high_minus_low")
    colors = ["#D55E00" if f == "candidate_ligand_panel" else "#0072B2" for f in e["feature"]]
    ax.barh(e["feature"], e["mean_delta_high_minus_low"], color=colors)
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_xlabel("SRRS-high minus SRRS-low delta")
    ax.set_title("GSE223964 ligand-excluded SRRS")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[2]
    s = spatial_tests.copy()
    y = np.arange(len(s))
    ax.errorbar(
        s["mean_delta"],
        y,
        xerr=[s["mean_delta"] - s["bootstrap_ci_low"], s["bootstrap_ci_high"] - s["mean_delta"]],
        fmt="o",
        color="#009E73",
        ecolor="#777777",
        capsize=2,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(s["comparison"])
    ax.set_xlabel("Spatial delta vs skin")
    ax.set_title("GSE241124 no ligand panel")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(FIG_DIR / "ligand_excluded_SRRS_sensitivity_summary.png", dpi=300)
    fig.savefig(FIG_DIR / "ligand_excluded_SRRS_sensitivity_summary.pdf")
    plt.close(fig)


def main():
    discovery_tests = discovery_ligand_excluded()
    external_tests = external_gse223964_ligand_excluded()
    spatial_tests = spatial_ligand_excluded()
    figures(discovery_tests, external_tests, spatial_tests)
    print("Wrote ligand-excluded SRRS sensitivity outputs")


if __name__ == "__main__":
    main()
