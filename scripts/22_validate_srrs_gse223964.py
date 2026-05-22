from pathlib import Path
import math

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
import matplotlib.pyplot as plt


PROJECT = Path.cwd()
H5AD = PROJECT / "data" / "raw" / "GSE223964" / "GSE223964_Integrated_all_cells.h5ad"
GENE_SET_FILE = PROJECT / "results" / "srrs_framework" / "SRRS_locked_gene_sets.tsv"
OUT_DIR = PROJECT / "results" / "srrs_framework" / "external_gse223964"
FIG_DIR = PROJECT / "figures" / "srrs_framework" / "external_gse223964"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


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

SRRS_COMPONENTS = [
    "fibroblast_repair_activation",
    "stromal_ligand_panel",
    "vascular_perivascular_receiver_coupling",
    "d7_acute_wound_alignment",
]
FIB_COMPONENTS = [
    "fibroblast_repair_activation",
    "stromal_ligand_panel",
    "d7_acute_wound_alignment",
]
RECEIVER_TYPES = {"endothelial", "pericyte_smc"}
KEY_LIGANDS = ["TNC", "IL11", "IL6", "INHBA", "THBS1", "SERPINE1", "CCL20", "WNT5A", "ADAM12", "PTGS2"]


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
    adj = np.clip(adj, 0, 1)
    tmp = np.empty_like(adj)
    tmp[order] = adj
    out[ok] = tmp
    return out


def get_matrix(adata):
    x = adata.X
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(np.asarray(x))


def score_gene_set(x, var_names, genes):
    present = [g for g in genes if g in var_names]
    if len(present) < 3:
        return np.full(x.shape[0], np.nan), present
    idx = [var_names.get_loc(g) for g in present]
    arr = x[:, idx].toarray().astype(float)
    means = arr.mean(axis=0)
    stds = arr.std(axis=0)
    stds[stds == 0] = 1.0
    return ((arr - means) / stds).mean(axis=1), present


def z_sample(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.full_like(x, np.nan)
    return (x - np.nanmean(x)) / sd


def mannwhitney(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return math.nan
    return stats.mannwhitneyu(a, b, alternative="two-sided").pvalue


def read_gene_sets():
    tbl = pd.read_csv(GENE_SET_FILE, sep="\t")
    return {
        mod: sub["gene"].dropna().astype(str).tolist()
        for mod, sub in tbl.groupby("module", sort=False)
    }


def annotate_cell_types(obs):
    marker_cols = [f"score_marker_{x}" for x in MARKER_SETS]
    cluster_scores = obs.groupby("leiden", observed=True)[marker_cols].mean()
    cluster_scores["assigned_broad_cell_type"] = [
        col.replace("score_marker_", "") for col in cluster_scores.idxmax(axis=1)
    ]
    cluster_scores["assigned_score"] = cluster_scores[marker_cols].max(axis=1)
    mapping = cluster_scores["assigned_broad_cell_type"].to_dict()
    obs["broad_cell_type"] = obs["leiden"].map(mapping).astype(str)
    cluster_scores.reset_index().to_csv(
        OUT_DIR / "GSE223964_SRRS_leiden_marker_annotation.tsv",
        sep="\t",
        index=False,
    )
    return obs


def test_condition(sample_scores, feature):
    dat = sample_scores[["condition", feature]].dropna()
    diabetic = dat.loc[dat["condition"] == "diabetic", feature]
    nondiabetic = dat.loc[dat["condition"] == "non diabetic", feature]
    return {
        "feature": feature,
        "n_diabetic": len(diabetic),
        "n_non_diabetic": len(nondiabetic),
        "mean_diabetic": float(np.mean(diabetic)) if len(diabetic) else math.nan,
        "mean_non_diabetic": float(np.mean(nondiabetic)) if len(nondiabetic) else math.nan,
        "delta_diabetic_minus_non_diabetic": (
            float(np.mean(diabetic) - np.mean(nondiabetic)) if len(diabetic) and len(nondiabetic) else math.nan
        ),
        "p_mannwhitney": mannwhitney(diabetic, nondiabetic),
    }


def main():
    print(f"Reading {H5AD}")
    gene_sets = read_gene_sets()
    adata = ad.read_h5ad(H5AD)
    adata.var_names_make_unique()
    x = get_matrix(adata)
    obs = adata.obs.copy()

    presence_rows = []
    for name, genes in MARKER_SETS.items():
        score, present = score_gene_set(x, adata.var_names, genes)
        obs[f"score_marker_{name}"] = score
        presence_rows.append({
            "set_name": f"marker_{name}",
            "n_requested": len(genes),
            "n_present": len(present),
            "present_genes": ";".join(present),
        })

    obs = annotate_cell_types(obs)

    for component in SRRS_COMPONENTS + ["negative_housekeeping_control"]:
        score, present = score_gene_set(x, adata.var_names, gene_sets.get(component, []))
        obs[f"score_{component}"] = score
        presence_rows.append({
            "set_name": component,
            "n_requested": len(gene_sets.get(component, [])),
            "n_present": len(present),
            "present_genes": ";".join(present),
        })
    pd.DataFrame(presence_rows).to_csv(OUT_DIR / "GSE223964_SRRS_gene_presence.tsv", sep="\t", index=False)

    obs["fibroblast_cell_SRRS"] = np.nan
    fib_mask = obs["broad_cell_type"].eq("fibroblast_stromal").to_numpy()
    if fib_mask.sum() > 0:
        fib_mat = np.vstack([z_sample(obs.loc[fib_mask, f"score_{c}"].to_numpy()) for c in FIB_COMPONENTS]).T
        obs.loc[fib_mask, "fibroblast_cell_SRRS"] = np.nanmean(fib_mat, axis=1)

    sample_rows = []
    for (condition, sample), sub in obs.groupby(["condition", "sample"], observed=True):
        row = {"condition": condition, "sample": sample, "total_cells": len(sub)}
        fib = sub[sub["broad_cell_type"] == "fibroblast_stromal"]
        recv = sub[sub["broad_cell_type"].isin(RECEIVER_TYPES)]
        row["fibroblast_cell_n"] = len(fib)
        row["receiver_cell_n"] = len(recv)
        for component in FIB_COMPONENTS:
            row[component] = float(np.nanmean(fib[f"score_{component}"])) if len(fib) else math.nan
        row["vascular_perivascular_receiver_coupling"] = (
            float(np.nanmean(recv["score_vascular_perivascular_receiver_coupling"])) if len(recv) else math.nan
        )
        row["negative_housekeeping_control"] = float(np.nanmean(sub["score_negative_housekeeping_control"]))
        row["fibroblast_cell_SRRS_mean"] = float(np.nanmean(fib["fibroblast_cell_SRRS"])) if len(fib) else math.nan
        sample_rows.append(row)
    sample_scores = pd.DataFrame(sample_rows).sort_values(["condition", "sample"])

    for component in SRRS_COMPONENTS:
        sample_scores[f"z_{component}"] = z_sample(sample_scores[component])
    sample_scores["SRRS_projected"] = sample_scores[[f"z_{c}" for c in SRRS_COMPONENTS]].mean(axis=1, skipna=True)
    sample_scores["negative_housekeeping_control_z"] = z_sample(sample_scores["negative_housekeeping_control"])
    sample_scores.to_csv(OUT_DIR / "GSE223964_SRRS_projected_sample_scores.tsv", sep="\t", index=False)

    test_features = SRRS_COMPONENTS + ["SRRS_projected", "fibroblast_cell_SRRS_mean", "negative_housekeeping_control_z"]
    condition_tests = pd.DataFrame([test_condition(sample_scores, f) for f in test_features])
    condition_tests["padj_bh"] = bh_adjust(condition_tests["p_mannwhitney"])
    condition_tests.to_csv(OUT_DIR / "GSE223964_SRRS_diabetic_vs_nondiabetic_tests.tsv", sep="\t", index=False)

    fib_obs = obs[obs["broad_cell_type"] == "fibroblast_stromal"].copy()
    q25 = np.nanquantile(fib_obs["fibroblast_cell_SRRS"], 0.25)
    q75 = np.nanquantile(fib_obs["fibroblast_cell_SRRS"], 0.75)
    fib_obs["SRRS_state"] = np.where(
        fib_obs["fibroblast_cell_SRRS"] >= q75,
        "SRRS_high",
        np.where(fib_obs["fibroblast_cell_SRRS"] <= q25, "SRRS_low", "intermediate"),
    )

    gene_rows = []
    for gene in KEY_LIGANDS:
        if gene not in adata.var_names:
            continue
        idx = adata.var_names.get_loc(gene)
        vals = x[fib_mask, idx].toarray().ravel().astype(float)
        tmp = fib_obs[["SRRS_state"]].copy()
        tmp["expr"] = vals
        hi = tmp.loc[tmp["SRRS_state"] == "SRRS_high", "expr"].to_numpy()
        lo = tmp.loc[tmp["SRRS_state"] == "SRRS_low", "expr"].to_numpy()
        gene_rows.append({
            "gene": gene,
            "mean_SRRS_high": float(np.mean(hi)) if len(hi) else math.nan,
            "mean_SRRS_low": float(np.mean(lo)) if len(lo) else math.nan,
            "delta_high_minus_low": float(np.mean(hi) - np.mean(lo)) if len(hi) and len(lo) else math.nan,
            "p_mannwhitney": mannwhitney(hi, lo),
            "pct_high_expr": float(np.mean(hi > 0)) if len(hi) else math.nan,
            "pct_low_expr": float(np.mean(lo > 0)) if len(lo) else math.nan,
        })
    ligand_tests = pd.DataFrame(gene_rows)
    ligand_tests["padj_bh"] = bh_adjust(ligand_tests["p_mannwhitney"])
    ligand_tests.to_csv(OUT_DIR / "GSE223964_SRRS_high_vs_low_fibroblast_ligand_tests.tsv", sep="\t", index=False)

    rng = np.random.default_rng(1)
    expressed = np.asarray(adata.var_names[(np.asarray((x[fib_mask, :] > 0).mean(axis=0)).ravel()) >= 0.02])
    key_present = [g for g in KEY_LIGANDS if g in expressed]
    observed = ligand_tests["delta_high_minus_low"].mean()
    random_deltas = []
    if len(key_present) >= 3 and len(expressed) > len(key_present):
        high_mask = fib_obs["SRRS_state"].eq("SRRS_high").to_numpy()
        low_mask = fib_obs["SRRS_state"].eq("SRRS_low").to_numpy()
        fib_x = x[fib_mask, :]
        for _ in range(1000):
            genes = rng.choice(expressed, size=len(key_present), replace=False)
            idx = [adata.var_names.get_loc(g) for g in genes]
            arr = fib_x[:, idx].toarray().astype(float)
            random_deltas.append(float(arr[high_mask, :].mean() - arr[low_mask, :].mean()))
    random_deltas = np.asarray(random_deltas, dtype=float)
    empirical_p = float((np.sum(random_deltas >= observed) + 1) / (len(random_deltas) + 1)) if len(random_deltas) else math.nan
    pd.DataFrame({
        "observed_candidate_ligand_delta": [observed],
        "random_panel_n": [len(random_deltas)],
        "random_mean_delta": [float(np.nanmean(random_deltas)) if len(random_deltas) else math.nan],
        "random_sd_delta": [float(np.nanstd(random_deltas, ddof=1)) if len(random_deltas) > 1 else math.nan],
        "empirical_p_random_ge_observed": [empirical_p],
    }).to_csv(OUT_DIR / "GSE223964_SRRS_random_ligand_panel_control.tsv", sep="\t", index=False)

    corr_rows = []
    if sample_scores["SRRS_projected"].notna().sum() >= 4:
        for feature in ["vascular_perivascular_receiver_coupling", "fibroblast_cell_SRRS_mean"]:
            dat = sample_scores[["SRRS_projected", feature]].dropna()
            if len(dat) >= 4:
                rho, p = stats.spearmanr(dat["SRRS_projected"], dat[feature])
                corr_rows.append({"comparison": f"SRRS_projected_vs_{feature}", "spearman_rho": rho, "p_spearman": p, "n_samples": len(dat)})
    pd.DataFrame(corr_rows).to_csv(OUT_DIR / "GSE223964_SRRS_receiver_coupling_correlations.tsv", sep="\t", index=False)

    colors = {"diabetic": "#D55E00", "non diabetic": "#0072B2"}
    plot_features = ["SRRS_projected"] + SRRS_COMPONENTS + ["negative_housekeeping_control_z"]
    fig, axes = plt.subplots(2, 3, figsize=(10, 6), constrained_layout=True)
    axes = axes.ravel()
    for ax, feature in zip(axes, plot_features):
        vals = []
        labels = []
        cols = []
        for cond in ["non diabetic", "diabetic"]:
            y = sample_scores.loc[sample_scores["condition"] == cond, feature].dropna().to_numpy()
            vals.append(y)
            labels.append(cond)
            cols.append(colors[cond])
        bp = ax.boxplot(vals, patch_artist=True, widths=0.5, showfliers=False)
        for patch, c in zip(bp["boxes"], cols):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        for i, (y, c) in enumerate(zip(vals, cols), start=1):
            jitter = np.linspace(-0.07, 0.07, len(y)) if len(y) > 1 else np.array([0])
            ax.scatter(np.full(len(y), i) + jitter, y, color=c, s=26, alpha=0.9)
        ax.set_title(feature, fontsize=9)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
        ax.axhline(0, color="#999999", linewidth=0.6)
    fig.savefig(FIG_DIR / "GSE223964_SRRS_projected_sample_scores.png", dpi=240)
    plt.close(fig)

    if not ligand_tests.empty:
        ligand_tests_sorted = ligand_tests.sort_values("delta_high_minus_low")
        fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
        ax.barh(ligand_tests_sorted["gene"], ligand_tests_sorted["delta_high_minus_low"], color="#0072B2")
        ax.axvline(0, color="#777777", linewidth=0.7)
        ax.set_xlabel("Expression delta, SRRS-high minus SRRS-low fibroblasts")
        ax.set_ylabel("")
        fig.savefig(FIG_DIR / "GSE223964_SRRS_high_fibroblast_candidate_ligands.png", dpi=240)
        plt.close(fig)

    print("Wrote GSE223964 SRRS validation outputs")


if __name__ == "__main__":
    main()
