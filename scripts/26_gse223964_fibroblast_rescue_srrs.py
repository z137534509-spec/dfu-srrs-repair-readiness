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
OUT_DIR = PROJECT / "results" / "srrs_framework" / "external_gse223964_rescue"
FIG_DIR = PROJECT / "figures" / "srrs_framework" / "external_gse223964_rescue"
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
FIB_COMPONENTS = ["fibroblast_repair_activation", "stromal_ligand_panel", "d7_acute_wound_alignment"]
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


def z(x):
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


def wilcoxon_paired_delta(deltas):
    deltas = np.asarray(deltas, dtype=float)
    deltas = deltas[np.isfinite(deltas)]
    if len(deltas) < 3:
        return math.nan
    try:
        return stats.wilcoxon(deltas).pvalue
    except Exception:
        return math.nan


def annotate_cell_types(obs):
    marker_cols = [f"score_marker_{x}" for x in MARKER_SETS]
    cluster_scores = obs.groupby("leiden", observed=True)[marker_cols].mean()
    cluster_scores["assigned_broad_cell_type"] = [
        col.replace("score_marker_", "") for col in cluster_scores.idxmax(axis=1)
    ]
    obs["broad_cell_type"] = obs["leiden"].map(cluster_scores["assigned_broad_cell_type"].to_dict()).astype(str)
    cluster_scores.reset_index().to_csv(OUT_DIR / "GSE223964_rescue_leiden_annotation.tsv", sep="\t", index=False)
    return obs


def main():
    print("Reading GSE223964 h5ad")
    gene_sets = pd.read_csv(GENE_SET_FILE, sep="\t")
    sets = {
        module: sub["gene"].dropna().astype(str).tolist()
        for module, sub in gene_sets.groupby("module", sort=False)
    }
    adata = ad.read_h5ad(H5AD)
    adata.var_names_make_unique()
    x = get_matrix(adata)
    obs = adata.obs.copy()

    for name, genes in MARKER_SETS.items():
        obs[f"score_marker_{name}"], _ = score_gene_set(x, adata.var_names, genes)
    obs = annotate_cell_types(obs)

    for comp in FIB_COMPONENTS + ["negative_housekeeping_control"]:
        obs[f"score_{comp}"], _ = score_gene_set(x, adata.var_names, sets.get(comp, []))

    fib_mask = obs["broad_cell_type"].eq("fibroblast_stromal").to_numpy()
    fib_obs = obs.loc[fib_mask].copy()
    fib_mat = np.vstack([z(fib_obs[f"score_{c}"].to_numpy()) for c in FIB_COMPONENTS]).T
    fib_obs["fibroblast_SRRS"] = np.nanmean(fib_mat, axis=1)
    fib_obs["negative_housekeeping_control_z"] = z(fib_obs["score_negative_housekeeping_control"])
    q25 = np.nanquantile(fib_obs["fibroblast_SRRS"], 0.25)
    q75 = np.nanquantile(fib_obs["fibroblast_SRRS"], 0.75)
    fib_obs["SRRS_state"] = np.where(
        fib_obs["fibroblast_SRRS"] >= q75,
        "SRRS_high",
        np.where(fib_obs["fibroblast_SRRS"] <= q25, "SRRS_low", "intermediate"),
    )

    # Candidate-ligand expression matrix in fibroblasts.
    present_ligands = [g for g in KEY_LIGANDS if g in adata.var_names]
    lig_idx = [adata.var_names.get_loc(g) for g in present_ligands]
    lig_expr = pd.DataFrame(
        x[fib_mask, :][:, lig_idx].toarray().astype(float),
        columns=present_ligands,
        index=fib_obs.index,
    )
    lig_expr["candidate_ligand_panel"] = lig_expr.mean(axis=1)
    for col in lig_expr.columns:
        fib_obs[f"expr_{col}"] = lig_expr[col].to_numpy()

    sample_rows = []
    for (condition, sample), sub in fib_obs.groupby(["condition", "sample"], observed=True):
        row = {
            "condition": condition,
            "sample": sample,
            "fibroblast_cell_n": len(sub),
            "fibroblast_SRRS_mean": sub["fibroblast_SRRS"].mean(),
            "fibroblast_SRRS_high_fraction": (sub["SRRS_state"] == "SRRS_high").mean(),
            "fibroblast_SRRS_low_fraction": (sub["SRRS_state"] == "SRRS_low").mean(),
            "candidate_ligand_panel_mean": sub["expr_candidate_ligand_panel"].mean(),
            "negative_housekeeping_control_z": sub["negative_housekeeping_control_z"].mean(),
        }
        for comp in FIB_COMPONENTS:
            row[f"{comp}_mean"] = sub[f"score_{comp}"].mean()
        sample_rows.append(row)
    sample_scores = pd.DataFrame(sample_rows).sort_values(["condition", "sample"])
    sample_scores.to_csv(OUT_DIR / "GSE223964_fibroblast_specific_SRRS_sample_scores.tsv", sep="\t", index=False)

    test_features = [
        "fibroblast_SRRS_mean",
        "fibroblast_SRRS_high_fraction",
        "fibroblast_SRRS_low_fraction",
        "candidate_ligand_panel_mean",
        "negative_housekeeping_control_z",
    ] + [f"{c}_mean" for c in FIB_COMPONENTS]
    condition_rows = []
    for feat in test_features:
        d = sample_scores.loc[sample_scores["condition"] == "diabetic", feat]
        nd = sample_scores.loc[sample_scores["condition"] == "non diabetic", feat]
        condition_rows.append({
            "feature": feat,
            "n_diabetic": len(d),
            "n_non_diabetic": len(nd),
            "mean_diabetic": d.mean(),
            "mean_non_diabetic": nd.mean(),
            "delta_diabetic_minus_non_diabetic": d.mean() - nd.mean(),
            "p_mannwhitney": mannwhitney(d, nd),
        })
    condition_tests = pd.DataFrame(condition_rows)
    condition_tests["padj_bh"] = bh_adjust(condition_tests["p_mannwhitney"])
    condition_tests.to_csv(OUT_DIR / "GSE223964_fibroblast_specific_condition_tests.tsv", sep="\t", index=False)

    # Within-sample paired high-vs-low evidence avoids treating cells as independent samples.
    paired_rows = []
    for (condition, sample), sub in fib_obs.groupby(["condition", "sample"], observed=True):
        hi = sub[sub["SRRS_state"] == "SRRS_high"]
        lo = sub[sub["SRRS_state"] == "SRRS_low"]
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
    paired.to_csv(OUT_DIR / "GSE223964_SRRS_high_low_within_sample_ligand_deltas.tsv", sep="\t", index=False)
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
    paired_tests.to_csv(OUT_DIR / "GSE223964_SRRS_high_low_within_sample_ligand_tests.tsv", sep="\t", index=False)

    # Random panel control using the same within-sample high-low contrast.
    expressed_frac = np.asarray((x[fib_mask, :] > 0).mean(axis=0)).ravel()
    pool = np.asarray(adata.var_names[expressed_frac >= 0.02])
    pool = np.asarray([g for g in pool if not (g.startswith("MT-") or g.startswith("RPL") or g.startswith("RPS"))])
    observed = paired_tests.loc[paired_tests["feature"] == "candidate_ligand_panel", "mean_delta_high_minus_low"].iloc[0]
    rng = np.random.default_rng(20260522)
    null = []
    state = fib_obs[["condition", "sample", "SRRS_state"]].copy()
    for _ in range(1000):
        genes = rng.choice(pool, size=len(present_ligands), replace=False)
        idx = [adata.var_names.get_loc(g) for g in genes]
        arr = x[fib_mask, :][:, idx].toarray().astype(float).mean(axis=1)
        tmp = state.copy()
        tmp["panel"] = arr
        deltas = []
        for _, sub in tmp.groupby(["condition", "sample"], observed=True):
            hi = sub.loc[sub["SRRS_state"] == "SRRS_high", "panel"]
            lo = sub.loc[sub["SRRS_state"] == "SRRS_low", "panel"]
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
    }).to_csv(OUT_DIR / "GSE223964_SRRS_within_sample_random_panel_control.tsv", sep="\t", index=False)

    # Figures.
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    for ax, feat in zip(axes, ["fibroblast_SRRS_mean", "fibroblast_SRRS_high_fraction", "candidate_ligand_panel_mean"]):
        vals = [sample_scores.loc[sample_scores["condition"] == c, feat].dropna().to_numpy() for c in ["non diabetic", "diabetic"]]
        bp = ax.boxplot(vals, patch_artist=True, widths=0.5, showfliers=False)
        for patch, color in zip(bp["boxes"], ["#0072B2", "#D55E00"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.55)
        for i, (y, color) in enumerate(zip(vals, ["#0072B2", "#D55E00"]), start=1):
            jitter = np.linspace(-0.06, 0.06, len(y)) if len(y) > 1 else np.array([0])
            ax.scatter(np.full(len(y), i) + jitter, y, color=color, s=30, alpha=0.9)
        ax.set_title(feat, fontsize=9)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["non diabetic", "diabetic"], rotation=25, ha="right", fontsize=8)
    fig.savefig(FIG_DIR / "GSE223964_fibroblast_specific_sample_metrics.png", dpi=240)
    plt.close(fig)

    plot_paired = paired_tests.sort_values("mean_delta_high_minus_low")
    fig, ax = plt.subplots(figsize=(6.8, 4.5), constrained_layout=True)
    ax.barh(plot_paired["feature"], plot_paired["mean_delta_high_minus_low"], color="#0072B2")
    ax.axvline(0, color="#777777", linewidth=0.7)
    ax.set_xlabel("Mean within-sample delta, SRRS-high minus SRRS-low fibroblasts")
    fig.savefig(FIG_DIR / "GSE223964_SRRS_high_low_paired_ligand_tests.png", dpi=240)
    plt.close(fig)

    if "X_umap" in adata.obsm:
        um = pd.DataFrame(adata.obsm["X_umap"][fib_mask, :], columns=["UMAP_1", "UMAP_2"], index=fib_obs.index)
        um = pd.concat([um, fib_obs[["condition", "SRRS_state", "fibroblast_SRRS"]]], axis=1)
        fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.8), constrained_layout=True)
        sc = axes[0].scatter(um["UMAP_1"], um["UMAP_2"], c=um["fibroblast_SRRS"], s=2, cmap="viridis")
        axes[0].set_title("Fibroblast SRRS")
        fig.colorbar(sc, ax=axes[0], shrink=0.75)
        colors = {"SRRS_high": "#D55E00", "intermediate": "#BDBDBD", "SRRS_low": "#0072B2"}
        axes[1].scatter(um["UMAP_1"], um["UMAP_2"], c=um["SRRS_state"].map(colors), s=2, alpha=0.85)
        axes[1].set_title("SRRS state")
        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])
        fig.savefig(FIG_DIR / "GSE223964_fibroblast_SRRS_umap.png", dpi=240)
        plt.close(fig)

    print("Wrote GSE223964 fibroblast-specific rescue outputs")


if __name__ == "__main__":
    main()
