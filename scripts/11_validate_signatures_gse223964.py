from pathlib import Path
import math

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
import matplotlib.pyplot as plt


PROJECT = Path.cwd()
H5AD = PROJECT / "data" / "raw" / "GSE223964" / "GSE223964_Integrated_all_cells.h5ad"
OUT_DIR = PROJECT / "results" / "validation_gse223964"
FIG_DIR = PROJECT / "figures" / "validation_gse223964"
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


SIGNATURES = {
    "fibroblast_sender_ligand_shortlist": [
        "IL11", "CCL20", "INHBA", "SERPINE1", "IL6", "THBS1",
        "TNC", "WNT5A", "ADAM12", "PTGS2",
    ],
    "fibroblast_repair_activation": [
        "THBS1", "WISP1", "IL11", "INHBA", "HIF1A", "TNFRSF12A", "CCL20",
        "MMP10", "COL10A1", "MMP3", "SERPINE1", "IL6", "PTGS2", "WNT5A",
        "TNFAIP6", "PDPN", "ITGA5", "ADAM12", "PLOD2",
    ],
    "ecm_remodeling_migration": [
        "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
        "COL6A1", "COL6A2", "COL10A1", "THBS1", "TNC", "MMP1", "MMP3",
        "MMP10", "SERPINE1", "ITGA5", "ITGAV", "PDPN", "TNFAIP6", "CTHRC1",
    ],
    "repair_inflammatory_signaling": [
        "IL6", "IL11", "INHBA", "AIM2", "CCL20", "CXCL5", "CXCL8",
        "TNFAIP3", "TNFRSF12A", "PTGS2", "FPR2", "TNIP3", "SLC39A8",
    ],
    "resolution_metabolic_myeloid": [
        "FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8",
        "CES1", "CXCL5", "TNC", "COL4A1", "COL4A2",
    ],
}


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


def gene_score(adata, genes, score_name):
    present = [g for g in genes if g in adata.var_names]
    if len(present) == 0:
        return np.zeros(adata.n_obs, dtype=float), present
    idx = [adata.var_names.get_loc(g) for g in present]
    x = adata.X[:, idx]
    if sparse.issparse(x):
        arr = x.toarray()
    else:
        arr = np.asarray(x)
    arr = arr.astype(float, copy=False)
    means = arr.mean(axis=0)
    stds = arr.std(axis=0)
    stds[stds == 0] = 1.0
    z = (arr - means) / stds
    score = z.mean(axis=1)
    return score, present


def mannwhitney(diabetic, nondiabetic):
    if len(diabetic) < 2 or len(nondiabetic) < 2:
        return math.nan
    try:
        return stats.mannwhitneyu(diabetic, nondiabetic, alternative="two-sided").pvalue
    except Exception:
        return math.nan


def main():
    print(f"Reading {H5AD}")
    adata = ad.read_h5ad(H5AD)
    adata.var_names_make_unique()
    obs = adata.obs.copy()

    presence_rows = []
    for name, genes in {**MARKER_SETS, **SIGNATURES}.items():
        score, present = gene_score(adata, genes, name)
        obs[f"score_{name}"] = score
        presence_rows.append({
            "set_name": name,
            "n_requested": len(genes),
            "n_present": len(present),
            "present_genes": ";".join(present),
            "missing_genes": ";".join([g for g in genes if g not in present]),
        })

    presence = pd.DataFrame(presence_rows)
    presence.to_csv(OUT_DIR / "GSE223964_signature_gene_presence.tsv", sep="\t", index=False)

    marker_cols = [f"score_{x}" for x in MARKER_SETS]
    cluster_scores = obs.groupby("leiden", observed=True)[marker_cols].mean()
    cluster_scores["assigned_broad_cell_type"] = [
        col.replace("score_", "") for col in cluster_scores.idxmax(axis=1)
    ]
    cluster_scores["assigned_score"] = cluster_scores[marker_cols].max(axis=1)
    cluster_scores.reset_index().to_csv(
        OUT_DIR / "GSE223964_leiden_marker_annotation.tsv",
        sep="\t",
        index=False,
    )

    mapping = cluster_scores["assigned_broad_cell_type"].to_dict()
    obs["broad_cell_type"] = obs["leiden"].map(mapping).astype(str)
    obs["condition_clean"] = obs["condition"].astype(str).str.replace(" ", "_", regex=False)

    cell_counts = (
        obs.groupby(["condition", "sample", "leiden", "broad_cell_type"], observed=True)
        .size()
        .reset_index(name="cells")
    )
    cell_counts.to_csv(OUT_DIR / "GSE223964_cell_counts_by_sample_cluster_celltype.tsv", sep="\t", index=False)

    score_cols = [f"score_{x}" for x in SIGNATURES]
    sample_scores = (
        obs.groupby(["condition", "sample", "broad_cell_type"], observed=True)[score_cols]
        .mean()
        .reset_index()
    )
    sample_counts = (
        obs.groupby(["condition", "sample", "broad_cell_type"], observed=True)
        .size()
        .reset_index(name="cells")
    )
    sample_scores = sample_scores.merge(sample_counts, on=["condition", "sample", "broad_cell_type"], how="left")
    sample_scores.to_csv(OUT_DIR / "GSE223964_signature_scores_by_sample_celltype.tsv", sep="\t", index=False)

    test_rows = []
    for ct in sorted(sample_scores["broad_cell_type"].unique()):
        zct = sample_scores[sample_scores["broad_cell_type"] == ct]
        for sig in SIGNATURES:
            col = f"score_{sig}"
            diabetic = zct.loc[zct["condition"] == "diabetic", col].dropna().to_numpy()
            nondiabetic = zct.loc[zct["condition"] == "non diabetic", col].dropna().to_numpy()
            test_rows.append({
                "broad_cell_type": ct,
                "signature": sig,
                "n_diabetic": len(diabetic),
                "n_non_diabetic": len(nondiabetic),
                "mean_diabetic": float(np.mean(diabetic)) if len(diabetic) else math.nan,
                "mean_non_diabetic": float(np.mean(nondiabetic)) if len(nondiabetic) else math.nan,
                "delta_diabetic_minus_non_diabetic": (
                    float(np.mean(diabetic) - np.mean(nondiabetic))
                    if len(diabetic) and len(nondiabetic)
                    else math.nan
                ),
                "p_mannwhitney": mannwhitney(diabetic, nondiabetic),
            })
    tests = pd.DataFrame(test_rows)
    tests["padj_bh"] = bh_adjust(tests["p_mannwhitney"].to_numpy())
    tests = tests.sort_values(["p_mannwhitney", "broad_cell_type", "signature"], na_position="last")
    tests.to_csv(OUT_DIR / "GSE223964_signature_diabetic_vs_nondiabetic_tests.tsv", sep="\t", index=False)

    focus_cts = ["fibroblast_stromal", "myeloid", "endothelial", "keratinocyte", "pericyte_smc"]
    focus = sample_scores[sample_scores["broad_cell_type"].isin(focus_cts)].copy()
    if not focus.empty:
        sig_order = list(SIGNATURES.keys())
        fig, axes = plt.subplots(len(sig_order), 1, figsize=(8, 2.5 * len(sig_order)), constrained_layout=True)
        if len(sig_order) == 1:
            axes = [axes]
        colors = {"diabetic": "#D55E00", "non diabetic": "#0072B2"}
        for ax, sig in zip(axes, sig_order):
            col = f"score_{sig}"
            positions = []
            labels = []
            values = []
            box_colors = []
            pos = 1
            for ct in focus_cts:
                for cond in ["non diabetic", "diabetic"]:
                    vals = focus.loc[(focus["broad_cell_type"] == ct) & (focus["condition"] == cond), col].dropna().to_numpy()
                    if len(vals) > 0:
                        values.append(vals)
                        positions.append(pos)
                        labels.append(f"{ct}\n{cond}")
                        box_colors.append(colors[cond])
                        pos += 1
                pos += 0.5
            bp = ax.boxplot(values, positions=positions, patch_artist=True, widths=0.55, showfliers=False)
            for patch, c in zip(bp["boxes"], box_colors):
                patch.set_facecolor(c)
                patch.set_alpha(0.65)
            for p, vals, c in zip(positions, values, box_colors):
                jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) > 1 else np.array([0])
                ax.scatter(np.full(len(vals), p) + jitter, vals, color=c, s=18, alpha=0.9)
            ax.set_title(sig)
            ax.set_ylabel("mean z-score")
            ax.set_xticks(positions)
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
            ax.axhline(0, color="#999999", linewidth=0.6)
        fig.savefig(FIG_DIR / "GSE223964_signature_scores_by_sample_celltype.png", dpi=220)
        plt.close(fig)

    print("Wrote GSE223964 validation outputs")


if __name__ == "__main__":
    main()
