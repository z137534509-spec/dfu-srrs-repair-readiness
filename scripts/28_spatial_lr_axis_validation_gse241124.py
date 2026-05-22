from pathlib import Path
import math
import re
import zipfile

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt


PROJECT = Path.cwd()
RAW_DIR = PROJECT / "data" / "raw" / "GSE241124"
EXTRACTED = RAW_DIR / "extracted"
META_FILE = RAW_DIR / "GSE241124_spatialseq_metadata_acutewound.txt.gz"
OUT_DIR = PROJECT / "results" / "srrs_framework" / "spatial_lr_axis_gse241124"
FIG_DIR = PROJECT / "figures" / "srrs_framework" / "spatial_lr_axis_gse241124"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

AXES = {
    "TNC_integrin_syndecan": {
        "ligands": ["TNC"],
        "receptors": ["ITGB1", "SDC4", "ITGA5", "ITGAV", "ITGA8"],
    },
    "IL11_gp130": {
        "ligands": ["IL11"],
        "receptors": ["IL11RA", "IL6ST"],
    },
    "IL6_gp130": {
        "ligands": ["IL6"],
        "receptors": ["IL6R", "IL6ST"],
    },
    "INHBA_activin_TGF_family": {
        "ligands": ["INHBA"],
        "receptors": ["ENG", "ACVR1", "ACVR1B", "ACVR2A", "ACVR2B", "TGFBR3", "BAMBI"],
    },
    "THBS1_matrix_receptors": {
        "ligands": ["THBS1"],
        "receptors": ["ITGB1", "SDC4", "CD47", "LRP1"],
    },
    "SERPINE1_LRP1": {
        "ligands": ["SERPINE1"],
        "receptors": ["LRP1"],
    },
    "CCL20_CCR6": {
        "ligands": ["CCL20"],
        "receptors": ["CCR6"],
    },
}
CONDITIONS = ["Skin", "Wound1", "Wound7", "Wound30"]


def key_from_file(path, suffix):
    name = re.sub(r"^GSM\d+_", "", path.name)
    return name.replace(suffix, "")


def read_10x_h5_gene_subset(path, wanted_genes):
    with h5py.File(path, "r") as f:
        grp = f["matrix"]
        barcodes = np.array([x.decode("utf-8") for x in grp["barcodes"][:]])
        genes = np.array([x.decode("utf-8") for x in grp["features/name"][:]])
        shape = tuple(grp["shape"][:])
        mat = sparse.csc_matrix(
            (grp["data"][:], grp["indices"][:], grp["indptr"][:]),
            shape=shape,
        ).tocsr()
    gene_to_indices = {}
    for i, gene in enumerate(genes):
        if gene in wanted_genes:
            gene_to_indices.setdefault(gene, []).append(i)
    present = sorted(gene_to_indices)
    if not present:
        return barcodes, present, sparse.csr_matrix((len(barcodes), 0))
    cols = []
    for gene in present:
        cols.append(np.asarray(mat[gene_to_indices[gene], :].sum(axis=0)).ravel())
    return barcodes, present, sparse.csr_matrix(np.vstack(cols).T)


def read_positions(zip_path):
    with zipfile.ZipFile(zip_path) as zf:
        candidates = [n for n in zf.namelist() if n.endswith("tissue_positions_list.csv") or n.endswith("tissue_positions.csv")]
        if not candidates:
            return pd.DataFrame()
        with zf.open(candidates[0]) as fh:
            pos = pd.read_csv(fh, header=None)
    pos = pos.iloc[:, :6]
    pos.columns = ["raw_barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    return pos


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


def bootstrap_mean(vals, n_boot=5000, seed=1):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return math.nan, math.nan, math.nan
    rng = np.random.default_rng(seed)
    boots = np.array([np.mean(rng.choice(vals, size=len(vals), replace=True)) for _ in range(n_boot)])
    return float(np.mean(vals)), float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))


def main():
    wanted_genes = sorted(set(g for ax in AXES.values() for side in ["ligands", "receptors"] for g in ax[side]))
    meta = pd.read_csv(META_FILE, sep="\t")
    meta["raw_barcode"] = meta["barcode"].str.extract(r"([ACGT]+-1)$")

    h5_files = list(EXTRACTED.glob("*_filtered_feature_bc_matrix.h5"))
    zip_files = list(EXTRACTED.glob("*_spatial_images.zip"))
    h5_by_key = {key_from_file(p, "_filtered_feature_bc_matrix.h5"): p for p in h5_files}
    zip_by_key = {key_from_file(p, "_spatial_images.zip"): p for p in zip_files}

    blocks = []
    presence_rows = []
    for (orig_ident, condition, sample_name, donor), sub_meta in meta.groupby(["orig.ident", "Condition", "Sample_name", "Donor"], observed=True):
        h5_path = h5_by_key.get(orig_ident) or h5_by_key.get(condition)
        zip_path = zip_by_key.get(orig_ident) or zip_by_key.get(condition)
        if h5_path is None:
            continue
        barcodes, present, mat = read_10x_h5_gene_subset(h5_path, wanted_genes)
        expr = pd.DataFrame(mat.toarray(), index=barcodes, columns=present)
        sub = sub_meta.copy().set_index("raw_barcode")
        common = expr.index.intersection(sub.index)
        if len(common) == 0:
            continue
        expr = expr.loc[common]
        sub = sub.loc[common].copy()
        sub.index.name = "raw_barcode"
        sub = sub.reset_index()
        lib = expr.sum(axis=1).replace(0, np.nan)
        log_expr = np.log1p(expr.div(lib, axis=0) * 10000.0).fillna(0.0)
        block = pd.concat([sub.reset_index(drop=True), log_expr.add_prefix("gene__").reset_index(drop=True)], axis=1)
        if zip_path is not None:
            pos = read_positions(zip_path)
            if not pos.empty:
                block = block.merge(pos, on="raw_barcode", how="left")
        blocks.append(block)
        presence_rows.append({
            "Sample_name": sample_name,
            "orig.ident": orig_ident,
            "Condition": condition,
            "Donor": donor,
            "spots": len(block),
            "n_wanted_genes": len(wanted_genes),
            "n_present_genes": len(present),
            "present_genes": ";".join(present),
        })

    spots = pd.concat(blocks, ignore_index=True)
    gene_cols = [c for c in spots.columns if c.startswith("gene__")]
    gene_names = [c.replace("gene__", "", 1) for c in gene_cols]
    expr = spots[gene_cols].copy()
    expr.columns = gene_names
    expr_z = expr.apply(z, axis=0)

    for axis, spec in AXES.items():
        ligs = [g for g in spec["ligands"] if g in expr_z.columns]
        recs = [g for g in spec["receptors"] if g in expr_z.columns]
        spots[f"{axis}_ligand_score"] = expr_z[ligs].mean(axis=1) if ligs else np.nan
        spots[f"{axis}_receptor_score"] = expr_z[recs].mean(axis=1) if recs else np.nan
        spots[f"{axis}_axis_score"] = np.nanmean(
            np.vstack([spots[f"{axis}_ligand_score"], spots[f"{axis}_receptor_score"]]).T,
            axis=1,
        )
        spots[f"{axis}_ligand_n"] = len(ligs)
        spots[f"{axis}_receptor_n"] = len(recs)

    keep_cols = [
        "barcode", "raw_barcode", "orig.ident", "Condition", "Sample_name", "Donor", "AnnoType",
        "nCount_Spatial", "nFeature_Spatial", "pxl_row", "pxl_col",
    ]
    score_cols = [c for c in spots.columns if c.endswith("_ligand_score") or c.endswith("_receptor_score") or c.endswith("_axis_score")]
    spots[keep_cols + score_cols].to_csv(OUT_DIR / "GSE241124_spatial_LR_axis_spot_scores.tsv", sep="\t", index=False)
    pd.DataFrame(presence_rows).to_csv(OUT_DIR / "GSE241124_spatial_LR_axis_gene_presence.tsv", sep="\t", index=False)

    sample_rows = []
    corr_rows = []
    for axis in AXES:
        axis_col = f"{axis}_axis_score"
        lig_col = f"{axis}_ligand_score"
        rec_col = f"{axis}_receptor_score"
        sample = (
            spots.groupby(["Donor", "Condition", "Sample_name"], observed=True)[[axis_col, lig_col, rec_col]]
            .mean()
            .reset_index()
        )
        sample["axis"] = axis
        sample_rows.append(sample)
        for (condition, donor), sub in spots.groupby(["Condition", "Donor"], observed=True):
            dat = sub[[lig_col, rec_col, axis_col]].dropna()
            if len(dat) < 20:
                continue
            rho, p = stats.spearmanr(dat[lig_col], dat[rec_col])
            corr_rows.append({
                "axis": axis,
                "Condition": condition,
                "Donor": donor,
                "ligand_receptor_spearman_rho": rho,
                "p_spearman": p,
                "n_spots": len(dat),
                "mean_axis_score": dat[axis_col].mean(),
            })
    sample_summary = pd.concat(sample_rows, ignore_index=True)
    sample_summary.to_csv(OUT_DIR / "GSE241124_spatial_LR_axis_sample_summary.tsv", sep="\t", index=False)
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(OUT_DIR / "GSE241124_spatial_LR_axis_ligand_receptor_correlations_by_donor.tsv", sep="\t", index=False)

    corr_summary_rows = []
    for (axis, condition), sub in corr.groupby(["axis", "Condition"], observed=True):
        mean_rho, lo, hi = bootstrap_mean(sub["ligand_receptor_spearman_rho"])
        corr_summary_rows.append({
            "axis": axis,
            "Condition": condition,
            "mean_donor_rho": mean_rho,
            "bootstrap_ci_low": lo,
            "bootstrap_ci_high": hi,
            "positive_donors": int((sub["ligand_receptor_spearman_rho"] > 0).sum()),
            "n_donors": sub["Donor"].nunique(),
        })
    corr_summary = pd.DataFrame(corr_summary_rows)
    corr_summary.to_csv(OUT_DIR / "GSE241124_spatial_LR_axis_correlation_donor_bootstrap_summary.tsv", sep="\t", index=False)

    # Cluster-robust condition models for axis scores.
    model_rows = []
    spots["log_nCount"] = np.log1p(spots["nCount_Spatial"])
    spots["log_nFeature"] = np.log1p(spots["nFeature_Spatial"])
    for axis in AXES:
        col = f"{axis}_axis_score"
        dat = spots.dropna(subset=[col, "Condition", "Sample_name", "log_nCount", "log_nFeature"]).copy()
        fit = smf.ols(f"{col} ~ C(Condition, Treatment(reference='Skin')) + log_nCount + log_nFeature", data=dat).fit(
            cov_type="cluster",
            cov_kwds={"groups": dat["Sample_name"]},
        )
        for term in fit.params.index:
            if "Condition" not in term:
                continue
            model_rows.append({
                "axis": axis,
                "term": term,
                "estimate": fit.params[term],
                "std_error_cluster_sample": fit.bse[term],
                "p_cluster_sample": fit.pvalues[term],
                "n_spots": int(fit.nobs),
                "n_clusters": dat["Sample_name"].nunique(),
            })
    models = pd.DataFrame(model_rows)
    models["padj_bh"] = bh_adjust(models["p_cluster_sample"])
    models.to_csv(OUT_DIR / "GSE241124_spatial_LR_axis_condition_clustered_models.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
    plot = corr_summary[corr_summary["Condition"].astype(str).isin(["Wound1", "Wound7", "Wound30"])].copy()
    plot = plot.sort_values("mean_donor_rho")
    labels = plot["Condition"].astype(str) + " " + plot["axis"]
    y = np.arange(len(plot))
    ax.errorbar(
        plot["mean_donor_rho"],
        y,
        xerr=[plot["mean_donor_rho"] - plot["bootstrap_ci_low"], plot["bootstrap_ci_high"] - plot["mean_donor_rho"]],
        fmt="o",
        color="#0072B2",
        ecolor="#777777",
        capsize=2,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Donor-bootstrapped ligand-receptor spatial Spearman rho")
    fig.savefig(FIG_DIR / "GSE241124_spatial_LR_axis_correlations.png", dpi=240)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.4), constrained_layout=True)
    w7 = models[models["term"].str.contains("Wound7")].sort_values("estimate")
    ax.barh(w7["axis"], w7["estimate"], color="#D55E00")
    ax.axvline(0, color="#777777", linewidth=0.7)
    ax.set_xlabel("Adjusted Wound7 vs Skin axis-score coefficient")
    fig.savefig(FIG_DIR / "GSE241124_spatial_LR_axis_Wound7_coefficients.png", dpi=240)
    plt.close(fig)

    print("Wrote GSE241124 spatial LR axis validation outputs")


if __name__ == "__main__":
    main()
