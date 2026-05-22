from pathlib import Path
import math
import re
import zipfile

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
import matplotlib.pyplot as plt


PROJECT = Path.cwd()
RAW_DIR = PROJECT / "data" / "raw" / "GSE241124"
EXTRACTED = RAW_DIR / "extracted"
META_FILE = RAW_DIR / "GSE241124_spatialseq_metadata_acutewound.txt.gz"
GENE_SET_FILE = PROJECT / "results" / "srrs_framework" / "SRRS_locked_gene_sets.tsv"
OUT_DIR = PROJECT / "results" / "srrs_framework" / "spatial_gse241124"
FIG_DIR = PROJECT / "figures" / "srrs_framework" / "spatial_gse241124"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


SRRS_COMPONENTS = [
    "fibroblast_repair_activation",
    "stromal_ligand_panel",
    "vascular_perivascular_receiver_coupling",
    "d7_acute_wound_alignment",
]
MARKER_SETS = {
    "fibroblast_stromal_marker": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A1"],
    "endothelial_marker": ["PECAM1", "VWF", "KDR", "RAMP2", "CLDN5", "ESAM"],
    "pericyte_smc_marker": ["RGS5", "PDGFRB", "MCAM", "CSPG4", "ACTA2", "TAGLN", "MYH11"],
    "ecm_remodeling_marker": ["TNC", "THBS1", "SERPINE1", "MMP1", "MMP3", "MMP10", "ITGA5", "ITGAV", "PDPN"],
}


def read_gene_sets():
    tbl = pd.read_csv(GENE_SET_FILE, sep="\t")
    return {
        mod: sub["gene"].dropna().astype(str).tolist()
        for mod, sub in tbl.groupby("module", sort=False)
    }


def key_from_file(path, suffix):
    name = path.name
    name = re.sub(r"^GSM\d+_", "", name)
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
        sub = mat[gene_to_indices[gene], :]
        summed = np.asarray(sub.sum(axis=0)).ravel()
        cols.append(summed)
    spot_by_gene = sparse.csr_matrix(np.vstack(cols).T)
    return barcodes, present, spot_by_gene


def read_positions(zip_path):
    with zipfile.ZipFile(zip_path) as zf:
        candidates = [n for n in zf.namelist() if n.endswith("tissue_positions_list.csv") or n.endswith("tissue_positions.csv")]
        if not candidates:
            return pd.DataFrame()
        name = candidates[0]
        with zf.open(name) as fh:
            pos = pd.read_csv(fh, header=None)
    if pos.shape[1] >= 6:
        pos = pos.iloc[:, :6]
        pos.columns = ["raw_barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    return pos


def score_modules(log_expr, genes, module_sets):
    # log_expr: spots x genes dataframe, genes already globally z-scored outside.
    out = {}
    for module, module_genes in module_sets.items():
        present = [g for g in module_genes if g in log_expr.columns]
        if len(present) < 3:
            out[module] = np.full(log_expr.shape[0], np.nan)
        else:
            out[module] = log_expr[present].mean(axis=1).to_numpy()
    return pd.DataFrame(out)


def z(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.full_like(x, np.nan)
    return (x - np.nanmean(x)) / sd


def paired_condition_test(sample_summary, feature, condition):
    wide = sample_summary.pivot_table(index="Donor", columns="Condition", values=feature, aggfunc="mean")
    if "Skin" not in wide.columns or condition not in wide.columns:
        return None
    dat = wide[["Skin", condition]].dropna()
    if len(dat) < 3:
        p = math.nan
    else:
        try:
            p = stats.wilcoxon(dat[condition], dat["Skin"]).pvalue
        except Exception:
            p = math.nan
    return {
        "feature": feature,
        "comparison": f"{condition}_vs_Skin_paired_by_donor",
        "n_donors": len(dat),
        "mean_skin": float(dat["Skin"].mean()) if len(dat) else math.nan,
        f"mean_{condition}": float(dat[condition].mean()) if len(dat) else math.nan,
        "delta_condition_minus_skin": float((dat[condition] - dat["Skin"]).mean()) if len(dat) else math.nan,
        "p_paired_wilcoxon": p,
    }


def main():
    gene_sets = read_gene_sets()
    module_sets = {k: gene_sets[k] for k in SRRS_COMPONENTS}
    module_sets["negative_housekeeping_control"] = gene_sets.get("negative_housekeeping_control", [])
    module_sets.update(MARKER_SETS)
    wanted_genes = sorted(set(g for genes in module_sets.values() for g in genes))

    meta = pd.read_csv(META_FILE, sep="\t")
    meta["raw_barcode"] = meta["barcode"].str.extract(r"([ACGT]+-1)$")

    h5_files = list(EXTRACTED.glob("*_filtered_feature_bc_matrix.h5"))
    zip_files = list(EXTRACTED.glob("*_spatial_images.zip"))
    h5_by_key = {key_from_file(p, "_filtered_feature_bc_matrix.h5"): p for p in h5_files}
    zip_by_key = {key_from_file(p, "_spatial_images.zip"): p for p in zip_files}

    all_blocks = []
    presence_rows = []
    for (orig_ident, condition, sample_name, donor), sub_meta in meta.groupby(
        ["orig.ident", "Condition", "Sample_name", "Donor"], observed=True
    ):
        h5_path = h5_by_key.get(orig_ident) or h5_by_key.get(condition)
        zip_path = zip_by_key.get(orig_ident) or zip_by_key.get(condition)
        if h5_path is None:
            print(f"Skipping {sample_name}: no H5 file for {orig_ident}")
            continue
        barcodes, present, mat = read_10x_h5_gene_subset(h5_path, wanted_genes)
        if mat.shape[1] == 0:
            continue
        expr = pd.DataFrame(mat.toarray(), index=barcodes, columns=present)
        expr = expr.loc[:, ~expr.columns.duplicated()]
        sub = sub_meta.copy().set_index("raw_barcode")
        common = expr.index.intersection(sub.index)
        if len(common) == 0:
            print(f"Skipping {sample_name}: no metadata/H5 barcode overlap")
            continue
        expr = expr.loc[common]
        sub = sub.loc[common].copy()
        sub.index.name = "raw_barcode"
        sub = sub.reset_index()
        lib = expr.sum(axis=1).replace(0, np.nan)
        log_expr = np.log1p(expr.div(lib, axis=0) * 10000.0).fillna(0.0)
        block = sub.reset_index(drop=True)
        gene_block = log_expr.add_prefix("gene__").reset_index(drop=True)
        block = pd.concat([block, gene_block], axis=1)
        if zip_path is not None:
            pos = read_positions(zip_path)
            if not pos.empty:
                block = block.merge(pos, on="raw_barcode", how="left")
        all_blocks.append(block)
        presence_rows.append({
            "sample_name": sample_name,
            "orig_ident": orig_ident,
            "condition": condition,
            "spots_used": len(block),
            "n_wanted_genes": len(wanted_genes),
            "n_present_genes": len(present),
            "present_gene_fraction": len(present) / len(wanted_genes),
        })

    all_spots = pd.concat(all_blocks, ignore_index=True)
    gene_cols = [c for c in all_spots.columns if c.startswith("gene__")]
    gene_names = [c.replace("gene__", "", 1) for c in gene_cols]
    expr_df = all_spots[gene_cols].copy()
    expr_df.columns = gene_names
    expr_z = expr_df.apply(z, axis=0)

    scores = score_modules(expr_z, expr_z.columns, module_sets)
    for col in scores.columns:
        all_spots[col] = scores[col]
    for component in SRRS_COMPONENTS:
        all_spots[f"z_{component}"] = z(all_spots[component])
    all_spots["spatial_SRRS"] = all_spots[[f"z_{c}" for c in SRRS_COMPONENTS]].mean(axis=1, skipna=True)
    all_spots["negative_housekeeping_control_z"] = z(all_spots["negative_housekeeping_control"])

    keep_cols = [
        "barcode", "raw_barcode", "orig.ident", "Condition", "Sample_name", "Donor", "AnnoType",
        "nCount_Spatial", "nFeature_Spatial", "spatial_SRRS", "negative_housekeeping_control_z",
        "pxl_row", "pxl_col", "array_row", "array_col",
    ] + list(module_sets.keys())
    keep_cols = [c for c in keep_cols if c in all_spots.columns]
    all_spots[keep_cols].to_csv(OUT_DIR / "GSE241124_spatial_SRRS_spot_scores.tsv", sep="\t", index=False)
    pd.DataFrame(presence_rows).to_csv(OUT_DIR / "GSE241124_spatial_SRRS_gene_presence_by_sample.tsv", sep="\t", index=False)

    sample_summary = (
        all_spots.groupby(["Donor", "Condition", "Sample_name", "orig.ident"], observed=True)[
            ["spatial_SRRS", "negative_housekeeping_control_z"] + SRRS_COMPONENTS +
            ["fibroblast_stromal_marker", "endothelial_marker", "pericyte_smc_marker", "ecm_remodeling_marker"]
        ]
        .mean()
        .reset_index()
    )
    sample_summary["spots"] = all_spots.groupby(["Donor", "Condition", "Sample_name", "orig.ident"], observed=True).size().to_numpy()
    sample_summary.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_sample_summary.tsv", sep="\t", index=False)

    region_summary = (
        all_spots.groupby(["Condition", "AnnoType"], observed=True)
        .agg(
            spots=("spatial_SRRS", "size"),
            spatial_SRRS_mean=("spatial_SRRS", "mean"),
            fibroblast_marker_mean=("fibroblast_stromal_marker", "mean"),
            endothelial_marker_mean=("endothelial_marker", "mean"),
            pericyte_marker_mean=("pericyte_smc_marker", "mean"),
            ecm_marker_mean=("ecm_remodeling_marker", "mean"),
        )
        .reset_index()
        .sort_values(["Condition", "spatial_SRRS_mean"], ascending=[True, False])
    )
    region_summary.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_region_summary.tsv", sep="\t", index=False)

    test_features = ["spatial_SRRS", "negative_housekeeping_control_z"] + SRRS_COMPONENTS
    test_rows = []
    for feature in test_features:
        for condition in ["Wound1", "Wound7", "Wound30"]:
            res = paired_condition_test(sample_summary, feature, condition)
            if res:
                test_rows.append(res)
    tests = pd.DataFrame(test_rows)
    if not tests.empty:
        tests["padj_bh"] = np.nan
        ok = tests["p_paired_wilcoxon"].notna().to_numpy()
        if ok.sum():
            p = tests.loc[ok, "p_paired_wilcoxon"].to_numpy()
            order = np.argsort(p)
            adj = p[order] * len(p) / (np.arange(len(p)) + 1)
            adj = np.minimum.accumulate(adj[::-1])[::-1]
            tmp = np.empty_like(adj)
            tmp[order] = np.clip(adj, 0, 1)
            tests.loc[ok, "padj_bh"] = tmp
    tests.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_condition_tests.tsv", sep="\t", index=False)

    corr_features = ["fibroblast_stromal_marker", "endothelial_marker", "pericyte_smc_marker", "ecm_remodeling_marker"]
    corr_rows = []
    for condition, sub in all_spots.groupby("Condition", observed=True):
        for feature in corr_features:
            dat = sub[["spatial_SRRS", feature]].dropna()
            if len(dat) >= 20:
                rho, p = stats.spearmanr(dat["spatial_SRRS"], dat[feature])
                corr_rows.append({
                    "condition": condition,
                    "feature": feature,
                    "spearman_rho": rho,
                    "p_spearman": p,
                    "n_spots": len(dat),
                })
    pd.DataFrame(corr_rows).to_csv(OUT_DIR / "GSE241124_spatial_SRRS_marker_correlations.tsv", sep="\t", index=False)

    colors = {"Skin": "#009E73", "Wound1": "#E69F00", "Wound7": "#0072B2", "Wound30": "#CC79A7"}
    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    order = ["Skin", "Wound1", "Wound7", "Wound30"]
    for donor, sub in sample_summary.groupby("Donor", observed=True):
        sub = sub.set_index("Condition").reindex(order)
        ax.plot(range(len(order)), sub["spatial_SRRS"], color="#777777", alpha=0.45, linewidth=0.9)
    vals = [sample_summary.loc[sample_summary["Condition"] == c, "spatial_SRRS"].dropna().to_numpy() for c in order]
    bp = ax.boxplot(vals, positions=range(len(order)), widths=0.45, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], [colors[o] for o in order]):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    for i, cnd in enumerate(order):
        y = vals[i]
        jitter = np.linspace(-0.06, 0.06, len(y)) if len(y) > 1 else np.array([0])
        ax.scatter(np.full(len(y), i) + jitter, y, color=colors[cnd], s=30, alpha=0.95, zorder=3)
    ax.axhline(0, color="#999999", linewidth=0.6)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=20, ha="right")
    ax.set_ylabel("Spatial SRRS")
    fig.savefig(FIG_DIR / "GSE241124_spatial_SRRS_condition_paired_samples.png", dpi=240)
    plt.close(fig)

    top_regions = region_summary[region_summary["spots"] >= 30].copy()
    top_regions = top_regions.sort_values("spatial_SRRS_mean", ascending=False).head(18)
    fig, ax = plt.subplots(figsize=(7.5, 5), constrained_layout=True)
    labels = top_regions["Condition"] + ": " + top_regions["AnnoType"]
    ax.barh(labels[::-1], top_regions["spatial_SRRS_mean"].to_numpy()[::-1], color="#0072B2")
    ax.axvline(0, color="#777777", linewidth=0.7)
    ax.set_xlabel("Mean spatial SRRS")
    ax.set_ylabel("")
    fig.savefig(FIG_DIR / "GSE241124_spatial_SRRS_top_regions.png", dpi=240)
    plt.close(fig)

    mapped = all_spots.dropna(subset=["pxl_row", "pxl_col"]).copy()
    if not mapped.empty:
        reps = []
        for condition in order:
            sub = mapped[mapped["Condition"] == condition]
            if sub.empty:
                continue
            sample = sub.groupby("Sample_name")["spatial_SRRS"].mean().sort_values(ascending=False).index[0]
            reps.append((condition, sample))
        fig, axes = plt.subplots(1, len(reps), figsize=(4.0 * len(reps), 3.8), constrained_layout=True)
        if len(reps) == 1:
            axes = [axes]
        vmin = np.nanquantile(mapped["spatial_SRRS"], 0.02)
        vmax = np.nanquantile(mapped["spatial_SRRS"], 0.98)
        for ax, (condition, sample) in zip(axes, reps):
            sub = mapped[mapped["Sample_name"] == sample]
            sc = ax.scatter(sub["pxl_col"], -sub["pxl_row"], c=sub["spatial_SRRS"], s=7, cmap="viridis", vmin=vmin, vmax=vmax)
            ax.set_title(sample)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
        fig.colorbar(sc, ax=axes, shrink=0.7, label="Spatial SRRS")
        fig.savefig(FIG_DIR / "GSE241124_spatial_SRRS_representative_maps.png", dpi=240)
        plt.close(fig)

    print("Wrote GSE241124 spatial SRRS outputs")


if __name__ == "__main__":
    main()
