from pathlib import Path
import re
import zipfile

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt


PROJECT = Path.cwd()
RAW_DIR = PROJECT / "data" / "raw" / "GSE241124"
EXTRACTED = RAW_DIR / "extracted"
META_FILE = RAW_DIR / "GSE241124_spatialseq_metadata_acutewound.txt.gz"
OUT_DIR = PROJECT / "results" / "srrs_framework" / "spatial_receptor_expression_gse241124"
FIG_DIR = PROJECT / "figures" / "srrs_framework" / "spatial_receptor_expression_gse241124"
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


def main():
    wanted_genes = sorted(set(g for spec in AXES.values() for side in ["ligands", "receptors"] for g in spec[side]))
    receptor_genes = sorted(set(g for spec in AXES.values() for g in spec["receptors"]))
    ligand_genes = sorted(set(g for spec in AXES.values() for g in spec["ligands"]))

    meta = pd.read_csv(META_FILE, sep="\t")
    meta["raw_barcode"] = meta["barcode"].str.extract(r"([ACGT]+-1)$")

    h5_by_key = {key_from_file(p, "_filtered_feature_bc_matrix.h5"): p for p in EXTRACTED.glob("*_filtered_feature_bc_matrix.h5")}
    zip_by_key = {key_from_file(p, "_spatial_images.zip"): p for p in EXTRACTED.glob("*_spatial_images.zip")}

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
            "missing_genes": ";".join([g for g in wanted_genes if g not in present]),
        })

    spots = pd.concat(blocks, ignore_index=True)
    gene_cols = [c for c in spots.columns if c.startswith("gene__")]
    gene_names = [c.replace("gene__", "", 1) for c in gene_cols]
    expr = spots[gene_cols].copy()
    expr.columns = gene_names
    for gene in wanted_genes:
        if gene not in expr.columns:
            expr[gene] = np.nan
    for gene in wanted_genes:
        spots[f"logexpr__{gene}"] = expr[gene]

    keep_cols = [
        "barcode", "raw_barcode", "orig.ident", "Condition", "Sample_name", "Donor", "AnnoType",
        "nCount_Spatial", "nFeature_Spatial", "pxl_row", "pxl_col",
    ]
    expr_cols = [f"logexpr__{g}" for g in wanted_genes]
    spots[keep_cols + expr_cols].to_csv(OUT_DIR / "GSE241124_spatial_LR_gene_logexpr_spots.tsv", sep="\t", index=False)
    pd.DataFrame(presence_rows).to_csv(OUT_DIR / "GSE241124_spatial_LR_gene_presence.tsv", sep="\t", index=False)

    sample_rows = []
    for gene in wanted_genes:
        col = f"logexpr__{gene}"
        sample = (
            spots.groupby(["Donor", "Condition", "Sample_name"], observed=True)[col]
            .mean()
            .reset_index()
            .rename(columns={col: "mean_logexpr"})
        )
        sample["gene"] = gene
        sample["role"] = "ligand" if gene in ligand_genes else "receptor"
        axes = [axis for axis, spec in AXES.items() if gene in spec["ligands"] or gene in spec["receptors"]]
        sample["axis"] = ";".join(axes)
        sample_rows.append(sample)
    sample_summary = pd.concat(sample_rows, ignore_index=True)
    sample_summary.to_csv(OUT_DIR / "GSE241124_spatial_LR_gene_sample_summary.tsv", sep="\t", index=False)

    spots["log_nCount"] = np.log1p(spots["nCount_Spatial"])
    spots["log_nFeature"] = np.log1p(spots["nFeature_Spatial"])
    model_rows = []
    for gene in wanted_genes:
        col = f"logexpr__{gene}"
        dat = spots.dropna(subset=[col, "Condition", "Sample_name", "log_nCount", "log_nFeature"]).copy()
        fit = smf.ols(f"{col} ~ C(Condition, Treatment(reference='Skin')) + log_nCount + log_nFeature", data=dat).fit(
            cov_type="cluster",
            cov_kwds={"groups": dat["Sample_name"]},
        )
        for term in fit.params.index:
            if "Condition" not in term:
                continue
            axes = [axis for axis, spec in AXES.items() if gene in spec["ligands"] or gene in spec["receptors"]]
            model_rows.append({
                "gene": gene,
                "role": "ligand" if gene in ligand_genes else "receptor",
                "axis": ";".join(axes),
                "term": term,
                "estimate": fit.params[term],
                "std_error_cluster_sample": fit.bse[term],
                "p_cluster_sample": fit.pvalues[term],
                "n_spots": int(fit.nobs),
                "n_clusters": dat["Sample_name"].nunique(),
            })
    models = pd.DataFrame(model_rows)
    models["padj_bh"] = bh_adjust(models["p_cluster_sample"])
    models.to_csv(OUT_DIR / "GSE241124_spatial_LR_gene_condition_clustered_models.tsv", sep="\t", index=False)

    # Compact Wound7 coefficient figure for receptor-side supplement.
    w7 = models[models["term"].str.contains("Wound7")].copy()
    w7 = w7.sort_values(["role", "estimate"])
    fig, ax = plt.subplots(figsize=(8.2, 6.2), constrained_layout=True)
    colors = np.where(w7["role"].eq("ligand"), "#D55E00", "#0072B2")
    ax.barh(w7["gene"] + " (" + w7["role"] + ")", w7["estimate"], color=colors)
    ax.axvline(0, color="#777777", linewidth=0.8)
    ax.set_xlabel("Adjusted Wound7 vs Skin log-expression coefficient")
    fig.savefig(FIG_DIR / "GSE241124_spatial_LR_gene_Wound7_coefficients.png", dpi=240)
    plt.close(fig)

    # Spatial side-by-side maps for one representative Wound7 sample with coordinates.
    selected = ["TNC", "ITGB1", "SDC4", "INHBA", "ACVR1", "ACVR2A", "IL6ST", "IL11RA", "CCR6"]
    w7_spots = spots[(spots["Condition"] == "Wound7") & spots["pxl_row"].notna() & spots["pxl_col"].notna()].copy()
    if not w7_spots.empty:
        chosen_sample = w7_spots.groupby("Sample_name", observed=True).size().sort_values(ascending=False).index[0]
        sub = w7_spots[w7_spots["Sample_name"] == chosen_sample].copy()
        ncols = 3
        nrows = int(np.ceil(len(selected) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(9.0, 8.5), constrained_layout=True)
        axes = np.ravel(axes)
        for ax, gene in zip(axes, selected):
            col = f"logexpr__{gene}"
            vals = sub[col] if col in sub else np.nan
            sc = ax.scatter(sub["pxl_col"], -sub["pxl_row"], c=vals, s=8, cmap="viridis", linewidths=0)
            ax.set_title(gene, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
        for ax in axes[len(selected):]:
            ax.axis("off")
        fig.suptitle(f"Wound7 representative spatial expression: {chosen_sample}", fontsize=10)
        fig.savefig(FIG_DIR / "GSE241124_spatial_LR_representative_Wound7_gene_maps.png", dpi=240)
        plt.close(fig)

    print("Wrote GSE241124 spatial receptor expression outputs")


if __name__ == "__main__":
    main()
