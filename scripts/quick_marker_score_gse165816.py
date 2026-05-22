#!/usr/bin/env python3
"""Quick marker-score audit for GSE165816 foot-skin cells.

This is a lightweight pre-analysis for Go/No-Go only. It reads marker genes from
the per-sample dense CSV count matrices without building a full Seurat/AnnData
object. Formal clustering, normalization, integration, and cell annotation must
be redone in the production pipeline.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
COUNT_DIR = ROOT / "data" / "raw" / "GSE165816" / "counts_csv"
QC = ROOT / "data" / "metadata" / "GSE165816_file_qc.tsv"
OUT_DIR = ROOT / "results" / "quick_marker_score"
OUT_DIR.mkdir(parents=True, exist_ok=True)


CELL_TYPE_MARKERS = {
    "keratinocyte": ["KRT14", "KRT5", "KRT1", "KRT10", "KRT6A", "KRT16", "KRT17"],
    "fibroblast_stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"],
    "endothelial": ["PECAM1", "VWF", "KDR", "FLT1", "CLDN5", "RAMP2", "ESAM"],
    "myeloid": ["LYZ", "LST1", "S100A8", "S100A9", "FCGR3A", "TYROBP", "CTSS"],
    "t_nk": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"],
    "b_plasma": ["MS4A1", "CD79A", "CD74", "MZB1", "JCHAIN"],
    "pericyte_smc": ["RGS5", "PDGFRB", "MCAM", "ACTA2", "TAGLN", "MYH11"],
    "melanocyte": ["MLANA", "PMEL", "TYR", "DCT"],
    "schwann": ["SOX10", "MPZ", "PLP1", "S100B"],
}


PROGRAM_MARKERS = {
    "inflammatory_arrest": [
        "IL1B",
        "TNF",
        "CXCL8",
        "CXCL2",
        "CCL2",
        "CCL3",
        "CCL4",
        "S100A8",
        "S100A9",
        "NFKBIA",
        "PTGS2",
    ],
    "angiogenesis": ["PECAM1", "VWF", "KDR", "FLT1", "ESAM", "EMCN", "ANGPT2", "PLVAP"],
    "ecm_remodeling": [
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "FN1",
        "POSTN",
        "MMP1",
        "MMP2",
        "MMP3",
        "MMP9",
        "MMP11",
        "TIMP1",
    ],
    "epithelial_migration": ["KRT6A", "KRT6B", "KRT16", "KRT17", "ITGA3", "ITGB1", "LAMC2", "MMP9", "AREG", "HBEGF"],
    "hypoxia_oxidative_stress": ["HIF1A", "VEGFA", "SOD2", "HMOX1", "NQO1", "TXN", "JUN", "FOS"],
    "senescence_sasp": ["CDKN1A", "CDKN2A", "SERPINE1", "IGFBP7", "MMP3", "CXCL8", "CCL2"],
}


MARKER_GENES = sorted(
    {gene for genes in CELL_TYPE_MARKERS.values() for gene in genes}
    | {gene for genes in PROGRAM_MARKERS.values() for gene in genes}
)


def read_marker_counts(path: Path) -> tuple[list[str], dict[str, np.ndarray]]:
    marker_set = set(MARKER_GENES)
    counts: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        barcodes = next(reader)
        for row in reader:
            if not row:
                continue
            gene = row[0]
            if gene not in marker_set:
                continue
            values = np.fromiter((int(x) if x else 0 for x in row[1:]), dtype=np.float32)
            counts[gene] = values
    return barcodes, counts


def score_gene_set(counts: dict[str, np.ndarray], genes: list[str], n_cells: int) -> np.ndarray:
    present = [counts[gene] for gene in genes if gene in counts]
    if not present:
        return np.zeros(n_cells, dtype=np.float32)
    return np.log1p(np.vstack(present).mean(axis=0))


def process_sample(row: pd.Series) -> pd.DataFrame:
    path = COUNT_DIR / row["file"]
    barcodes, counts = read_marker_counts(path)
    n_cells = len(barcodes)
    cell_type_scores = {
        name: score_gene_set(counts, genes, n_cells) for name, genes in CELL_TYPE_MARKERS.items()
    }
    program_scores = {
        name: score_gene_set(counts, genes, n_cells) for name, genes in PROGRAM_MARKERS.items()
    }
    score_matrix = np.vstack([cell_type_scores[name] for name in CELL_TYPE_MARKERS])
    best_index = score_matrix.argmax(axis=0)
    best_score = score_matrix.max(axis=0)
    cell_types = np.array(list(CELL_TYPE_MARKERS.keys()), dtype=object)[best_index]
    cell_types[best_score <= 0] = "unknown"

    out = pd.DataFrame(
        {
            "sample": row["sample_code"],
            "geo_accession": row["geo_accession"],
            "barcode": barcodes,
            "disease": row["char_disease"],
            "tissue": row["char_tissue"],
            "quick_cell_type": cell_types,
            "quick_cell_type_score": best_score,
        }
    )
    for name, values in program_scores.items():
        out[name] = values
    return out


def main() -> None:
    qc = pd.read_csv(QC, sep="\t", dtype=str)
    qc["cells_from_header"] = qc["cells_from_header"].astype(int)
    foot = qc[qc["char_tissue"] == "Foot skin"].copy()
    cell_tables = []
    for _, row in foot.iterrows():
        print(f"processing {row['geo_accession']} {row['sample_code']} {row['char_disease']}")
        cell_tables.append(process_sample(row))
    cells = pd.concat(cell_tables, ignore_index=True)
    cells.to_csv(OUT_DIR / "GSE165816_foot_skin_quick_cell_scores.tsv.gz", sep="\t", index=False)

    summary = (
        cells.groupby(["disease", "quick_cell_type"], dropna=False)
        .agg(
            cells=("barcode", "size"),
            inflammatory_arrest_mean=("inflammatory_arrest", "mean"),
            angiogenesis_mean=("angiogenesis", "mean"),
            ecm_remodeling_mean=("ecm_remodeling", "mean"),
            epithelial_migration_mean=("epithelial_migration", "mean"),
            hypoxia_oxidative_stress_mean=("hypoxia_oxidative_stress", "mean"),
            senescence_sasp_mean=("senescence_sasp", "mean"),
        )
        .reset_index()
    )
    totals = cells.groupby("disease")["barcode"].size().rename("disease_total_cells")
    summary = summary.merge(totals, on="disease")
    summary["fraction_within_disease"] = summary["cells"] / summary["disease_total_cells"]
    summary.to_csv(OUT_DIR / "GSE165816_foot_skin_quick_summary_by_disease_celltype.tsv", sep="\t", index=False)

    sample_summary = (
        cells.groupby(["geo_accession", "sample", "disease", "quick_cell_type"], dropna=False)
        .size()
        .rename("cells")
        .reset_index()
    )
    sample_summary.to_csv(OUT_DIR / "GSE165816_foot_skin_quick_sample_celltype_counts.tsv", sep="\t", index=False)

    print(f"cells: {len(cells)}")
    print(f"wrote {OUT_DIR}")
    print(summary.sort_values(["disease", "cells"], ascending=[True, False]).head(30).to_string(index=False))


if __name__ == "__main__":
    main()
