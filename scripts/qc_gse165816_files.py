#!/usr/bin/env python3
"""Quick file-level QC for GSE165816 per-sample count CSVs."""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
COUNT_DIR = ROOT / "data" / "raw" / "GSE165816" / "counts_csv"
SAMPLE_META = ROOT / "data" / "metadata" / "dfu_geo_sample_metadata.tsv"
OUT = ROOT / "data" / "metadata" / "GSE165816_file_qc.tsv"


def inspect_count_file(path: Path) -> dict[str, str | int]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n")
        cell_count = len(header.split(","))
        gene_count = sum(1 for _ in handle)
    match = re.match(r"(GSM\d+)_(.+?)counts\.csv\.gz$", path.name)
    return {
        "file": path.name,
        "geo_accession": match.group(1) if match else "",
        "sample_code": match.group(2) if match else "",
        "cells_from_header": cell_count,
        "genes_from_rows": gene_count,
        "file_size_bytes": path.stat().st_size,
    }


def main() -> None:
    rows = [inspect_count_file(path) for path in sorted(COUNT_DIR.glob("*.csv.gz"))]
    qc = pd.DataFrame(rows)
    sample_meta = pd.read_csv(SAMPLE_META, sep="\t", dtype=str).query("series == 'GSE165816'")
    keep_cols = [
        "geo_accession",
        "title",
        "char_disease",
        "char_tissue",
        "char_age",
        "char_sex",
    ]
    keep_cols = [col for col in keep_cols if col in sample_meta.columns]
    merged = qc.merge(sample_meta[keep_cols], on="geo_accession", how="left")
    merged.to_csv(OUT, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"Wrote {OUT}")
    print(f"Files: {len(merged)}")
    print(f"Cells total: {int(merged['cells_from_header'].sum())}")
    print(merged.groupby(["char_tissue", "char_disease"], dropna=False)["cells_from_header"].sum())


if __name__ == "__main__":
    main()
