#!/usr/bin/env python3
"""Audit GEO accessions for the DFU virtual perturbation project.

The script checks NCBI GEO FTP folders, records available series matrix and
supplementary files, and downloads the small series-matrix metadata files for
sample-level inspection.
"""

from __future__ import annotations

import csv
import gzip
import html.parser
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import requests


ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = ROOT / "data" / "metadata"
METADATA_DIR.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class Dataset:
    accession: str
    planned_role: str
    analysis_type: str
    priority: str
    notes: str


DATASETS = [
    Dataset(
        "GSE165816",
        "Primary DFU scRNA-seq dataset; cell states, repair scores, GRN, composition",
        "scRNA-seq",
        "primary",
        "Main DFU single-cell landscape; avoid MT1X-style composition-only story.",
    ),
    Dataset(
        "GSE223964",
        "Independent DFU/NDFU scRNA-seq validation; candidate localization",
        "scRNA-seq",
        "validation",
        "Processed h5ad/h5seurat files available; group labels must be rechecked.",
    ),
    Dataset(
        "GSE248247",
        "Small DFU vs non-diabetic foot ulcer scRNA-seq descriptive validation",
        "scRNA-seq",
        "supportive",
        "Very small sample/cell number; do not use as main statistical validation.",
    ),
    Dataset(
        "GSE80178",
        "Bulk array external validation; repair/inflammation/DNA repair signatures",
        "bulk array",
        "validation",
        "Previously used in MT1X work; avoid duplicate figures and conclusions.",
    ),
    Dataset(
        "GSE134431",
        "Bulk RNA-seq validation; healing vs non-healing immune recruitment logic",
        "bulk RNA-seq",
        "validation",
        "Useful for FOXM1/immune-stalling background and healer/non-healer contrasts.",
    ),
    Dataset(
        "GSE199939",
        "Second bulk RNA-seq validation; DFU vs non-DFU skin",
        "bulk RNA-seq",
        "validation",
        "Has gene TPM matrix; confirm sample labels before differential testing.",
    ),
    Dataset(
        "GSE166120",
        "GeoMx spatial-context validation of modules/ROI background",
        "GeoMx spatial transcriptomics",
        "spatial_context",
        "Panel-limited; validate modules/ROI context, not necessarily individual candidates.",
    ),
]


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.links.append(value)


def geo_bucket(accession: str) -> str:
    match = re.fullmatch(r"GSE(\d+)", accession)
    if not match:
        raise ValueError(f"Unsupported GEO accession: {accession}")
    number = match.group(1)
    return f"GSE{number[:-3]}nnn"


def base_url(accession: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_bucket(accession)}/{accession}"


def get_links(url: str) -> list[str]:
    try:
        response = requests.get(url, timeout=30)
    except requests.RequestException:
        return []
    if response.status_code != 200 or "Object not found" in response.text:
        return []
    parser = LinkParser()
    parser.feed(response.text)
    return [link for link in parser.links if link and not link.startswith("/geo/series")]


def download_small_file(url: str, output: Path, max_bytes: int = 10_000_000) -> bool:
    if output.exists() and output.stat().st_size > 0:
        return True
    try:
        with requests.get(url, stream=True, timeout=60) as response:
            if response.status_code != 200:
                return False
            total = int(response.headers.get("Content-Length", "0") or "0")
            if total and total > max_bytes:
                return False
            output.parent.mkdir(parents=True, exist_ok=True)
            written = 0
            with output.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 256):
                    if not chunk:
                        continue
                    written += len(chunk)
                    if written > max_bytes:
                        output.unlink(missing_ok=True)
                        return False
                    handle.write(chunk)
    except requests.RequestException:
        return False
    return output.exists() and output.stat().st_size > 0


def parse_series_matrix(path: Path) -> dict[str, str]:
    result: dict[str, str] = {
        "series_title": "",
        "series_status": "",
        "platforms": "",
        "sample_count": "0",
        "sample_titles": "",
        "sample_geo_accessions": "",
        "sample_sources": "",
        "sample_characteristics": "",
    }
    if not path.exists():
        return result
    opener = gzip.open if path.suffix == ".gz" else open
    sample_titles: list[str] = []
    sample_accessions: list[str] = []
    sample_sources: list[str] = []
    sample_characteristics: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!Series_title"):
                result["series_title"] = first_value(line)
            elif line.startswith("!Series_status"):
                result["series_status"] = first_value(line)
            elif line.startswith("!Series_platform_id"):
                result["platforms"] = values_joined(line)
            elif line.startswith("!Sample_title"):
                sample_titles = split_values(line)
            elif line.startswith("!Sample_geo_accession"):
                sample_accessions = split_values(line)
            elif line.startswith("!Sample_source_name_ch1"):
                sample_sources = split_values(line)
            elif line.startswith("!Sample_characteristics_ch1"):
                values = split_values(line)
                sample_characteristics.extend(values)
            elif line.startswith("!series_matrix_table_begin"):
                break
    result["sample_count"] = str(len(sample_accessions) or len(sample_titles))
    result["sample_titles"] = " || ".join(sample_titles[:20])
    result["sample_geo_accessions"] = " || ".join(sample_accessions[:20])
    result["sample_sources"] = " || ".join(sample_sources[:20])
    result["sample_characteristics"] = " || ".join(sorted(set(sample_characteristics))[:30])
    return result


def split_values(line: str) -> list[str]:
    parts = line.split("\t")[1:]
    return [part.strip().strip('"') for part in parts if part.strip()]


def values_joined(line: str) -> str:
    return " || ".join(split_values(line))


def first_value(line: str) -> str:
    values = split_values(line)
    return values[0] if values else ""


def classify_download_status(supplementary_files: Iterable[str]) -> str:
    files = list(supplementary_files)
    if any(name.endswith((".h5ad.gz", ".h5seurat")) for name in files):
        return "processed_single_cell_object_available"
    if any("matrix" in name.lower() or name.endswith((".txt.gz", ".xlsx")) for name in files):
        return "processed_matrix_available"
    if any(name.endswith("_RAW.tar") for name in files):
        return "raw_or_custom_matrix_available"
    return "metadata_only_or_unconfirmed"


def main() -> None:
    rows: list[dict[str, str]] = []
    for dataset in DATASETS:
        url = base_url(dataset.accession)
        supplementary = sorted(
            link for link in get_links(f"{url}/suppl/") if link != "../"
        )
        matrix_files = sorted(link for link in get_links(f"{url}/matrix/") if link.endswith(".gz"))
        soft_files = sorted(link for link in get_links(f"{url}/soft/") if link.endswith(".gz"))
        matrix_path = METADATA_DIR / f"{dataset.accession}_series_matrix.txt.gz"
        matrix_downloaded = False
        if matrix_files:
            matrix_downloaded = download_small_file(f"{url}/matrix/{matrix_files[0]}", matrix_path)
        matrix_meta = parse_series_matrix(matrix_path if matrix_downloaded else Path())
        rows.append(
            {
                "accession": dataset.accession,
                "analysis_type": dataset.analysis_type,
                "priority": dataset.priority,
                "planned_role": dataset.planned_role,
                "ftp_base": url,
                "series_title": matrix_meta["series_title"],
                "series_status": matrix_meta["series_status"],
                "platforms": matrix_meta["platforms"],
                "sample_count_from_series_matrix": matrix_meta["sample_count"],
                "supplementary_files": " || ".join(supplementary),
                "matrix_files": " || ".join(matrix_files),
                "soft_files": " || ".join(soft_files),
                "download_status": classify_download_status(supplementary),
                "sample_titles_preview": matrix_meta["sample_titles"],
                "sample_sources_preview": matrix_meta["sample_sources"],
                "sample_characteristics_preview": matrix_meta["sample_characteristics"],
                "notes": dataset.notes,
            }
        )

    out_csv = METADATA_DIR / "dfu_geo_accession_audit.csv"
    out_tsv = METADATA_DIR / "dfu_geo_accession_audit.tsv"
    fieldnames = list(rows[0].keys())
    with out_csv.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with out_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_tsv}")


if __name__ == "__main__":
    main()
