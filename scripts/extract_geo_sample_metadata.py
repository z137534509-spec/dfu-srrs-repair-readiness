#!/usr/bin/env python3
"""Extract sample-level metadata from downloaded GEO series matrix files."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = ROOT / "data" / "metadata"
OUT = METADATA_DIR / "dfu_geo_sample_metadata.tsv"


def split_values(line: str) -> list[str]:
    return [part.strip().strip('"') for part in line.rstrip("\n").split("\t")[1:]]


def parse_characteristic(value: str) -> tuple[str, str] | None:
    if ":" not in value:
        return None
    key, val = value.split(":", 1)
    key = key.strip().lower().replace(" ", "_")
    val = val.strip()
    if not key:
        return None
    return key, val


def parse_matrix(path: Path) -> list[dict[str, str]]:
    table: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                key = line.split("\t", 1)[0].replace("!Sample_", "").lower()
                values = split_values(line)
                if key == "characteristics_ch1":
                    characteristics.append(values)
                else:
                    table[key] = values
    count = max((len(values) for values in table.values()), default=0)
    rows: list[dict[str, str]] = []
    for index in range(count):
        row: dict[str, str] = {
            "series": path.name.split("_series_matrix", 1)[0],
            "sample_index": str(index + 1),
        }
        for key, values in table.items():
            if index < len(values):
                row[key] = values[index]
        parsed_chars: dict[str, list[str]] = defaultdict(list)
        for char_values in characteristics:
            if index >= len(char_values):
                continue
            parsed = parse_characteristic(char_values[index])
            if parsed:
                key, val = parsed
                parsed_chars[key].append(val)
        for key, values in parsed_chars.items():
            row[f"char_{key}"] = " | ".join(dict.fromkeys(values))
        rows.append(row)
    return rows


def main() -> None:
    matrix_files = sorted(METADATA_DIR.glob("GSE*_series_matrix.txt.gz"))
    rows: list[dict[str, str]] = []
    for matrix in matrix_files:
        rows.extend(parse_matrix(matrix))
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT}")
    print(f"Samples: {len(rows)}")


if __name__ == "__main__":
    main()
