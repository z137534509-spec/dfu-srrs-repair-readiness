from pathlib import Path
import re
import time

import pandas as pd
import requests


PROJECT = Path.cwd()
GENE_SET_FILE = PROJECT / "results" / "srrs_framework" / "SRRS_locked_gene_sets.tsv"
OUT_DIR = PROJECT / "results" / "srrs_framework" / "perturbation_enrichr"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ENRICHR = "https://maayanlab.cloud/Enrichr"
LIBRARIES = [
    "Ligand_Perturbations_from_GEO_up",
    "Ligand_Perturbations_from_GEO_down",
    "LINCS_L1000_Ligand_Perturbations_up",
    "LINCS_L1000_Ligand_Perturbations_down",
    "Gene_Perturbations_from_GEO_up",
    "Gene_Perturbations_from_GEO_down",
    "LINCS_L1000_Chem_Pert_up",
    "LINCS_L1000_Chem_Pert_down",
    "Diabetes_Perturbations_GEO_2022",
]
FOCUS_RE = re.compile(
    r"IL[-_ ]?6|IL[-_ ]?11|TGFB|TGF[-_ ]?B|ACTIVIN|INHBA|TNF|STAT3|JAK|"
    r"HIF|HYPOX|TNC|TENASCIN|THBS1|SERPINE1|PAI[-_ ]?1|PTGS2|PROSTAGLANDIN|"
    r"FIBROBLAST|WOUND|MATRIX|INTEGRIN",
    re.IGNORECASE,
)


def add_list(name, genes):
    payload = {
        "list": "\n".join(genes),
        "description": name,
    }
    r = requests.post(f"{ENRICHR}/addList", files=payload, timeout=30)
    r.raise_for_status()
    return r.json()["userListId"]


def enrich(user_list_id, library):
    r = requests.get(
        f"{ENRICHR}/enrich",
        params={"userListId": user_list_id, "backgroundType": library},
        timeout=60,
    )
    r.raise_for_status()
    rows = r.json().get(library, [])
    out = []
    for row in rows:
        out.append({
            "rank": row[0],
            "term": row[1],
            "p_value": row[2],
            "z_score": row[3],
            "combined_score": row[4],
            "overlap_genes": ";".join(row[5]) if isinstance(row[5], list) else str(row[5]),
            "adjusted_p_value": row[6],
        })
    return pd.DataFrame(out)


def unique_keep_order(items):
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


def main():
    gene_sets = pd.read_csv(GENE_SET_FILE, sep="\t")
    sets = {
        module: sub["gene"].dropna().astype(str).tolist()
        for module, sub in gene_sets.groupby("module", sort=False)
    }
    query_sets = {
        "SRRS_repair_ligand_receiver_core": unique_keep_order(
            sets["fibroblast_repair_activation"]
            + sets["stromal_ligand_panel"]
            + sets["vascular_perivascular_receiver_coupling"]
        ),
        "SRRS_d7_wound_alignment": unique_keep_order(sets["d7_acute_wound_alignment"][:100]),
        "SRRS_all_modules": unique_keep_order(
            sets["fibroblast_repair_activation"]
            + sets["stromal_ligand_panel"]
            + sets["vascular_perivascular_receiver_coupling"]
            + sets["d7_acute_wound_alignment"][:100]
        ),
    }

    all_rows = []
    for query_name, genes in query_sets.items():
        genes = [g for g in genes if isinstance(g, str) and g]
        user_list_id = add_list(query_name, genes)
        for library in LIBRARIES:
            try:
                res = enrich(user_list_id, library)
                if res.empty:
                    continue
                res.insert(0, "library", library)
                res.insert(0, "query_gene_n", len(genes))
                res.insert(0, "query", query_name)
                all_rows.append(res.head(100))
                time.sleep(0.2)
            except Exception as exc:
                all_rows.append(pd.DataFrame([{
                    "query": query_name,
                    "query_gene_n": len(genes),
                    "library": library,
                    "rank": None,
                    "term": f"ERROR: {exc}",
                    "p_value": None,
                    "z_score": None,
                    "combined_score": None,
                    "overlap_genes": "",
                    "adjusted_p_value": None,
                }]))
    full = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    full.to_csv(OUT_DIR / "SRRS_Enrichr_perturbation_top100_by_library.tsv", sep="\t", index=False)

    if not full.empty:
        focus = full[full["term"].fillna("").str.contains(FOCUS_RE)].copy()
        focus = focus.sort_values(["query", "adjusted_p_value", "p_value", "rank"], na_position="last")
        focus.to_csv(OUT_DIR / "SRRS_Enrichr_focused_perturbation_hits.tsv", sep="\t", index=False)
        top = (
            full.sort_values(["query", "adjusted_p_value", "p_value", "rank"], na_position="last")
            .groupby(["query", "library"], as_index=False)
            .head(10)
        )
        top.to_csv(OUT_DIR / "SRRS_Enrichr_top10_hits.tsv", sep="\t", index=False)
    print("Wrote SRRS Enrichr perturbation outputs")


if __name__ == "__main__":
    main()
