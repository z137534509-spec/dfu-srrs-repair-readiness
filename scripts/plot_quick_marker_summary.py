#!/usr/bin/env python3
"""Plot quick marker-score summaries for GSE165816."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "quick_marker_score" / "GSE165816_foot_skin_quick_summary_by_disease_celltype.tsv"
FIG_DIR = ROOT / "figures" / "quick_marker_score"
FIG_DIR.mkdir(parents=True, exist_ok=True)


DISEASE_ORDER = ["Non-diabetic", "Non-DFU Diabetic", "DFU-healer", "DFU-nonhealer"]


def main() -> None:
    df = pd.read_csv(IN, sep="\t")
    df["disease"] = pd.Categorical(df["disease"], categories=DISEASE_ORDER, ordered=True)
    frac = df.pivot_table(
        index="disease",
        columns="quick_cell_type",
        values="fraction_within_disease",
        aggfunc="sum",
        observed=False,
    ).fillna(0)
    preferred_cols = [
        "keratinocyte",
        "fibroblast_stromal",
        "endothelial",
        "myeloid",
        "t_nk",
        "b_plasma",
        "pericyte_smc",
        "melanocyte",
        "schwann",
        "unknown",
    ]
    frac = frac[[col for col in preferred_cols if col in frac.columns]]
    ax = frac.plot(kind="bar", stacked=True, figsize=(10, 5), width=0.8)
    ax.set_ylabel("Fraction of foot-skin cells")
    ax.set_xlabel("")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "GSE165816_quick_celltype_fraction_by_disease.png", dpi=200)
    plt.close()

    modules = [
        "inflammatory_arrest_mean",
        "angiogenesis_mean",
        "ecm_remodeling_mean",
        "epithelial_migration_mean",
        "hypoxia_oxidative_stress_mean",
        "senescence_sasp_mean",
    ]
    focus = df[df["quick_cell_type"].isin(["keratinocyte", "fibroblast_stromal", "endothelial", "myeloid"])]
    for module in modules:
        mat = focus.pivot_table(
            index="quick_cell_type",
            columns="disease",
            values=module,
            aggfunc="mean",
            observed=False,
        ).reindex(columns=DISEASE_ORDER)
        plt.figure(figsize=(7, 3.5))
        sns.heatmap(mat, annot=True, fmt=".2f", cmap="viridis")
        plt.title(module.replace("_mean", ""))
        plt.xlabel("")
        plt.ylabel("")
        plt.tight_layout()
        plt.savefig(FIG_DIR / f"GSE165816_{module}_heatmap.png", dpi=200)
        plt.close()
    print(f"Wrote figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
