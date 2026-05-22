from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from openpyxl import load_workbook


PROJECT = Path.cwd()
RESULTS = PROJECT / "results" / "srrs_framework"
FIG_DIR = PROJECT / "figures" / "ccs_submission"
MAN_DIR = PROJECT / "manuscript" / "Cell_Communication_Signaling"
SUPP_DIR = MAN_DIR / "supplementary_tables"
OUT_DIR = RESULTS / "evidence_hierarchy"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7E57C2"
GREY = "#6E6E6E"


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def savefig(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def letter(ax, label):
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")


def donor_first_figure5():
    paired = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_paired_donor_bootstrap_tests.tsv")
    corr = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_marker_correlation_donor_bootstrap_summary.tsv")
    sample = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_strengthened_sample_summary.tsv")
    models = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_spot_level_clustered_models.tsv")

    order = ["Skin", "Wound1", "Wound7", "Wound30"]
    colors = {"Skin": "#8A8A8A", "Wound1": GREEN, "Wound7": ORANGE, "Wound30": BLUE}

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.3))

    ax = axes[0, 0]
    for donor, sub in sample.groupby("Donor", observed=True):
        sub = sub.set_index("Condition").reindex(order)
        x = np.arange(len(order))
        ax.plot(x, sub["spatial_SRRS"], color="#B6B6B6", linewidth=1.2, zorder=1)
        ax.scatter(x, sub["spatial_SRRS"], s=42, color=[colors[c] for c in order], edgecolor="white", linewidth=0.6, zorder=2)
        if np.all(np.isfinite(sub["spatial_SRRS"].to_numpy())):
            ax.text(x[-1] + 0.04, sub["spatial_SRRS"].iloc[-1], str(donor), fontsize=7, va="center", color="#555555")
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order)
    ax.set_ylabel("Sample-level spatial SRRS")
    ax.set_title("Donor-first spatial SRRS trajectories")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "A")

    ax = axes[0, 1]
    sub = paired[(paired["feature"].isin(["spatial_SRRS", "spatial_SRRS_resid"])) & (paired["comparison"].isin(["Wound1_vs_Skin", "Wound7_vs_Skin", "Wound30_vs_Skin"]))].copy()
    sub["label"] = sub["feature"].str.replace("spatial_SRRS", "SRRS", regex=False) + " " + sub["comparison"].str.replace("_vs_Skin", "", regex=False)
    sub = sub.sort_values("mean_delta")
    y = np.arange(len(sub))
    ax.errorbar(
        sub["mean_delta"],
        y,
        xerr=[sub["mean_delta"] - sub["bootstrap_ci_low"], sub["bootstrap_ci_high"] - sub["mean_delta"]],
        fmt="o",
        color=GREEN,
        ecolor="#777777",
        capsize=2,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["label"], fontsize=8)
    ax.set_xlabel("Paired donor delta vs Skin")
    ax.set_title("Bootstrap donor-level effects")
    ax.text(0.02, 0.05, "Wound7 raw SRRS: 4/4 donors positive\nPaired Wilcoxon P = 0.125", transform=ax.transAxes, fontsize=8.5, va="bottom")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "B")

    ax = axes[1, 0]
    sub = corr[corr["Condition"] == "Wound7"].copy()
    sub["label"] = sub["feature"].str.replace("_marker", "", regex=False).str.replace("_", " ", regex=False)
    sub = sub.sort_values("mean_donor_rho")
    y = np.arange(len(sub))
    ax.errorbar(
        sub["mean_donor_rho"],
        y,
        xerr=[sub["mean_donor_rho"] - sub["bootstrap_ci_low"], sub["bootstrap_ci_high"] - sub["mean_donor_rho"]],
        fmt="o",
        color=PURPLE,
        ecolor="#777777",
        capsize=2,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["label"])
    ax.set_xlabel("Mean donor Spearman rho")
    ax.set_title("Wound7 donor-level marker correlations")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "C")

    ax = axes[1, 1]
    sub = models[models["term"].str.contains("Wound7")].copy()
    sub["model"] = pd.Categorical(sub["model"], categories=["unadjusted", "adjusted", "residual"], ordered=True)
    sub = sub.sort_values("model")
    ax.bar(sub["model"].astype(str), sub["estimate"], color=[BLUE, ORANGE, GREEN])
    ax.set_ylabel("Wound7 vs Skin coefficient")
    ax.set_title("Spot-level sensitivity model")
    ax.text(0.02, 0.95, "Clustered by spatial sample;\nshown as sensitivity support", transform=ax.transAxes, va="top", fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "D")

    fig.tight_layout()
    savefig(fig, "Figure5_Spatial_SRRS_localisation")


def evidence_hierarchy():
    rows = [
        {
            "Evidence layer": "DFU outcome discovery",
            "Dataset/resource": "GSE165816",
            "What it supports": "SRRS is associated with 12-week DFU healing status in the discovery cohort.",
            "What it cannot prove": "Independent generalisability, causality or clinical readiness.",
            "How handled in manuscript": "Primary discovery result; explicitly described as one small outcome cohort.",
        },
        {
            "Evidence layer": "Permutation, LODO and bootstrap robustness",
            "Dataset/resource": "GSE165816",
            "What it supports": "Within-cohort stability of the SRRS-healing association.",
            "What it cannot prove": "External validation in another DFU outcome cohort.",
            "How handled in manuscript": "Reported as robustness rather than validation.",
        },
        {
            "Evidence layer": "Random gene-panel and housekeeping controls",
            "Dataset/resource": "GSE165816",
            "What it supports": "SRRS modules outperform random expressed-gene panels and are not explained by housekeeping activity.",
            "What it cannot prove": "Absence of all confounding or causal mechanism.",
            "How handled in manuscript": "Used to reduce arbitrary-gene-set concern.",
        },
        {
            "Evidence layer": "Ligand-excluded SRRS sensitivity",
            "Dataset/resource": "GSE165816, GSE223964, GSE241124",
            "What it supports": "Candidate ligand findings are not a trivial consequence of including the ligand panel in SRRS.",
            "What it cannot prove": "Ligands are causal therapeutic targets.",
            "How handled in manuscript": "Presented as a non-circular sensitivity analysis.",
        },
        {
            "Evidence layer": "Acute wound reference alignment",
            "Dataset/resource": "GSE241132",
            "What it supports": "Healing-associated DFU fibroblasts align with a published human acute wound repair-phase programme.",
            "What it cannot prove": "DFU fully recapitulates normal acute wound healing or that the reference is a DFU validation cohort.",
            "How handled in manuscript": "Described as reference alignment, not validation.",
        },
        {
            "Evidence layer": "External chronic wound state conservation",
            "Dataset/resource": "GSE223964",
            "What it supports": "SRRS-high fibroblasts in an independent chronic wound dataset retain the candidate ligand programme.",
            "What it cannot prove": "Independent prediction of DFU healing outcome.",
            "How handled in manuscript": "Called state conservation rather than outcome validation.",
        },
        {
            "Evidence layer": "Acute wound spatial tissue context",
            "Dataset/resource": "GSE241124",
            "What it supports": "SRRS localises to repair-phase acute wound tissue and correlates with stromal/vascular marker programmes.",
            "What it cannot prove": "DFU spatial validation or direct cell-cell signalling.",
            "How handled in manuscript": "Donor-first spatial contextualisation; spot-level models treated as sensitivity support.",
        },
        {
            "Evidence layer": "Spatial ligand-receptor and receptor-side analyses",
            "Dataset/resource": "GSE241124",
            "What it supports": "Tiered candidate communication modules, strongest for TNC-integrin/syndecan and ligand-led INHBA/TGF-family.",
            "What it cannot prove": "Protein-level receptor activation or causal signalling.",
            "How handled in manuscript": "Candidate communication axes with explicit tiering.",
        },
        {
            "Evidence layer": "Perturbation-informed enrichment",
            "Dataset/resource": "GEO/LINCS perturbation libraries",
            "What it supports": "Biological plausibility of inducible TGF, TNF, IL6/FGF and fibroblast-related programmes.",
            "What it cannot prove": "DFU-specific functional perturbation or therapeutic efficacy.",
            "How handled in manuscript": "Framed as plausibility support, not functional validation.",
        },
    ]
    tbl = pd.DataFrame(rows)
    tbl.to_csv(OUT_DIR / "SRRS_evidence_hierarchy_table.tsv", sep="\t", index=False)
    tbl.to_excel(SUPP_DIR / "Supplementary_Table_Evidence_Hierarchy.xlsx", index=False)
    md = ["# Evidence Hierarchy Table", ""]
    md.append("| Evidence layer | Dataset/resource | What it supports | What it cannot prove |")
    md.append("|---|---|---|---|")
    for _, row in tbl.iterrows():
        md.append(f"| {row['Evidence layer']} | {row['Dataset/resource']} | {row['What it supports']} | {row['What it cannot prove']} |")
    (MAN_DIR / "05_Evidence_Hierarchy_Table.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return tbl


def append_supplementary_sheet():
    xlsx = SUPP_DIR / "Supplementary_Tables_SRRS_Framework.xlsx"
    if not xlsx.exists():
        return
    with pd.ExcelWriter(xlsx, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        read_tsv(RESULTS / "ligand_excluded_srrs" / "GSE165816_ligand_panel_excluded_existing_tests.tsv").to_excel(writer, sheet_name="S16_ligand_excluded_DFU", index=False)
        read_tsv(RESULTS / "ligand_excluded_srrs" / "GSE223964_ligand_excluded_SRRS_high_low_ligand_tests.tsv").to_excel(writer, sheet_name="S17_ligand_excluded_ext", index=False)
        read_tsv(RESULTS / "ligand_excluded_srrs" / "GSE241124_spatial_ligand_panel_excluded_donor_tests.tsv").to_excel(writer, sheet_name="S18_ligand_excluded_spatial", index=False)
        read_tsv(OUT_DIR / "SRRS_evidence_hierarchy_table.tsv").to_excel(writer, sheet_name="S19_evidence_hierarchy", index=False)

    index_path = SUPP_DIR / "Supplementary_Table_Index.tsv"
    index = pd.read_csv(index_path, sep="\t") if index_path.exists() else pd.DataFrame(columns=["sheet", "source_file"])
    extra = pd.DataFrame([
        {"sheet": "S16_ligand_excluded_DFU", "source_file": "results/srrs_framework/ligand_excluded_srrs/GSE165816_ligand_panel_excluded_existing_tests.tsv"},
        {"sheet": "S17_ligand_excluded_ext", "source_file": "results/srrs_framework/ligand_excluded_srrs/GSE223964_ligand_excluded_SRRS_high_low_ligand_tests.tsv"},
        {"sheet": "S18_ligand_excluded_spatial", "source_file": "results/srrs_framework/ligand_excluded_srrs/GSE241124_spatial_ligand_panel_excluded_donor_tests.tsv"},
        {"sheet": "S19_evidence_hierarchy", "source_file": "results/srrs_framework/evidence_hierarchy/SRRS_evidence_hierarchy_table.tsv"},
    ])
    index = pd.concat([index[~index["sheet"].isin(extra["sheet"])], extra], ignore_index=True)
    index.to_csv(index_path, sep="\t", index=False)


def main():
    donor_first_figure5()
    evidence_hierarchy()
    append_supplementary_sheet()
    print("Wrote donor-first Figure 5 and evidence hierarchy table")


if __name__ == "__main__":
    main()
