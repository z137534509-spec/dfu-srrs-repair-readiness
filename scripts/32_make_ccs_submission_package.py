from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


PROJECT = Path.cwd()
RESULTS = PROJECT / "results" / "srrs_framework"
FIG_DIR = PROJECT / "figures" / "ccs_submission"
MAN_DIR = PROJECT / "manuscript" / "Cell_Communication_Signaling"
SUPP_DIR = MAN_DIR / "supplementary_tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
SUPP_DIR.mkdir(parents=True, exist_ok=True)

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7E57C2"
GREY = "#6E6E6E"
LIGHT = "#F5F6F8"


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def savefig(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def letter(ax, label):
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")


def nice_feature(name):
    mapping = {
        "SRRS": "SRRS",
        "fibroblast_repair_activation": "Fibroblast repair",
        "stromal_ligand_panel": "Stromal ligands",
        "vascular_perivascular_receiver_coupling": "Receiver coupling",
        "d7_acute_wound_alignment": "D7 alignment",
        "negative_housekeeping_control_z": "Housekeeping control",
        "SRRS_4_equal": "4-component SRRS",
        "SRRS_3_no_receiver": "No receiver",
        "SRRS_receiver_10pct": "Receiver 10%",
        "SRRS_receiver_20pct": "Receiver 20%",
        "SRRS_2_fibroblast_ligand": "Fibroblast+ligand",
        "SRRS_3_no_d7": "No D7",
        "SRRS_3_no_ligand": "No ligand",
    }
    return mapping.get(name, name.replace("_", " "))


def boxplot_with_points(ax, df, feature, ylabel=None):
    order = ["Non-healer", "Healer"]
    vals = [df.loc[df["healing_status"] == g, feature].dropna().to_numpy() for g in order]
    bp = ax.boxplot(vals, tick_labels=order, patch_artist=True, widths=0.55, showfliers=False)
    for patch, color in zip(bp["boxes"], [ORANGE, BLUE]):
        patch.set_facecolor(color)
        patch.set_alpha(0.18)
        patch.set_edgecolor(color)
    for med in bp["medians"]:
        med.set_color("#222222")
    rng = np.random.default_rng(2)
    for i, (g, color) in enumerate(zip(order, [ORANGE, BLUE]), start=1):
        y = df.loc[df["healing_status"] == g, feature].dropna().to_numpy()
        x = rng.normal(i, 0.045, size=len(y))
        ax.scatter(x, y, s=32, color=color, edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel or nice_feature(feature))
    ax.spines[["top", "right"]].set_visible(False)


def errorbar_effects(ax, df, label_col, effect_col, lo_col, hi_col, color=BLUE, title=None, xlab="Delta"):
    plot = df.copy()
    plot = plot.sort_values(effect_col)
    y = np.arange(len(plot))
    ax.errorbar(
        plot[effect_col],
        y,
        xerr=[plot[effect_col] - plot[lo_col], plot[hi_col] - plot[effect_col]],
        fmt="o",
        color=color,
        ecolor="#808080",
        capsize=2,
        markersize=5,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(plot[label_col])
    ax.set_xlabel(xlab)
    if title:
        ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)


def figure1_framework():
    fig, ax = plt.subplots(figsize=(11.5, 6.6))
    ax.set_axis_off()

    def box(x, y, w, h, title, body, color):
        patch = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.02,rounding_size=0.025",
            linewidth=1.2,
            edgecolor=color,
            facecolor="white",
        )
        ax.add_patch(patch)
        ax.text(x + 0.03, y + h - 0.08, title, fontsize=11, fontweight="bold", color=color, va="top")
        ax.text(x + 0.03, y + h - 0.18, textwrap.fill(body, 28), fontsize=9.2, va="top", color="#222222")

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13, linewidth=1.2, color="#555555"))

    ax.text(0.04, 0.94, "Fibroblast-centred stromal repair-readiness framework", fontsize=16, fontweight="bold")
    ax.text(0.04, 0.89, "Locked SRRS formula projected across public human wound single-cell and spatial datasets", fontsize=10.5, color="#333333")

    box(0.05, 0.58, 0.21, 0.22, "Discovery outcome", "GSE165816 DFU scRNA-seq\n9 healers vs 5 non-healers", BLUE)
    box(0.39, 0.58, 0.23, 0.22, "Locked SRRS", "Fibroblast repair + stromal ligands + receiver context + D7 acute wound alignment", PURPLE)
    box(0.74, 0.58, 0.21, 0.22, "Primary output", "Repair-readiness state and tiered communication axes", GREEN)
    arrow(0.26, 0.69, 0.39, 0.69)
    arrow(0.62, 0.69, 0.74, 0.69)

    box(0.05, 0.23, 0.21, 0.23, "External chronic wound", "GSE223964 scRNA-seq\nSRRS-high vs low fibroblast conservation", ORANGE)
    box(0.31, 0.23, 0.21, 0.23, "Spatial context", "GSE241124 Visium\nWound7 tissue localisation and marker correlation", GREEN)
    box(0.57, 0.23, 0.21, 0.23, "Communication evidence", "Spatial LR axes and individual receptor-side models", BLUE)
    box(0.82, 0.23, 0.14, 0.23, "Controls", "Random panels\nLODO\nPermutation\nPerturbation", GREY)
    arrow(0.50, 0.58, 0.16, 0.46)
    arrow(0.50, 0.58, 0.42, 0.46)
    arrow(0.50, 0.58, 0.68, 0.46)
    arrow(0.50, 0.58, 0.89, 0.46)

    ax.text(0.05, 0.09, "Claim boundary", fontsize=11, fontweight="bold", color="#222222")
    ax.text(0.18, 0.09, "Candidate communication framework; not independent DFU outcome validation or causal perturbation.", fontsize=10, color="#222222")
    savefig(fig, "Figure1_Multicohort_SRRS_framework")


def figure2_discovery():
    sample = read_tsv(RESULTS / "GSE165816_SRRS_discovery_sample_scores.tsv")
    tests = read_tsv(RESULTS / "GSE165816_SRRS_discovery_healer_vs_nonhealer_tests.tsv")
    sens = read_tsv(RESULTS / "srrs_component_sensitivity" / "GSE165816_SRRS_variant_healer_nonhealer_tests.tsv")
    rand = read_tsv(RESULTS / "GSE165816_SRRS_random_gene_panel_controls.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.0))
    ax = axes[0, 0]
    boxplot_with_points(ax, sample, "SRRS", "Sample-level SRRS")
    ax.set_title("Discovery DFU outcome cohort")
    ax.text(0.55, 0.97, "Delta = 1.319\nExact P = 0.0015", transform=ax.transAxes, ha="left", va="top", fontsize=9)
    letter(ax, "A")

    ax = axes[0, 1]
    sub = tests[tests["feature"].isin(["fibroblast_repair_activation", "stromal_ligand_panel", "d7_acute_wound_alignment", "vascular_perivascular_receiver_coupling", "negative_housekeeping_control_z"])].copy()
    sub["label"] = sub["feature"].map(nice_feature)
    errorbar_effects(ax, sub, "label", "delta_healer_minus_nonhealer", "bootstrap_ci_low", "bootstrap_ci_high", color=BLUE, title="SRRS components", xlab="Healer minus non-healer")
    letter(ax, "B")

    ax = axes[1, 0]
    sens["label"] = sens["feature"].map(nice_feature)
    keep = ["SRRS_4_equal", "SRRS_3_no_receiver", "SRRS_receiver_10pct", "SRRS_receiver_20pct", "SRRS_2_fibroblast_ligand"]
    sub = sens[sens["feature"].isin(keep)].copy()
    errorbar_effects(ax, sub, "label", "delta_healer_minus_nonhealer", "bootstrap_ci_low", "bootstrap_ci_high", color=PURPLE, title="Component-weight sensitivity", xlab="Healer minus non-healer")
    letter(ax, "C")

    ax = axes[1, 1]
    rand = rand.copy()
    rand["label"] = rand["module"].map(nice_feature)
    y = np.arange(len(rand))
    observed_col = "observed_delta" if "observed_delta" in rand.columns else "observed_delta_healer_minus_nonhealer"
    ax.scatter(rand["random_mean_delta"], y, color="#888888", s=48, label="Random mean")
    ax.scatter(rand[observed_col], y, color=ORANGE, s=58, label="Observed")
    for _, row in rand.iterrows():
        yi = list(rand.index).index(row.name)
        ax.plot([row["random_mean_delta"], row[observed_col]], [yi, yi], color="#BBBBBB", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(rand["label"])
    ax.set_xlabel("Healer minus non-healer delta")
    ax.set_title("Random gene-panel controls")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "D")
    fig.tight_layout()
    savefig(fig, "Figure2_Discovery_SRRS_and_controls")


def figure3_alignment():
    align = read_tsv(RESULTS / "d7_alignment_robustness" / "GSE165816_alternative_alignment_sample_scores.tsv")
    tests = read_tsv(RESULTS / "d7_alignment_robustness" / "GSE165816_alternative_alignment_healer_nonhealer_tests.tsv")
    old = read_tsv(PROJECT / "results" / "human_wound_roadmap_gse241132" / "GSE165816_DFU_healer_vs_nonhealer_GSE241132_signature_tests.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))
    ax = axes[0]
    boxplot_with_points(ax, align, "D7_vs_D0_top100_original", "D7 vs D0 alignment score")
    ax.text(0.53, 0.96, "Delta = 0.251\nExact P = 0.0225", transform=ax.transAxes, ha="left", va="top", fontsize=9)
    ax.set_title("Published D7 acute wound alignment")
    letter(ax, "A")

    ax = axes[1]
    tests["label"] = tests["feature"].str.replace("_top100", "", regex=False).str.replace("_", " ", regex=False)
    errorbar_effects(ax, tests, "label", "delta_healer_minus_nonhealer", "bootstrap_ci_low", "bootstrap_ci_high", color=GREEN, title="Alternative acute wound signatures", xlab="Healer minus non-healer")
    letter(ax, "B")

    ax = axes[2]
    heat = old[old["feature"].isin([
        "GSE241132_D1_wound_vs_D0_skin_fibroblast_up",
        "GSE241132_D7_wound_vs_D0_skin_fibroblast_up",
        "GSE241132_D30_wound_vs_D0_skin_fibroblast_up",
        "mean_repair_axis",
    ])].copy()
    heat["label"] = heat["feature"].str.replace("GSE241132_", "", regex=False).str.replace("_fibroblast_up", "", regex=False).str.replace("_", " ", regex=False)
    y = np.arange(len(heat))
    ax.barh(y, heat["delta_healer_minus_nonhealer"], color=[BLUE if x != "mean repair axis" else PURPLE for x in heat["label"]])
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(heat["label"], fontsize=8)
    ax.set_xlabel("Healer minus non-healer")
    ax.set_title("Reference wound-programme scores")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "C")
    fig.tight_layout()
    savefig(fig, "Figure3_Acute_wound_alignment")


def figure4_external():
    tests = read_tsv(RESULTS / "external_gse223964_rescue" / "GSE223964_SRRS_high_low_within_sample_ligand_tests.tsv")
    deltas = read_tsv(RESULTS / "external_gse223964_rescue" / "GSE223964_SRRS_high_low_within_sample_ligand_deltas.tsv")
    cond = read_tsv(RESULTS / "external_gse223964_rescue" / "GSE223964_fibroblast_specific_condition_tests.tsv")
    rand = read_tsv(RESULTS / "external_gse223964_rescue" / "GSE223964_SRRS_within_sample_random_panel_control.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    ax = axes[0]
    plot = tests.sort_values("mean_delta_high_minus_low")
    colors = [PURPLE if f == "candidate_ligand_panel" else BLUE for f in plot["feature"]]
    ax.barh(plot["feature"], plot["mean_delta_high_minus_low"], color=colors)
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_xlabel("SRRS-high minus SRRS-low fibroblast delta")
    ax.set_title("Within-sample ligand conservation")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "A")

    ax = axes[1]
    panel = deltas[deltas["feature"].isin(["candidate_ligand_panel", "IL11", "IL6", "INHBA", "THBS1", "SERPINE1", "TNC"])].copy()
    samples = sorted(panel["sample_id"].unique()) if "sample_id" in panel.columns else sorted(panel.iloc[:, 0].unique())
    sample_col = "sample_id" if "sample_id" in panel.columns else panel.columns[0]
    pivot = panel.pivot_table(index="feature", columns=sample_col, values="delta_high_minus_low", aggfunc="mean")
    im = ax.imshow(pivot, aspect="auto", cmap="RdBu_r", vmin=-np.nanmax(abs(pivot.to_numpy())), vmax=np.nanmax(abs(pivot.to_numpy())))
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_title("Paired deltas by sample")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    letter(ax, "B")

    ax = axes[2]
    observed_col = "observed_candidate_ligand_delta" if "observed_candidate_ligand_delta" in rand.columns else "observed_candidate_ligand_within_sample_delta"
    observed = float(rand[observed_col].iloc[0])
    random_mean = float(rand["random_mean_delta"].iloc[0])
    ax.bar(["Observed", "Random mean"], [observed, random_mean], color=[ORANGE, "#A0A0A0"])
    ax.set_ylabel("Candidate-panel delta")
    ax.set_title("External random-panel control")
    ax.text(0.02, 0.95, "Empirical P = 0.001", transform=ax.transAxes, va="top", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "C")
    fig.tight_layout()
    savefig(fig, "Figure4_External_chronic_wound_conservation")


def figure5_spatial():
    paired = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_paired_donor_bootstrap_tests.tsv")
    models = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_spot_level_clustered_models.tsv")
    corr = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_marker_correlation_donor_bootstrap_summary.tsv")
    sample = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_strengthened_sample_summary.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2))
    ax = axes[0, 0]
    order = ["Skin", "Wound1", "Wound7", "Wound30"]
    vals = [sample.loc[sample["Condition"] == c, "spatial_SRRS"].to_numpy() for c in order]
    bp = ax.boxplot(vals, tick_labels=order, patch_artist=True, widths=0.55, showfliers=False)
    for patch, color in zip(bp["boxes"], ["#999999", GREEN, ORANGE, BLUE]):
        patch.set_facecolor(color)
        patch.set_alpha(0.18)
        patch.set_edgecolor(color)
    for donor, sub in sample.groupby("Donor"):
        sub = sub.set_index("Condition").reindex(order)
        ax.plot(np.arange(1, len(order) + 1), sub["spatial_SRRS"], color="#BBBBBB", linewidth=0.9, zorder=1)
        ax.scatter(np.arange(1, len(order) + 1), sub["spatial_SRRS"], color="#333333", s=18, zorder=2)
    ax.set_ylabel("Sample-level spatial SRRS")
    ax.set_title("Paired spatial samples")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "A")

    ax = axes[0, 1]
    sub = paired[(paired["feature"].isin(["spatial_SRRS", "spatial_SRRS_resid"])) & (paired["comparison"].isin(["Wound1_vs_Skin", "Wound7_vs_Skin", "Wound30_vs_Skin"]))].copy()
    sub["label"] = sub["feature"].str.replace("spatial_SRRS", "SRRS", regex=False) + " " + sub["comparison"].str.replace("_vs_Skin", "", regex=False)
    errorbar_effects(ax, sub, "label", "mean_delta", "bootstrap_ci_low", "bootstrap_ci_high", color=GREEN, title="Paired donor bootstrap", xlab="Delta vs skin")
    letter(ax, "B")

    ax = axes[1, 0]
    sub = models[models["term"].str.contains("Wound7")].copy()
    sub["model"] = pd.Categorical(sub["model"], categories=["unadjusted", "adjusted", "residual"], ordered=True)
    sub = sub.sort_values("model")
    ax.bar(sub["model"].astype(str), sub["estimate"], color=[BLUE, ORANGE, GREEN])
    ax.set_ylabel("Wound7 vs Skin coefficient")
    ax.set_title("Clustered spot-level models")
    ax.text(0.02, 0.95, "Adjusted beta = 0.441\nP = 3.08e-11", transform=ax.transAxes, va="top", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "C")

    ax = axes[1, 1]
    sub = corr[corr["Condition"] == "Wound7"].copy()
    sub["label"] = sub["feature"].str.replace("_marker", "", regex=False).str.replace("_", " ", regex=False)
    errorbar_effects(ax, sub, "label", "mean_donor_rho", "bootstrap_ci_low", "bootstrap_ci_high", color=PURPLE, title="Wound7 marker correlations", xlab="Mean donor Spearman rho")
    letter(ax, "D")
    fig.tight_layout()
    savefig(fig, "Figure5_Spatial_SRRS_localisation")


def figure6_communication():
    axes_df = read_tsv(RESULTS / "spatial_lr_axis_gse241124" / "GSE241124_spatial_LR_axis_condition_clustered_models.tsv")
    corr = read_tsv(RESULTS / "spatial_lr_axis_gse241124" / "GSE241124_spatial_LR_axis_correlation_donor_bootstrap_summary.tsv")
    gene = read_tsv(RESULTS / "spatial_receptor_expression_gse241124" / "GSE241124_spatial_LR_gene_condition_clustered_models.tsv")
    pert = read_tsv(RESULTS / "perturbation_enrichr" / "SRRS_Enrichr_focused_perturbation_hits.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.6))
    ax = axes[0, 0]
    sub = axes_df[axes_df["term"].str.contains("Wound7")].sort_values("estimate")
    ax.barh(sub["axis"], sub["estimate"], color=[ORANGE if p < 0.1 else "#9E9E9E" for p in sub["padj_bh"]])
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Adjusted Wound7 vs Skin axis-score coefficient")
    ax.set_title("Spatial LR axis scores")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "A")

    ax = axes[0, 1]
    sub = corr[corr["Condition"] == "Wound7"].sort_values("mean_donor_rho")
    errorbar_effects(ax, sub, "axis", "mean_donor_rho", "bootstrap_ci_low", "bootstrap_ci_high", color=BLUE, title="Wound7 ligand-receptor coexpression", xlab="Mean donor Spearman rho")
    letter(ax, "B")

    ax = axes[1, 0]
    w7 = gene[gene["term"].str.contains("Wound7")].copy()
    keep = ["TNC", "ITGB1", "ITGA5", "ITGAV", "INHBA", "BAMBI", "ENG", "IL11", "IL11RA", "IL6", "IL6R", "IL6ST", "CCL20", "CCR6"]
    w7 = w7[w7["gene"].isin(keep)].copy()
    w7 = w7.sort_values("estimate")
    colors = np.where(w7["role"].eq("ligand"), ORANGE, BLUE)
    ax.barh(w7["gene"] + " (" + w7["role"] + ")", w7["estimate"], color=colors)
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Adjusted Wound7 vs Skin gene coefficient")
    ax.set_title("Individual ligand/receptor expression")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "C")

    ax = axes[1, 1]
    p = pert.copy()
    # Prefer focused biologically named hits with adjusted P values.
    sort_col = "Adjusted P-value" if "Adjusted P-value" in p.columns else ("adjusted_p_value" if "adjusted_p_value" in p.columns else [c for c in p.columns if "adjusted" in c.lower()][0])
    term_col = "Term" if "Term" in p.columns else [c for c in p.columns if c.lower() == "term"][0]
    p = p.sort_values(sort_col).head(10).copy()
    p["neglog10"] = -np.log10(p[sort_col].astype(float).clip(lower=1e-300))
    p["short_term"] = p[term_col].astype(str).str.replace("_", " ", regex=False).str.slice(0, 44)
    ax.barh(p["short_term"][::-1], p["neglog10"][::-1], color=GREEN)
    ax.set_xlabel("-log10 adjusted P")
    ax.set_title("Perturbation-informed enrichment")
    ax.spines[["top", "right"]].set_visible(False)
    letter(ax, "D")
    fig.tight_layout()
    savefig(fig, "Figure6_Communication_and_perturbation")


def figure7_clinical_model():
    fig, ax = plt.subplots(figsize=(11.5, 6.8))
    ax.set_axis_off()
    ax.text(0.04, 0.94, "Working model: stromal repair readiness in chronic wound healing", fontsize=16, fontweight="bold")
    ax.text(0.04, 0.89, "Conceptual summary for prospective biopsy-based validation", fontsize=10.5, color="#333333")

    def panel(x, y, w, h, title, color, lines):
        patch = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.02,rounding_size=0.025",
            linewidth=1.4,
            edgecolor=color,
            facecolor="white",
        )
        ax.add_patch(patch)
        ax.text(x + 0.03, y + h - 0.06, title, fontsize=12, fontweight="bold", color=color, va="top")
        yy = y + h - 0.15
        for line in lines:
            ax.text(x + 0.04, yy, f"- {line}", fontsize=9.4, va="top", color="#222222")
            yy -= 0.075

    panel(
        0.05, 0.48, 0.40, 0.31,
        "Repair-ready wound bed",
        GREEN,
        [
            "High fibroblast repair activation",
            "TNC-integrin/syndecan matrix module",
            "INHBA/TGF-family ligand programme",
            "IL11/IL6-gp130 candidate activity",
            "ECM-vascular/perivascular coupling",
        ],
    )
    panel(
        0.55, 0.48, 0.40, 0.31,
        "Stalled wound bed",
        ORANGE,
        [
            "Low SRRS or incomplete repair mobilisation",
            "Weaker acute wound D7 alignment",
            "Reduced stromal ligand programme",
            "Uncoupled inflammation and repair",
            "Candidate non-healing molecular state",
        ],
    )

    # Stylised wound-bed signal arrows.
    for x, color, label in [(0.25, GREEN, "fibroblast -> matrix -> vessel"), (0.75, ORANGE, "inflammation without repair coupling")]:
        ax.plot([x - 0.14, x + 0.14], [0.40, 0.40], color=color, linewidth=3, solid_capstyle="round")
        ax.scatter([x - 0.12, x, x + 0.12], [0.40, 0.40, 0.40], s=[130, 180, 130], color=[color, "#FFFFFF", color], edgecolor=color, linewidth=2)
        ax.text(x, 0.35, label, ha="center", fontsize=9.2, color="#333333")

    # Workflow band.
    y = 0.18
    steps = [
        ("Wound biopsy", BLUE),
        ("Single-cell/spatial\nor targeted panel", PURPLE),
        ("SRRS scoring", GREEN),
        ("Repair-readiness\nstratification", ORANGE),
        ("Prospective and\nfunctional validation", GREY),
    ]
    x_positions = np.linspace(0.08, 0.84, len(steps))
    for i, ((label, color), x) in enumerate(zip(steps, x_positions)):
        patch = FancyBboxPatch(
            (x, y), 0.13, 0.105,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            linewidth=1.2,
            edgecolor=color,
            facecolor=LIGHT,
        )
        ax.add_patch(patch)
        ax.text(x + 0.065, y + 0.074, label, ha="center", va="top", fontsize=8.8, color="#222222")
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + 0.135, y + 0.052), (x_positions[i + 1] - 0.008, y + 0.052), arrowstyle="-|>", mutation_scale=12, color="#555555", linewidth=1.1))

    ax.text(0.05, 0.065, "Current evidence: computational prioritisation across public datasets.", fontsize=9.5, color="#333333")
    ax.text(0.05, 0.035, "Required next step: prospective DFU outcome cohort with protein-level and functional validation.", fontsize=9.5, color="#333333")
    savefig(fig, "Figure7_Clinical_repair_readiness_model")


def supplementary_tables():
    xlsx = SUPP_DIR / "Supplementary_Tables_SRRS_Framework.xlsx"
    tables = {
        "S1_locked_gene_sets": RESULTS / "SRRS_locked_gene_sets.tsv",
        "S2_discovery_scores": RESULTS / "GSE165816_SRRS_discovery_sample_scores.tsv",
        "S3_discovery_tests": RESULTS / "GSE165816_SRRS_discovery_healer_vs_nonhealer_tests.tsv",
        "S4_component_sensitivity": RESULTS / "srrs_component_sensitivity" / "GSE165816_SRRS_variant_healer_nonhealer_tests.tsv",
        "S5_D7_robustness": RESULTS / "d7_alignment_robustness" / "GSE165816_alternative_alignment_healer_nonhealer_tests.tsv",
        "S6_random_controls": RESULTS / "GSE165816_SRRS_random_gene_panel_controls.tsv",
        "S7_GSE223964_ligands": RESULTS / "external_gse223964_rescue" / "GSE223964_SRRS_high_low_within_sample_ligand_tests.tsv",
        "S8_GSE223964_condition": RESULTS / "external_gse223964_rescue" / "GSE223964_fibroblast_specific_condition_tests.tsv",
        "S9_spatial_sample": RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_strengthened_sample_summary.tsv",
        "S10_spatial_models": RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_spot_level_clustered_models.tsv",
        "S11_spatial_marker_corr": RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_marker_correlation_donor_bootstrap_summary.tsv",
        "S12_spatial_LR_axes": RESULTS / "spatial_lr_axis_gse241124" / "GSE241124_spatial_LR_axis_condition_clustered_models.tsv",
        "S13_spatial_LR_corr": RESULTS / "spatial_lr_axis_gse241124" / "GSE241124_spatial_LR_axis_correlation_donor_bootstrap_summary.tsv",
        "S14_gene_receptor_models": RESULTS / "spatial_receptor_expression_gse241124" / "GSE241124_spatial_LR_gene_condition_clustered_models.tsv",
        "S15_perturbation_hits": RESULTS / "perturbation_enrichr" / "SRRS_Enrichr_focused_perturbation_hits.tsv",
    }
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        for sheet, path in tables.items():
            df = read_tsv(path)
            if len(df) > 50000:
                df = df.head(50000)
            df.to_excel(writer, sheet_name=sheet[:31], index=False)

    index_rows = [{"sheet": sheet[:31], "source_file": str(path.relative_to(PROJECT))} for sheet, path in tables.items()]
    pd.DataFrame(index_rows).to_csv(SUPP_DIR / "Supplementary_Table_Index.tsv", sep="\t", index=False)


def figure_legends():
    text = """# Figure Legends

## Figure 1. Multi-cohort SRRS framework

Overview of the study design. GSE165816 was used as the DFU outcome discovery cohort to define a locked Stromal Repair Readiness Score (SRRS). The framework was then evaluated by independent chronic wound single-cell projection (GSE223964), acute wound spatial transcriptomics (GSE241124), spatial ligand-receptor analysis, random-panel controls and perturbation-informed enrichment. The intended claim is a candidate communication framework, not independent DFU outcome validation or causal perturbation.

## Figure 2. SRRS separates healing and non-healing DFU samples in the discovery cohort

(A) Sample-level SRRS in healer and non-healer DFU samples from GSE165816. (B) Healer-minus-non-healer effect sizes for SRRS components with bootstrap confidence intervals. (C) Component-weight sensitivity showing that the healer/non-healer signal persists when the receiver-coupling component is removed or downweighted. (D) Observed module effects compared with random expressed-gene panels.

## Figure 3. DFU healing aligns with a published acute wound repair programme

(A) Projection of the original GSE241132 D7-vs-D0 fibroblast alignment signature into GSE165816 DFU samples. (B) Robustness across alternative acute wound alignment definitions, including D7-vs-D1, D7-vs-D30, D7/D30-vs-D0 and D7 consensus signatures. (C) Reference wound-programme score deltas in healer versus non-healer DFU samples.

## Figure 4. SRRS-high fibroblast ligand programmes are conserved in independent chronic wound scRNA-seq

(A) Paired within-sample SRRS-high minus SRRS-low fibroblast ligand differences in GSE223964. (B) Heatmap of sample-level paired ligand deltas. (C) Candidate ligand-panel difference compared with random fibroblast-expressed panels. GSE223964 was used to test fibroblast-state conservation rather than outcome prediction.

## Figure 5. Spatial transcriptomics localises SRRS to repair-active wound tissue

(A) Paired donor-level spatial SRRS across intact skin and wound time points in GSE241124. (B) Bootstrap paired donor effects for raw and residualised SRRS. (C) Clustered spot-level Wound7-vs-skin model coefficients. (D) Wound7 donor-level correlations between spatial SRRS and ECM, fibroblast, endothelial and pericyte/smooth-muscle marker programmes.

## Figure 6. Spatial communication and perturbation-informed prioritisation

(A) Adjusted Wound7-vs-skin spatial ligand-receptor axis coefficients. (B) Wound7 donor-level ligand-receptor spatial coexpression. (C) Individual ligand and receptor Wound7-vs-skin expression coefficients, showing strongest receptor-side support for TNC-integrin/syndecan and limited spatial gp130 receptor-side support for IL11/IL6. (D) Perturbation-informed enrichment of SRRS genes in GEO/LINCS signatures.

## Figure 7. Clinical repair-readiness model

Conceptual model contrasting a repair-ready wound bed with coordinated fibroblast, matrix and vascular/perivascular programmes against a stalled wound bed with incomplete stromal repair mobilisation. The proposed translational workflow links wound biopsy, single-cell/spatial or targeted profiling, SRRS scoring, repair-readiness stratification and prospective functional validation.
"""
    path = MAN_DIR / "02_Figure_Legends.md"
    path.write_text(text, encoding="utf-8")


def package_readme():
    figures = sorted(FIG_DIR.glob("Figure*.png"))
    lines = [
        "# CCS Submission Package",
        "",
        "Generated from existing SRRS framework results.",
        "",
        "## Main Figures",
        "",
    ]
    for fig in figures:
        lines.append(f"- `{fig.relative_to(PROJECT)}`")
    lines += [
        "",
        "## Supplementary Tables",
        "",
        f"- `{(SUPP_DIR / 'Supplementary_Tables_SRRS_Framework.xlsx').relative_to(PROJECT)}`",
        f"- `{(SUPP_DIR / 'Supplementary_Table_Index.tsv').relative_to(PROJECT)}`",
        "",
        "## Notes",
        "",
        "- GSE223964 is used for fibroblast-state conservation, not independent DFU outcome validation.",
        "- Spatial LR evidence is strongest for TNC-integrin/syndecan; INHBA is best framed as a ligand-led TGF-family module with mixed receptor-side behaviour.",
        "- IL11/IL6-gp130 should remain candidate-level because spatial receptor-side support is limited.",
    ]
    (MAN_DIR / "03_CCS_Submission_Package_Readme.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    figure1_framework()
    figure2_discovery()
    figure3_alignment()
    figure4_external()
    figure5_spatial()
    figure6_communication()
    figure7_clinical_model()
    supplementary_tables()
    figure_legends()
    package_readme()
    print(f"Wrote CCS figures to {FIG_DIR}")
    print(f"Wrote supplementary tables to {SUPP_DIR}")


if __name__ == "__main__":
    main()
