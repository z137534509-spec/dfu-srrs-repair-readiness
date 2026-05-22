from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle


PROJECT = Path.cwd()
RESULTS = PROJECT / "results" / "srrs_framework"
FIG_DIR = PROJECT / "figures" / "ccs_submission"
FIG_DIR.mkdir(parents=True, exist_ok=True)

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7E57C2"
GREY = "#6E6E6E"
LIGHT = "#F6F7F9"
INK = "#222222"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8,
        "figure.titlesize": 13,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def savefig(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_label(ax, label):
    ax.text(
        -0.06,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def clean_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.6, zorder=0)


def wrap_label(text, width=18):
    return "\n".join(textwrap.wrap(str(text), width=width, break_long_words=False))


def trim_white_image(img, pad=25, threshold=0.985):
    """Trim empty white margins from imported raster panels."""
    arr = np.asarray(img)
    if arr.ndim == 2:
        rgb = arr
    else:
        rgb = arr[..., :3]
    mask = np.any(rgb < threshold, axis=-1) if rgb.ndim == 3 else rgb < threshold
    if not np.any(mask):
        return img
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    r0 = max(rows.min() - pad, 0)
    r1 = min(rows.max() + pad + 1, arr.shape[0])
    c0 = max(cols.min() - pad, 0)
    c1 = min(cols.max() + pad + 1, arr.shape[1])
    return arr[r0:r1, c0:c1]


def nice_feature(name):
    mapping = {
        "SRRS": "SRRS",
        "fibroblast_repair_activation": "Fibroblast\nrepair",
        "stromal_ligand_panel": "Stromal\nligands",
        "vascular_perivascular_receiver_coupling": "Receiver\ncontext",
        "d7_acute_wound_alignment": "D7 wound\nalignment",
        "negative_housekeeping_control_z": "Housekeeping\ncontrol",
        "SRRS_4_equal": "4-component\nSRRS",
        "SRRS_3_no_receiver": "No receiver\ncontext",
        "SRRS_receiver_10pct": "Receiver\n10%",
        "SRRS_receiver_20pct": "Receiver\n20%",
        "SRRS_2_fibroblast_ligand": "Fibroblast +\nligands",
    }
    return mapping.get(str(name), str(name).replace("_", " "))


def boxplot_points(ax, df, feature, ylabel, title=None):
    order = ["Non-healer", "Healer"]
    colors = [ORANGE, BLUE]
    vals = [df.loc[df["healing_status"] == g, feature].dropna().to_numpy() for g in order]
    bp = ax.boxplot(vals, tick_labels=order, widths=0.55, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.13)
        patch.set_edgecolor(color)
        patch.set_linewidth(1.2)
    for med in bp["medians"]:
        med.set_color("#222222")
        med.set_linewidth(1.2)
    rng = np.random.default_rng(11)
    for i, (g, color) in enumerate(zip(order, colors), 1):
        y = df.loc[df["healing_status"] == g, feature].dropna().to_numpy()
        x = rng.normal(i, 0.035, size=len(y))
        ax.scatter(x, y, s=26, color=color, edgecolor="white", linewidth=0.45, zorder=3)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, pad=8)
    clean_axes(ax)
    ax.grid(axis="x", visible=False)


def forest(ax, df, label_col, effect_col, lo_col, hi_col, color=BLUE, xlabel="Effect"):
    plot = df.copy().sort_values(effect_col)
    y = np.arange(len(plot))
    lo = plot[effect_col] - plot[lo_col]
    hi = plot[hi_col] - plot[effect_col]
    ax.errorbar(
        plot[effect_col],
        y,
        xerr=[lo, hi],
        fmt="o",
        color=color,
        ecolor="#7A7A7A",
        capsize=2.5,
        markersize=5,
        linewidth=1,
        zorder=3,
    )
    ax.axvline(0, color="#9A9A9A", linestyle="--", linewidth=0.9, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(plot[label_col])
    ax.set_xlabel(xlabel)
    clean_axes(ax)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    ax.grid(axis="y", visible=False)


def draw_box(ax, xy, wh, title, lines, color, title_size=10.5, body_size=8.4):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.015,rounding_size=0.025",
        linewidth=1.4,
        edgecolor=color,
        facecolor="white",
    )
    ax.add_patch(patch)
    ax.text(x + 0.025, y + h - 0.045, title, color=color, fontsize=title_size, fontweight="bold", va="top")
    yy = y + h - 0.105
    for line in lines:
        ax.text(x + 0.025, yy, line, color=INK, fontsize=body_size, va="top")
        yy -= 0.045


def draw_arrow(ax, start, end):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.2,
            color="#555555",
            shrinkA=4,
            shrinkB=4,
            connectionstyle="arc3,rad=0.0",
        )
    )


def figure1():
    fig, ax = plt.subplots(figsize=(12.0, 6.6))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.04, 0.94, "Stromal repair readiness framework", fontsize=15, fontweight="bold", color=INK)
    ax.text(
        0.04,
        0.895,
        "A locked SRRS formula is projected across public human wound single-cell and spatial datasets.",
        fontsize=9.5,
        color="#444444",
    )

    draw_box(ax, (0.06, 0.61), (0.22, 0.18), "Discovery", ["GSE165816 DFU scRNA-seq", "9 healers / 5 non-healers"], BLUE)
    draw_box(ax, (0.39, 0.61), (0.22, 0.18), "Locked SRRS", ["repair + ligands", "receiver + D7 alignment"], PURPLE)
    draw_box(ax, (0.72, 0.61), (0.22, 0.18), "Output", ["Repair-readiness state", "Tiered communication axes"], GREEN)
    draw_arrow(ax, (0.28, 0.70), (0.39, 0.70))
    draw_arrow(ax, (0.61, 0.70), (0.72, 0.70))

    ax.text(0.06, 0.535, "Supporting evidence layers", fontsize=10.5, fontweight="bold", color=INK)
    ax.text(
        0.06,
        0.505,
        "Each layer strengthens biological plausibility without being treated as independent DFU outcome validation.",
        fontsize=8.6,
        color="#444444",
        va="baseline",
    )

    y2 = 0.235
    boxes = [
        ((0.04, y2), (0.21, 0.22), "External state", ["GSE223964", "state conservation", "ligand panel retained"], ORANGE),
        ((0.28, y2), (0.21, 0.22), "Spatial context", ["GSE241124", "donor-first spatial", "stromal-vascular context"], GREEN),
        ((0.52, y2), (0.21, 0.22), "Communication", ["spatial LR axes", "receptor-side models", "TNC strongest"], BLUE),
        ((0.76, y2), (0.20, 0.22), "Robustness", ["permutation / LODO", "random panels", "ligand-excluded score"], GREY),
    ]
    for xy, wh, title, lines, color in boxes:
        draw_box(ax, xy, wh, title, lines, color, title_size=9.8, body_size=8.0)

    ax.add_patch(Rectangle((0.05, 0.095), 0.90, 0.105, facecolor=LIGHT, edgecolor="#D5D8DE", linewidth=1.0))
    ax.text(0.075, 0.158, "Claim boundary", fontsize=10.5, fontweight="bold", color=INK, va="center")
    ax.text(
        0.23,
        0.158,
        "Candidate communication framework; not independent DFU outcome validation and not causal perturbation.",
        fontsize=9,
        color=INK,
        va="center",
    )
    savefig(fig, "Figure1_Multicohort_SRRS_framework")


def figure2():
    sample = read_tsv(RESULTS / "GSE165816_SRRS_discovery_sample_scores.tsv")
    tests = read_tsv(RESULTS / "GSE165816_SRRS_discovery_healer_vs_nonhealer_tests.tsv")
    sens = read_tsv(RESULTS / "srrs_component_sensitivity" / "GSE165816_SRRS_variant_healer_nonhealer_tests.tsv")
    rand = read_tsv(RESULTS / "GSE165816_SRRS_random_gene_panel_controls.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.5), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.08, h_pad=0.08, wspace=0.18, hspace=0.16)

    ax = axes[0, 0]
    boxplot_points(ax, sample, "SRRS", "Sample-level SRRS", "Discovery cohort\nDelta 1.319; exact P = 0.0015")
    panel_label(ax, "A")

    ax = axes[0, 1]
    keep = [
        "fibroblast_repair_activation",
        "stromal_ligand_panel",
        "d7_acute_wound_alignment",
        "vascular_perivascular_receiver_coupling",
        "negative_housekeeping_control_z",
    ]
    sub = tests[tests["feature"].isin(keep)].copy()
    sub["label"] = sub["feature"].map(nice_feature)
    forest(ax, sub, "label", "delta_healer_minus_nonhealer", "bootstrap_ci_low", "bootstrap_ci_high", color=BLUE, xlabel="Healer minus non-healer")
    ax.set_title("SRRS components", pad=8)
    panel_label(ax, "B")

    ax = axes[1, 0]
    keep = ["SRRS_4_equal", "SRRS_3_no_receiver", "SRRS_receiver_10pct", "SRRS_receiver_20pct", "SRRS_2_fibroblast_ligand"]
    sub = sens[sens["feature"].isin(keep)].copy()
    sub["label"] = sub["feature"].map(nice_feature)
    forest(ax, sub, "label", "delta_healer_minus_nonhealer", "bootstrap_ci_low", "bootstrap_ci_high", color=PURPLE, xlabel="Healer minus non-healer")
    ax.set_title("Component-weight sensitivity", pad=8)
    panel_label(ax, "C")

    ax = axes[1, 1]
    rand = rand.copy()
    rand["label"] = rand["module"].map(nice_feature)
    obs_col = "observed_delta_healer_minus_nonhealer" if "observed_delta_healer_minus_nonhealer" in rand.columns else "observed_delta"
    y = np.arange(len(rand))
    ax.barh(y - 0.16, rand["random_mean_delta"], height=0.28, color="#B4B4B4", label="Random mean")
    ax.barh(y + 0.16, rand[obs_col], height=0.28, color=ORANGE, label="Observed")
    ax.set_yticks(y)
    ax.set_yticklabels(rand["label"])
    ax.set_xlabel("Healer minus non-healer delta")
    ax.set_title("Random gene-panel controls", pad=8)
    ax.legend(frameon=False, loc="lower right")
    clean_axes(ax)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    panel_label(ax, "D")
    savefig(fig, "Figure2_Discovery_SRRS_and_controls")


def figure3():
    align = read_tsv(RESULTS / "d7_alignment_robustness" / "GSE165816_alternative_alignment_sample_scores.tsv")
    tests = read_tsv(RESULTS / "d7_alignment_robustness" / "GSE165816_alternative_alignment_healer_nonhealer_tests.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.7), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.08, h_pad=0.08, wspace=0.20)

    ax = axes[0]
    boxplot_points(ax, align, "D7_vs_D0_top100_original", "Alignment score", "D7 acute wound alignment\nDelta 0.251; exact P = 0.0225")
    panel_label(ax, "A")

    ax = axes[1]
    sub = tests.copy()
    label_map = {
        "D1_vs_D0_top100": "D1 vs skin",
        "D7_vs_D0_top100_original": "D7 vs skin",
        "D30_vs_D0_top100": "D30 vs skin",
        "D7_vs_D1_top100": "D7 vs D1",
        "D7_vs_D30_top100": "D7 vs D30",
        "D7D30_vs_D0_top100": "D7/D30 vs skin",
        "D7_consensus_vs_D0_D1_D30_top100": "D7 consensus",
    }
    sub["label"] = sub["feature"].map(label_map).fillna(sub["feature"])
    forest(ax, sub, "label", "delta_healer_minus_nonhealer", "bootstrap_ci_low", "bootstrap_ci_high", color=GREEN, xlabel="Healer minus non-healer")
    ax.set_title("Alternative signatures", pad=8)
    panel_label(ax, "B")

    ax = axes[2]
    sub = tests.copy()
    sub["label"] = sub["feature"].map(label_map).fillna(sub["feature"])
    order = sub.sort_values("exact_permutation_p")
    colors = [PURPLE if x == "D7 vs skin" else BLUE for x in order["label"]]
    ax.barh(np.arange(len(order)), -np.log10(order["exact_permutation_p"]), color=colors)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order["label"])
    ax.set_xlabel("-log10 exact permutation P")
    ax.set_title("Permutation support", pad=8)
    clean_axes(ax)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    panel_label(ax, "C")
    savefig(fig, "Figure3_Acute_wound_alignment")


def figure4():
    tests = read_tsv(RESULTS / "external_gse223964_rescue" / "GSE223964_SRRS_high_low_within_sample_ligand_tests.tsv")
    deltas = read_tsv(RESULTS / "external_gse223964_rescue" / "GSE223964_SRRS_high_low_within_sample_ligand_deltas.tsv")
    rand = read_tsv(RESULTS / "external_gse223964_rescue" / "GSE223964_SRRS_within_sample_random_panel_control.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.8), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.08, h_pad=0.08, wspace=0.20)

    ax = axes[0]
    panel = deltas[deltas["feature"].eq("candidate_ligand_panel")].copy()
    panel = panel.sort_values("sample")
    ax.axhline(0, color="#AAAAAA", linestyle="--", linewidth=0.8)
    ax.scatter(np.arange(len(panel)), panel["delta_high_minus_low"], color=PURPLE, s=40, zorder=3)
    ax.plot(np.arange(len(panel)), panel["delta_high_minus_low"], color=PURPLE, linewidth=1.2, alpha=0.75)
    ax.set_xticks(np.arange(len(panel)))
    ax.set_xticklabels(panel["sample"], rotation=35, ha="right")
    ax.set_ylabel("SRRS-high minus low")
    ax.set_title("Candidate ligand panel\n7/7 samples positive; P = 0.0156", pad=8)
    clean_axes(ax)
    panel_label(ax, "A")

    ax = axes[1]
    keep = ["TNC", "IL11", "IL6", "INHBA", "THBS1", "SERPINE1", "CCL20", "candidate_ligand_panel"]
    sub = tests[tests["feature"].isin(keep)].copy()
    sub["label"] = sub["feature"].replace({"candidate_ligand_panel": "Ligand panel"})
    sub = sub.sort_values("mean_delta_high_minus_low")
    colors = [PURPLE if x == "Ligand panel" else BLUE for x in sub["label"]]
    ax.barh(np.arange(len(sub)), sub["mean_delta_high_minus_low"], color=colors)
    ax.set_yticks(np.arange(len(sub)))
    ax.set_yticklabels(sub["label"])
    ax.axvline(0, color="#AAAAAA", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Mean paired delta")
    ax.set_title("Individual ligand support", pad=8)
    clean_axes(ax)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    panel_label(ax, "B")

    ax = axes[2]
    observed = float(rand["observed_candidate_ligand_within_sample_delta"].iloc[0])
    random = float(rand["random_mean_delta"].iloc[0])
    sd = float(rand["random_sd_delta"].iloc[0])
    ax.bar([0, 1], [random, observed], color=["#B4B4B4", ORANGE], width=0.55)
    ax.errorbar([0], [random], yerr=[sd], color="#555555", capsize=4, fmt="none")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Random\npanels", "Observed\npanel"])
    ax.set_ylabel("Mean paired delta")
    ax.set_title("Random-panel control\nempirical P = 0.001", pad=8)
    clean_axes(ax)
    ax.grid(axis="x", visible=False)
    panel_label(ax, "C")
    savefig(fig, "Figure4_External_chronic_wound_conservation")


def figure5():
    sample = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_strengthened_sample_summary.tsv")
    paired = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_paired_donor_bootstrap_tests.tsv")
    corr = read_tsv(RESULTS / "spatial_gse241124_strengthened" / "GSE241124_spatial_SRRS_marker_correlation_donor_bootstrap_summary.tsv")
    maps_path = PROJECT / "figures" / "srrs_framework" / "spatial_gse241124" / "GSE241124_spatial_SRRS_representative_maps.png"

    fig = plt.figure(figsize=(13.5, 8.0), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.08, h_pad=0.08, wspace=0.16, hspace=0.08)
    gs = fig.add_gridspec(2, 3, height_ratios=[0.78, 1.0])

    ax = fig.add_subplot(gs[0, :])
    img = trim_white_image(plt.imread(maps_path), pad=18)
    ax.imshow(img)
    ax.axis("off")
    ax.set_title("Representative spatial SRRS maps", pad=6)
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[1, 0])
    cond_order = ["Skin", "Wound1", "Wound7", "Wound30"]
    xmap = {c: i for i, c in enumerate(cond_order)}
    for donor, g in sample.groupby("Donor"):
        g = g[g["Condition"].isin(cond_order)].sort_values("Condition", key=lambda s: s.map(xmap))
        ax.plot([xmap[c] for c in g["Condition"]], g["spatial_SRRS"], marker="o", linewidth=1.2, alpha=0.85)
    ax.set_xticks(range(len(cond_order)))
    ax.set_xticklabels(["Skin", "D1", "D7", "D30"])
    ax.set_ylabel("Spatial SRRS")
    ax.set_title("Donor-level trajectories", pad=8)
    clean_axes(ax)
    panel_label(ax, "B")

    ax = fig.add_subplot(gs[1, 1])
    sub = paired[paired["comparison"].isin(["Wound1_vs_Skin", "Wound7_vs_Skin", "Wound30_vs_Skin"]) & paired["feature"].eq("spatial_SRRS")].copy()
    sub["label"] = sub["comparison"].replace({"Wound1_vs_Skin": "D1 vs skin", "Wound7_vs_Skin": "D7 vs skin", "Wound30_vs_Skin": "D30 vs skin"})
    forest(ax, sub, "label", "mean_delta", "bootstrap_ci_low", "bootstrap_ci_high", color=GREEN, xlabel="Mean donor delta")
    ax.set_title("Paired donor effects", pad=8)
    panel_label(ax, "C")

    ax = fig.add_subplot(gs[1, 2])
    sub = corr[corr["Condition"].eq("Wound7")].copy()
    sub["label"] = sub["feature"].replace(
        {
            "ecm_remodeling_marker": "ECM",
            "fibroblast_stromal_marker": "Fibroblast",
            "endothelial_marker": "Endothelial",
            "pericyte_smc_marker": "Pericyte/SMC",
        }
    )
    forest(ax, sub, "label", "mean_donor_rho", "bootstrap_ci_low", "bootstrap_ci_high", color=PURPLE, xlabel="Mean donor Spearman rho")
    ax.set_title("Wound day 7 marker correlations", pad=8)
    panel_label(ax, "D")
    savefig(fig, "Figure5_Spatial_SRRS_localisation")


def figure6():
    axes_df = read_tsv(RESULTS / "spatial_lr_axis_gse241124" / "GSE241124_spatial_LR_axis_condition_clustered_models.tsv")
    corr = read_tsv(RESULTS / "spatial_lr_axis_gse241124" / "GSE241124_spatial_LR_axis_correlation_donor_bootstrap_summary.tsv")
    gene = read_tsv(RESULTS / "spatial_receptor_expression_gse241124" / "GSE241124_spatial_LR_gene_condition_clustered_models.tsv")
    pert = read_tsv(RESULTS / "perturbation_enrichr" / "SRRS_Enrichr_focused_perturbation_hits.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.8), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.08, h_pad=0.08, wspace=0.18, hspace=0.14)

    ax = axes[0, 0]
    sub = axes_df[axes_df["term"].str.contains("Wound7")].copy().sort_values("estimate")
    sub["label"] = sub["axis"].str.replace("_", "\n", regex=False)
    colors = [ORANGE if p < 0.1 else "#A8A8A8" for p in sub["padj_bh"]]
    ax.barh(np.arange(len(sub)), sub["estimate"], color=colors)
    ax.set_yticks(np.arange(len(sub)))
    ax.set_yticklabels(sub["label"])
    ax.axvline(0, color="#AAAAAA", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Adjusted Wound7 vs skin coefficient")
    ax.set_title("Spatial LR axis scores", pad=8)
    clean_axes(ax)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    panel_label(ax, "A")

    ax = axes[0, 1]
    sub = corr[corr["Condition"].eq("Wound7")].copy()
    sub["label"] = sub["axis"].str.replace("_", "\n", regex=False)
    forest(ax, sub, "label", "mean_donor_rho", "bootstrap_ci_low", "bootstrap_ci_high", color=BLUE, xlabel="Mean donor Spearman rho")
    ax.set_title("Wound7 LR coexpression", pad=8)
    panel_label(ax, "B")

    ax = axes[1, 0]
    w7 = gene[gene["term"].str.contains("Wound7")].copy()
    keep = ["TNC", "ITGB1", "ITGA5", "ITGAV", "INHBA", "BAMBI", "ENG", "IL11", "IL11RA", "IL6", "IL6R", "IL6ST", "CCL20", "CCR6"]
    w7 = w7[w7["gene"].isin(keep)].sort_values("estimate")
    label = [f"{g} ({r[0].upper()})" for g, r in zip(w7["gene"], w7["role"])]
    colors = np.where(w7["role"].eq("ligand"), ORANGE, BLUE)
    ax.barh(np.arange(len(w7)), w7["estimate"], color=colors)
    ax.set_yticks(np.arange(len(w7)))
    ax.set_yticklabels(label)
    ax.axvline(0, color="#AAAAAA", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Adjusted Wound7 vs skin coefficient")
    ax.set_title("Ligand and receptor genes", pad=8)
    clean_axes(ax)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    panel_label(ax, "C")

    ax = axes[1, 1]
    p = pert.sort_values("adjusted_p_value").head(8).copy()
    p["score"] = -np.log10(p["adjusted_p_value"].clip(lower=1e-300))
    p["label"] = [wrap_label(t.replace("_", " "), 30) for t in p["term"]]
    ax.barh(np.arange(len(p))[::-1], p["score"], color=GREEN)
    ax.set_yticks(np.arange(len(p))[::-1])
    ax.set_yticklabels(p["label"], fontsize=7.5)
    ax.set_xlabel("-log10 adjusted P")
    ax.set_title("Perturbation-informed enrichment", pad=8)
    clean_axes(ax)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    panel_label(ax, "D")
    savefig(fig, "Figure6_Communication_and_perturbation")


def figure7():
    fig, ax = plt.subplots(figsize=(12.0, 6.8))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.04, 0.94, "Clinical model: wound-bed repair readiness", fontsize=15, fontweight="bold", color=INK)
    ax.text(0.04, 0.895, "A research framework for prospective biopsy-based validation.", fontsize=9.5, color="#444444")

    draw_box(
        ax,
        (0.06, 0.55),
        (0.39, 0.28),
        "Repair-ready wound bed",
        ["High fibroblast repair activation", "Matrix-remodelling ligand programme", "TNC-integrin/syndecan strongest", "Stromal-vascular tissue context"],
        GREEN,
        title_size=11,
        body_size=8.7,
    )
    draw_box(
        ax,
        (0.55, 0.55),
        (0.39, 0.28),
        "Stalled wound bed",
        ["Low or incomplete SRRS programme", "Weaker D7 acute wound alignment", "Uncoupled inflammation and repair", "Candidate non-healing state"],
        ORANGE,
        title_size=11,
        body_size=8.7,
    )

    ax.add_patch(Rectangle((0.08, 0.38), 0.34, 0.035, color=GREEN, alpha=0.18))
    ax.text(0.25, 0.397, "fibroblast -> matrix -> vessel", ha="center", va="center", fontsize=8.5, color=INK)
    ax.add_patch(Rectangle((0.58, 0.38), 0.32, 0.035, color=ORANGE, alpha=0.18))
    ax.text(0.74, 0.397, "inflammation without repair coupling", ha="center", va="center", fontsize=8.5, color=INK)

    steps = [
        ("Wound\nbiopsy", BLUE),
        ("Molecular\nprofiling", PURPLE),
        ("SRRS\nscore", GREEN),
        ("Repair-readiness\nstrata", ORANGE),
        ("Prospective\nvalidation", GREY),
    ]
    y = 0.18
    xs = np.linspace(0.06, 0.80, len(steps))
    for i, ((label, color), x) in enumerate(zip(steps, xs)):
        patch = FancyBboxPatch((x, y), 0.13, 0.105, boxstyle="round,pad=0.012,rounding_size=0.018", linewidth=1.1, edgecolor=color, facecolor=LIGHT)
        ax.add_patch(patch)
        ax.text(x + 0.065, y + 0.074, label, ha="center", va="top", fontsize=8.2, color=INK)
        if i < len(steps) - 1:
            draw_arrow(ax, (x + 0.13, y + 0.052), (xs[i + 1], y + 0.052))

    ax.add_patch(Rectangle((0.06, 0.055), 0.88, 0.07, facecolor=LIGHT, edgecolor="#D5D8DE", linewidth=1))
    ax.text(0.08, 0.092, "Current status:", fontsize=9, fontweight="bold", va="center", color=INK)
    ax.text(
        0.205,
        0.092,
        "computational prioritisation across public datasets; not clinical deployment.",
        fontsize=8.7,
        va="center",
        color=INK,
    )
    savefig(fig, "Figure7_Clinical_repair_readiness_model")


def main():
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    figure6()
    figure7()
    print(f"Clean CCS figures written to {FIG_DIR}")


if __name__ == "__main__":
    main()
