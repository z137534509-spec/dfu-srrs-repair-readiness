from pathlib import Path
import itertools
import math

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt


PROJECT = Path.cwd()
SPOT_FILE = PROJECT / "results" / "srrs_framework" / "spatial_gse241124" / "GSE241124_spatial_SRRS_spot_scores.tsv"
OUT_DIR = PROJECT / "results" / "srrs_framework" / "spatial_gse241124_strengthened"
FIG_DIR = PROJECT / "figures" / "srrs_framework" / "spatial_gse241124_strengthened"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

CONDITIONS = ["Skin", "Wound1", "Wound7", "Wound30"]
MARKERS = ["fibroblast_stromal_marker", "endothelial_marker", "pericyte_smc_marker", "ecm_remodeling_marker"]


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
    adj = np.clip(adj, 0, 1)
    tmp = np.empty_like(adj)
    tmp[order] = adj
    out[ok] = tmp
    return out


def paired_bootstrap(wide, condition, feature, n_boot=10000, seed=1):
    dat = wide[["Skin", condition]].dropna()
    diffs = (dat[condition] - dat["Skin"]).to_numpy(dtype=float)
    obs = float(np.mean(diffs)) if len(diffs) else math.nan
    rng = np.random.default_rng(seed)
    boots = np.array([np.mean(rng.choice(diffs, size=len(diffs), replace=True)) for _ in range(n_boot)]) if len(diffs) else np.array([])
    ci = np.quantile(boots, [0.025, 0.975]) if len(boots) else [math.nan, math.nan]
    # Exact sign-flip permutation for paired donor differences.
    if len(diffs):
        perms = []
        for signs in itertools.product([-1, 1], repeat=len(diffs)):
            perms.append(np.mean(diffs * np.asarray(signs)))
        perms = np.asarray(perms)
        p_perm = (np.sum(np.abs(perms) >= abs(obs)) + 1) / (len(perms) + 1)
        sign_same = int(np.sum(diffs > 0))
    else:
        p_perm = math.nan
        sign_same = 0
    try:
        p_wil = stats.wilcoxon(diffs).pvalue if len(diffs) >= 3 else math.nan
    except Exception:
        p_wil = math.nan
    return {
        "feature": feature,
        "comparison": f"{condition}_vs_Skin",
        "n_donors": len(diffs),
        "mean_delta": obs,
        "bootstrap_ci_low": float(ci[0]),
        "bootstrap_ci_high": float(ci[1]),
        "positive_donors": sign_same,
        "p_exact_signflip": p_perm,
        "p_paired_wilcoxon": p_wil,
    }


def donor_level_correlations(spots, feature):
    rows = []
    for (condition, donor), sub in spots.groupby(["Condition", "Donor"], observed=True):
        dat = sub[["spatial_SRRS", feature]].dropna()
        if len(dat) < 20:
            continue
        rho, p = stats.spearmanr(dat["spatial_SRRS"], dat[feature])
        rows.append({"Condition": condition, "Donor": donor, "feature": feature, "rho": rho, "p": p, "n_spots": len(dat)})
    return pd.DataFrame(rows)


def bootstrap_donor_rho(donor_rhos, seed=1, n_boot=10000):
    vals = donor_rhos[np.isfinite(donor_rhos)]
    if len(vals) == 0:
        return math.nan, math.nan, math.nan, 0, 0
    rng = np.random.default_rng(seed)
    boots = np.array([np.mean(rng.choice(vals, size=len(vals), replace=True)) for _ in range(n_boot)])
    return float(np.mean(vals)), float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975)), int(np.sum(vals > 0)), len(vals)


def main():
    spots = pd.read_csv(SPOT_FILE, sep="\t")
    spots = spots[spots["Condition"].isin(CONDITIONS)].copy()
    spots["Condition"] = pd.Categorical(spots["Condition"], categories=CONDITIONS, ordered=True)
    spots["log_nCount"] = np.log1p(spots["nCount_Spatial"])
    spots["log_nFeature"] = np.log1p(spots["nFeature_Spatial"])

    # Residualised SRRS controls generic library/housekeeping activity.
    adj_model = smf.ols(
        "spatial_SRRS ~ negative_housekeeping_control_z + log_nCount + log_nFeature",
        data=spots.dropna(subset=["spatial_SRRS", "negative_housekeeping_control_z", "log_nCount", "log_nFeature"]),
    ).fit()
    spots.loc[adj_model.model.data.row_labels, "spatial_SRRS_resid"] = adj_model.resid
    spots.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_spot_scores_with_residual.tsv", sep="\t", index=False)

    sample_summary = (
        spots.groupby(["Donor", "Condition", "Sample_name"], observed=True)[
            ["spatial_SRRS", "spatial_SRRS_resid", "negative_housekeeping_control_z"] + MARKERS
        ]
        .mean()
        .reset_index()
    )
    sample_summary["spots"] = spots.groupby(["Donor", "Condition", "Sample_name"], observed=True).size().to_numpy()
    sample_summary.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_strengthened_sample_summary.tsv", sep="\t", index=False)

    paired_rows = []
    for feature in ["spatial_SRRS", "spatial_SRRS_resid", "negative_housekeeping_control_z"] + MARKERS:
        wide = sample_summary.pivot_table(index="Donor", columns="Condition", values=feature, aggfunc="mean")
        for condition in ["Wound1", "Wound7", "Wound30"]:
            paired_rows.append(paired_bootstrap(wide, condition, feature))
    paired = pd.DataFrame(paired_rows)
    paired["padj_signflip_bh"] = bh_adjust(paired["p_exact_signflip"])
    paired.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_paired_donor_bootstrap_tests.tsv", sep="\t", index=False)

    # Spot-level models are supporting tissue-localisation analyses; cluster-robust SEs are by spatial sample.
    model_rows = []
    formulas = {
        "unadjusted": "spatial_SRRS ~ C(Condition, Treatment(reference='Skin'))",
        "adjusted": "spatial_SRRS ~ C(Condition, Treatment(reference='Skin')) + negative_housekeeping_control_z + log_nCount + log_nFeature",
        "residual": "spatial_SRRS_resid ~ C(Condition, Treatment(reference='Skin'))",
    }
    for name, formula in formulas.items():
        dat = spots.dropna(subset=["spatial_SRRS", "spatial_SRRS_resid", "negative_housekeeping_control_z", "log_nCount", "log_nFeature", "Sample_name"])
        fit = smf.ols(formula, data=dat).fit(cov_type="cluster", cov_kwds={"groups": dat["Sample_name"]})
        for term in fit.params.index:
            if "Condition" not in term:
                continue
            model_rows.append({
                "model": name,
                "term": term,
                "estimate": fit.params[term],
                "std_error_cluster_sample": fit.bse[term],
                "p_cluster_sample": fit.pvalues[term],
                "n_spots": int(fit.nobs),
                "n_clusters": dat["Sample_name"].nunique(),
            })
    models = pd.DataFrame(model_rows)
    models["padj_bh"] = bh_adjust(models["p_cluster_sample"])
    models.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_spot_level_clustered_models.tsv", sep="\t", index=False)

    corr_all = pd.concat([donor_level_correlations(spots, m) for m in MARKERS], ignore_index=True)
    corr_all.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_marker_correlations_by_donor.tsv", sep="\t", index=False)
    corr_summary_rows = []
    for (condition, feature), sub in corr_all.groupby(["Condition", "feature"], observed=True):
        mean_rho, lo, hi, positive, n = bootstrap_donor_rho(sub["rho"].to_numpy())
        corr_summary_rows.append({
            "Condition": condition,
            "feature": feature,
            "mean_donor_rho": mean_rho,
            "bootstrap_ci_low": lo,
            "bootstrap_ci_high": hi,
            "positive_donors": positive,
            "n_donors": n,
        })
    corr_summary = pd.DataFrame(corr_summary_rows)
    corr_summary.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_marker_correlation_donor_bootstrap_summary.tsv", sep="\t", index=False)

    # Region enrichment for top-quartile SRRS spots within each condition.
    region_rows = []
    for condition, sub in spots.groupby("Condition", observed=True):
        cutoff = np.nanquantile(sub["spatial_SRRS"], 0.75)
        sub = sub.copy()
        sub["SRRS_high_spot"] = sub["spatial_SRRS"] >= cutoff
        total_high = int(sub["SRRS_high_spot"].sum())
        total_low = int((~sub["SRRS_high_spot"]).sum())
        for anno, z in sub.groupby("AnnoType", observed=True):
            high = int(z["SRRS_high_spot"].sum())
            low = int((~z["SRRS_high_spot"]).sum())
            other_high = total_high - high
            other_low = total_low - low
            if high + low < 20:
                continue
            odds, p = stats.fisher_exact([[high, low], [other_high, other_low]], alternative="greater")
            region_rows.append({
                "Condition": condition,
                "AnnoType": anno,
                "spots": int(high + low),
                "high_spots": high,
                "high_fraction": high / (high + low),
                "odds_ratio_high_spot_enrichment": odds,
                "p_fisher_greater": p,
                "mean_SRRS": float(z["spatial_SRRS"].mean()),
            })
    region = pd.DataFrame(region_rows)
    region["padj_bh"] = bh_adjust(region["p_fisher_greater"])
    region = region.sort_values(["Condition", "p_fisher_greater", "mean_SRRS"], ascending=[True, True, False])
    region.to_csv(OUT_DIR / "GSE241124_spatial_SRRS_high_spot_region_enrichment.tsv", sep="\t", index=False)

    # Figures.
    fig, ax = plt.subplots(figsize=(7.2, 4.3), constrained_layout=True)
    plot = paired[(paired["feature"].isin(["spatial_SRRS", "spatial_SRRS_resid", "negative_housekeeping_control_z"])) &
                  (paired["comparison"].isin(["Wound1_vs_Skin", "Wound7_vs_Skin", "Wound30_vs_Skin"]))].copy()
    ylabels = plot["feature"] + " " + plot["comparison"]
    y = np.arange(len(plot))
    ax.errorbar(
        plot["mean_delta"],
        y,
        xerr=[plot["mean_delta"] - plot["bootstrap_ci_low"], plot["bootstrap_ci_high"] - plot["mean_delta"]],
        fmt="o",
        color="#0072B2",
        ecolor="#777777",
        capsize=2,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Paired donor delta vs Skin, bootstrap 95% CI")
    fig.savefig(FIG_DIR / "GSE241124_spatial_SRRS_paired_bootstrap_effects.png", dpi=240)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
    cdat = corr_summary[corr_summary["Condition"].astype(str).isin(["Wound1", "Wound7", "Wound30"])].copy()
    cdat["label"] = cdat["Condition"].astype(str) + " " + cdat["feature"].str.replace("_marker", "", regex=False)
    cdat = cdat.sort_values("mean_donor_rho")
    y = np.arange(len(cdat))
    ax.errorbar(
        cdat["mean_donor_rho"],
        y,
        xerr=[cdat["mean_donor_rho"] - cdat["bootstrap_ci_low"], cdat["bootstrap_ci_high"] - cdat["mean_donor_rho"]],
        fmt="o",
        color="#D55E00",
        ecolor="#777777",
        capsize=2,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(cdat["label"], fontsize=8)
    ax.set_xlabel("Mean donor-level Spearman rho with spatial SRRS")
    fig.savefig(FIG_DIR / "GSE241124_spatial_SRRS_marker_correlation_bootstrap.png", dpi=240)
    plt.close(fig)

    print("Wrote strengthened GSE241124 spatial SRRS outputs")


if __name__ == "__main__":
    main()
