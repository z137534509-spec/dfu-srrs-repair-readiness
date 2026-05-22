from pathlib import Path
import itertools
import math

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


PROJECT = Path.cwd()
IN_FILE = PROJECT / "results" / "srrs_framework" / "GSE165816_SRRS_discovery_sample_scores.tsv"
OUT_DIR = PROJECT / "results" / "srrs_framework" / "srrs_component_sensitivity"
FIG_DIR = PROJECT / "figures" / "srrs_framework" / "srrs_component_sensitivity"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


Z = {
    "fibroblast": "z_fibroblast_repair_activation",
    "ligand": "z_stromal_ligand_panel",
    "receiver": "z_vascular_perivascular_receiver_coupling",
    "d7": "z_d7_acute_wound_alignment",
}


def exact_perm_p(values, labels):
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    ok = np.isfinite(values)
    values = values[ok]
    labels = labels[ok]
    n_h = int(np.sum(labels == "Healer"))
    obs = values[labels == "Healer"].mean() - values[labels == "Non-healer"].mean()
    deltas = []
    for idx in itertools.combinations(range(len(values)), n_h):
        idx = np.asarray(idx)
        mask = np.zeros(len(values), dtype=bool)
        mask[idx] = True
        deltas.append(values[mask].mean() - values[~mask].mean())
    deltas = np.asarray(deltas)
    return float((np.sum(np.abs(deltas) >= abs(obs)) + 1) / (len(deltas) + 1))


def bootstrap_ci(values, labels, n_boot=10000, seed=1):
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    h = values[(labels == "Healer") & np.isfinite(values)]
    n = values[(labels == "Non-healer") & np.isfinite(values)]
    if len(h) < 2 or len(n) < 2:
        return math.nan, math.nan
    boot = np.array([rng.choice(h, size=len(h), replace=True).mean() - rng.choice(n, size=len(n), replace=True).mean() for _ in range(n_boot)])
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def test_feature(df, feature):
    values = df[feature].to_numpy(dtype=float)
    labels = df["healing_status"].to_numpy()
    h = values[labels == "Healer"]
    n = values[labels == "Non-healer"]
    delta = float(np.nanmean(h) - np.nanmean(n))
    p = stats.mannwhitneyu(h, n, alternative="two-sided").pvalue
    lo, hi = bootstrap_ci(values, labels)
    return {
        "feature": feature,
        "healer_n": int(np.sum(labels == "Healer")),
        "nonhealer_n": int(np.sum(labels == "Non-healer")),
        "healer_mean": float(np.nanmean(h)),
        "nonhealer_mean": float(np.nanmean(n)),
        "delta_healer_minus_nonhealer": delta,
        "p_mannwhitney": float(p),
        "exact_permutation_p": exact_perm_p(values, labels),
        "bootstrap_ci_low": lo,
        "bootstrap_ci_high": hi,
    }


def lodo(df, feature):
    full_delta = test_feature(df, feature)["delta_healer_minus_nonhealer"]
    rows = []
    for i, row in df.iterrows():
        sub = df.drop(index=i)
        delta = test_feature(sub, feature)["delta_healer_minus_nonhealer"]
        rows.append({
            "feature": feature,
            "dropped_sample": row["sample_code"],
            "delta_healer_minus_nonhealer": delta,
            "same_direction_as_full": np.sign(delta) == np.sign(full_delta),
        })
    return pd.DataFrame(rows)


def main():
    df = pd.read_csv(IN_FILE, sep="\t")
    df["SRRS_4_equal"] = df["SRRS"]
    df["SRRS_3_no_receiver"] = df[[Z["fibroblast"], Z["ligand"], Z["d7"]]].mean(axis=1)
    df["SRRS_3_no_d7"] = df[[Z["fibroblast"], Z["ligand"], Z["receiver"]]].mean(axis=1)
    df["SRRS_3_no_ligand"] = df[[Z["fibroblast"], Z["receiver"], Z["d7"]]].mean(axis=1)
    df["SRRS_2_fibroblast_ligand"] = df[[Z["fibroblast"], Z["ligand"]]].mean(axis=1)
    df["SRRS_receiver_10pct"] = (
        0.30 * df[Z["fibroblast"]] + 0.30 * df[Z["ligand"]] + 0.30 * df[Z["d7"]] + 0.10 * df[Z["receiver"]]
    )
    df["SRRS_receiver_20pct"] = (
        (0.80 / 3.0) * df[Z["fibroblast"]] + (0.80 / 3.0) * df[Z["ligand"]] +
        (0.80 / 3.0) * df[Z["d7"]] + 0.20 * df[Z["receiver"]]
    )
    variants = [
        "SRRS_4_equal",
        "SRRS_3_no_receiver",
        "SRRS_receiver_10pct",
        "SRRS_receiver_20pct",
        "SRRS_2_fibroblast_ligand",
        "SRRS_3_no_d7",
        "SRRS_3_no_ligand",
    ]
    tests = pd.DataFrame([test_feature(df, v) for v in variants])
    tests["BH_mannwhitney"] = multipletests_bh(tests["p_mannwhitney"].to_numpy())
    lod = pd.concat([lodo(df, v) for v in variants], ignore_index=True)
    lod_summary = (
        lod.groupby("feature")
        .agg(
            lodo_min_delta=("delta_healer_minus_nonhealer", "min"),
            lodo_max_delta=("delta_healer_minus_nonhealer", "max"),
            lodo_same_direction=("same_direction_as_full", lambda x: f"{int(x.sum())}/{len(x)}"),
        )
        .reset_index()
    )
    tests = tests.merge(lod_summary, on="feature", how="left")

    corr = df[variants].corr(method="spearman").reset_index().rename(columns={"index": "feature"})
    df[["sample_code", "healing_status"] + variants].to_csv(OUT_DIR / "GSE165816_SRRS_variant_sample_scores.tsv", sep="\t", index=False)
    tests.to_csv(OUT_DIR / "GSE165816_SRRS_variant_healer_nonhealer_tests.tsv", sep="\t", index=False)
    lod.to_csv(OUT_DIR / "GSE165816_SRRS_variant_leave_one_sample_out.tsv", sep="\t", index=False)
    corr.to_csv(OUT_DIR / "GSE165816_SRRS_variant_spearman_correlations.tsv", sep="\t", index=False)

    plot = tests.sort_values("delta_healer_minus_nonhealer")
    fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
    y = np.arange(len(plot))
    ax.errorbar(
        plot["delta_healer_minus_nonhealer"],
        y,
        xerr=[
            plot["delta_healer_minus_nonhealer"] - plot["bootstrap_ci_low"],
            plot["bootstrap_ci_high"] - plot["delta_healer_minus_nonhealer"],
        ],
        fmt="o",
        color="#0072B2",
        ecolor="#777777",
        capsize=2,
    )
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["feature"], fontsize=8)
    ax.set_xlabel("Healer minus non-healer delta, bootstrap 95% CI")
    fig.savefig(FIG_DIR / "GSE165816_SRRS_component_sensitivity_effects.png", dpi=240)
    plt.close(fig)
    print("Wrote SRRS component sensitivity outputs")


def multipletests_bh(pvalues):
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
    tmp = np.empty_like(adj)
    tmp[order] = np.clip(adj, 0, 1)
    out[ok] = tmp
    return out


if __name__ == "__main__":
    main()
