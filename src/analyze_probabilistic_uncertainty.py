#!/usr/bin/env python3
"""Evaluate probabilistic intensity predictions with empirical and site uncertainty.

The deterministic thinning experiment gives one predicted instrumental
intensity per held-out station.  This script turns those predictions into
calibrated predictive distributions and separates the empirical residual
spread from uncertainty in the site-amplification-to-intensity conversion.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


DEFAULT_PREDICTIONS = Path(
    "outputs/csv/station_thinning_interpolation_6lower_plus_class/thinning_prediction_errors.csv.gz"
)
DEFAULT_CSV_DIR = Path("outputs/csv/probabilistic_uncertainty_6lower_plus")
DEFAULT_PNG_DIR = Path("outputs/png/probabilistic_uncertainty_6lower_plus")

INTENSITY_CLASS_LABELS = np.array(["0", "1", "2", "3", "4", "5-", "5+", "6-", "6+", "7"], dtype=object)
CLASS_LOWER = np.array([-np.inf, 0.5, 1.5, 2.5, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5], dtype=float)
CLASS_UPPER = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5, np.inf], dtype=float)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Add probabilistic uncertainty diagnostics to thinning validation.")
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--amp-intensity-coef", type=float, default=1.72)
    parser.add_argument(
        "--sigma-amp-intensity-coef",
        type=float,
        default=0.20,
        help="Sensitivity spread of the intensity-amplification coefficient c in c*log10(ARV).",
    )
    parser.add_argument(
        "--sigma-log10-amplification",
        type=float,
        default=0.08,
        help="Assumed one-sigma uncertainty of log10(ARV) for the site amplification grid.",
    )
    parser.add_argument("--min-sigma", type=float, default=0.05)
    parser.add_argument("--no-png", action="store_true")
    return parser.parse_args()


def density_bins(series: pd.Series) -> pd.Categorical:
    bins = [0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000, np.inf]
    return pd.cut(series, bins=bins, right=False)


def robust_sigma(values: pd.Series, min_sigma: float) -> float:
    arr = values.to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 8:
        return float("nan")
    med = np.nanmedian(arr)
    abs_dev = np.abs(arr - med)
    sigma_mad = 1.4826 * np.nanmedian(abs_dev)
    sigma_q68 = np.nanquantile(abs_dev, 0.682689492)
    sigma_rmse = math.sqrt(float(np.nanmean((arr - med) ** 2)))
    return float(max(min_sigma, sigma_mad, sigma_q68, min(sigma_rmse, 2.0 * sigma_q68 + min_sigma)))


def finite_median(values: pd.Series) -> float:
    arr = values.to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan")
    return float(np.median(arr))


def read_predictions(path: Path) -> pd.DataFrame:
    cols = [
        "event_id",
        "method",
        "keep_fraction",
        "random_index",
        "observed_intensity",
        "predicted_intensity",
        "prediction_error",
        "abs_error",
        "class_match",
        "observed_intensity_class_index",
        "observed_intensity_class",
        "effective_train_density_per_10000km2",
        "effective_train_median_nn_distance_km",
        "amplification_vs400",
        "local_site_delta_std_20km",
    ]
    df = pd.read_csv(path, usecols=cols)
    df["density_bin"] = density_bins(df["effective_train_density_per_10000km2"]).astype(str)
    df["observed_intensity_group"] = np.select(
        [
            df["observed_intensity_class_index"] >= 7,
            df["observed_intensity_class_index"] >= 5,
        ],
        ["I>=6-", "I>=5-"],
        default="all",
    )
    return df


def site_conversion_sigma(
    amplification: pd.Series,
    amp_coef: float,
    sigma_amp_coef: float,
    sigma_log10_amp: float,
) -> np.ndarray:
    log10_amp = np.log10(np.clip(amplification.to_numpy(dtype=float), 0.05, None))
    return np.sqrt((log10_amp * sigma_amp_coef) ** 2 + (amp_coef * sigma_log10_amp) ** 2)


def empirical_uncertainty_table(df: pd.DataFrame, min_sigma: float) -> pd.DataFrame:
    rows = []
    for (method, density_bin), group in df.groupby(["method", "density_bin"], observed=True):
        centered = group["prediction_error"] - group["prediction_error"].median()
        rows.append(
            {
                "method": method,
                "density_bin": density_bin,
                "n_predictions": len(group),
                "median_density_per_10000km2": group["effective_train_density_per_10000km2"].median(),
                "median_train_nn_km": group["effective_train_median_nn_distance_km"].median(),
                "bias_median": group["prediction_error"].median(),
                "bias_mean": group["prediction_error"].mean(),
                "sigma_empirical": robust_sigma(centered, min_sigma),
                "rmse": math.sqrt(float(np.mean(group["prediction_error"].to_numpy(dtype=float) ** 2))),
                "mae": group["abs_error"].mean(),
                "deterministic_class_accuracy": group["class_match"].mean(),
            }
        )
    return pd.DataFrame(rows)


def attach_uncertainty(df: pd.DataFrame, unc: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    keyed = unc[["method", "density_bin", "bias_median", "sigma_empirical"]].copy()
    out = df.merge(keyed, on=["method", "density_bin"], how="left")
    fallback_sigma = robust_sigma(out["prediction_error"] - out["prediction_error"].median(), args.min_sigma)
    out["sigma_empirical"] = out["sigma_empirical"].fillna(fallback_sigma).clip(lower=args.min_sigma)
    out["bias_median"] = out["bias_median"].fillna(0.0)
    out["mu_bias_corrected"] = out["predicted_intensity"] - out["bias_median"]
    out["sigma_site_conversion"] = site_conversion_sigma(
        out["amplification_vs400"],
        args.amp_intensity_coef,
        args.sigma_amp_intensity_coef,
        args.sigma_log10_amplification,
    )
    # The empirical residual sigma is the main predictive uncertainty.  The
    # conservative variant adds site-conversion uncertainty to show sensitivity.
    out["sigma_predictive"] = out["sigma_empirical"]
    out["sigma_predictive_with_site"] = np.sqrt(out["sigma_empirical"] ** 2 + out["sigma_site_conversion"] ** 2)
    return out


def true_class_probability(mu: np.ndarray, sigma: np.ndarray, class_index: np.ndarray) -> np.ndarray:
    idx = np.asarray(class_index, dtype=int)
    lower = np.full(len(idx), np.nan, dtype=float)
    upper = np.full(len(idx), np.nan, dtype=float)
    valid = (0 <= idx) & (idx < len(INTENSITY_CLASS_LABELS))
    lower[valid] = CLASS_LOWER[idx[valid]]
    upper[valid] = CLASS_UPPER[idx[valid]]
    z_hi = (upper - mu) / sigma
    z_lo = (lower - mu) / sigma
    return np.clip(norm.cdf(z_hi) - norm.cdf(z_lo), 0.0, 1.0)


def interval_coverage(error: pd.Series, sigma: pd.Series, level: float) -> float:
    z = float(norm.ppf((1.0 + level) / 2.0))
    return float(np.mean(np.abs(error.to_numpy(dtype=float)) <= z * sigma.to_numpy(dtype=float)))


def aggregate_probabilistic(out: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    corrected_error = out["mu_bias_corrected"] - out["observed_intensity"]
    out["corrected_error"] = corrected_error
    for suffix, sigma_col in [("", "sigma_predictive"), ("_with_site", "sigma_predictive_with_site")]:
        out[f"prob_true_class{suffix}"] = true_class_probability(
            out["mu_bias_corrected"].to_numpy(dtype=float),
            out[sigma_col].to_numpy(dtype=float),
            out["observed_intensity_class_index"].to_numpy(dtype=int),
        )
        out[f"log_score_class{suffix}"] = -np.log(np.clip(out[f"prob_true_class{suffix}"], 1e-12, 1.0))

    station_summary_rows = []
    group_cols = ["method", "density_bin", "observed_intensity_group"]
    grouped_source = pd.concat(
        [
            out.assign(observed_intensity_group="all"),
            out[out["observed_intensity_class_index"] >= 5].assign(observed_intensity_group="I>=5-"),
            out[out["observed_intensity_class_index"] >= 7].assign(observed_intensity_group="I>=6-"),
        ],
        ignore_index=True,
    )
    for keys, group in grouped_source.groupby(group_cols, observed=True):
        method, density_bin, intensity_group = keys
        row = {
            "method": method,
            "density_bin": density_bin,
            "observed_intensity_group": intensity_group,
            "n_predictions": len(group),
            "median_density_per_10000km2": group["effective_train_density_per_10000km2"].median(),
            "median_train_nn_km": finite_median(group["effective_train_median_nn_distance_km"]),
            "deterministic_class_accuracy": group["class_match"].mean(),
            "mae_bias_corrected": np.mean(np.abs(group["corrected_error"])),
            "rmse_bias_corrected": math.sqrt(float(np.mean(group["corrected_error"] ** 2))),
            "sigma_empirical_median": group["sigma_empirical"].median(),
            "sigma_site_conversion_median": group["sigma_site_conversion"].median(),
            "mean_prob_true_class": group["prob_true_class"].mean(),
            "median_prob_true_class": group["prob_true_class"].median(),
            "mean_prob_true_class_with_site": group["prob_true_class_with_site"].mean(),
            "median_prob_true_class_with_site": group["prob_true_class_with_site"].median(),
            "mean_log_score_class": group["log_score_class"].mean(),
            "coverage_50": interval_coverage(group["corrected_error"], group["sigma_predictive"], 0.50),
            "coverage_80": interval_coverage(group["corrected_error"], group["sigma_predictive"], 0.80),
            "coverage_90": interval_coverage(group["corrected_error"], group["sigma_predictive"], 0.90),
            "coverage_80_with_site": interval_coverage(
                group["corrected_error"], group["sigma_predictive_with_site"], 0.80
            ),
            "coverage_90_with_site": interval_coverage(
                group["corrected_error"], group["sigma_predictive_with_site"], 0.90
            ),
        }
        station_summary_rows.append(row)
    station_summary = pd.DataFrame(station_summary_rows)

    trial_keys = ["event_id", "method", "keep_fraction", "random_index", "density_bin", "observed_intensity_group"]
    trial_rows = []
    for keys, group in grouped_source.groupby(trial_keys, observed=True):
        event_id, method, keep_fraction, random_index, density_bin, intensity_group = keys
        trial_rows.append(
            {
                "event_id": event_id,
                "method": method,
                "keep_fraction": keep_fraction,
                "random_index": random_index,
                "density_bin": density_bin,
                "observed_intensity_group": intensity_group,
                "n_predictions": len(group),
                "median_density_per_10000km2": group["effective_train_density_per_10000km2"].median(),
                "median_train_nn_km": finite_median(group["effective_train_median_nn_distance_km"]),
                "deterministic_class_accuracy": group["class_match"].mean(),
                "mae_bias_corrected": np.mean(np.abs(group["corrected_error"])),
                "rmse_bias_corrected": math.sqrt(float(np.mean(group["corrected_error"] ** 2))),
                "sigma_empirical_median": group["sigma_empirical"].median(),
                "sigma_site_conversion_median": group["sigma_site_conversion"].median(),
                "mean_prob_true_class": group["prob_true_class"].mean(),
                "median_prob_true_class": group["prob_true_class"].median(),
                "mean_prob_true_class_with_site": group["prob_true_class_with_site"].mean(),
                "median_prob_true_class_with_site": group["prob_true_class_with_site"].median(),
                "mean_log_score_class": group["log_score_class"].mean(),
                "coverage_50": interval_coverage(group["corrected_error"], group["sigma_predictive"], 0.50),
                "coverage_80": interval_coverage(group["corrected_error"], group["sigma_predictive"], 0.80),
                "coverage_90": interval_coverage(group["corrected_error"], group["sigma_predictive"], 0.90),
                "coverage_80_with_site": interval_coverage(
                    group["corrected_error"], group["sigma_predictive_with_site"], 0.80
                ),
                "coverage_90_with_site": interval_coverage(
                    group["corrected_error"], group["sigma_predictive_with_site"], 0.90
                ),
            }
        )
    trial_metrics = pd.DataFrame(trial_rows)
    prob_summary = (
        trial_metrics.groupby(["method", "density_bin", "observed_intensity_group"], observed=True)
        .agg(
            n_trials=("event_id", "size"),
            n_predictions_median=("n_predictions", "median"),
            median_density_per_10000km2=("median_density_per_10000km2", "median"),
            median_train_nn_km=("median_train_nn_km", "median"),
            deterministic_class_accuracy=("deterministic_class_accuracy", "median"),
            mae_bias_corrected=("mae_bias_corrected", "median"),
            rmse_bias_corrected=("rmse_bias_corrected", "median"),
            sigma_empirical_median=("sigma_empirical_median", "median"),
            sigma_site_conversion_median=("sigma_site_conversion_median", "median"),
            mean_prob_true_class=("mean_prob_true_class", "median"),
            median_prob_true_class=("median_prob_true_class", "median"),
            mean_prob_true_class_with_site=("mean_prob_true_class_with_site", "median"),
            median_prob_true_class_with_site=("median_prob_true_class_with_site", "median"),
            mean_log_score_class=("mean_log_score_class", "median"),
            coverage_50=("coverage_50", "median"),
            coverage_80=("coverage_80", "median"),
            coverage_90=("coverage_90", "median"),
            coverage_80_with_site=("coverage_80_with_site", "median"),
            coverage_90_with_site=("coverage_90_with_site", "median"),
        )
        .reset_index()
    )

    class_summary = (
        out.groupby(["method", "density_bin", "observed_intensity_class"], observed=True)
        .agg(
            n_predictions=("class_match", "size"),
            median_density_per_10000km2=("effective_train_density_per_10000km2", "median"),
            deterministic_class_accuracy=("class_match", "mean"),
            mae_bias_corrected=("corrected_error", lambda x: float(np.mean(np.abs(x)))),
            mean_prob_true_class=("prob_true_class", "mean"),
            mean_prob_true_class_with_site=("prob_true_class_with_site", "mean"),
            sigma_empirical_median=("sigma_empirical", "median"),
            sigma_site_conversion_median=("sigma_site_conversion", "median"),
        )
        .reset_index()
    )

    component_summary = (
        out.groupby(["method", "density_bin"], observed=True)
        .agg(
            n_predictions=("class_match", "size"),
            median_density_per_10000km2=("effective_train_density_per_10000km2", "median"),
            sigma_empirical_median=("sigma_empirical", "median"),
            sigma_site_conversion_median=("sigma_site_conversion", "median"),
            local_site_delta_std_20km_median=("local_site_delta_std_20km", "median"),
            residual_sigma_after_site_median=(
                "sigma_empirical",
                lambda x: float("nan"),
            ),
        )
        .reset_index()
    )
    component_summary["residual_sigma_after_site_median"] = np.sqrt(
        np.maximum(
            component_summary["sigma_empirical_median"] ** 2
            - component_summary["sigma_site_conversion_median"] ** 2,
            0.0,
        )
    )
    component_summary["site_variance_fraction"] = (
        component_summary["sigma_site_conversion_median"] ** 2
        / component_summary["sigma_empirical_median"].clip(lower=1e-9) ** 2
    ).clip(upper=1.0)
    return prob_summary, station_summary, class_summary, component_summary


def sufficient_density(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = [
        ("deterministic_class_accuracy", [0.70, 0.80, 0.90]),
        ("mean_prob_true_class", [0.40, 0.50, 0.60]),
        ("coverage_80", [0.75, 0.80]),
        ("coverage_90", [0.85, 0.90]),
    ]
    for (method, intensity_group), group in summary.groupby(["method", "observed_intensity_group"], observed=True):
        group = group.sort_values("median_density_per_10000km2")
        for metric, targets in metrics:
            for target in targets:
                ok = group[group[metric] >= target]
                if ok.empty:
                    rows.append(
                        {
                            "method": method,
                            "observed_intensity_group": intensity_group,
                            "metric": metric,
                            "target": target,
                            "sufficient_density_per_10000km2": np.nan,
                            "sufficient_median_train_nn_km": np.nan,
                            "achieved_value": np.nan,
                            "status": "not_reached",
                        }
                    )
                else:
                    first = ok.iloc[0]
                    rows.append(
                        {
                            "method": method,
                            "observed_intensity_group": intensity_group,
                            "metric": metric,
                            "target": target,
                            "sufficient_density_per_10000km2": first["median_density_per_10000km2"],
                            "sufficient_median_train_nn_km": first["median_train_nn_km"],
                            "achieved_value": first[metric],
                            "status": "reached",
                        }
                    )
    return pd.DataFrame(rows)


def setup_plot_style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 320,
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#222222",
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "grid.color": "#d8d8d8",
            "grid.linewidth": 0.55,
            "grid.alpha": 0.85,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
        }
    )


def savefig(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(fig)


def format_density_axis(ax, values: pd.Series) -> None:
    max_density = float(np.nanmax(values)) if len(values) else 1.0
    upper = max(1.0, max_density * 1.15)
    ticks = [0, 1, 3, 10, 30, 100, 300, 1000]
    ticks = [tick for tick in ticks if tick <= upper]
    ax.set_xscale("symlog", linthresh=1.0, linscale=0.75)
    ax.set_xlim(0, upper)
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0" if tick == 0 else f"{tick:g}" for tick in ticks])


def make_plots(prob_summary: pd.DataFrame, component_summary: pd.DataFrame, png_dir: Path) -> None:
    import matplotlib.pyplot as plt

    setup_plot_style()
    png_dir.mkdir(parents=True, exist_ok=True)
    colors = {
        "gmpe_raw": "#666666",
        "gmpe_calibrated": "#A6761D",
        "idw": "#D55E00",
        "kriging": "#CC79A7",
        "gmpe_kriging": "#7A3B00",
    }

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), sharey=True)
    for ax, intensity_group in zip(axes, ["all", "I>=5-"]):
        sub = prob_summary[prob_summary["observed_intensity_group"] == intensity_group]
        for method, group in sub.groupby("method", observed=True):
            group = group.sort_values("median_density_per_10000km2")
            ax.plot(
                group["median_density_per_10000km2"],
                group["mean_prob_true_class"],
                marker="o",
                linewidth=1.7,
                color=colors.get(method),
                label=method,
            )
        ax.set_title(f"Observed class group: {intensity_group}")
        ax.set_xlabel("Effective station density (/10,000 km$^2$)")
        format_density_axis(ax, sub["median_density_per_10000km2"])
    axes[0].set_ylabel("Mean probability assigned to true class")
    axes[1].legend(loc="lower right")
    savefig(fig, png_dir / "prob_true_class_vs_density.png")

    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    high = prob_summary[prob_summary["observed_intensity_group"] == "I>=5-"]
    for method, group in high.groupby("method", observed=True):
        group = group.sort_values("median_density_per_10000km2")
        ax.plot(
            group["median_density_per_10000km2"],
            group["deterministic_class_accuracy"],
            marker="o",
            linewidth=1.7,
            color=colors.get(method),
            label=f"{method}: exact class",
        )
        ax.plot(
            group["median_density_per_10000km2"],
            group["mean_prob_true_class"],
            marker="s",
            linewidth=1.2,
            linestyle="--",
            color=colors.get(method),
            alpha=0.75,
        )
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Effective station density (/10,000 km$^2$)")
    ax.set_ylabel("Accuracy / probability")
    ax.set_title("High-intensity stations: deterministic and probabilistic scores")
    format_density_axis(ax, high["median_density_per_10000km2"])
    ax.legend(ncol=2)
    savefig(fig, png_dir / "high_intensity_probabilistic_scores.png")

    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    all_group = prob_summary[prob_summary["observed_intensity_group"] == "all"]
    for method, group in all_group.groupby("method", observed=True):
        group = group.sort_values("median_density_per_10000km2")
        ax.plot(
            group["median_density_per_10000km2"],
            group["coverage_80"],
            marker="o",
            linewidth=1.7,
            color=colors.get(method),
            label=method,
        )
    ax.axhline(0.80, color="#222222", linestyle="--", linewidth=1.0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Effective station density (/10,000 km$^2$)")
    ax.set_ylabel("Observed coverage of 80% predictive interval")
    ax.set_title("Predictive interval calibration")
    format_density_axis(ax, all_group["median_density_per_10000km2"])
    ax.legend(ncol=2)
    savefig(fig, png_dir / "predictive_interval_coverage80.png")

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), sharey=True)
    for ax, method in zip(axes, ["idw", "gmpe_kriging"]):
        sub = component_summary[component_summary["method"] == method].sort_values("median_density_per_10000km2")
        ax.plot(
            sub["median_density_per_10000km2"],
            sub["sigma_empirical_median"],
            marker="o",
            linewidth=1.8,
            color="#222222",
            label="Empirical predictive sigma",
        )
        ax.plot(
            sub["median_density_per_10000km2"],
            sub["sigma_site_conversion_median"],
            marker="s",
            linewidth=1.6,
            color="#0072B2",
            label="Site-conversion sigma",
        )
        ax.plot(
            sub["median_density_per_10000km2"],
            sub["residual_sigma_after_site_median"],
            marker="^",
            linewidth=1.6,
            color="#D55E00",
            label="Residual after site term",
        )
        ax.set_title(method)
        ax.set_xlabel("Effective station density (/10,000 km$^2$)")
        format_density_axis(ax, sub["median_density_per_10000km2"])
    axes[0].set_ylabel("Intensity sigma")
    axes[1].legend(loc="upper right")
    savefig(fig, png_dir / "uncertainty_components_vs_density.png")


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)

    df = read_predictions(args.predictions)
    unc = empirical_uncertainty_table(df, args.min_sigma)
    out = attach_uncertainty(df, unc, args)
    prob_summary, station_summary, class_summary, component_summary = aggregate_probabilistic(out)
    thresholds = sufficient_density(prob_summary)

    unc.to_csv(args.csv_dir / "empirical_uncertainty_by_method_density.csv", index=False)
    prob_summary.to_csv(args.csv_dir / "probabilistic_validation_summary.csv", index=False)
    station_summary.to_csv(args.csv_dir / "probabilistic_validation_station_weighted_summary.csv", index=False)
    class_summary.to_csv(args.csv_dir / "probabilistic_validation_by_class.csv", index=False)
    component_summary.to_csv(args.csv_dir / "uncertainty_components_by_method_density.csv", index=False)
    thresholds.to_csv(args.csv_dir / "probabilistic_sufficient_density_thresholds.csv", index=False)

    if not args.no_png:
        make_plots(prob_summary, component_summary, args.png_dir)

    print(f"Saved CSV directory: {args.csv_dir}")
    if not args.no_png:
        print(f"Saved PNG directory: {args.png_dir}")


if __name__ == "__main__":
    main()
