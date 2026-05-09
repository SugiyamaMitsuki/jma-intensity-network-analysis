#!/usr/bin/env python3
"""Diagnose intensity-dependent predictability and high-intensity nonlinearity."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


DEFAULT_PREDICTIONS = Path(
    "outputs/csv/station_thinning_interpolation_6lower_plus_class/thinning_prediction_errors.csv.gz"
)
DEFAULT_CSV_DIR = Path("outputs/csv/intensity_dependent_predictability")
DEFAULT_PNG_DIR = Path("outputs/png/intensity_dependent_predictability")
METHOD_ORDER = ["gmpe_raw", "gmpe_calibrated", "idw", "kriging", "gmpe_kriging"]
CLASS_ORDER = ["1", "2", "3", "4", "5-", "5+", "6-", "6+", "7"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize prediction errors by observed intensity class.")
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--no-png", action="store_true")
    return parser.parse_args()


def density_bins(series: pd.Series) -> pd.Categorical:
    bins = [0, 5, 10, 20, 30, 50, 75, 100, np.inf]
    return pd.cut(series, bins=bins, right=False)


def read_predictions(path: Path) -> pd.DataFrame:
    cols = [
        "method",
        "observed_intensity",
        "observed_intensity_class",
        "observed_intensity_class_index",
        "predicted_intensity",
        "prediction_error",
        "abs_error",
        "class_match",
        "class_error",
        "effective_train_density_per_10000km2",
        "effective_train_median_nn_distance_km",
    ]
    df = pd.read_csv(path, usecols=cols)
    df["density_bin"] = density_bins(df["effective_train_density_per_10000km2"]).astype(str)
    df["class_abs_error"] = np.abs(df["class_error"].to_numpy(dtype=float))
    df["under_class"] = df["class_error"] < 0
    df["over_class"] = df["class_error"] > 0
    df["high_intensity_group"] = np.select(
        [
            df["observed_intensity_class_index"] >= 7,
            df["observed_intensity_class_index"] >= 5,
        ],
        ["I>=6-", "I>=5-"],
        default="I<5-",
    )
    return df


def rmse(values: pd.Series) -> float:
    arr = values.to_numpy(dtype=float)
    return float(math.sqrt(np.mean(arr**2)))


def mad_sigma(values: pd.Series) -> float:
    arr = values.to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan")
    med = np.median(arr)
    return float(1.4826 * np.median(np.abs(arr - med)))


def summarize_by_class(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby(["method", "density_bin", "observed_intensity_class"], observed=True)
        .agg(
            n_predictions=("class_match", "size"),
            median_density_per_10000km2=("effective_train_density_per_10000km2", "median"),
            median_train_nn_km=("effective_train_median_nn_distance_km", "median"),
            observed_intensity_median=("observed_intensity", "median"),
            predicted_intensity_median=("predicted_intensity", "median"),
            bias=("prediction_error", "mean"),
            median_error=("prediction_error", "median"),
            mae=("abs_error", "mean"),
            rmse=("prediction_error", rmse),
            residual_sigma_mad=("prediction_error", mad_sigma),
            class_accuracy=("class_match", "mean"),
            class_mae=("class_abs_error", "mean"),
            under_class_rate=("under_class", "mean"),
            over_class_rate=("over_class", "mean"),
        )
        .reset_index()
    )
    summary["method"] = pd.Categorical(summary["method"], categories=METHOD_ORDER, ordered=True)
    summary["observed_intensity_class"] = pd.Categorical(
        summary["observed_intensity_class"], categories=CLASS_ORDER, ordered=True
    )
    return summary.sort_values(["method", "density_bin", "observed_intensity_class"]).reset_index(drop=True)


def summarize_high_groups(df: pd.DataFrame) -> pd.DataFrame:
    frames = [
        df.assign(intensity_group="all"),
        df[df["observed_intensity_class_index"] < 5].assign(intensity_group="I<5-"),
        df[df["observed_intensity_class_index"] >= 5].assign(intensity_group="I>=5-"),
        df[df["observed_intensity_class_index"] >= 7].assign(intensity_group="I>=6-"),
    ]
    long = pd.concat(frames, ignore_index=True)
    summary = (
        long.groupby(["method", "density_bin", "intensity_group"], observed=True)
        .agg(
            n_predictions=("class_match", "size"),
            median_density_per_10000km2=("effective_train_density_per_10000km2", "median"),
            median_train_nn_km=("effective_train_median_nn_distance_km", "median"),
            bias=("prediction_error", "mean"),
            median_error=("prediction_error", "median"),
            mae=("abs_error", "mean"),
            rmse=("prediction_error", rmse),
            residual_sigma_mad=("prediction_error", mad_sigma),
            class_accuracy=("class_match", "mean"),
            class_mae=("class_abs_error", "mean"),
            under_class_rate=("under_class", "mean"),
            over_class_rate=("over_class", "mean"),
        )
        .reset_index()
    )
    summary["method"] = pd.Categorical(summary["method"], categories=METHOD_ORDER, ordered=True)
    return summary.sort_values(["method", "density_bin", "intensity_group"]).reset_index(drop=True)


def trend_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (method, density_bin), group in df.groupby(["method", "density_bin"], observed=True):
        if len(group) < 50:
            continue
        for metric, values in [
            ("abs_error", group["abs_error"]),
            ("class_abs_error", group["class_abs_error"]),
            ("under_class", group["under_class"].astype(float)),
        ]:
            rho, p_value = stats.spearmanr(group["observed_intensity"], values)
            rows.append(
                {
                    "method": method,
                    "density_bin": density_bin,
                    "metric": metric,
                    "n_predictions": len(group),
                    "spearman_rho_vs_observed_intensity": float(rho),
                    "p_value": float(p_value),
                    "median_density_per_10000km2": float(group["effective_train_density_per_10000km2"].median()),
                }
            )
    out = pd.DataFrame(rows)
    out["method"] = pd.Categorical(out["method"], categories=METHOD_ORDER, ordered=True)
    return out.sort_values(["metric", "method", "median_density_per_10000km2"]).reset_index(drop=True)


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


def make_plots(class_summary: pd.DataFrame, group_summary: pd.DataFrame, png_dir: Path) -> None:
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

    high_density = class_summary[class_summary["density_bin"] == "[50.0, 75.0)"].copy()
    methods = ["gmpe_raw", "idw", "kriging", "gmpe_kriging"]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharex=True)
    for method in methods:
        sub = high_density[high_density["method"].astype(str) == method].sort_values("observed_intensity_class")
        axes[0].plot(
            sub["observed_intensity_class"].astype(str),
            sub["mae"],
            marker="o",
            linewidth=1.7,
            color=colors.get(method),
            label=method,
        )
        axes[1].plot(
            sub["observed_intensity_class"].astype(str),
            sub["class_accuracy"],
            marker="o",
            linewidth=1.7,
            color=colors.get(method),
            label=method,
        )
    axes[0].set_ylabel("Mean absolute error")
    axes[0].set_title("Continuous error by observed class")
    axes[1].set_ylabel("Exact class hit rate")
    axes[1].set_ylim(0, 1.02)
    axes[1].set_title("Class accuracy by observed class")
    axes[0].set_xlabel("Observed intensity class")
    axes[1].set_xlabel("Observed intensity class")
    axes[1].legend(loc="upper right")
    savefig(fig, png_dir / "highest_density_error_by_observed_class.png")

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    sub = group_summary[
        group_summary["intensity_group"].isin(["I<5-", "I>=5-", "I>=6-"])
        & group_summary["method"].astype(str).isin(["idw", "kriging", "gmpe_kriging"])
    ].copy()
    sub = sub[sub["density_bin"] == "[50.0, 75.0)"]
    x_labels = ["I<5-", "I>=5-", "I>=6-"]
    x = np.arange(len(x_labels))
    offsets = {"idw": -0.22, "kriging": 0.0, "gmpe_kriging": 0.22}
    for method, group in sub.groupby("method", observed=True):
        y = group.set_index("intensity_group").reindex(x_labels)["under_class_rate"]
        ax.bar(x + offsets.get(str(method), 0.0), y, width=0.2, color=colors.get(str(method)), label=str(method))
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Under-predicted class fraction")
    ax.set_title("High observed intensities are more often underpredicted")
    ax.legend()
    savefig(fig, png_dir / "highest_density_underprediction_by_intensity_group.png")


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)

    df = read_predictions(args.predictions)
    class_summary = summarize_by_class(df)
    group_summary = summarize_high_groups(df)
    trend = trend_tests(df)

    class_summary.to_csv(args.csv_dir / "error_by_observed_intensity_class.csv", index=False)
    group_summary.to_csv(args.csv_dir / "error_by_intensity_group.csv", index=False)
    trend.to_csv(args.csv_dir / "intensity_error_trend_tests.csv", index=False)
    if not args.no_png:
        make_plots(class_summary, group_summary, args.png_dir)
    print(f"Saved CSV directory: {args.csv_dir}")
    if not args.no_png:
        print(f"Saved PNG directory: {args.png_dir}")


if __name__ == "__main__":
    main()
