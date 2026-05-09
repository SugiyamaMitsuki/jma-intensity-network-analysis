#!/usr/bin/env python3
"""Plot observed versus predicted intensity for station-thinning validation."""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_PREDICTIONS = Path(
    "outputs/csv/station_thinning_interpolation_6lower_plus_class/thinning_prediction_errors.csv.gz"
)
DEFAULT_CSV_DIR = Path("outputs/csv/observed_predicted_intensity")
DEFAULT_PNG_DIR = Path("outputs/png/observed_predicted_intensity")
DEFAULT_SUMMARY = Path("data/derived/observed_predicted_intensity_scatter_summary.csv")
DEFAULT_ASSET_JA = Path("paper/assets_ja/fig05b_observed_predicted_intensity_scatter.png")
DEFAULT_ASSET_EN = Path("paper/assets_en/fig05b_observed_predicted_intensity_scatter.png")

PANEL_SPECS = [
    ("gmpe_raw", "zero", "a  Simplified attenuation baseline", "Zero retained stations"),
    ("idw", "high_density", "b  IDW", "Highest density bin: 50-75 / 10,000 km2"),
    ("kriging", "high_density", "c  Smoothing kriging", "Highest density bin: 50-75 / 10,000 km2"),
    ("gmpe_kriging", "high_density", "d  Attenuation-residual kriging", "Highest density bin: 50-75 / 10,000 km2"),
]
CLASS_BOUNDS = [0.5, 1.5, 2.5, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5, 7.25]
CLASS_TICKS = [1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 6.5, 7.0]
CLASS_LABELS = ["1", "2", "3", "4", "5-", "5+", "6-", "6+", "7"]
AXIS_LIMITS = (0.8, 7.15)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--derived-summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--assets-ja", type=Path, default=DEFAULT_ASSET_JA)
    parser.add_argument("--assets-en", type=Path, default=DEFAULT_ASSET_EN)
    parser.add_argument("--high-density-min", type=float, default=50.0)
    parser.add_argument("--high-density-max", type=float, default=75.0)
    return parser.parse_args()


def read_predictions(path: Path) -> pd.DataFrame:
    columns = [
        "event_id",
        "station_code",
        "method",
        "observed_intensity",
        "predicted_intensity",
        "prediction_error",
        "abs_error",
        "class_error",
        "class_match",
        "observed_intensity_class_index",
        "effective_train_density_per_10000km2",
        "effective_train_median_nn_distance_km",
    ]
    df = pd.read_csv(path, usecols=columns)
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["observed_intensity", "predicted_intensity"])
    df["class_abs_error"] = np.abs(df["class_error"].to_numpy(dtype=float))
    return df


def panel_data(df: pd.DataFrame, method: str, scope: str, density_min: float, density_max: float) -> pd.DataFrame:
    sub = df[df["method"] == method].copy()
    if scope == "zero":
        sub = sub[sub["effective_train_density_per_10000km2"] == 0.0]
    elif scope == "high_density":
        density = sub["effective_train_density_per_10000km2"]
        sub = sub[(density >= density_min) & (density < density_max)]
    else:
        raise ValueError(f"Unknown panel scope: {scope}")
    return sub.reset_index(drop=True)


def rmse(values: pd.Series) -> float:
    arr = values.to_numpy(dtype=float)
    return float(math.sqrt(np.mean(arr * arr)))


def finite_median(values: pd.Series) -> float:
    arr = values.to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan")
    return float(np.median(arr))


def summarize_panel(data: pd.DataFrame, method: str, scope: str, label: str) -> dict[str, object]:
    high = data[data["observed_intensity_class_index"] >= 5]
    return {
        "method": method,
        "scope": scope,
        "label": label,
        "n_predictions": int(len(data)),
        "median_density_per_10000km2": finite_median(data["effective_train_density_per_10000km2"]),
        "median_train_nn_km": finite_median(data["effective_train_median_nn_distance_km"]),
        "bias": float(data["prediction_error"].mean()),
        "mae": float(data["abs_error"].mean()),
        "rmse": rmse(data["prediction_error"]),
        "exact_class_accuracy": float(data["class_match"].mean()),
        "within_one_class": float((data["class_abs_error"] <= 1).mean()),
        "underprediction_rate": float((data["class_error"] < 0).mean()),
        "overprediction_rate": float((data["class_error"] > 0).mean()),
        "n_predictions_i_ge_5lower": int(len(high)),
        "exact_class_accuracy_i_ge_5lower": float(high["class_match"].mean()) if len(high) else np.nan,
        "underprediction_rate_i_ge_5lower": float((high["class_error"] < 0).mean()) if len(high) else np.nan,
    }


def setup_plot_style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 360,
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#222222",
            "axes.linewidth": 0.8,
            "axes.grid": False,
            "axes.labelsize": 8.6,
            "axes.titlesize": 9.2,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 7.8,
            "legend.fontsize": 7.4,
        }
    )


def draw_class_guides(ax) -> None:
    for lo, hi in zip(CLASS_BOUNDS[:-1], CLASS_BOUNDS[1:]):
        ax.add_patch(
            plt_rectangle(
                (lo, lo),
                hi - lo,
                hi - lo,
                facecolor="#7fcdbb",
                alpha=0.12,
                edgecolor="none",
                zorder=0,
            )
        )
    for value in CLASS_BOUNDS[:-1]:
        ax.axvline(value, color="#d0d0d0", linewidth=0.45, zorder=1)
        ax.axhline(value, color="#d0d0d0", linewidth=0.45, zorder=1)
    x0, x1 = AXIS_LIMITS
    ax.plot([x0, x1], [x0, x1], color="#111111", linewidth=1.0, zorder=3)
    ax.plot([x0, x1], [x0 - 0.5, x1 - 0.5], color="#555555", linewidth=0.65, linestyle=":", zorder=3)
    ax.plot([x0, x1], [x0 + 0.5, x1 + 0.5], color="#555555", linewidth=0.65, linestyle=":", zorder=3)


def plt_rectangle(*args, **kwargs):
    from matplotlib.patches import Rectangle

    return Rectangle(*args, **kwargs)


def draw_panel(ax, data: pd.DataFrame, summary: dict[str, object], title: str, subtitle: str):
    from matplotlib.colors import LogNorm

    draw_class_guides(ax)
    hb = ax.hexbin(
        data["observed_intensity"],
        data["predicted_intensity"],
        gridsize=46,
        extent=(AXIS_LIMITS[0], AXIS_LIMITS[1], AXIS_LIMITS[0], AXIS_LIMITS[1]),
        mincnt=1,
        cmap="Greys",
        norm=LogNorm(vmin=1),
        linewidths=0,
        alpha=0.90,
        zorder=2,
    )

    high = data[data["observed_intensity_class_index"] >= 5].copy()
    if len(high) > 6000:
        high = high.sample(6000, random_state=20260509)
    if len(high):
        under = high["class_error"] < 0
        exact = high["class_error"] == 0
        over = high["class_error"] > 0
        ax.scatter(
            high.loc[under, "observed_intensity"],
            high.loc[under, "predicted_intensity"],
            s=10,
            c="#b2182b",
            alpha=0.58,
            linewidths=0,
            label="I>=5-, under",
            zorder=4,
        )
        ax.scatter(
            high.loc[exact, "observed_intensity"],
            high.loc[exact, "predicted_intensity"],
            s=9,
            c="#2166ac",
            alpha=0.52,
            linewidths=0,
            label="I>=5-, exact",
            zorder=4,
        )
        ax.scatter(
            high.loc[over, "observed_intensity"],
            high.loc[over, "predicted_intensity"],
            s=9,
            c="#f4a582",
            alpha=0.45,
            linewidths=0,
            label="I>=5-, over",
            zorder=4,
        )

    text = (
        f"n={summary['n_predictions']:,}\n"
        f"RMSE={summary['rmse']:.2f}, MAE={summary['mae']:.2f}\n"
        f"Exact={summary['exact_class_accuracy']:.2f}, under={summary['underprediction_rate']:.2f}"
    )
    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=7.4,
        color="#111111",
        bbox={"boxstyle": "round,pad=0.24", "facecolor": "white", "edgecolor": "#777777", "linewidth": 0.45, "alpha": 0.88},
        zorder=5,
    )
    ax.set_title(f"{title}\n{subtitle}", loc="left", pad=5)
    ax.set_xlim(*AXIS_LIMITS)
    ax.set_ylim(*AXIS_LIMITS)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(CLASS_TICKS)
    ax.set_yticks(CLASS_TICKS)
    ax.set_xticklabels(CLASS_LABELS)
    ax.set_yticklabels(CLASS_LABELS)
    return hb


def draw_scatter_figure(panel_frames: list[pd.DataFrame], summaries: list[dict[str, object]], output_path: Path) -> None:
    import matplotlib.pyplot as plt

    setup_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.9, 9.7), constrained_layout=False)
    fig.subplots_adjust(left=0.075, right=0.875, bottom=0.105, top=0.955, wspace=0.17, hspace=0.34)
    last_hb = None
    for ax, (_, _, title, subtitle), frame, summary in zip(axes.ravel(), PANEL_SPECS, panel_frames, summaries):
        last_hb = draw_panel(ax, frame, summary, title, subtitle)

    for ax in axes[:, 0]:
        ax.set_ylabel("Predicted JMA intensity")
    for ax in axes[-1, :]:
        ax.set_xlabel("Observed JMA intensity")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.475, 0.025))
    if last_hb is not None:
        cbar_ax = fig.add_axes([0.902, 0.235, 0.028, 0.52])
        cbar = fig.colorbar(last_hb, cax=cbar_ax)
        cbar.set_label("Withheld-station count per hexagon")
        cbar.ax.tick_params(labelsize=7.2)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)

    df = read_predictions(args.predictions)
    panel_frames: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []
    for method, scope, label, _subtitle in PANEL_SPECS:
        frame = panel_data(df, method, scope, args.high_density_min, args.high_density_max)
        panel_frames.append(frame)
        summaries.append(summarize_panel(frame, method, scope, label))

    summary = pd.DataFrame(summaries)
    summary_path = args.csv_dir / "observed_predicted_intensity_scatter_summary.csv"
    summary.to_csv(summary_path, index=False)
    args.derived_summary.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.derived_summary, index=False)

    figure_path = args.png_dir / "observed_predicted_intensity_scatter.png"
    draw_scatter_figure(panel_frames, summaries, figure_path)
    for asset in (args.assets_ja, args.assets_en):
        asset.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(figure_path, asset)

    print(f"Summary CSV -> {summary_path}")
    print(f"Derived summary -> {args.derived_summary}")
    print(f"PNG -> {figure_path}")


if __name__ == "__main__":
    main()
