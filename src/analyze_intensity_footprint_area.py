#!/usr/bin/env python3
"""Summarize areal footprints from estimated seismic-intensity grids."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_METHOD_SUMMARY = Path("outputs/csv/estimated_intensity_distribution_6lower_plus/estimated_intensity_method_summary.csv")
DEFAULT_TARGET_EVENTS = Path("outputs/csv/hypocenter_catalog/target_events_intensity_6lower_plus_with_hypocenter.csv")
DEFAULT_CSV_DIR = Path("outputs/csv/estimated_intensity_distribution_6lower_plus/synthesis")
DEFAULT_PNG_DIR = Path("outputs/png/estimated_intensity_distribution_6lower_plus/synthesis")
THRESHOLDS = (5.5, 6.0, 6.5, 7.0)
METHOD_ORDER = ["gmpe_raw", "gmpe_calibrated", "idw", "kriging", "gmpe_kriging"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute intensity-footprint areas from estimated intensity grids.")
    parser.add_argument("--method-summary", type=Path, default=DEFAULT_METHOD_SUMMARY)
    parser.add_argument("--target-events", type=Path, default=DEFAULT_TARGET_EVENTS)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--no-png", action="store_true")
    return parser.parse_args()


def inferred_spacing(values: pd.Series) -> float:
    unique = np.sort(values.dropna().astype(float).unique())
    if len(unique) < 2:
        return 0.02
    diffs = np.diff(unique)
    diffs = diffs[diffs > 1e-8]
    return float(np.nanmedian(diffs)) if len(diffs) else 0.02


def grid_cell_area_km2(grid: pd.DataFrame) -> np.ndarray:
    lon_spacing = inferred_spacing(grid["longitude"])
    lat_spacing = inferred_spacing(grid["latitude"])
    lat = grid["latitude"].to_numpy(dtype=float)
    return (lon_spacing * 111.32 * np.cos(np.radians(lat))) * (lat_spacing * 110.57)


def intensity_class_bin(value: float) -> str:
    if value >= 7.0:
        return "7"
    if value >= 6.5:
        return "6_upper"
    if value >= 6.0:
        return "6_lower"
    if value >= 5.5:
        return "5_upper"
    return "below_5_upper"


def period_bin(year: int) -> str:
    if year <= 1995:
        return "pre_1996"
    if year <= 2010:
        return "1996_2010"
    return "2011_2022"


def summarize_grid(row: pd.Series) -> dict[str, object]:
    grid_path = Path(str(row["grid_csv"]))
    grid = pd.read_csv(grid_path, usecols=["longitude", "latitude", "estimated_surface_intensity"])
    area = grid_cell_area_km2(grid)
    intensity = grid["estimated_surface_intensity"].to_numpy(dtype=float)
    total_area = float(np.nansum(area[np.isfinite(intensity)]))
    out: dict[str, object] = {
        "event_id": row["event_id"],
        "method": row["method"],
        "n_station_used": int(row["n_station_used"]) if "n_station_used" in row and pd.notna(row["n_station_used"]) else np.nan,
        "n_grid_cells": int(len(grid)),
        "grid_total_area_km2": total_area,
        "max_estimated_surface_intensity": float(np.nanmax(intensity)),
        "mean_estimated_surface_intensity": float(np.nanmean(intensity)),
        "grid_csv": str(grid_path),
    }
    for threshold in THRESHOLDS:
        key = str(threshold).replace(".", "p")
        mask = intensity >= threshold
        footprint = float(np.nansum(area[mask]))
        out[f"area_ge_{key}_km2"] = footprint
        out[f"area_ge_{key}_fraction"] = footprint / total_area if total_area > 0 else np.nan
    return out


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
            "grid.color": "#d7d7d7",
            "grid.linewidth": 0.55,
            "grid.alpha": 0.85,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
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


def make_plots(footprints: pd.DataFrame, method_summary: pd.DataFrame, png_dir: Path) -> None:
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
    labels = {
        "gmpe_raw": "GMPE only",
        "gmpe_calibrated": "GMPE calibrated",
        "idw": "IDW",
        "kriging": "Kriging",
        "gmpe_kriging": "GMPE-residual kriging",
    }

    long_rows = []
    for _, row in footprints.iterrows():
        for threshold in THRESHOLDS:
            key = str(threshold).replace(".", "p")
            long_rows.append(
                {
                    "event_id": row["event_id"],
                    "method": row["method"],
                    "threshold": threshold,
                    "area_km2": row[f"area_ge_{key}_km2"],
                    "max_observed_intensity": row["max_intensity_value"],
                    "period": row["period"],
                }
            )
    long = pd.DataFrame(long_rows)

    metric_rows = []
    for method, group in long.groupby("method"):
        for threshold, sub in group.groupby("threshold"):
            values = sub["area_km2"].to_numpy(dtype=float)
            metric_rows.append(
                {
                    "method": method,
                    "threshold": threshold,
                    "positive_rate": float(np.nanmean(values > 0.0)),
                    "p90_area_km2": float(np.nanquantile(values, 0.90)),
                }
            )
    metrics = pd.DataFrame(metric_rows)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    for method in METHOD_ORDER:
        sub = metrics[metrics["method"] == method].sort_values("threshold")
        axes[0].plot(
            sub["threshold"],
            sub["positive_rate"],
            marker="o",
            linewidth=1.7,
            color=colors[method],
            label=labels[method],
        )
        axes[1].plot(
            sub["threshold"],
            sub["p90_area_km2"],
            marker="o",
            linewidth=1.7,
            color=colors[method],
            label=labels[method],
        )
    axes[0].set_ylim(0, 1.02)
    axes[0].set_xlabel("Intensity threshold")
    axes[0].set_ylabel("Fraction of events with nonzero footprint")
    axes[0].set_title("Detection of high-intensity footprint")
    axes[1].set_yscale("symlog", linthresh=1.0)
    axes[1].set_xlabel("Intensity threshold")
    axes[1].set_ylabel("P90 footprint area (km$^2$)")
    axes[1].set_title("Upper-tail footprint size")
    axes[0].set_xticks(THRESHOLDS)
    axes[1].set_xticks(THRESHOLDS)
    axes[1].legend(loc="upper right")
    savefig(fig, png_dir / "intensity_footprint_area_by_method_threshold.png")

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.9), sharey=False)
    for ax, threshold in zip(axes, [5.5, 6.0, 6.5]):
        sub = long[(long["threshold"] == threshold) & (long["method"].isin(["gmpe_raw", "idw", "kriging", "gmpe_kriging"]))]
        for method, group in sub.groupby("method"):
            ax.scatter(
                group["max_observed_intensity"],
                group["area_km2"],
                s=18,
                alpha=0.62,
                color=colors.get(method),
                label=labels.get(method, method),
                edgecolor="white",
                linewidth=0.3,
            )
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_xlabel("Observed maximum intensity")
        ax.set_title(f"I >= {threshold:g} footprint")
    axes[0].set_ylabel("Estimated area (km$^2$)")
    axes[-1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))
    savefig(fig, png_dir / "footprint_area_vs_observed_max_intensity.png")

    raw = footprints[footprints["method"] == "gmpe_raw"].set_index("event_id")
    ratio_rows = []
    for _, row in footprints[footprints["method"].isin(["idw", "kriging", "gmpe_kriging"])].iterrows():
        if row["event_id"] not in raw.index:
            continue
        for threshold in [5.5, 6.0, 6.5]:
            key = str(threshold).replace(".", "p")
            denom = float(raw.loc[row["event_id"], f"area_ge_{key}_km2"])
            ratio_rows.append(
                {
                    "event_id": row["event_id"],
                    "method": row["method"],
                    "threshold": threshold,
                    "area_ratio_vs_gmpe_raw": (float(row[f"area_ge_{key}_km2"]) + 1.0) / (denom + 1.0),
                    "period": row["period"],
                }
            )
    ratio = pd.DataFrame(ratio_rows)
    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    offset = {"idw": -0.18, "kriging": 0.0, "gmpe_kriging": 0.18}
    for method, group in ratio.groupby("method"):
        xs = group["threshold"] + offset.get(method, 0.0)
        ax.scatter(xs, group["area_ratio_vs_gmpe_raw"], s=18, alpha=0.45, color=colors.get(method), label=method)
        med = group.groupby("threshold")["area_ratio_vs_gmpe_raw"].median()
        ax.plot(med.index + offset.get(method, 0.0), med.values, color=colors.get(method), linewidth=2.0)
    ax.axhline(1.0, color="#222222", linestyle="--", linewidth=1.0)
    ax.set_yscale("log")
    ax.set_xticks([5.5, 6.0, 6.5])
    ax.set_xlabel("Footprint threshold")
    ax.set_ylabel("(spatial method area + 1) / (GMPE-only area + 1)")
    ax.set_title("Observation-constrained footprint area relative to GMPE-only prediction")
    ax.legend(ncol=3)
    savefig(fig, png_dir / "footprint_area_ratio_vs_gmpe_raw.png")

    fig, ax = plt.subplots(figsize=(9.4, 4.6))
    sub = footprints[footprints["method"].isin(["idw", "kriging", "gmpe_kriging"])].copy()
    sub["year"] = sub["year"].astype(int)
    for method, group in sub.groupby("method"):
        yearly = group.groupby("year")["area_ge_5p5_km2"].median()
        ax.plot(yearly.index, yearly.values, marker="o", markersize=3.2, linewidth=1.6, color=colors.get(method), label=labels.get(method, method))
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("Year")
    ax.set_ylabel("Median I >= 5.5 footprint area (km$^2$)")
    ax.set_title("Annual median high-intensity footprint for Imax >= 6-lower target events")
    ax.legend(ncol=3)
    savefig(fig, png_dir / "annual_median_i55_footprint_area.png")


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    if not args.no_png:
        args.png_dir.mkdir(parents=True, exist_ok=True)

    method_summary = pd.read_csv(args.method_summary)
    target = pd.read_csv(args.target_events)
    event_cols = [
        "event_id",
        "intensity_origin_time",
        "year",
        "analysis_magnitude",
        "analysis_depth_km",
        "max_intensity_value",
        "hypocenter_region",
        "n_valid_intensity",
    ]
    target_meta = target[[col for col in event_cols if col in target.columns]].copy()
    target_meta["max_intensity_class_bin"] = target_meta["max_intensity_value"].map(intensity_class_bin)
    target_meta["period"] = target_meta["year"].astype(int).map(period_bin)

    rows = []
    for i, (_, row) in enumerate(method_summary.iterrows(), start=1):
        rows.append(summarize_grid(row))
        if i % 50 == 0:
            print(f"Processed {i}/{len(method_summary)} method grids")
    footprints = pd.DataFrame(rows).merge(target_meta, on="event_id", how="left")
    footprints.to_csv(args.csv_dir / "intensity_footprint_area_summary.csv", index=False)

    summary_rows = []
    for threshold in THRESHOLDS:
        key = str(threshold).replace(".", "p")
        col = f"area_ge_{key}_km2"
        grouped = (
            footprints.groupby(["method", "max_intensity_class_bin", "period"], dropna=False)
            .agg(
                n_events=("event_id", "nunique"),
                median_area_km2=(col, "median"),
                mean_area_km2=(col, "mean"),
                p90_area_km2=(col, lambda x: float(np.nanquantile(x, 0.90))),
                median_max_estimated_intensity=("max_estimated_surface_intensity", "median"),
                median_n_station_used=("n_station_used", "median") if "n_station_used" in footprints else ("event_id", "size"),
            )
            .reset_index()
        )
        grouped["threshold"] = threshold
        summary_rows.append(grouped)
    grouped_summary = pd.concat(summary_rows, ignore_index=True)
    grouped_summary.to_csv(args.csv_dir / "intensity_footprint_area_group_summary.csv", index=False)

    method_rows = []
    for method, group in footprints.groupby("method"):
        row: dict[str, object] = {"method": method, "n_events": int(group["event_id"].nunique())}
        for threshold in THRESHOLDS:
            key = str(threshold).replace(".", "p")
            col = f"area_ge_{key}_km2"
            values = group[col].to_numpy(dtype=float)
            row[f"area_ge_{key}_positive_event_rate"] = float(np.nanmean(values > 0.0))
            row[f"area_ge_{key}_median_km2"] = float(np.nanmedian(values))
            row[f"area_ge_{key}_mean_km2"] = float(np.nanmean(values))
            row[f"area_ge_{key}_p75_km2"] = float(np.nanquantile(values, 0.75))
            row[f"area_ge_{key}_p90_km2"] = float(np.nanquantile(values, 0.90))
            row[f"area_ge_{key}_max_km2"] = float(np.nanmax(values))
        row["median_max_estimated_surface_intensity"] = float(np.nanmedian(group["max_estimated_surface_intensity"]))
        row["p90_max_estimated_surface_intensity"] = float(np.nanquantile(group["max_estimated_surface_intensity"], 0.90))
        method_rows.append(row)
    method_area_summary = pd.DataFrame(method_rows)
    method_area_summary.to_csv(args.csv_dir / "intensity_footprint_area_method_summary.csv", index=False)

    event_wide = footprints.pivot_table(
        index=["event_id", "intensity_origin_time", "year", "analysis_magnitude", "analysis_depth_km", "max_intensity_value", "hypocenter_region"],
        columns="method",
        values=["max_estimated_surface_intensity", "area_ge_5p5_km2", "area_ge_6p0_km2", "area_ge_6p5_km2"],
        aggfunc="first",
    )
    event_wide.columns = ["__".join([str(v) for v in col if v]) for col in event_wide.columns]
    event_wide = event_wide.reset_index()
    event_wide.to_csv(args.csv_dir / "intensity_footprint_area_event_wide.csv", index=False)

    if not args.no_png:
        make_plots(footprints, method_summary, args.png_dir)

    print(f"Saved: {args.csv_dir / 'intensity_footprint_area_summary.csv'}")
    print(f"Saved: {args.csv_dir / 'intensity_footprint_area_group_summary.csv'}")
    print(f"Saved: {args.csv_dir / 'intensity_footprint_area_method_summary.csv'}")
    print(f"Saved: {args.csv_dir / 'intensity_footprint_area_event_wide.csv'}")
    if not args.no_png:
        print(f"Saved PNG directory: {args.png_dir}")


if __name__ == "__main__":
    main()
