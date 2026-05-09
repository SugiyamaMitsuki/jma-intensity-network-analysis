#!/usr/bin/env python3
"""Counterfactual network-density analysis for the 2018 northern Osaka earthquake."""

from __future__ import annotations

import argparse
import math
import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from analyze_jma_intensity import active_station_records_at_year_end, haversine_km, load_station_index
from analyze_station_thinning_interpolation import error_metrics, intensity_class_index, intensity_class_labels
from estimate_jma_intensity_distribution import (
    DEFAULT_GROUND_GRID,
    DEFAULT_TARGET_CSV,
    EventContext,
    coarsen_ground_grid,
    crop_ground_grid,
    draw_map,
    estimate_event_method,
    event_context,
    inferred_grid_spacing,
    load_event_station_observations,
    load_ground_grid,
    nearest_ground_values,
    project_km,
    site_intensity_delta,
)
from map_prefecture_boundaries import DEFAULT_PREFECTURE_BOUNDARY, read_gmt_multisegment


DEFAULT_TARGET_6LOWER_CSV = Path("outputs/csv/hypocenter_catalog/target_events_intensity_6lower_plus_with_hypocenter.csv")
DEFAULT_CSV_DIR = Path("outputs/csv/osaka_2018_network_counterfactual")
DEFAULT_PNG_DIR = Path("outputs/png/osaka_2018_network_counterfactual")
DEFAULT_MANUSCRIPT_ASSETS_JA = Path("outputs/manuscript/assets_ja")
DEFAULT_MANUSCRIPT_ASSETS_EN = Path("outputs/manuscript/assets_en")
INTENSITY_THRESHOLDS = [4.5, 5.0, 5.5, 6.0, 6.5]
RADII_KM = [10.0, 20.0, 50.0, 100.0]
SCENARIO_LABELS = {
    "current_2018_observed": "2018 observed network",
    "counterfactual_1994_active": "1994 active-site geometry",
    "strict_1994_code_overlap": "Strict 1994 code overlap",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare estimated intensity maps for the 2018 northern Osaka earthquake "
            "using the actual 2018 network and a counterfactual pre-expansion network."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--target-events-csv", type=Path, default=DEFAULT_TARGET_6LOWER_CSV)
    parser.add_argument("--ground-grid", type=Path, default=DEFAULT_GROUND_GRID)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--manuscript-assets-ja", type=Path, default=DEFAULT_MANUSCRIPT_ASSETS_JA)
    parser.add_argument("--manuscript-assets-en", type=Path, default=DEFAULT_MANUSCRIPT_ASSETS_EN)
    parser.add_argument("--event-id", default="i2018_000842")
    parser.add_argument("--counterfactual-year", type=int, default=1994)
    parser.add_argument("--radius-km", type=float, default=100.0)
    parser.add_argument("--map-radius-km", type=float, default=115.0)
    parser.add_argument("--match-radius-km", type=float, default=10.0)
    parser.add_argument("--grid-spacing", type=float, default=0.02)
    parser.add_argument("--methods", default="gmpe_raw,idw,kriging,gmpe_kriging")
    parser.add_argument("--map-method", default="gmpe_kriging")
    parser.add_argument("--min-station-intensity", type=float, default=1.0)
    parser.add_argument("--amp-intensity-coef", type=float, default=1.72)
    parser.add_argument("--idw-neighbors", type=int, default=16)
    parser.add_argument("--idw-power", type=float, default=2.0)
    parser.add_argument("--rbf-neighbors", type=int, default=80)
    parser.add_argument("--rbf-smoothing", type=float, default=0.08)
    parser.add_argument("--kriging-neighbors", type=int, default=18)
    parser.add_argument("--kriging-range-km", type=float, default=None)
    parser.add_argument("--kriging-nugget", type=float, default=0.02)
    parser.add_argument("--gmpe-fault-type", default="crustal", choices=["crustal", "interplate", "intraslab"])
    parser.add_argument("--gmpe-min-distance-km", type=float, default=3.0)
    parser.add_argument("--gmpe-bias-stat", default="median", choices=["median", "mean", "none"])
    parser.add_argument("--pgv-intensity-intercept", type=float, default=2.68)
    parser.add_argument("--kriging-max-grid-cells", type=int, default=15_000)
    parser.add_argument("--clip-min", type=float, default=0.0)
    parser.add_argument("--clip-max", type=float, default=7.2)
    parser.add_argument("--allow-overshoot", action="store_true")
    parser.add_argument("--prefecture-boundary", type=Path, default=DEFAULT_PREFECTURE_BOUNDARY)
    parser.add_argument("--no-pygmt-maps", action="store_true")
    return parser.parse_args()


def estimator_args(args: argparse.Namespace) -> SimpleNamespace:
    return SimpleNamespace(
        mode="auto",
        shallow_depth_threshold_km=150.0,
        amp_intensity_coef=args.amp_intensity_coef,
        idw_neighbors=args.idw_neighbors,
        idw_power=args.idw_power,
        rbf_neighbors=args.rbf_neighbors,
        rbf_smoothing=args.rbf_smoothing,
        kriging_neighbors=args.kriging_neighbors,
        kriging_range_km=args.kriging_range_km,
        kriging_nugget=args.kriging_nugget,
        gmpe_fault_type=args.gmpe_fault_type,
        gmpe_min_distance_km=args.gmpe_min_distance_km,
        gmpe_bias_stat=args.gmpe_bias_stat,
        pgv_intensity_intercept=args.pgv_intensity_intercept,
        kriging_max_grid_cells=args.kriging_max_grid_cells,
        clip_min=args.clip_min,
        clip_max=args.clip_max,
        allow_overshoot=args.allow_overshoot,
    )


def load_event(args: argparse.Namespace) -> EventContext:
    candidates = [args.target_events_csv, DEFAULT_TARGET_CSV]
    for path in candidates:
        if not path.exists():
            continue
        df = pd.read_csv(path, low_memory=False)
        hit = df[df["event_id"].astype(str) == args.event_id]
        if not hit.empty:
            return event_context(hit.iloc[0])
    raise SystemExit(f"Event not found in target catalog: {args.event_id}")


def region_from_radius(event: EventContext, radius_km: float) -> list[float]:
    lat_delta = radius_km / 110.57
    lon_delta = radius_km / (111.32 * math.cos(math.radians(event.latitude)))
    return [
        max(122.0, event.longitude - lon_delta),
        min(146.5, event.longitude + lon_delta),
        max(24.0, event.latitude - lat_delta),
        min(46.5, event.latitude + lat_delta),
    ]


def in_radius(df: pd.DataFrame, event: EventContext, radius_km: float) -> pd.DataFrame:
    out = df.copy()
    if "epicentral_distance_km" not in out:
        out["epicentral_distance_km"] = [
            haversine_km(event.latitude, event.longitude, lat, lon)
            for lat, lon in zip(out["latitude"], out["longitude"])
        ]
    return out[out["epicentral_distance_km"].astype(float) <= radius_km].copy()


def load_current_stations(args: argparse.Namespace, event: EventContext) -> pd.DataFrame:
    stations = load_event_station_observations(event, args.data_dir)
    stations = stations.dropna(subset=["latitude", "longitude", "observed_intensity"]).copy()
    stations = stations[stations["observed_intensity"].astype(float) >= args.min_station_intensity].copy()
    stations = in_radius(stations, event, args.radius_km)
    if len(stations) < 3:
        raise SystemExit(f"{event.event_id}: too few current stations within {args.radius_km:g} km.")
    return stations.reset_index(drop=True)


def build_counterfactual_stations(
    args: argparse.Namespace,
    event: EventContext,
    current_all: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    station_index = load_station_index(args.data_dir)
    if station_index is None:
        raise SystemExit(f"Station index not found under {args.data_dir}")

    active_records = active_station_records_at_year_end(station_index, args.counterfactual_year)
    obs_xy = project_km(
        current_all["longitude"].to_numpy(dtype=float),
        current_all["latitude"].to_numpy(dtype=float),
        event.longitude,
        event.latitude,
    )
    tree = cKDTree(obs_xy)
    obs = current_all.reset_index(drop=True)
    rows: list[dict[str, object]] = []
    match_rows: list[dict[str, object]] = []
    for rec in active_records:
        if rec.latitude is None or rec.longitude is None:
            continue
        epicentral_distance = haversine_km(event.latitude, event.longitude, rec.latitude, rec.longitude)
        if epicentral_distance is None or epicentral_distance > args.radius_km:
            continue
        query_xy = project_km(np.array([rec.longitude]), np.array([rec.latitude]), event.longitude, event.latitude)
        nearest_distance, nearest_idx = tree.query(query_xy, k=1)
        nearest_distance = float(np.ravel(nearest_distance)[0])
        nearest_idx = int(np.ravel(nearest_idx)[0])
        src = obs.iloc[nearest_idx]
        keep = nearest_distance <= args.match_radius_km
        match_rows.append(
            {
                "counterfactual_year": args.counterfactual_year,
                "station_code": rec.code,
                "station_name": rec.name,
                "latitude": rec.latitude,
                "longitude": rec.longitude,
                "epicentral_distance_km": epicentral_distance,
                "matched": keep,
                "matched_distance_km": nearest_distance,
                "matched_current_station_code": src["station_code"],
                "matched_current_station_name": src["station_name"],
                "matched_current_observed_intensity": src["observed_intensity"],
                "matched_current_epicentral_distance_km": src["epicentral_distance_km"],
            }
        )
        if not keep:
            continue
        rows.append(
            {
                "event_id": event.event_id,
                "station_code": f"cf{args.counterfactual_year}_{rec.code}",
                "station_name": rec.name,
                "latitude": rec.latitude,
                "longitude": rec.longitude,
                "intensity_code": src.get("intensity_code", ""),
                "intensity_class_value": src.get("intensity_class_value", np.nan),
                "measured_intensity": src.get("measured_intensity", np.nan),
                "observed_intensity": float(src["observed_intensity"]),
                "pga_total_gal": src.get("pga_total_gal", np.nan),
                "epicentral_distance_km": epicentral_distance,
                "source_station_code": src["station_code"],
                "source_station_name": src["station_name"],
                "source_station_distance_km": float(src["epicentral_distance_km"]),
                "pseudo_match_distance_km": nearest_distance,
            }
        )
    counterfactual = pd.DataFrame(rows)
    matches = pd.DataFrame(match_rows)
    if len(counterfactual) < 3:
        raise SystemExit(
            f"Only {len(counterfactual)} matched {args.counterfactual_year} stations within "
            f"{args.radius_km:g} km. Increase --match-radius-km."
        )
    return counterfactual.reset_index(drop=True), matches


def build_strict_overlap(current: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    station_index = load_station_index(args.data_dir)
    if station_index is None:
        return pd.DataFrame()
    active_codes = {rec.code for rec in active_station_records_at_year_end(station_index, args.counterfactual_year)}
    strict = current[current["station_code"].astype(str).isin(active_codes)].copy()
    return strict.reset_index(drop=True)


def station_nearest_distance_summary(stations: pd.DataFrame, event: EventContext) -> dict[str, float]:
    out: dict[str, float] = {}
    if len(stations) >= 2:
        xy = project_km(
            stations["longitude"].to_numpy(dtype=float),
            stations["latitude"].to_numpy(dtype=float),
            event.longitude,
            event.latitude,
        )
        dist, _ = cKDTree(xy).query(xy, k=2)
        nn = dist[:, 1]
        out["station_nn_median_km"] = float(np.nanmedian(nn))
        out["station_nn_mean_km"] = float(np.nanmean(nn))
        out["station_nn_q90_km"] = float(np.nanquantile(nn, 0.90))
    else:
        out["station_nn_median_km"] = np.nan
        out["station_nn_mean_km"] = np.nan
        out["station_nn_q90_km"] = np.nan
    return out


def summarize_stations(
    scenario: str,
    stations: pd.DataFrame,
    event: EventContext,
    radius_km: float,
) -> dict[str, object]:
    distances = stations["epicentral_distance_km"].astype(float).to_numpy()
    summary: dict[str, object] = {
        "event_id": event.event_id,
        "scenario": scenario,
        "scenario_label": SCENARIO_LABELS.get(scenario, scenario),
        "n_station_used": int(len(stations)),
        "radius_km": radius_km,
        "density_per_10000km2": float(len(stations) / (math.pi * radius_km**2) * 10_000.0),
        "nearest_epicentral_distance_km": float(np.nanmin(distances)),
        "median_epicentral_distance_km": float(np.nanmedian(distances)),
        "q90_epicentral_distance_km": float(np.nanquantile(distances, 0.90)),
        "max_observed_intensity": float(np.nanmax(stations["observed_intensity"].astype(float))),
        "mean_observed_intensity": float(np.nanmean(stations["observed_intensity"].astype(float))),
    }
    for radius in RADII_KM:
        summary[f"station_count_within_{int(radius)}km"] = int(np.sum(distances <= radius))
    summary.update(station_nearest_distance_summary(stations, event))
    return summary


def grid_cell_area_km2(grid: pd.DataFrame) -> np.ndarray:
    lon_spacing = inferred_grid_spacing(grid[["longitude", "latitude"]])
    lat_values = np.sort(grid["latitude"].astype(float).unique())
    lat_diffs = np.diff(lat_values)
    lat_diffs = lat_diffs[lat_diffs > 1e-8]
    lat_spacing = float(np.nanmedian(lat_diffs)) if len(lat_diffs) else lon_spacing
    lat = grid["latitude"].to_numpy(dtype=float)
    return (lon_spacing * 111.32 * np.cos(np.radians(lat))) * (lat_spacing * 110.57)


def summarize_footprint(scenario: str, method: str, grid: pd.DataFrame) -> dict[str, object]:
    intensity = grid["estimated_surface_intensity"].to_numpy(dtype=float)
    area = grid_cell_area_km2(grid)
    total = float(np.nansum(area[np.isfinite(intensity)]))
    row: dict[str, object] = {
        "scenario": scenario,
        "method": method,
        "n_grid_cells": int(len(grid)),
        "grid_total_area_km2": total,
        "mean_estimated_surface_intensity": float(np.nanmean(intensity)),
        "max_estimated_surface_intensity": float(np.nanmax(intensity)),
    }
    for threshold in INTENSITY_THRESHOLDS:
        key = str(threshold).replace(".", "p")
        footprint = float(np.nansum(area[intensity >= threshold]))
        row[f"area_ge_{key}_km2"] = footprint
        row[f"area_ge_{key}_fraction"] = footprint / total if total > 0 else np.nan
    return row


def predict_grid_at_points(grid: pd.DataFrame, points: pd.DataFrame, event: EventContext) -> np.ndarray:
    grid_xy = project_km(
        grid["longitude"].to_numpy(dtype=float),
        grid["latitude"].to_numpy(dtype=float),
        event.longitude,
        event.latitude,
    )
    point_xy = project_km(
        points["longitude"].to_numpy(dtype=float),
        points["latitude"].to_numpy(dtype=float),
        event.longitude,
        event.latitude,
    )
    _, idx = cKDTree(grid_xy).query(point_xy, k=1)
    return grid["estimated_surface_intensity"].to_numpy(dtype=float)[idx]


def validation_rows(
    scenario: str,
    method: str,
    grid: pd.DataFrame,
    validation: pd.DataFrame,
    event: EventContext,
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    pred = predict_grid_at_points(grid, validation, event)
    truth = validation["observed_intensity"].to_numpy(dtype=float)
    errors = pred - truth
    prediction_df = validation[
        ["event_id", "station_code", "station_name", "latitude", "longitude", "observed_intensity", "epicentral_distance_km"]
    ].copy()
    prediction_df["scenario"] = scenario
    prediction_df["method"] = method
    prediction_df["predicted_intensity"] = pred
    prediction_df["prediction_error"] = errors
    obs_class = intensity_class_index(truth)
    pred_class = intensity_class_index(pred)
    prediction_df["observed_intensity_class"] = intensity_class_labels(obs_class)
    prediction_df["predicted_intensity_class"] = intensity_class_labels(pred_class)
    prediction_df["class_match"] = obs_class == pred_class

    rows: list[dict[str, object]] = []
    subsets = [
        ("all_I>=1", np.isfinite(truth)),
        ("I>=4.5", truth >= 4.5),
        ("I>=5.0", truth >= 5.0),
        ("I>=5.5", truth >= 5.5),
        ("I>=6.0", truth >= 6.0),
    ]
    for subset, mask in subsets:
        if not np.any(mask):
            continue
        metrics = error_metrics(errors[mask], truth[mask], pred[mask])
        metrics.update(
            {
                "scenario": scenario,
                "method": method,
                "validation_subset": subset,
                "n_validation": int(np.sum(mask)),
            }
        )
        rows.append(metrics)
    return rows, prediction_df


def configure_matplotlib() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 330,
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.labelsize": 7.6,
            "axes.titlesize": 8.2,
            "axes.edgecolor": "#222222",
            "axes.linewidth": 0.8,
            "axes.grid": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 6.8,
            "legend.frameon": False,
            "xtick.color": "#222222",
            "ytick.color": "#222222",
            "axes.labelcolor": "#222222",
            "axes.titleweight": "bold",
        }
    )


def pivot_grid(grid: pd.DataFrame, column: str = "estimated_surface_intensity") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lons = np.sort(grid["longitude"].astype(float).unique())
    lats = np.sort(grid["latitude"].astype(float).unique())
    arr = np.full((len(lats), len(lons)), np.nan)
    lon_map = {round(v, 8): i for i, v in enumerate(lons)}
    lat_map = {round(v, 8): i for i, v in enumerate(lats)}
    for lon, lat, value in zip(grid["longitude"], grid["latitude"], grid[column]):
        arr[lat_map[round(float(lat), 8)], lon_map[round(float(lon), 8)]] = float(value)
    return lons, lats, arr


def draw_boundaries(ax, region: list[float], color: str = "#595959", linewidth: float = 0.35) -> None:
    lon_min, lon_max, lat_min, lat_max = region
    for segment in read_gmt_multisegment(DEFAULT_PREFECTURE_BOUNDARY):
        mask = (
            (segment[:, 0] >= lon_min - 0.2)
            & (segment[:, 0] <= lon_max + 0.2)
            & (segment[:, 1] >= lat_min - 0.2)
            & (segment[:, 1] <= lat_max + 0.2)
        )
        if np.any(mask):
            ax.plot(segment[:, 0], segment[:, 1], color=color, linewidth=linewidth, alpha=0.85, zorder=4)


def setup_map_axis(ax, region: list[float], event: EventContext) -> None:
    ax.set_xlim(region[0], region[1])
    ax.set_ylim(region[2], region[3])
    ax.set_aspect(1.0 / math.cos(math.radians(event.latitude)))
    ax.set_xlabel("Longitude", fontsize=7.6)
    ax.set_ylabel("Latitude", fontsize=7.6)
    ax.tick_params(labelsize=7.5, length=2.5)


def plot_map_comparison(
    event: EventContext,
    region: list[float],
    grids: dict[tuple[str, str], pd.DataFrame],
    stations: dict[str, pd.DataFrame],
    method: str,
    png_dir: Path,
    intensity_to_log10_pgv_coef: float,
) -> Path:
    import matplotlib.pyplot as plt

    current = grids[("current_2018_observed", method)]
    counter = grids[("counterfactual_1994_active", method)]
    merged = current[["longitude", "latitude", "estimated_surface_intensity"]].merge(
        counter[["longitude", "latitude", "estimated_surface_intensity"]],
        on=["longitude", "latitude"],
        suffixes=("_current", "_1994"),
    )
    merged["difference"] = merged["estimated_surface_intensity_1994"] - merged["estimated_surface_intensity_current"]
    merged["log10_pgv_ratio"] = merged["difference"] / intensity_to_log10_pgv_coef
    lons, lats, current_arr = pivot_grid(current)
    _, _, counter_arr = pivot_grid(counter)
    _, _, diff_arr = pivot_grid(merged, "difference")
    _, _, log_ratio_arr = pivot_grid(merged, "log10_pgv_ratio")
    fig, axes = plt.subplots(1, 4, figsize=(12.8, 3.55), constrained_layout=True)
    cmap_intensity = "turbo"
    levels = [4.5, 5.0, 5.5, 6.0, 6.5]
    image0 = None
    for ax, arr, scenario, title in [
        (axes[0], current_arr, "current_2018_observed", "2018 observed network"),
        (axes[1], counter_arr, "counterfactual_1994_active", "1994 active-site geometry"),
    ]:
        image0 = ax.pcolormesh(lons, lats, arr, cmap=cmap_intensity, vmin=0.0, vmax=7.0, shading="nearest")
        ax.contour(lons, lats, arr, levels=levels, colors="white", linewidths=0.45, alpha=0.85)
        draw_boundaries(ax, region)
        sub = stations[scenario]
        ax.scatter(
            sub["longitude"],
            sub["latitude"],
            s=8 if scenario == "current_2018_observed" else 14,
            facecolors="white",
            edgecolors="#111111",
            linewidths=0.25,
            alpha=0.78,
            zorder=5,
        )
        ax.scatter([event.longitude], [event.latitude], marker="*", s=72, c="#111111", edgecolors="white", linewidths=0.45, zorder=6)
        ax.set_title(f"{title}\nN={len(sub)} within 100 km", fontsize=8.2)
        setup_map_axis(ax, region, event)
    image2 = axes[2].pcolormesh(lons, lats, diff_arr, cmap="RdBu_r", vmin=-1.2, vmax=1.2, shading="nearest")
    axes[2].contour(lons, lats, diff_arr, levels=[-0.5, 0.0, 0.5], colors=["#2b6cb0", "#222222", "#b83232"], linewidths=0.55)
    draw_boundaries(axes[2], region)
    axes[2].scatter([event.longitude], [event.latitude], marker="*", s=72, c="#111111", edgecolors="white", linewidths=0.45, zorder=6)
    axes[2].set_title("Counterfactual minus current\nestimated intensity", fontsize=8.2)
    setup_map_axis(axes[2], region, event)
    image3 = axes[3].pcolormesh(lons, lats, log_ratio_arr, cmap="RdBu_r", vmin=-0.65, vmax=0.65, shading="nearest")
    axes[3].contour(
        lons,
        lats,
        log_ratio_arr,
        levels=[math.log10(0.5), 0.0, math.log10(2.0)],
        colors=["#2b6cb0", "#222222", "#b83232"],
        linewidths=0.55,
    )
    draw_boundaries(axes[3], region)
    axes[3].scatter([event.longitude], [event.latitude], marker="*", s=72, c="#111111", edgecolors="white", linewidths=0.45, zorder=6)
    axes[3].set_title("Equivalent PGV ratio\n1994/current", fontsize=8.2)
    setup_map_axis(axes[3], region, event)
    fig.colorbar(image0, ax=axes[:2], location="bottom", fraction=0.05, pad=0.08, label="Estimated JMA instrumental intensity")
    fig.colorbar(image2, ax=axes[2], location="bottom", fraction=0.05, pad=0.08, label="Difference")
    ratio_bar = fig.colorbar(image3, ax=axes[3], location="bottom", fraction=0.05, pad=0.08, label="Equivalent PGV ratio")
    ratio_ticks = np.log10(np.array([0.25, 0.5, 1.0, 2.0, 4.0]))
    ratio_bar.set_ticks(ratio_ticks)
    ratio_bar.set_ticklabels(["0.25", "0.5", "1", "2", "4"])
    fig.suptitle("2018 Northern Osaka earthquake: network-density counterfactual", fontsize=9.5, fontweight="bold")
    out = png_dir / "osaka_2018_current_vs_1994_gmpe_kriging_maps.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_summary(
    station_summary: pd.DataFrame,
    validation: pd.DataFrame,
    footprint: pd.DataFrame,
    png_dir: Path,
) -> Path:
    import matplotlib.pyplot as plt

    method = "gmpe_kriging"
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.05), constrained_layout=True)
    colors = {
        "current_2018_observed": "#1f77b4",
        "counterfactual_1994_active": "#d62728",
        "strict_1994_code_overlap": "#6b6b6b",
    }
    for _, row in station_summary.iterrows():
        scenario = row["scenario"]
        counts = [row[f"station_count_within_{int(radius)}km"] for radius in RADII_KM]
        axes[0].plot(RADII_KM, counts, marker="o", linewidth=1.7, color=colors.get(scenario, "#333333"), label=SCENARIO_LABELS.get(scenario, scenario))
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Epicentral radius (km)", fontsize=7.6)
    axes[0].set_ylabel("Station count", fontsize=7.6)
    axes[0].grid(True, color="#dddddd", linewidth=0.55)
    axes[0].legend(fontsize=6.6)
    axes[0].set_title("Network support")

    val = validation[(validation["method"] == method) & (validation["validation_subset"].isin(["all_I>=1", "I>=5.0", "I>=5.5"]))]
    scenarios = ["current_2018_observed", "counterfactual_1994_active", "strict_1994_code_overlap"]
    x = np.arange(len(scenarios))
    width = 0.24
    for offset, subset in zip([-width, 0.0, width], ["all_I>=1", "I>=5.0", "I>=5.5"]):
        values = []
        for scenario in scenarios:
            hit = val[(val["scenario"] == scenario) & (val["validation_subset"] == subset)]
            values.append(float(hit["class_accuracy"].iloc[0]) if not hit.empty else np.nan)
        axes[1].bar(x + offset, values, width=width, label=subset)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(["2018", "1994 geom.", "strict"], rotation=0)
    axes[1].set_ylim(0.0, 1.0)
    axes[1].set_ylabel("Class accuracy", fontsize=7.6)
    axes[1].grid(True, axis="y", color="#dddddd", linewidth=0.55)
    axes[1].legend(fontsize=6.6)
    axes[1].set_title("Validation at 2018 stations")

    fp = footprint[footprint["method"] == method].copy()
    threshold_cols = [("area_ge_5p0_km2", "I>=5.0"), ("area_ge_5p5_km2", "I>=5.5"), ("area_ge_6p0_km2", "I>=6.0")]
    x2 = np.arange(len(threshold_cols))
    width2 = 0.24
    for idx, scenario in enumerate(scenarios):
        sub = fp[fp["scenario"] == scenario]
        if sub.empty:
            continue
        values = [float(sub[col].iloc[0]) for col, _ in threshold_cols]
        axes[2].bar(x2 + (idx - 1) * width2, values, width=width2, color=colors.get(scenario, "#333333"), label=SCENARIO_LABELS.get(scenario, scenario))
    axes[2].set_xticks(x2)
    axes[2].set_xticklabels([label for _, label in threshold_cols])
    axes[2].set_ylabel("Footprint area (km$^2$)", fontsize=7.6)
    axes[2].grid(True, axis="y", color="#dddddd", linewidth=0.55)
    axes[2].legend(fontsize=6.4)
    axes[2].set_title("Estimated footprint")

    out = png_dir / "osaka_2018_counterfactual_summary.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def write_report(
    args: argparse.Namespace,
    event: EventContext,
    station_summary: pd.DataFrame,
    method_summary: pd.DataFrame,
    validation: pd.DataFrame,
    footprint: pd.DataFrame,
    map_png: Path,
    summary_png: Path,
) -> Path:
    def md_table(df: pd.DataFrame, cols: list[str], float_fmt: str = ".3f") -> str:
        sub = df[cols].copy()
        def fmt(value: object) -> str:
            if pd.isna(value):
                return ""
            if isinstance(value, (float, np.floating)):
                return f"{float(value):{float_fmt}}"
            if isinstance(value, (int, np.integer)):
                return str(int(value))
            return str(value)

        header = "| " + " | ".join(cols) + " |"
        separator = "| " + " | ".join(["---"] * len(cols)) + " |"
        rows = ["| " + " | ".join(fmt(value) for value in row) + " |" for row in sub.to_numpy(dtype=object)]
        return "\n".join([header, separator, *rows])

    station_cols = [
        "scenario",
        "n_station_used",
        "density_per_10000km2",
        "station_count_within_10km",
        "station_count_within_20km",
        "station_count_within_50km",
        "station_count_within_100km",
        "nearest_epicentral_distance_km",
        "station_nn_median_km",
        "max_observed_intensity",
    ]
    validation_focus = validation[
        (validation["method"] == args.map_method)
        & (validation["validation_subset"].isin(["all_I>=1", "I>=5.0", "I>=5.5", "I>=6.0"]))
    ].copy()
    validation_cols = [
        "scenario",
        "validation_subset",
        "n_validation",
        "class_accuracy",
        "class_within_1",
        "mae",
        "rmse",
        "bias",
        "class_under_rate",
        "class_over_rate",
    ]
    footprint_focus = footprint[footprint["method"] == args.map_method].copy()
    footprint_cols = [
        "scenario",
        "max_estimated_surface_intensity",
        "area_ge_5p0_km2",
        "area_ge_5p5_km2",
        "area_ge_6p0_km2",
        "area_ge_6p5_km2",
    ]
    out = args.csv_dir.parent.parent / "osaka_2018_network_counterfactual_report.md"
    map_rel = Path(os.path.relpath(map_png, out.parent))
    summary_rel = Path(os.path.relpath(summary_png, out.parent))
    lines = [
        "# 2018年大阪府北部地震に対する観測網密度の反実仮想分析",
        "",
        f"- 対象地震: `{event.event_id}` ({event.origin_time}, M{event.magnitude:.1f}, depth {event.depth_km:.1f} km)",
        f"- 比較範囲: 震央から {args.radius_km:g} km以内",
        f"- 反実仮想: {args.counterfactual_year}年末時点で稼働していた観測点位置に，2018年実観測の最近傍値を {args.match_radius_km:g} km以内で割り当てた。",
        "",
        "この分析は，1990年代相当の観測点密度・幾何が同じ地震の推計震度分布をどの程度変えるかを評価するための再サンプリングである。",
        "したがって，当時の観測点固有の設置条件・機器特性・地盤応答を完全再現するものではない。",
        "",
        "## 観測点支持",
        "",
        md_table(station_summary, station_cols),
        "",
        "## 推計震度分布の比較",
        "",
        f"![2018 current network versus 1994 counterfactual]({map_rel})",
        "",
        "震度は線形振幅ではないため，主比較には震度差ΔIを用いる。",
        "右端の等価PGV比は `I = 2.68 + 1.72 log10(PGV)` に基づく参考換算であり，震度値そのものの比ではない。",
        "",
        f"![Counterfactual summary]({summary_rel})",
        "",
        "## 2018年実観測点での検証",
        "",
        md_table(validation_focus, validation_cols),
        "",
        "## 面積指標",
        "",
        md_table(footprint_focus, footprint_cols),
        "",
        "## 原稿に入れるべき解釈",
        "",
        "2018年の実観測網では，震央100 km以内に多数の観測点があり，震源近傍から大阪平野・京都盆地・奈良盆地にかけての局所的な震度勾配を直接拘束できる。",
        f"一方，{args.counterfactual_year}年末の配置へ落とすと，近傍観測点数は大幅に減少し，推計震度分布は距離減衰式事前場と少数の残差観測に強く依存する。",
        "この差は最大震度域そのものだけでなく，震度5弱以上・5強以上の面積評価にも現れる。",
        "したがって，近年の地震被害域の面的把握を過去観測網へ外挿する場合，観測点密度の増加前後で同じ推計震度分布品質を仮定してはいけない。",
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def copy_for_manuscript(args: argparse.Namespace, map_png: Path, summary_png: Path) -> None:
    args.manuscript_assets_ja.mkdir(parents=True, exist_ok=True)
    args.manuscript_assets_en.mkdir(parents=True, exist_ok=True)
    targets = [
        (map_png, args.manuscript_assets_ja / "fig11_osaka_2018_counterfactual_maps.png"),
        (summary_png, args.manuscript_assets_ja / "fig12_osaka_2018_counterfactual_summary.png"),
        (map_png, args.manuscript_assets_en / "fig11_osaka_2018_counterfactual_maps.png"),
        (summary_png, args.manuscript_assets_en / "fig12_osaka_2018_counterfactual_summary.png"),
    ]
    for src, dst in targets:
        shutil.copy2(src, dst)


def main() -> None:
    args = parse_args()
    configure_matplotlib()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)
    methods = [token.strip() for token in args.methods.split(",") if token.strip()]
    event = load_event(args)
    est_args = estimator_args(args)

    current = load_current_stations(args, event)
    current_all = load_event_station_observations(event, args.data_dir).dropna(
        subset=["latitude", "longitude", "observed_intensity"]
    )
    counterfactual, matches = build_counterfactual_stations(args, event, current_all)
    strict = build_strict_overlap(current, args)
    scenarios = {
        "current_2018_observed": current,
        "counterfactual_1994_active": counterfactual,
    }
    if len(strict) >= 3:
        scenarios["strict_1994_code_overlap"] = strict

    region = region_from_radius(event, args.map_radius_km)
    ground_all = load_ground_grid(args.ground_grid)
    ground = coarsen_ground_grid(crop_ground_grid(ground_all, region), args.grid_spacing)
    validation = nearest_ground_values(current, ground_all, event.longitude, event.latitude)
    validation["site_intensity_delta"] = site_intensity_delta(validation["amplification_vs400"], args.amp_intensity_coef)
    validation["bedrock_intensity"] = validation["observed_intensity"] - validation["site_intensity_delta"]

    matches.to_csv(args.csv_dir / "counterfactual_station_matches.csv", index=False)
    station_rows: list[dict[str, object]] = []
    method_rows: list[dict[str, object]] = []
    footprint_rows: list[dict[str, object]] = []
    validation_metric_rows: list[dict[str, object]] = []
    prediction_rows: list[pd.DataFrame] = []
    grids: dict[tuple[str, str], pd.DataFrame] = {}
    stations_with_ground: dict[str, pd.DataFrame] = {}

    for scenario, raw_stations in scenarios.items():
        stations = nearest_ground_values(raw_stations, ground_all, event.longitude, event.latitude)
        stations["site_intensity_delta"] = site_intensity_delta(stations["amplification_vs400"], args.amp_intensity_coef)
        stations["bedrock_intensity"] = stations["observed_intensity"] - stations["site_intensity_delta"]
        stations_with_ground[scenario] = stations
        stations.to_csv(args.csv_dir / f"{scenario}_stations_with_ground.csv", index=False)
        station_rows.append(summarize_stations(scenario, stations, event, args.radius_km))
        for method in methods:
            method_ground = ground
            grid, summary = estimate_event_method(event, stations, method_ground, method, est_args)
            summary.update(
                {
                    "scenario": scenario,
                    "scenario_label": SCENARIO_LABELS.get(scenario, scenario),
                    "region_lon_min": region[0],
                    "region_lon_max": region[1],
                    "region_lat_min": region[2],
                    "region_lat_max": region[3],
                    "map_grid_spacing_deg": inferred_grid_spacing(method_ground),
                }
            )
            grid_path = args.csv_dir / f"{scenario}_estimated_intensity_{method}.csv"
            grid.to_csv(grid_path, index=False)
            summary["grid_csv"] = str(grid_path)
            method_rows.append(summary)
            footprint_rows.append(summarize_footprint(scenario, method, grid))
            rows, pred_df = validation_rows(scenario, method, grid, validation, event)
            validation_metric_rows.extend(rows)
            prediction_rows.append(pred_df)
            grids[(scenario, method)] = grid
            if not args.no_pygmt_maps and method == args.map_method and scenario in {
                "current_2018_observed",
                "counterfactual_1994_active",
            }:
                draw_map(
                    event,
                    stations,
                    grid,
                    f"{method}; {SCENARIO_LABELS.get(scenario, scenario)}",
                    region,
                    args.png_dir / f"{scenario}_estimated_intensity_{method}.png",
                    args.prefecture_boundary,
                )

    station_summary = pd.DataFrame(station_rows)
    method_summary = pd.DataFrame(method_rows)
    footprint = pd.DataFrame(footprint_rows)
    validation_metrics = pd.DataFrame(validation_metric_rows)
    predictions = pd.concat(prediction_rows, ignore_index=True)
    station_summary.to_csv(args.csv_dir / "scenario_station_summary.csv", index=False)
    method_summary.to_csv(args.csv_dir / "scenario_method_summary.csv", index=False)
    footprint.to_csv(args.csv_dir / "scenario_footprint_summary.csv", index=False)
    validation_metrics.to_csv(args.csv_dir / "scenario_validation_metrics.csv", index=False)
    predictions.to_csv(args.csv_dir / "scenario_validation_predictions.csv.gz", index=False, compression="gzip")

    map_png = plot_map_comparison(
        event,
        region,
        grids,
        stations_with_ground,
        args.map_method,
        args.png_dir,
        args.amp_intensity_coef,
    )
    summary_png = plot_summary(station_summary, validation_metrics, footprint, args.png_dir)
    copy_for_manuscript(args, map_png, summary_png)
    report = write_report(args, event, station_summary, method_summary, validation_metrics, footprint, map_png, summary_png)
    print(f"Saved CSV: {args.csv_dir}")
    print(f"Saved PNG: {args.png_dir}")
    print(f"Saved report: {report}")
    print(station_summary.to_string(index=False))


if __name__ == "__main__":
    main()
