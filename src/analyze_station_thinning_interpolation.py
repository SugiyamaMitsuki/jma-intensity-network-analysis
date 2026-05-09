#!/usr/bin/env python3
"""Analyze how station thinning affects estimated seismic intensity interpolation."""

from __future__ import annotations

import argparse
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

from estimate_jma_intensity_distribution import (
    DEFAULT_GROUND_GRID,
    DEFAULT_TARGET_CSV,
    EventContext,
    calibrate_reference_trend,
    crop_ground_grid,
    event_context,
    event_region,
    fit_source_trend,
    gmpe_event_bias,
    gmpe_reference_bedrock_intensity,
    hypocentral_distance_km,
    inferred_grid_spacing,
    interpolate_idw,
    interpolate_values,
    load_event_station_observations,
    load_ground_grid,
    nearest_ground_values,
    predict_source_trend,
    processing_mode,
    project_km,
    site_intensity_delta,
)


DEFAULT_CSV_DIR = Path("outputs/csv/station_thinning_interpolation")
DEFAULT_PNG_DIR = Path("outputs/png/station_thinning_interpolation")
INTENSITY_CLASS_LABELS = np.array(["0", "1", "2", "3", "4", "5-", "5+", "6-", "6+", "7"], dtype=object)
INTENSITY_CLASS_VALUES = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 6.5, 7.0], dtype=float)
INTENSITY_CLASS_BOUNDS = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5], dtype=float)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Randomly thin seismic intensity stations, predict held-out station "
            "intensity with multiple interpolation methods, and relate errors to "
            "station density and J-SHIS surface-ground complexity."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--target-events-csv", type=Path, default=DEFAULT_TARGET_CSV)
    parser.add_argument("--ground-grid", type=Path, default=DEFAULT_GROUND_GRID)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--event-id", default="all-targets")
    parser.add_argument("--methods", default="linear,idw,spline,kriging")
    parser.add_argument("--keep-fractions", default="0.1,0.2,0.3,0.5,0.7,0.9")
    parser.add_argument("--n-random", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260509)
    parser.add_argument("--min-station-intensity", type=float, default=1.0)
    parser.add_argument("--min-train-stations", type=int, default=8)
    parser.add_argument("--station-lat-min", type=float, default=20.0)
    parser.add_argument("--station-lat-max", type=float, default=50.0)
    parser.add_argument("--station-lon-min", type=float, default=120.0)
    parser.add_argument("--station-lon-max", type=float, default=155.0)
    parser.add_argument("--region-padding-deg", type=float, default=0.55)
    parser.add_argument("--amp-intensity-coef", type=float, default=1.72)
    parser.add_argument("--mode", choices=["auto", "observed", "source_residual"], default="auto")
    parser.add_argument("--shallow-depth-threshold-km", type=float, default=150.0)
    parser.add_argument("--idw-neighbors", type=int, default=16)
    parser.add_argument("--idw-power", type=float, default=2.0)
    parser.add_argument("--rbf-neighbors", type=int, default=80)
    parser.add_argument("--rbf-smoothing", type=float, default=0.08)
    parser.add_argument("--kriging-neighbors", type=int, default=10)
    parser.add_argument("--kriging-range-km", type=float, default=None)
    parser.add_argument("--kriging-nugget", type=float, default=0.02)
    parser.add_argument("--gmpe-fault-type", choices=["crustal", "interplate", "intraslab"], default="crustal")
    parser.add_argument("--gmpe-min-distance-km", type=float, default=3.0)
    parser.add_argument("--gmpe-bias-stat", choices=["median", "mean", "none"], default="median")
    parser.add_argument("--pgv-intensity-intercept", type=float, default=2.68)
    parser.add_argument("--allow-overshoot", action="store_true")
    parser.add_argument("--clip-min", type=float, default=0.0)
    parser.add_argument("--clip-max", type=float, default=7.2)
    parser.add_argument("--local-complexity-radius-km", type=float, default=20.0)
    parser.add_argument("--no-png", action="store_true")
    return parser.parse_args()


def parse_float_list(value: str) -> list[float]:
    out: list[float] = []
    for part in value.split(","):
        if part.strip():
            out.append(float(part))
    if not out:
        raise SystemExit("At least one numeric value is required.")
    return out


def parse_methods(value: str) -> list[str]:
    allowed = {
        "linear",
        "cubic",
        "spline",
        "idw",
        "kriging",
        "gmpe_raw",
        "gmpe_calibrated",
        "gmpe_kriging",
        "nearest",
    }
    methods = []
    for part in value.split(","):
        method = part.strip().lower()
        if not method:
            continue
        if method not in allowed:
            raise SystemExit(f"Unknown method {method}; allowed={sorted(allowed)}")
        if method not in methods:
            methods.append(method)
    if not methods:
        raise SystemExit("At least one interpolation method is required.")
    return methods


def selected_events(path: Path, event_id: str) -> list[EventContext]:
    df = pd.read_csv(path, low_memory=False)
    if event_id != "all-targets":
        df = df[df["event_id"] == event_id]
    if df.empty:
        raise SystemExit(f"No target events found for {event_id}")
    return [event_context(row) for _, row in df.iterrows()]


def approximate_ground_area_km2(ground: pd.DataFrame) -> float:
    if ground.empty:
        return float("nan")
    spacing = inferred_grid_spacing(ground)
    lat = ground["latitude"].to_numpy(dtype=float)
    cell_area = (spacing * 111.32 * np.cos(np.radians(lat))) * (spacing * 110.57)
    return float(np.nansum(np.maximum(cell_area, 0.0)))


def q_iqr(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return float("nan")
    return float(np.nanquantile(finite, 0.75) - np.nanquantile(finite, 0.25))


def event_ground_complexity(ground: pd.DataFrame, amp_coef: float) -> dict[str, float]:
    site_delta = site_intensity_delta(ground["amplification_vs400"], amp_coef)
    return {
        "ground_area_km2": approximate_ground_area_km2(ground),
        "event_amp_mean": float(np.nanmean(ground["amplification_vs400"])),
        "event_amp_std": float(np.nanstd(ground["amplification_vs400"])),
        "event_amp_iqr": q_iqr(ground["amplification_vs400"].to_numpy(dtype=float)),
        "event_avs30_mean": float(np.nanmean(ground["avs30_m_s"])),
        "event_avs30_std": float(np.nanstd(ground["avs30_m_s"])),
        "event_avs30_iqr": q_iqr(ground["avs30_m_s"].to_numpy(dtype=float)),
        "event_site_delta_std": float(np.nanstd(site_delta)),
        "event_site_delta_iqr": q_iqr(site_delta),
    }


def valid_station_coordinate_mask(stations: pd.DataFrame, args: argparse.Namespace) -> pd.Series:
    """Reject aggregate or malformed station records that cannot be spatially validated."""
    lat = stations["latitude"].astype(float)
    lon = stations["longitude"].astype(float)
    return (
        lat.between(args.station_lat_min, args.station_lat_max)
        & lon.between(args.station_lon_min, args.station_lon_max)
        & np.isfinite(lat)
        & np.isfinite(lon)
        & ~((lat == 0.0) & (lon == 0.0))
    )


def add_local_ground_complexity(
    stations: pd.DataFrame,
    ground: pd.DataFrame,
    event: EventContext,
    radius_km: float,
    amp_coef: float,
) -> pd.DataFrame:
    out = stations.copy()
    ground_xy = project_km(ground["longitude"], ground["latitude"], event.longitude, event.latitude)
    station_xy = project_km(out["longitude"], out["latitude"], event.longitude, event.latitude)
    tree = cKDTree(ground_xy)
    amp = ground["amplification_vs400"].to_numpy(dtype=float)
    avs = ground["avs30_m_s"].to_numpy(dtype=float)
    delta = site_intensity_delta(amp, amp_coef)
    amp_std: list[float] = []
    amp_iqr: list[float] = []
    avs_std: list[float] = []
    avs_iqr: list[float] = []
    delta_std: list[float] = []
    n_cells: list[int] = []
    for xy in station_xy:
        idx = tree.query_ball_point(xy, radius_km)
        if not idx:
            amp_std.append(float("nan"))
            amp_iqr.append(float("nan"))
            avs_std.append(float("nan"))
            avs_iqr.append(float("nan"))
            delta_std.append(float("nan"))
            n_cells.append(0)
            continue
        idx_arr = np.asarray(idx, dtype=int)
        amp_vals = amp[idx_arr]
        avs_vals = avs[idx_arr]
        delta_vals = delta[idx_arr]
        amp_std.append(float(np.nanstd(amp_vals)))
        amp_iqr.append(q_iqr(amp_vals))
        avs_std.append(float(np.nanstd(avs_vals)))
        avs_iqr.append(q_iqr(avs_vals))
        delta_std.append(float(np.nanstd(delta_vals)))
        n_cells.append(int(len(idx_arr)))
    out[f"local_amp_std_{int(radius_km)}km"] = amp_std
    out[f"local_amp_iqr_{int(radius_km)}km"] = amp_iqr
    out[f"local_avs30_std_{int(radius_km)}km"] = avs_std
    out[f"local_avs30_iqr_{int(radius_km)}km"] = avs_iqr
    out[f"local_site_delta_std_{int(radius_km)}km"] = delta_std
    out[f"local_ground_cell_count_{int(radius_km)}km"] = n_cells
    return out


def nearest_neighbor_summary(train_xy: np.ndarray) -> tuple[float, float]:
    if len(train_xy) < 2:
        return float("nan"), float("nan")
    dist, _ = cKDTree(train_xy).query(train_xy, k=2)
    nn = dist[:, 1]
    return float(np.nanmean(nn)), float(np.nanmedian(nn))


def predict_holdout(
    method: str,
    event: EventContext,
    train: pd.DataFrame,
    test: pd.DataFrame,
    args: argparse.Namespace,
) -> np.ndarray:
    train_xy = project_km(train["longitude"], train["latitude"], event.longitude, event.latitude)
    test_xy = project_km(test["longitude"], test["latitude"], event.longitude, event.latitude)
    train_values = train["bedrock_intensity"].to_numpy(dtype=float)
    mode = processing_mode(args, event)
    if method == "gmpe_raw":
        bedrock_pred = gmpe_reference_bedrock_intensity(event, test_xy, args)
    elif method == "gmpe_calibrated":
        reference_train = gmpe_reference_bedrock_intensity(event, train_xy, args)
        reference_test = gmpe_reference_bedrock_intensity(event, test_xy, args)
        intercept, slope, fitted_train = calibrate_reference_trend(train_values, reference_train)
        bias = gmpe_event_bias(train_values, fitted_train, args)
        bedrock_pred = intercept + slope * reference_test + bias
    elif method == "gmpe_kriging":
        reference_train = gmpe_reference_bedrock_intensity(event, train_xy, args)
        reference_test = gmpe_reference_bedrock_intensity(event, test_xy, args)
        intercept, slope, fitted_train = calibrate_reference_trend(train_values, reference_train)
        fitted_test = intercept + slope * reference_test
        bias = gmpe_event_bias(train_values, fitted_train, args)
        residual_train = train_values - fitted_train - bias
        try:
            residual_test = interpolate_values(method, train_xy, residual_train, test_xy, args)
        except Exception:
            residual_test = interpolate_idw(
                train_xy,
                residual_train,
                test_xy,
                neighbors=min(args.idw_neighbors, len(train_values)),
                power=args.idw_power,
            )
        if not args.allow_overshoot:
            lo, hi = np.nanquantile(residual_train, [0.005, 0.995])
            residual_test = np.clip(residual_test, lo, hi)
        bedrock_pred = fitted_test + bias + residual_test
    elif mode == "source_residual" and len(train) >= 6:
        trend = fit_source_trend(train_xy, train_values, event.depth_km)
        residual_train = train_values - np.asarray(trend["fitted"], dtype=float)
        try:
            residual_test = interpolate_values(method, train_xy, residual_train, test_xy, args)
        except Exception:
            residual_test = interpolate_idw(
                train_xy,
                residual_train,
                test_xy,
                neighbors=min(args.idw_neighbors, len(train_values)),
                power=args.idw_power,
            )
        if not args.allow_overshoot:
            lo, hi = np.nanquantile(residual_train, [0.005, 0.995])
            residual_test = np.clip(residual_test, lo, hi)
        bedrock_pred = predict_source_trend(
            test_xy,
            event.depth_km,
            np.asarray(trend["coef"], dtype=float),
        ) + residual_test
    else:
        try:
            bedrock_pred = interpolate_values(method, train_xy, train_values, test_xy, args)
        except Exception:
            bedrock_pred = interpolate_idw(
                train_xy,
                train_values,
                test_xy,
                neighbors=min(args.idw_neighbors, len(train_values)),
                power=args.idw_power,
            )
        if not args.allow_overshoot:
            lo, hi = np.nanquantile(train_values, [0.005, 0.995])
            bedrock_pred = np.clip(bedrock_pred, lo, hi)
    pred = bedrock_pred + test["site_intensity_delta"].to_numpy(dtype=float)
    return np.clip(pred, args.clip_min, args.clip_max)


def intensity_class_index(values: np.ndarray | pd.Series) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    out = np.searchsorted(INTENSITY_CLASS_BOUNDS, arr, side="right")
    out[~np.isfinite(arr)] = -1
    return out.astype(int)


def intensity_class_labels(class_index: np.ndarray | pd.Series) -> np.ndarray:
    idx = np.asarray(class_index, dtype=int)
    out = np.full(idx.shape, "", dtype=object)
    valid = (0 <= idx) & (idx < len(INTENSITY_CLASS_LABELS))
    out[valid] = INTENSITY_CLASS_LABELS[idx[valid]]
    return out


def intensity_class_values(class_index: np.ndarray | pd.Series) -> np.ndarray:
    idx = np.asarray(class_index, dtype=int)
    out = np.full(idx.shape, np.nan, dtype=float)
    valid = (0 <= idx) & (idx < len(INTENSITY_CLASS_VALUES))
    out[valid] = INTENSITY_CLASS_VALUES[idx[valid]]
    return out


def error_metrics(errors: np.ndarray, truth: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    abs_error = np.abs(errors)
    truth_class = intensity_class_index(truth)
    pred_class = intensity_class_index(predicted)
    class_diff = pred_class - truth_class
    valid_class = (truth_class >= 0) & (pred_class >= 0)
    valid_diff = class_diff[valid_class]
    return {
        "bias": float(np.nanmean(errors)),
        "mae": float(np.nanmean(abs_error)),
        "rmse": float(np.sqrt(np.nanmean(errors**2))),
        "median_abs_error": float(np.nanmedian(abs_error)),
        "p90_abs_error": float(np.nanquantile(abs_error, 0.90)),
        "within_0p25": float(np.nanmean(abs_error <= 0.25)),
        "within_0p50": float(np.nanmean(abs_error <= 0.50)),
        "within_1p00": float(np.nanmean(abs_error <= 1.00)),
        "class_accuracy": float(np.nanmean(valid_diff == 0)) if len(valid_diff) else float("nan"),
        "class_within_1": float(np.nanmean(np.abs(valid_diff) <= 1)) if len(valid_diff) else float("nan"),
        "class_over_rate": float(np.nanmean(valid_diff > 0)) if len(valid_diff) else float("nan"),
        "class_under_rate": float(np.nanmean(valid_diff < 0)) if len(valid_diff) else float("nan"),
        "class_mean_error": float(np.nanmean(valid_diff)) if len(valid_diff) else float("nan"),
        "class_mae": float(np.nanmean(np.abs(valid_diff))) if len(valid_diff) else float("nan"),
    }


def process_event(
    event: EventContext,
    methods: list[str],
    keep_fractions: list[float],
    rng: np.random.Generator,
    ground_all: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[list[dict[str, object]], list[pd.DataFrame], dict[str, object]]:
    stations = load_event_station_observations(event, args.data_dir)
    n_loaded = len(stations)
    valid_coord = valid_station_coordinate_mask(stations, args)
    n_invalid_coord = int((~valid_coord).sum())
    if n_invalid_coord:
        print(f"  Excluding {n_invalid_coord} stations outside coordinate bounds for {event.event_id}")
    stations = stations[valid_coord].copy()
    stations = stations[stations["observed_intensity"] >= args.min_station_intensity].copy()
    if len(stations) < args.min_train_stations + 3:
        raise SystemExit(f"{event.event_id}: too few stations after filtering: {len(stations)}")

    region = event_region(stations, event, args.region_padding_deg)
    ground_region = crop_ground_grid(ground_all, region)
    complexity = event_ground_complexity(ground_region, args.amp_intensity_coef)
    stations = nearest_ground_values(stations, ground_all, event.longitude, event.latitude)
    stations["site_intensity_delta"] = site_intensity_delta(stations["amplification_vs400"], args.amp_intensity_coef)
    stations["bedrock_intensity"] = stations["observed_intensity"] - stations["site_intensity_delta"]
    stations = add_local_ground_complexity(
        stations,
        ground_region,
        event,
        args.local_complexity_radius_km,
        args.amp_intensity_coef,
    )

    station_out = args.csv_dir / "stations" / f"{event.event_id}_station_observations_ground_complexity.csv"
    station_out.parent.mkdir(parents=True, exist_ok=True)
    stations.to_csv(station_out, index=False)

    n_all = len(stations)
    station_xy_all = project_km(stations["longitude"], stations["latitude"], event.longitude, event.latitude)
    all_mean_nn, all_median_nn = nearest_neighbor_summary(station_xy_all)
    summary_rows: list[dict[str, object]] = []
    prediction_frames: list[pd.DataFrame] = []
    ground_area = complexity["ground_area_km2"]

    for keep_fraction in keep_fractions:
        n_train = int(round(n_all * keep_fraction))
        n_train = max(args.min_train_stations, min(n_train, n_all - 3))
        if n_train >= n_all:
            continue
        for random_index in range(args.n_random):
            train_idx = np.sort(rng.choice(n_all, size=n_train, replace=False))
            train_mask = np.zeros(n_all, dtype=bool)
            train_mask[train_idx] = True
            train = stations.iloc[train_mask].copy()
            test = stations.iloc[~train_mask].copy()
            train_xy = project_km(train["longitude"], train["latitude"], event.longitude, event.latitude)
            train_mean_nn, train_median_nn = nearest_neighbor_summary(train_xy)
            density = n_train / ground_area * 10000.0 if ground_area and math.isfinite(ground_area) else float("nan")
            for method in methods:
                uses_observations = method != "gmpe_raw"
                effective_n_train = n_train if uses_observations else 0
                effective_density = density if uses_observations else 0.0
                effective_train_mean_nn = train_mean_nn if uses_observations else float("nan")
                effective_train_median_nn = train_median_nn if uses_observations else float("nan")
                pred = predict_holdout(method, event, train, test, args)
                truth = test["observed_intensity"].to_numpy(dtype=float)
                error = pred - truth
                metrics = error_metrics(error, truth, pred)
                base = {
                    "event_id": event.event_id,
                    "origin_time": event.origin_time,
                    "year": event.year,
                    "magnitude": event.magnitude,
                    "depth_km": event.depth_km,
                    "method": method,
                    "keep_fraction": keep_fraction,
                    "random_index": random_index,
                    "n_station_all": n_all,
                    "n_station_loaded": n_loaded,
                    "n_station_invalid_coordinates": n_invalid_coord,
                    "n_train": n_train,
                    "n_test": len(test),
                    "train_density_per_10000km2": density,
                    "uses_observations": uses_observations,
                    "effective_n_train": effective_n_train,
                    "effective_train_density_per_10000km2": effective_density,
                    "ground_area_km2": ground_area,
                    "all_mean_nn_distance_km": all_mean_nn,
                    "all_median_nn_distance_km": all_median_nn,
                    "train_mean_nn_distance_km": train_mean_nn,
                    "train_median_nn_distance_km": train_median_nn,
                    "effective_train_mean_nn_distance_km": effective_train_mean_nn,
                    "effective_train_median_nn_distance_km": effective_train_median_nn,
                    **complexity,
                }
                summary_rows.append({**base, **metrics})
                pred_df = test[
                    [
                        "station_code",
                        "station_name",
                        "latitude",
                        "longitude",
                        "observed_intensity",
                        "amplification_vs400",
                        "avs30_m_s",
                        f"local_amp_std_{int(args.local_complexity_radius_km)}km",
                        f"local_amp_iqr_{int(args.local_complexity_radius_km)}km",
                        f"local_avs30_std_{int(args.local_complexity_radius_km)}km",
                        f"local_site_delta_std_{int(args.local_complexity_radius_km)}km",
                    ]
                ].copy()
                pred_df["event_id"] = event.event_id
                pred_df["method"] = method
                pred_df["keep_fraction"] = keep_fraction
                pred_df["random_index"] = random_index
                pred_df["n_train"] = n_train
                pred_df["train_density_per_10000km2"] = density
                pred_df["uses_observations"] = uses_observations
                pred_df["effective_n_train"] = effective_n_train
                pred_df["effective_train_density_per_10000km2"] = effective_density
                pred_df["effective_train_mean_nn_distance_km"] = effective_train_mean_nn
                pred_df["effective_train_median_nn_distance_km"] = effective_train_median_nn
                pred_df["predicted_intensity"] = pred
                pred_df["prediction_error"] = error
                pred_df["abs_error"] = np.abs(error)
                observed_class = intensity_class_index(truth)
                predicted_class = intensity_class_index(pred)
                pred_df["observed_intensity_class_index"] = observed_class
                pred_df["predicted_intensity_class_index"] = predicted_class
                pred_df["observed_intensity_class"] = intensity_class_labels(observed_class)
                pred_df["predicted_intensity_class"] = intensity_class_labels(predicted_class)
                pred_df["observed_intensity_class_value"] = intensity_class_values(observed_class)
                pred_df["predicted_intensity_class_value"] = intensity_class_values(predicted_class)
                pred_df["class_error"] = predicted_class - observed_class
                pred_df["class_match"] = pred_df["class_error"] == 0
                prediction_frames.append(pred_df)

    event_complexity = {
        "event_id": event.event_id,
        "origin_time": event.origin_time,
        "year": event.year,
        "magnitude": event.magnitude,
        "depth_km": event.depth_km,
        "n_station_all": n_all,
        "n_station_loaded": n_loaded,
        "n_station_invalid_coordinates": n_invalid_coord,
        "station_csv": str(station_out),
        **complexity,
        "all_mean_nn_distance_km": all_mean_nn,
        "all_median_nn_distance_km": all_median_nn,
    }
    return summary_rows, prediction_frames, event_complexity


def density_bins(series: pd.Series) -> pd.Categorical:
    bins = [0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000, np.inf]
    return pd.cut(series, bins=bins, right=False)


def aggregate_outputs(summary: pd.DataFrame, predictions: pd.DataFrame, out_dir: Path) -> dict[str, Path]:
    out_paths: dict[str, Path] = {}
    binned = summary.copy()
    density_col = (
        "effective_train_density_per_10000km2"
        if "effective_train_density_per_10000km2" in binned.columns
        else "train_density_per_10000km2"
    )
    n_train_col = "effective_n_train" if "effective_n_train" in binned.columns else "n_train"
    train_nn_col = (
        "effective_train_median_nn_distance_km"
        if "effective_train_median_nn_distance_km" in binned.columns
        else "train_median_nn_distance_km"
    )
    binned["density_bin"] = density_bins(binned[density_col]).astype(str)
    density_summary = (
        binned.groupby(["method", "density_bin"], observed=False)
        .agg(
            n_trials=("rmse", "size"),
            median_density=(density_col, "median"),
            median_n_train=(n_train_col, "median"),
            median_rmse=("rmse", "median"),
            median_mae=("mae", "median"),
            median_within_0p50=("within_0p50", "median"),
            median_class_accuracy=("class_accuracy", "median"),
            median_class_within_1=("class_within_1", "median"),
            median_class_over_rate=("class_over_rate", "median"),
            median_class_under_rate=("class_under_rate", "median"),
            median_class_mae=("class_mae", "median"),
            median_train_nn_km=(train_nn_col, "median"),
        )
        .reset_index()
    )
    density_summary = density_summary[density_summary["n_trials"] > 0]
    path = out_dir / "thinning_density_bin_summary.csv"
    density_summary.to_csv(path, index=False)
    out_paths["density_bin_summary"] = path

    threshold_rows = []
    class_accuracy_targets = [0.70, 0.80, 0.90]
    for method, group in density_summary.sort_values("median_density").groupby("method"):
        for target in class_accuracy_targets:
            ok = group[group["median_class_accuracy"] >= target]
            criterion = f"median_class_accuracy>={target:.2f}"
            if ok.empty:
                threshold_rows.append(
                    {
                        "method": method,
                        "criterion": criterion,
                        "target_class_accuracy": target,
                        "sufficient_density_per_10000km2": np.nan,
                        "sufficient_median_train_nn_km": np.nan,
                        "sufficient_class_accuracy": np.nan,
                        "status": "not_reached",
                    }
                )
            else:
                first = ok.iloc[0]
                threshold_rows.append(
                    {
                        "method": method,
                        "criterion": criterion,
                        "target_class_accuracy": target,
                        "sufficient_density_per_10000km2": first["median_density"],
                        "sufficient_median_train_nn_km": first["median_train_nn_km"],
                        "sufficient_class_accuracy": first["median_class_accuracy"],
                        "status": "reached",
                    }
                )
    thresholds = pd.DataFrame(threshold_rows)
    path = out_dir / "thinning_sufficient_density_thresholds.csv"
    thresholds.to_csv(path, index=False)
    out_paths["thresholds"] = path

    radius = [col for col in predictions.columns if col.startswith("local_site_delta_std_")]
    if radius:
        complexity_col = radius[0]
        pred = predictions.copy()
        if "class_match" not in pred.columns:
            pred["class_match"] = (
                intensity_class_index(pred["predicted_intensity"])
                == intensity_class_index(pred["observed_intensity"])
            )
        pred["complexity_quantile_bin"] = pd.qcut(
            pred[complexity_col].rank(method="first"),
            q=5,
            labels=["very_low", "low", "middle", "high", "very_high"],
        )
        complexity_summary = (
            pred.groupby(["method", "complexity_quantile_bin"], observed=True)
            .agg(
                n_predictions=("abs_error", "size"),
                median_local_site_delta_std=(complexity_col, "median"),
                mae=("abs_error", "mean"),
                rmse=("prediction_error", lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))),
                p90_abs_error=("abs_error", lambda x: float(np.quantile(x, 0.90))),
                class_accuracy=("class_match", "mean"),
                class_mae=("class_error", lambda x: float(np.mean(np.abs(np.asarray(x))))),
            )
            .reset_index()
        )
        path = out_dir / "thinning_error_by_local_ground_complexity.csv"
        complexity_summary.to_csv(path, index=False)
        out_paths["local_complexity_summary"] = path

    corr_rows = []
    metrics = ["event_amp_std", "event_amp_iqr", "event_avs30_std", "event_avs30_iqr", "event_site_delta_std"]
    event_method = summary.groupby(["event_id", "method"], as_index=False).agg(
        rmse=("rmse", "median"),
        mae=("mae", "median"),
        class_accuracy=("class_accuracy", "median"),
        **{metric: (metric, "first") for metric in metrics},
    )
    for method, group in event_method.groupby("method"):
        for metric in metrics:
            valid = group[[metric, "rmse", "mae", "class_accuracy"]].dropna()
            if len(valid) < 4:
                continue
            for err_col in ["rmse", "mae", "class_accuracy"]:
                rho, p = stats.spearmanr(valid[metric], valid[err_col])
                corr_rows.append(
                    {
                        "method": method,
                        "complexity_metric": metric,
                        "error_metric": err_col,
                        "spearman_rho": rho,
                        "p_value": p,
                        "n_events": len(valid),
                    }
                )
    corr = pd.DataFrame(corr_rows)
    path = out_dir / "event_ground_complexity_error_correlations.csv"
    corr.to_csv(path, index=False)
    out_paths["complexity_correlations"] = path
    return out_paths


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
            "grid.color": "#d9d9d9",
            "grid.linewidth": 0.55,
            "grid.alpha": 0.85,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
        }
    )


def savefig(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(fig)


def format_effective_density_axis(ax, density_summary: pd.DataFrame) -> None:
    max_density = float(np.nanmax(density_summary["median_density"])) if len(density_summary) else 1.0
    upper = max(1.0, max_density * 1.18)
    ticks = [0, 1, 3, 10, 30, 100, 300, 1000]
    ticks = [tick for tick in ticks if tick <= upper]
    rounded_upper = float(math.ceil(max_density / 10.0) * 10.0) if max_density > 10 else max_density
    if rounded_upper > 0 and all(abs(rounded_upper - tick) > 1e-9 for tick in ticks):
        ticks.append(rounded_upper)
    ticks = sorted(ticks)
    ax.set_xscale("symlog", linthresh=1.0, linscale=0.75)
    ax.set_xlim(0, upper)
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0" if tick == 0 else f"{tick:g}" for tick in ticks])


def plot_outputs(summary: pd.DataFrame, predictions: pd.DataFrame, density_summary: pd.DataFrame, thresholds: pd.DataFrame, png_dir: Path) -> None:
    import matplotlib.pyplot as plt

    setup_plot_style()
    colors = {
        "linear": "#0072B2",
        "idw": "#D55E00",
        "spline": "#009E73",
        "kriging": "#CC79A7",
        "gmpe_raw": "#666666",
        "gmpe_calibrated": "#A6761D",
        "gmpe_kriging": "#7A3B00",
        "cubic": "#E69F00",
        "nearest": "#666666",
    }

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for method, group in density_summary.groupby("method"):
        group = group.sort_values("median_density")
        ax.plot(
            group["median_density"],
            group["median_rmse"],
            marker="o",
            linewidth=1.8,
            color=colors.get(method, None),
            label=method,
        )
    ax.axhline(0.5, color="#333333", linestyle="--", linewidth=1.0)
    format_effective_density_axis(ax, density_summary)
    ax.set_xlabel("Effective training station density (/10,000 km$^2$)")
    ax.set_ylabel("Median RMSE at held-out stations")
    ax.set_title("Interpolation error decreases with station density")
    ax.legend(ncol=2)
    savefig(fig, png_dir / "thinning_rmse_vs_station_density.png")

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for method, group in density_summary.groupby("method"):
        group = group.sort_values("median_density")
        ax.plot(
            group["median_density"],
            group["median_within_0p50"],
            marker="o",
            linewidth=1.8,
            color=colors.get(method, None),
            label=method,
        )
    for target, alpha in [(0.70, 0.55), (0.80, 0.75), (0.90, 0.95)]:
        ax.axhline(target, color="#333333", linestyle="--", linewidth=0.9, alpha=alpha)
    format_effective_density_axis(ax, density_summary)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Effective training station density (/10,000 km$^2$)")
    ax.set_ylabel("Fraction within ±0.5 intensity")
    ax.set_title("Held-out accuracy versus station density")
    ax.legend(ncol=2)
    savefig(fig, png_dir / "thinning_within05_vs_station_density.png")

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for method, group in density_summary.groupby("method"):
        group = group.sort_values("median_density")
        ax.plot(
            group["median_density"],
            group["median_class_accuracy"],
            marker="o",
            linewidth=1.8,
            color=colors.get(method, None),
            label=method,
        )
    for target, alpha in [(0.70, 0.55), (0.80, 0.75), (0.90, 0.95)]:
        ax.axhline(target, color="#333333", linestyle="--", linewidth=0.9, alpha=alpha)
    format_effective_density_axis(ax, density_summary)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Effective training station density (/10,000 km$^2$)")
    ax.set_ylabel("Exact seismic intensity class hit rate")
    ax.set_title("Held-out class accuracy versus station density")
    ax.legend(ncol=2)
    savefig(fig, png_dir / "thinning_class_accuracy_vs_station_density.png")

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    reached = thresholds[thresholds["status"] == "reached"].copy()
    if not reached.empty:
        target_offsets = {0.70: -0.24, 0.80: 0.0, 0.90: 0.24}
        target_markers = {0.70: "o", 0.80: "s", 0.90: "^"}
        methods = sorted(reached["method"].unique())
        method_pos = {method: idx for idx, method in enumerate(methods)}
        for target, group in reached.groupby("target_class_accuracy"):
            y = [method_pos[m] + target_offsets.get(float(target), 0.0) for m in group["method"]]
            ax.scatter(
                group["sufficient_density_per_10000km2"],
                y,
                s=54,
                marker=target_markers.get(float(target), "o"),
                color=[colors.get(m, "#888888") for m in group["method"]],
                edgecolor="#222222",
                linewidth=0.45,
                label=f"{int(round(float(target) * 100))}%",
            )
        ax.set_yticks(list(method_pos.values()))
        ax.set_yticklabels(methods)
        ax.legend(title="Class hit target")
        ax.set_xscale("symlog", linthresh=1.0, linscale=0.75)
    ax.set_xlabel("Sufficient density (/10,000 km$^2$)")
    ax.set_title("Density meeting exact class hit-rate targets")
    savefig(fig, png_dir / "thinning_sufficient_density_thresholds.png")

    radius_cols = [col for col in predictions.columns if col.startswith("local_site_delta_std_")]
    if radius_cols:
        complexity_col = radius_cols[0]
        pred = predictions.copy()
        pred["complexity_bin"] = pd.qcut(
            pred[complexity_col].rank(method="first"),
            q=5,
            labels=["very low", "low", "middle", "high", "very high"],
        )
        grouped = pred.groupby(["method", "complexity_bin"], observed=True)["abs_error"].mean().reset_index()
        fig, ax = plt.subplots(figsize=(7.2, 4.5))
        labels = ["very low", "low", "middle", "high", "very high"]
        x = np.arange(len(labels))
        for method, group in grouped.groupby("method"):
            y = group.set_index("complexity_bin").reindex(labels)["abs_error"]
            ax.plot(x, y, marker="o", linewidth=1.8, color=colors.get(method, None), label=method)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_xlabel(f"Local ground complexity ({complexity_col})")
        ax.set_ylabel("Mean absolute error")
        ax.set_title("Ground-amplification variability and interpolation error")
        ax.legend(ncol=2)
        savefig(fig, png_dir / "thinning_error_vs_local_ground_complexity.png")

    event_method = summary.groupby(["event_id", "method"], as_index=False).agg(
        rmse=("rmse", "median"),
        event_site_delta_std=("event_site_delta_std", "first"),
        n_station_all=("n_station_all", "first"),
    )
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    for method, group in event_method.groupby("method"):
        ax.scatter(
            group["event_site_delta_std"],
            group["rmse"],
            s=np.clip(group["n_station_all"] / 15.0, 18, 130),
            alpha=0.72,
            color=colors.get(method, None),
            label=method,
            edgecolor="white",
            linewidth=0.4,
        )
    ax.set_xlabel("Event-region std. of site intensity correction")
    ax.set_ylabel("Median trial RMSE")
    ax.set_title("Event-scale ground complexity versus interpolation error")
    ax.legend(ncol=2)
    savefig(fig, png_dir / "event_complexity_vs_rmse.png")


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)
    methods = parse_methods(args.methods)
    keep_fractions = parse_float_list(args.keep_fractions)
    events = selected_events(args.target_events_csv, args.event_id)
    ground = load_ground_grid(args.ground_grid)
    rng = np.random.default_rng(args.seed)

    all_summary_rows: list[dict[str, object]] = []
    prediction_frames: list[pd.DataFrame] = []
    event_complexity_rows: list[dict[str, object]] = []
    for event in events:
        print(f"Processing {event.event_id}: {len(methods)} methods, {len(keep_fractions)} keep fractions, {args.n_random} random draws")
        summary_rows, pred_frames, complexity = process_event(event, methods, keep_fractions, rng, ground, args)
        all_summary_rows.extend(summary_rows)
        prediction_frames.extend(pred_frames)
        event_complexity_rows.append(complexity)

    summary = pd.DataFrame(all_summary_rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    event_complexity = pd.DataFrame(event_complexity_rows)

    summary_path = args.csv_dir / "thinning_trial_summary.csv"
    pred_path = args.csv_dir / "thinning_prediction_errors.csv.gz"
    complexity_path = args.csv_dir / "event_ground_complexity_summary.csv"
    summary.to_csv(summary_path, index=False)
    predictions.to_csv(pred_path, index=False, compression="gzip")
    event_complexity.to_csv(complexity_path, index=False)
    extra = aggregate_outputs(summary, predictions, args.csv_dir)
    if not args.no_png:
        density_summary = pd.read_csv(extra["density_bin_summary"])
        thresholds = pd.read_csv(extra["thresholds"])
        plot_outputs(summary, predictions, density_summary, thresholds, args.png_dir)

    print(f"Saved trial summary: {summary_path}")
    print(f"Saved prediction errors: {pred_path}")
    print(f"Saved event complexity: {complexity_path}")
    for name, path in extra.items():
        print(f"Saved {name}: {path}")


if __name__ == "__main__":
    main()
