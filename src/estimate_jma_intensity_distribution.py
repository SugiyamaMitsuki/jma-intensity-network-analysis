#!/usr/bin/env python3
"""Estimate JMA-style seismic intensity distribution with multiple interpolators.

This script traces the public JMA estimated-intensity workflow:

1. Convert observed surface intensity at stations to engineering-bedrock
   intensity by removing the J-SHIS site amplification term.
2. Estimate a spatially continuous engineering-bedrock intensity field.
3. Add the J-SHIS site amplification term at each grid cell to obtain surface
   estimated intensity.

For shallow events, the default mode fits a simple hypocentral-distance trend
and interpolates station residuals, mirroring JMA's concept of combining
source-based estimates with observed differences. It is not a reproduction of
the full EEW prediction engine.
"""

from __future__ import annotations

import argparse
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import CloughTocher2DInterpolator, LinearNDInterpolator, NearestNDInterpolator, RBFInterpolator
from scipy.spatial import cKDTree, distance_matrix

from analyze_jma_intensity import (
    SOURCE_RECORD_TYPES,
    as_ascii,
    haversine_km,
    load_station_index,
    parse_intensity_code,
    parse_measured_intensity,
    parse_scaled_int,
)
from plot_station_maps_pygmt import configure_conda_gmt
from map_prefecture_boundaries import DEFAULT_PREFECTURE_BOUNDARY, plot_prefecture_boundaries_pygmt


DEFAULT_EVENTS_CSV = Path("outputs/csv/hypocenter_catalog/jma_intensity_events_with_hypocenter.csv")
DEFAULT_TARGET_CSV = Path("outputs/csv/hypocenter_catalog/target_events_intensity_6upper_plus_with_hypocenter.csv")
DEFAULT_GROUND_GRID = Path("outputs/csv/jshis_surface_ground/jshis_surface_ground_grid_0p02deg.csv")
DEFAULT_CSV_DIR = Path("outputs/csv/estimated_intensity_distribution")
DEFAULT_PNG_DIR = Path("outputs/png/estimated_intensity_distribution")

METHOD_ALIASES = {
    "linear": "linear",
    "nearest": "nearest",
    "idw": "idw",
    "spline": "spline",
    "spline_rbf": "spline",
    "rbf": "spline",
    "cubic": "cubic",
    "clough_tocher": "cubic",
    "kriging": "kriging",
    "ordinary_kriging": "kriging",
    "gmpe": "gmpe_raw",
    "gmpe_raw": "gmpe_raw",
    "attenuation": "gmpe_raw",
    "attenuation_only": "gmpe_raw",
    "gmpe_calibrated": "gmpe_calibrated",
    "calibrated_gmpe": "gmpe_calibrated",
    "gmpe_kriging": "gmpe_kriging",
    "reference_kriging": "gmpe_kriging",
    "residual_kriging": "gmpe_kriging",
    "gmpe_residual_kriging": "gmpe_kriging",
}


@dataclass(frozen=True)
class EventContext:
    event_id: str
    year: int
    time_key: int
    latitude: float
    longitude: float
    depth_km: float
    magnitude: float | None
    max_intensity_value: float | None
    origin_time: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create JMA-style estimated intensity maps using several spatial interpolation methods."
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--events-csv", type=Path, default=DEFAULT_EVENTS_CSV)
    parser.add_argument("--target-events-csv", type=Path, default=DEFAULT_TARGET_CSV)
    parser.add_argument("--ground-grid", type=Path, default=DEFAULT_GROUND_GRID)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument(
        "--event-id",
        default="i2016_000828",
        help="Event id to process, or 'all-targets' for all events with observed intensity >= 6 upper.",
    )
    parser.add_argument(
        "--methods",
        default="linear,cubic,spline,idw,kriging,nearest",
        help="Comma-separated methods: linear,cubic,spline,idw,kriging,gmpe_raw,gmpe_calibrated,gmpe_kriging,nearest.",
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "observed", "source_residual"],
        default="auto",
        help="auto uses source_residual for depth < 150 km and observed otherwise.",
    )
    parser.add_argument("--shallow-depth-threshold-km", type=float, default=150.0)
    parser.add_argument(
        "--grid-spacing",
        type=float,
        default=0.02,
        help="Target map grid spacing in degrees. Values coarser than the input J-SHIS grid are aggregated.",
    )
    parser.add_argument("--region-padding-deg", type=float, default=0.55)
    parser.add_argument("--min-station-intensity", type=float, default=1.0)
    parser.add_argument("--amp-intensity-coef", type=float, default=1.72)
    parser.add_argument("--idw-neighbors", type=int, default=16)
    parser.add_argument("--idw-power", type=float, default=2.0)
    parser.add_argument("--rbf-neighbors", type=int, default=80)
    parser.add_argument("--rbf-smoothing", type=float, default=0.08)
    parser.add_argument("--kriging-neighbors", type=int, default=18)
    parser.add_argument("--kriging-range-km", type=float, default=None)
    parser.add_argument("--kriging-nugget", type=float, default=0.02)
    parser.add_argument(
        "--gmpe-fault-type",
        choices=["crustal", "interplate", "intraslab"],
        default="crustal",
        help="Fault-type term for the Si and Midorikawa (1999) PGV attenuation reference field.",
    )
    parser.add_argument(
        "--gmpe-min-distance-km",
        type=float,
        default=3.0,
        help="Minimum proxy fault distance used by the GMPE reference field.",
    )
    parser.add_argument(
        "--gmpe-bias-stat",
        choices=["median", "mean", "none"],
        default="median",
        help="Event-term calibration applied to the GMPE reference before residual kriging.",
    )
    parser.add_argument("--pgv-intensity-intercept", type=float, default=2.68)
    parser.add_argument(
        "--kriging-max-grid-cells",
        type=int,
        default=15_000,
        help="Coarsen kriging output automatically when the event grid is larger than this.",
    )
    parser.add_argument("--clip-min", type=float, default=0.0)
    parser.add_argument("--clip-max", type=float, default=7.2)
    parser.add_argument(
        "--allow-overshoot",
        action="store_true",
        help="Allow interpolated bedrock/residual values to exceed the station value range.",
    )
    parser.add_argument("--prefecture-boundary", type=Path, default=DEFAULT_PREFECTURE_BOUNDARY)
    parser.add_argument("--no-png", action="store_true")
    return parser.parse_args()


def normalize_methods(value: str) -> list[str]:
    methods: list[str] = []
    for token in value.split(","):
        key = token.strip().lower()
        if not key:
            continue
        if key not in METHOD_ALIASES:
            raise SystemExit(f"Unknown interpolation method: {token}")
        method = METHOD_ALIASES[key]
        if method not in methods:
            methods.append(method)
    if not methods:
        raise SystemExit("At least one interpolation method is required.")
    return methods


def event_time_key(row: pd.Series) -> int:
    return int(
        f"{int(row['year']):04d}{int(row['month']):02d}{int(row['day']):02d}"
        f"{int(row['hour']):02d}{int(row['minute']):02d}"
    )


def event_context(row: pd.Series) -> EventContext:
    lat = row.get("analysis_latitude", row.get("latitude"))
    lon = row.get("analysis_longitude", row.get("longitude"))
    dep = row.get("analysis_depth_km", row.get("depth_km"))
    mag = row.get("analysis_magnitude", row.get("magnitude"))
    origin = row.get("intensity_origin_time")
    if pd.isna(origin):
        origin = (
            f"{int(row['year']):04d}-{int(row['month']):02d}-{int(row['day']):02d} "
            f"{int(row['hour']):02d}:{int(row['minute']):02d}:{float(row['second']):05.2f}"
        )
    return EventContext(
        event_id=str(row["event_id"]),
        year=int(row["year"]),
        time_key=event_time_key(row),
        latitude=float(lat),
        longitude=float(lon),
        depth_km=float(dep) if not pd.isna(dep) else 0.0,
        magnitude=float(mag) if not pd.isna(mag) else None,
        max_intensity_value=float(row["max_intensity_value"]) if not pd.isna(row.get("max_intensity_value")) else None,
        origin_time=str(origin),
    )


def selected_events(args: argparse.Namespace) -> list[EventContext]:
    if args.event_id == "all-targets":
        df = pd.read_csv(args.target_events_csv, low_memory=False)
    else:
        df = pd.read_csv(args.events_csv, low_memory=False)
        df = df[df["event_id"] == args.event_id]
    if df.empty:
        raise SystemExit(f"No event rows found for {args.event_id}")
    return [event_context(row) for _, row in df.iterrows()]


def load_event_station_observations(event: EventContext, data_dir: Path):
    station_index = load_station_index(data_dir)
    if station_index is None:
        raise SystemExit(f"Station index not found under {data_dir}")

    path = data_dir / f"i{event.year}.zip"
    if not path.exists():
        raise SystemExit(f"Intensity zip not found: {path}")

    rows: dict[str, dict[str, object]] = {}
    current_id: str | None = None
    current_has_observations = False
    event_counter = 0

    with zipfile.ZipFile(path) as zf:
        names = [name for name in zf.namelist() if not name.endswith("/")]
        if not names:
            return pd.DataFrame()
        data = zf.read(names[0])

    for raw_line in data.splitlines():
        if not raw_line:
            continue
        first = raw_line[:1]
        if first in SOURCE_RECORD_TYPES:
            if current_id is None or current_has_observations:
                event_counter += 1
                current_id = f"{path.stem}_{event_counter:06d}"
                current_has_observations = False
            continue

        current_has_observations = True
        if current_id != event.event_id:
            continue
        padded = raw_line.ljust(96, b" ")
        station_code = as_ascii(padded[0:7]).strip()
        if not station_code:
            continue
        intensity_code, intensity_value = parse_intensity_code(padded[18:19])
        measured_intensity = parse_measured_intensity(padded[20:22])
        observed = measured_intensity if measured_intensity is not None else intensity_value
        if observed is None:
            continue
        rec = station_index.get(station_code, event.time_key)
        if rec is None or rec.latitude is None or rec.longitude is None:
            continue
        pga_total = parse_scaled_int(padded[29:34], 10.0)
        distance = haversine_km(event.latitude, event.longitude, rec.latitude, rec.longitude)
        new_row = {
            "event_id": event.event_id,
            "station_code": station_code,
            "station_name": rec.name,
            "latitude": rec.latitude,
            "longitude": rec.longitude,
            "intensity_code": intensity_code,
            "intensity_class_value": intensity_value,
            "measured_intensity": measured_intensity,
            "observed_intensity": observed,
            "pga_total_gal": pga_total,
            "epicentral_distance_km": distance,
        }
        old_row = rows.get(station_code)
        if old_row is None or float(new_row["observed_intensity"]) > float(old_row["observed_intensity"]):
            rows[station_code] = new_row

    return pd.DataFrame(rows.values())


def project_km(lon: np.ndarray, lat: np.ndarray, lon0: float, lat0: float) -> np.ndarray:
    x = (np.asarray(lon, dtype=float) - lon0) * 111.32 * math.cos(math.radians(lat0))
    y = (np.asarray(lat, dtype=float) - lat0) * 110.57
    return np.column_stack([x, y])


def load_ground_grid(path: Path) -> pd.DataFrame:
    cols = ["longitude", "latitude", "avs30_m_s", "amplification_vs400", "n_250m_mesh"]
    df = pd.read_csv(path, usecols=lambda c: c in cols)
    required = {"longitude", "latitude", "avs30_m_s", "amplification_vs400"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Ground grid is missing columns: {sorted(missing)}")
    return df.dropna(subset=["longitude", "latitude", "avs30_m_s", "amplification_vs400"]).copy()


def event_region(stations: pd.DataFrame, event: EventContext, padding: float) -> list[float]:
    lon_values = list(stations["longitude"].astype(float)) + [event.longitude]
    lat_values = list(stations["latitude"].astype(float)) + [event.latitude]
    lon_min = max(122.0, min(lon_values) - padding)
    lon_max = min(146.5, max(lon_values) + padding)
    lat_min = max(24.0, min(lat_values) - padding)
    lat_max = min(46.5, max(lat_values) + padding)
    return [lon_min, lon_max, lat_min, lat_max]


def crop_ground_grid(ground: pd.DataFrame, region: list[float]) -> pd.DataFrame:
    lon_min, lon_max, lat_min, lat_max = region
    out = ground[
        ground["longitude"].between(lon_min, lon_max)
        & ground["latitude"].between(lat_min, lat_max)
    ].copy()
    if out.empty:
        raise SystemExit(f"No J-SHIS ground grid cells in region {region}")
    return out


def inferred_grid_spacing(ground: pd.DataFrame) -> float:
    lon_unique = np.sort(ground["longitude"].unique())
    lat_unique = np.sort(ground["latitude"].unique())
    diffs: list[float] = []
    if len(lon_unique) > 1:
        diffs.extend(np.diff(lon_unique)[:2000])
    if len(lat_unique) > 1:
        diffs.extend(np.diff(lat_unique)[:2000])
    diffs = [float(v) for v in diffs if v > 1e-8]
    return float(np.nanmedian(diffs)) if diffs else 0.02


def coarsen_ground_grid(ground: pd.DataFrame, spacing: float) -> pd.DataFrame:
    base_spacing = inferred_grid_spacing(ground)
    if spacing <= base_spacing * 1.25:
        return ground.copy()
    lon_min = math.floor(float(ground["longitude"].min()) / spacing) * spacing
    lat_min = math.floor(float(ground["latitude"].min()) / spacing) * spacing
    out = ground.copy()
    out["_lon_bin"] = np.floor((out["longitude"] - lon_min) / spacing).astype(int)
    out["_lat_bin"] = np.floor((out["latitude"] - lat_min) / spacing).astype(int)
    grouped = out.groupby(["_lon_bin", "_lat_bin"], as_index=False).agg(
        avs30_m_s=("avs30_m_s", "mean"),
        amplification_vs400=("amplification_vs400", "mean"),
        n_250m_mesh=("n_250m_mesh", "sum") if "n_250m_mesh" in out else ("avs30_m_s", "size"),
    )
    grouped["longitude"] = lon_min + (grouped["_lon_bin"] + 0.5) * spacing
    grouped["latitude"] = lat_min + (grouped["_lat_bin"] + 0.5) * spacing
    grouped = grouped.drop(columns=["_lon_bin", "_lat_bin"])
    grouped = grouped[["longitude", "latitude", "avs30_m_s", "amplification_vs400", "n_250m_mesh"]]
    return grouped


def coarsen_ground_to_limit(ground: pd.DataFrame, max_cells: int) -> pd.DataFrame:
    if len(ground) <= max_cells:
        return ground.copy()
    base_spacing = inferred_grid_spacing(ground)
    factor = math.ceil(math.sqrt(len(ground) / max(max_cells, 1)))
    return coarsen_ground_grid(ground, base_spacing * factor)


def nearest_ground_values(points: pd.DataFrame, ground: pd.DataFrame, lon0: float, lat0: float) -> pd.DataFrame:
    tree_xy = project_km(ground["longitude"].to_numpy(), ground["latitude"].to_numpy(), lon0, lat0)
    point_xy = project_km(points["longitude"].to_numpy(), points["latitude"].to_numpy(), lon0, lat0)
    dist, idx = cKDTree(tree_xy).query(point_xy, k=1)
    out = points.copy()
    out["ground_nearest_distance_km"] = dist
    out["avs30_m_s"] = ground["avs30_m_s"].to_numpy()[idx]
    out["amplification_vs400"] = ground["amplification_vs400"].to_numpy()[idx]
    return out


def site_intensity_delta(amplification: np.ndarray | pd.Series, coefficient: float) -> np.ndarray:
    amp = np.clip(np.asarray(amplification, dtype=float), 0.05, None)
    return coefficient * np.log10(amp)


def hypocentral_distance_km(xy: np.ndarray, depth_km: float) -> np.ndarray:
    return np.sqrt(np.sum(xy**2, axis=1) + depth_km**2)


def fit_source_trend(train_xy: np.ndarray, bedrock: np.ndarray, depth_km: float) -> dict[str, object]:
    r = hypocentral_distance_km(train_xy, depth_km)
    xmat = np.column_stack([np.ones_like(r), np.log10(r + 1.0), r / 100.0])
    lam = 0.05
    penalty = np.diag([0.0, lam, lam])
    coef = np.linalg.solve(xmat.T @ xmat + penalty, xmat.T @ bedrock)
    fitted = xmat @ coef
    return {"coef": coef, "fitted": fitted}


def predict_source_trend(xy: np.ndarray, depth_km: float, coef: np.ndarray) -> np.ndarray:
    r = hypocentral_distance_km(xy, depth_km)
    xmat = np.column_stack([np.ones_like(r), np.log10(r + 1.0), r / 100.0])
    return xmat @ coef


def rupture_length_utsu_km(magnitude: float | None) -> float:
    if magnitude is None or not np.isfinite(magnitude):
        return 0.0
    return float(10.0 ** (0.5 * float(magnitude) - 1.85))


def proxy_fault_distance_km(xy: np.ndarray, event: EventContext, min_distance_km: float = 3.0) -> np.ndarray:
    """Approximate fault distance from hypocentral distance when no fault plane is available."""

    hypo = hypocentral_distance_km(xy, event.depth_km)
    half_length = 0.5 * rupture_length_utsu_km(event.magnitude)
    return np.maximum(float(min_distance_km), hypo - half_length)


def si_midorikawa_fault_type_term(fault_type: str) -> float:
    if fault_type == "interplate":
        return -0.02
    if fault_type == "intraslab":
        return 0.12
    return 0.0


def si_midorikawa_1999_log10_pgv600(event: EventContext, xy: np.ndarray, args) -> np.ndarray:
    """Return log10(PGV600[cm/s]) using the Si and Midorikawa (1999) PGV form."""

    magnitude = event.magnitude if event.magnitude is not None and np.isfinite(event.magnitude) else 6.0
    depth = max(0.0, float(event.depth_km))
    fault_distance = proxy_fault_distance_km(xy, event, args.gmpe_min_distance_km)
    source_saturation = 0.0028 * 10.0 ** (0.5 * float(magnitude))
    fault_term = si_midorikawa_fault_type_term(args.gmpe_fault_type)
    return (
        0.58 * float(magnitude)
        + 0.0038 * depth
        + fault_term
        - 1.29
        - np.log10(fault_distance + source_saturation)
        - 0.002 * fault_distance
    )


def gmpe_reference_bedrock_intensity(event: EventContext, xy: np.ndarray, args) -> np.ndarray:
    log10_pgv = si_midorikawa_1999_log10_pgv600(event, xy, args)
    return args.pgv_intensity_intercept + args.amp_intensity_coef * log10_pgv


def gmpe_event_bias(observed_bedrock: np.ndarray, reference_bedrock: np.ndarray, args) -> float:
    residual = observed_bedrock - reference_bedrock
    finite = residual[np.isfinite(residual)]
    if len(finite) == 0 or args.gmpe_bias_stat == "none":
        return 0.0
    if args.gmpe_bias_stat == "mean":
        return float(np.nanmean(finite))
    return float(np.nanmedian(finite))


def calibrate_reference_trend(
    observed_bedrock: np.ndarray,
    reference_bedrock: np.ndarray,
    ridge: float = 0.08,
) -> tuple[float, float, np.ndarray]:
    """Fit a stable event-specific linear calibration of a reference intensity field."""

    observed = np.asarray(observed_bedrock, dtype=float)
    reference = np.asarray(reference_bedrock, dtype=float)
    valid = np.isfinite(observed) & np.isfinite(reference)
    if valid.sum() < 4 or float(np.nanstd(reference[valid])) < 1e-6:
        intercept = float(np.nanmedian(observed[valid] - reference[valid])) if valid.any() else 0.0
        slope = 1.0
    else:
        x = reference[valid]
        y = observed[valid]
        x_mean = float(np.nanmean(x))
        y_mean = float(np.nanmean(y))
        xc = x - x_mean
        yc = y - y_mean
        slope = float(np.sum(xc * yc) / (np.sum(xc * xc) + ridge * len(x)))
        slope = float(np.clip(slope, 0.15, 2.50))
        intercept = y_mean - slope * x_mean
    fitted = intercept + slope * reference
    return intercept, slope, fitted


def fill_with_nearest(train_xy: np.ndarray, train_values: np.ndarray, eval_xy: np.ndarray, values: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=float).copy()
    missing = ~np.isfinite(out)
    if missing.any():
        nearest = NearestNDInterpolator(train_xy, train_values)
        out[missing] = nearest(eval_xy[missing])
    return out


def interpolate_idw(
    train_xy: np.ndarray,
    train_values: np.ndarray,
    eval_xy: np.ndarray,
    neighbors: int,
    power: float,
    chunk_size: int = 20_000,
) -> np.ndarray:
    tree = cKDTree(train_xy)
    k = min(max(1, neighbors), len(train_values))
    out = np.empty(len(eval_xy), dtype=float)
    for start in range(0, len(eval_xy), chunk_size):
        stop = min(start + chunk_size, len(eval_xy))
        dist, idx = tree.query(eval_xy[start:stop], k=k)
        if k == 1:
            out[start:stop] = train_values[idx]
            continue
        exact = dist[:, 0] <= 1e-10
        weights = 1.0 / np.maximum(dist, 1e-6) ** power
        pred = np.sum(weights * train_values[idx], axis=1) / np.sum(weights, axis=1)
        pred[exact] = train_values[idx[exact, 0]]
        out[start:stop] = pred
    return out


def estimate_kriging_range(train_xy: np.ndarray) -> float:
    if len(train_xy) < 3:
        return 30.0
    dist, _ = cKDTree(train_xy).query(train_xy, k=2)
    median_nn = float(np.nanmedian(dist[:, 1]))
    extent = float(np.nanmax(np.ptp(train_xy, axis=0)))
    return max(10.0, min(extent / 2.5, median_nn * 10.0))


def variogram_exponential(h: np.ndarray, sill: float, range_km: float, nugget: float) -> np.ndarray:
    h = np.asarray(h, dtype=float)
    gamma = nugget + sill * (1.0 - np.exp(-h / max(range_km, 1e-6)))
    gamma[h <= 1e-10] = 0.0
    return gamma


def covariance_exponential(h: np.ndarray, sill: float, range_km: float) -> np.ndarray:
    h = np.asarray(h, dtype=float)
    return sill * np.exp(-h / max(range_km, 1e-6))


def kriging_nugget_variance(sill: float, nugget: float) -> float:
    if nugget < 1.0:
        return max(0.0, nugget * sill)
    return max(0.0, nugget)


def interpolate_ordinary_kriging(
    train_xy: np.ndarray,
    train_values: np.ndarray,
    eval_xy: np.ndarray,
    neighbors: int,
    range_km: float | None,
    nugget: float,
) -> np.ndarray:
    if len(train_values) < 3:
        return interpolate_idw(train_xy, train_values, eval_xy, neighbors=1, power=2.0)
    tree = cKDTree(train_xy)
    k = min(max(3, neighbors), len(train_values))
    sill = float(np.nanvar(train_values))
    if sill <= 1e-10:
        return np.full(len(eval_xy), float(np.nanmean(train_values)))
    vrange = float(range_km) if range_km is not None else estimate_kriging_range(train_xy)
    out = np.empty(len(eval_xy), dtype=float)

    dist_all, idx_all = tree.query(eval_xy, k=k)
    for i in range(len(eval_xy)):
        dist = np.atleast_1d(dist_all[i])
        idx = np.atleast_1d(idx_all[i])
        if dist[0] <= 1e-10:
            out[i] = train_values[idx[0]]
            continue
        local_xy = train_xy[idx]
        local_values = train_values[idx]
        dmat = distance_matrix(local_xy, local_xy)
        cov = covariance_exponential(dmat, sill=sill, range_km=vrange)
        cov.flat[:: k + 1] += kriging_nugget_variance(sill, nugget)
        amat = np.ones((k + 1, k + 1), dtype=float)
        amat[:k, :k] = cov
        amat[k, k] = 0.0
        rhs = np.ones(k + 1, dtype=float)
        rhs[:k] = covariance_exponential(dist, sill=sill, range_km=vrange)
        rhs[k] = 1.0
        try:
            sol = np.linalg.solve(amat, rhs)
            out[i] = float(np.dot(sol[:k], local_values))
        except np.linalg.LinAlgError:
            out[i] = interpolate_idw(local_xy, local_values, eval_xy[i : i + 1], neighbors=k, power=2.0)[0]
    return out


def interpolate_values(method: str, train_xy: np.ndarray, train_values: np.ndarray, eval_xy: np.ndarray, args) -> np.ndarray:
    if method == "nearest":
        return NearestNDInterpolator(train_xy, train_values)(eval_xy)
    if method == "linear":
        pred = LinearNDInterpolator(train_xy, train_values)(eval_xy)
        return fill_with_nearest(train_xy, train_values, eval_xy, pred)
    if method == "cubic":
        pred = CloughTocher2DInterpolator(train_xy, train_values)(eval_xy)
        return fill_with_nearest(train_xy, train_values, eval_xy, pred)
    if method == "spline":
        neighbors = min(max(3, args.rbf_neighbors), len(train_values))
        model = RBFInterpolator(
            train_xy,
            train_values,
            kernel="thin_plate_spline",
            smoothing=args.rbf_smoothing,
            neighbors=neighbors,
        )
        return model(eval_xy)
    if method == "idw":
        return interpolate_idw(train_xy, train_values, eval_xy, args.idw_neighbors, args.idw_power)
    if method == "kriging":
        return interpolate_ordinary_kriging(
            train_xy,
            train_values,
            eval_xy,
            args.kriging_neighbors,
            args.kriging_range_km,
            args.kriging_nugget,
        )
    if method == "gmpe_kriging":
        return interpolate_ordinary_kriging(
            train_xy,
            train_values,
            eval_xy,
            args.kriging_neighbors,
            args.kriging_range_km,
            args.kriging_nugget,
        )
    raise ValueError(method)


def processing_mode(args: argparse.Namespace, event: EventContext) -> str:
    if args.mode != "auto":
        return args.mode
    return "source_residual" if event.depth_km < args.shallow_depth_threshold_km else "observed"


def estimate_event_method(
    event: EventContext,
    stations: pd.DataFrame,
    ground: pd.DataFrame,
    method: str,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, dict[str, object]]:
    lon0 = event.longitude
    lat0 = event.latitude
    station_xy = project_km(stations["longitude"].to_numpy(), stations["latitude"].to_numpy(), lon0, lat0)
    grid_xy = project_km(ground["longitude"].to_numpy(), ground["latitude"].to_numpy(), lon0, lat0)

    stations = stations.copy()
    stations["site_intensity_delta"] = site_intensity_delta(stations["amplification_vs400"], args.amp_intensity_coef)
    stations["bedrock_intensity"] = stations["observed_intensity"] - stations["site_intensity_delta"]

    ground_delta = site_intensity_delta(ground["amplification_vs400"], args.amp_intensity_coef)
    bedrock = stations["bedrock_intensity"].to_numpy(dtype=float)
    mode = processing_mode(args, event)

    gmpe_bias = np.nan
    gmpe_intercept = np.nan
    gmpe_slope = np.nan
    if method == "gmpe_raw":
        grid_bedrock = gmpe_reference_bedrock_intensity(event, grid_xy, args)
        residual = np.full(len(grid_bedrock), np.nan)
        trend_coef = np.array([np.nan, np.nan, np.nan])
        gmpe_bias = 0.0
        gmpe_intercept = 0.0
        gmpe_slope = 1.0
    elif method == "gmpe_calibrated":
        reference_station = gmpe_reference_bedrock_intensity(event, station_xy, args)
        reference_grid = gmpe_reference_bedrock_intensity(event, grid_xy, args)
        gmpe_intercept, gmpe_slope, fitted_station = calibrate_reference_trend(bedrock, reference_station)
        gmpe_bias = gmpe_event_bias(bedrock, fitted_station, args)
        grid_bedrock = gmpe_intercept + gmpe_slope * reference_grid + gmpe_bias
        residual = np.full(len(grid_bedrock), np.nan)
        trend_coef = np.array([np.nan, np.nan, np.nan])
    elif method == "gmpe_kriging":
        reference_station = gmpe_reference_bedrock_intensity(event, station_xy, args)
        reference_grid = gmpe_reference_bedrock_intensity(event, grid_xy, args)
        gmpe_intercept, gmpe_slope, fitted_station = calibrate_reference_trend(bedrock, reference_station)
        fitted_grid = gmpe_intercept + gmpe_slope * reference_grid
        gmpe_bias = gmpe_event_bias(bedrock, fitted_station, args)
        train_values = bedrock - fitted_station - gmpe_bias
        residual = interpolate_values(method, station_xy, train_values, grid_xy, args)
        if not args.allow_overshoot:
            lo, hi = np.nanquantile(train_values, [0.005, 0.995])
            residual = np.clip(residual, lo, hi)
        grid_bedrock = fitted_grid + gmpe_bias + residual
        trend_coef = np.array([np.nan, np.nan, np.nan])
    elif mode == "source_residual" and len(stations) >= 6:
        trend = fit_source_trend(station_xy, bedrock, event.depth_km)
        train_values = bedrock - np.asarray(trend["fitted"], dtype=float)
        residual = interpolate_values(method, station_xy, train_values, grid_xy, args)
        if method in {"linear", "cubic"}:
            residual = np.where(np.isfinite(residual), residual, 0.0)
        if not args.allow_overshoot:
            lo, hi = np.nanquantile(train_values, [0.005, 0.995])
            residual = np.clip(residual, lo, hi)
        grid_bedrock = predict_source_trend(grid_xy, event.depth_km, np.asarray(trend["coef"], dtype=float)) + residual
        trend_coef = np.asarray(trend["coef"], dtype=float)
    else:
        train_values = bedrock
        grid_bedrock = interpolate_values(method, station_xy, train_values, grid_xy, args)
        if not args.allow_overshoot:
            lo, hi = np.nanquantile(train_values, [0.005, 0.995])
            grid_bedrock = np.clip(grid_bedrock, lo, hi)
        residual = np.full(len(grid_bedrock), np.nan)
        trend_coef = np.array([np.nan, np.nan, np.nan])

    surface = np.clip(grid_bedrock + ground_delta, args.clip_min, args.clip_max)
    out = ground.copy()
    out["estimated_bedrock_intensity"] = grid_bedrock
    out["site_intensity_delta"] = ground_delta
    out["estimated_surface_intensity"] = surface
    out["interpolated_residual"] = residual
    out["method"] = method
    out["mode"] = mode
    out["event_id"] = event.event_id
    summary = {
        "event_id": event.event_id,
        "method": method,
        "mode": mode,
        "n_station_used": int(len(stations)),
        "n_grid_cells": int(len(out)),
        "min_estimated_surface_intensity": float(np.nanmin(surface)),
        "mean_estimated_surface_intensity": float(np.nanmean(surface)),
        "max_estimated_surface_intensity": float(np.nanmax(surface)),
        "trend_intercept": float(trend_coef[0]),
        "trend_log10_r_coef": float(trend_coef[1]),
        "trend_r100_coef": float(trend_coef[2]),
        "gmpe_bias": float(gmpe_bias) if np.isfinite(gmpe_bias) else np.nan,
        "gmpe_intercept": float(gmpe_intercept) if np.isfinite(gmpe_intercept) else np.nan,
        "gmpe_slope": float(gmpe_slope) if np.isfinite(gmpe_slope) else np.nan,
        "gmpe_fault_type": args.gmpe_fault_type if method.startswith("gmpe") else "",
        "gmpe_min_distance_km": args.gmpe_min_distance_km if method.startswith("gmpe") else np.nan,
    }
    return out, summary


def grid_to_dataarray(df: pd.DataFrame, value_col: str) -> xr.DataArray:
    def regular_axis(values: pd.Series) -> tuple[np.ndarray, float]:
        unique = np.sort(values.astype(float).unique())
        diffs = np.diff(unique)
        diffs = diffs[diffs > 1e-8]
        spacing = float(np.nanmin(diffs)) if len(diffs) else 0.02
        n = int(round((float(unique[-1]) - float(unique[0])) / spacing)) + 1
        axis = float(unique[0]) + np.arange(n) * spacing
        return axis, spacing

    lons, lon_spacing = regular_axis(df["longitude"])
    lats, lat_spacing = regular_axis(df["latitude"])
    data = np.full((len(lats), len(lons)), np.nan, dtype=np.float32)
    for lon, lat, value in zip(df["longitude"], df["latitude"], df[value_col]):
        ix = int(round((float(lon) - float(lons[0])) / lon_spacing))
        iy = int(round((float(lat) - float(lats[0])) / lat_spacing))
        if 0 <= ix < len(lons) and 0 <= iy < len(lats):
            data[iy, ix] = value
    return xr.DataArray(data, coords={"lat": lats, "lon": lons}, dims=("lat", "lon"), name=value_col)


def draw_map(
    event: EventContext,
    stations: pd.DataFrame,
    grid: pd.DataFrame,
    method: str,
    region: list[float],
    output_path: Path,
    prefecture_boundary: Path | None = DEFAULT_PREFECTURE_BOUNDARY,
) -> None:
    configure_conda_gmt()
    import pygmt

    da = grid_to_dataarray(grid, "estimated_surface_intensity")
    lon_span = region[1] - region[0]
    lat_span = region[3] - region[2]
    width = 13.5
    projection = f"M{width}c"
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="turbo", series=[0, 7, 0.5], continuous=True)
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN="0.7p,#222222",
        MAP_GRID_PEN_PRIMARY="0.12p,#d6d6d6",
        FONT_ANNOT_PRIMARY="7.5p,Helvetica,#222222",
        FONT_LABEL="8.5p,Helvetica,#222222",
        FONT_TITLE="11p,Helvetica-Bold,#111111",
        FORMAT_GEO_MAP="dddF",
        COLOR_NAN="white",
    ):
        fig.basemap(
            region=region,
            projection=projection,
            frame=[
                "xafg1+lLongitude" if lon_span <= 6 else "xafg2+lLongitude",
                "yafg1+lLatitude" if lat_span <= 6 else "yafg2+lLatitude",
                f"+tEstimated seismic intensity: {event.event_id} ({method})",
            ],
        )
        fig.coast(water="#edf4f7", land="#fbfbf8", shorelines="0.35p,#444444", resolution="i")
        fig.grdimage(grid=da, region=region, projection=projection, cmap=True, nan_transparent=True)
        fig.grdcontour(
            grid=da,
            region=region,
            projection=projection,
            levels=0.5,
            pen="0.18p,white@30",
        )
        fig.coast(shorelines="0.35p,#222222", borders="1/0.15p,#777777", resolution="i")
        plot_prefecture_boundaries_pygmt(fig, prefecture_boundary, pen="0.16p,#5f5f5f")
        subset = stations[stations["observed_intensity"] >= 3.0]
        fig.plot(
            x=subset["longitude"],
            y=subset["latitude"],
            style="c0.08c",
            fill="white",
            pen="0.2p,#111111",
        )
        fig.plot(x=[event.longitude], y=[event.latitude], style="a0.28c", fill="#111111", pen="0.25p,white")
        fig.text(
            x=region[0] + 0.04 * lon_span,
            y=region[3] - 0.06 * lat_span,
            text=f"{event.origin_time}  M{event.magnitude if event.magnitude is not None else 'NA'}  depth {event.depth_km:.1f} km",
            font="7p,Helvetica,#111111",
            justify="TL",
            fill="white@25",
            clearance="2p/2p",
        )
        fig.colorbar(position="JBC+w9.5c/0.28c+o0c/0.65c+h", frame=["x+lEstimated instrumental intensity", "y"])
    fig.savefig(output_path, crop=True, dpi=300)


def write_method_note(out_dir: Path) -> Path:
    path = out_dir / "method_note.md"
    path.write_text(
        "\n".join(
            [
                "# 気象庁方式を踏まえた推計震度分布の実装メモ",
                "",
                "このプログラムは、気象庁が公開している推計震度分布図の計算概念をトレースする研究用実装です。",
                "緊急地震速報の運用システムや係数を完全再現するものではありません。",
                "",
                "- 観測点の地表震度からJ-SHISの地盤増幅率を引き、工学的基盤上の震度に戻します:",
                "  `I_bedrock = I_surface - c log10(ARV)`。",
                "- `ARV` は Vs=400 m/s 工学的基盤から地表までの地盤増幅率で、`c` の既定値は 1.72 です。",
                "- 浅い地震では、基盤震度に対して震源距離トレンドを当て、観測点の差分を空間補間します。",
                "- 深い地震、または `--mode observed` では、観測点の基盤震度を直接補間します。",
                "- 格子点の地表推計震度は `I_surface_grid = I_bedrock_grid + c log10(ARV_grid)` です。",
                "- 実装済み補間法は、線形Delaunay補間、Clough-Tocher三次補間、thin-plate spline RBF、",
                "  IDW、局所通常クリギング、GMPEのみ、GMPE校正のみ、GMPE参照場付き残差クリギング、",
                "  最近傍補間です。",
                "- `gmpe_raw` は司・翠川(1999)型のPGV距離減衰式と `I = 2.68 + 1.72 log10(PGV)`",
                "  から基準場だけを作ります。",
                "- `gmpe_calibrated` は観測点で `I_bedrock = a + b I_GMPE + residual` の a,b を校正しますが、",
                "  residual は補間しません。",
                "- `gmpe_kriging` は校正後の観測点残差を通常クリギングします。",
                "  断層面が未設定のため、宇津式 `log10(L)=0.5M-1.85` で断層長を近似し、",
                "  震源距離から L/2 を差し引いた距離を断層最短距離の代理値にしています。",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def process_event(event: EventContext, methods: list[str], ground_all: pd.DataFrame, args: argparse.Namespace) -> list[dict[str, object]]:
    stations = load_event_station_observations(event, args.data_dir)
    stations = stations[stations["observed_intensity"] >= args.min_station_intensity].copy()
    if len(stations) < 3:
        raise SystemExit(f"{event.event_id}: need at least 3 stations, got {len(stations)}")

    region = event_region(stations, event, args.region_padding_deg)
    ground = coarsen_ground_grid(crop_ground_grid(ground_all, region), args.grid_spacing)
    stations = nearest_ground_values(stations, ground_all, event.longitude, event.latitude)
    stations["site_intensity_delta"] = site_intensity_delta(stations["amplification_vs400"], args.amp_intensity_coef)
    stations["bedrock_intensity"] = stations["observed_intensity"] - stations["site_intensity_delta"]

    event_csv_dir = args.csv_dir / event.event_id
    event_png_dir = args.png_dir / event.event_id
    event_csv_dir.mkdir(parents=True, exist_ok=True)
    event_png_dir.mkdir(parents=True, exist_ok=True)
    stations.to_csv(event_csv_dir / f"{event.event_id}_station_observations_with_ground.csv", index=False)

    summaries: list[dict[str, object]] = []
    for method in methods:
        method_ground = ground
        if method == "kriging":
            method_ground = coarsen_ground_to_limit(ground, args.kriging_max_grid_cells)
        grid, summary = estimate_event_method(event, stations, method_ground, method, args)
        summary["region_lon_min"] = region[0]
        summary["region_lon_max"] = region[1]
        summary["region_lat_min"] = region[2]
        summary["region_lat_max"] = region[3]
        summary["map_grid_spacing_deg"] = inferred_grid_spacing(method_ground)
        grid_path = event_csv_dir / f"{event.event_id}_estimated_intensity_{method}.csv"
        grid.to_csv(grid_path, index=False)
        summary["grid_csv"] = str(grid_path)
        if not args.no_png:
            png_path = event_png_dir / f"{event.event_id}_estimated_intensity_{method}.png"
            draw_map(event, stations, grid, method, region, png_path, args.prefecture_boundary)
            summary["png"] = str(png_path)
        summaries.append(summary)
    pd.DataFrame(summaries).to_csv(event_csv_dir / f"{event.event_id}_method_summary.csv", index=False)
    return summaries


def main() -> None:
    args = parse_args()
    methods = normalize_methods(args.methods)
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)
    write_method_note(args.csv_dir)

    events = selected_events(args)
    ground = load_ground_grid(args.ground_grid)
    all_summaries: list[dict[str, object]] = []
    for event in events:
        print(f"Processing {event.event_id} with methods: {','.join(methods)}")
        all_summaries.extend(process_event(event, methods, ground, args))
    out_summary = args.csv_dir / "estimated_intensity_method_summary.csv"
    pd.DataFrame(all_summaries).to_csv(out_summary, index=False)
    print(f"Saved summary: {out_summary}")


if __name__ == "__main__":
    main()
