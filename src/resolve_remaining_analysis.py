#!/usr/bin/env python3
"""Close remaining JMA intensity-analysis issues.

This script adds three sensitivity analyses that were intentionally left as
limitations in the previous report:

1. Fixed station-code subsets for average-intensity comparisons.
2. Epicentral-distance bands so low-intensity far-field observations do not
   dominate one undifferentiated event mean.
3. A cautious IDW-based intensity footprint proxy with station-support metrics.

The footprint is a local circular, non-land-clipped proxy, not an official
isoseismal polygon.  It is meant to test whether area-style metrics can be used
and what support information must accompany them.
"""

from __future__ import annotations

import argparse
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from analyze_jma_intensity import (
    INTENSITY_VALUE,
    MAG_BINS,
    SOURCE_RECORD_TYPES,
    StationIndex,
    active_station_records_at_year_end,
    as_ascii,
    configure_journal_matplotlib,
    haversine_km,
    is_land_region_approx,
    load_station_index,
    mag_bin,
    mean,
    parse_intensity_code,
    parse_measured_intensity,
    parse_source_record,
    period_label,
    qtile,
    save_journal_png,
    stdev,
    style_journal_axis,
)


DISTANCE_BANDS = [
    (0.0, 20.0, "0-20"),
    (20.0, 50.0, "20-50"),
    (50.0, 100.0, "50-100"),
    (100.0, 200.0, "100-200"),
]
FOOTPRINT_THRESHOLDS = [2.0, 3.0, 4.0, 5.0]
PERIOD_ORDER = ["1980-1995", "1996-2003", "2004-2010", "2011-2022"]
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
EARTH_RADIUS_KM = 6371.0088


@dataclass
class Observation:
    code: str
    latitude: float | None
    longitude: float | None
    distance_km: float | None
    class_intensity: float | None
    measured_intensity: float | None
    best_intensity: float | None


@dataclass
class EventObservations:
    event_id: str
    year: int | None
    month: int | None
    day: int | None
    hour: int | None
    minute: int | None
    latitude: float | None
    longitude: float | None
    depth_km: float | None
    magnitude: float | None
    max_intensity_value: float | None
    hypocenter_region: str
    observations: list[Observation]

    @property
    def period(self) -> str | None:
        return period_label(self.year) if self.year is not None else None

    @property
    def mag_bin(self) -> str | None:
        return mag_bin(self.magnitude)

    @property
    def is_land_region_approx(self) -> bool:
        return is_land_region_approx(self.hypocenter_region)

    @property
    def is_shallow_20km(self) -> bool:
        return self.depth_km is not None and self.depth_km <= 20.0


def parse_years(value: str | None, data_dir: Path) -> list[int]:
    if value is None or value.lower() == "auto":
        years = sorted(
            int(path.stem[1:])
            for path in data_dir.glob("i[0-9][0-9][0-9][0-9].zip")
            if path.stem[1:].isdigit()
        )
        if not years:
            raise SystemExit(f"No iYYYY.zip files found in {data_dir}")
        return years
    out: set[int] = set()
    for part in value.split(","):
        item = part.strip()
        if not item:
            continue
        if ":" in item:
            start, end = [int(token) for token in item.split(":", 1)]
            out.update(range(start, end + 1))
        elif "-" in item and len(item) == 9:
            start, end = [int(token) for token in item.split("-", 1)]
            out.update(range(start, end + 1))
        else:
            out.add(int(item))
    return sorted(out)


def parse_observation(line: bytes, event, station_index: StationIndex) -> Observation | None:
    padded = line.ljust(96, b" ")
    code = as_ascii(padded[0:7]).strip()
    if not code:
        return None
    _, class_value = parse_intensity_code(padded[18:19])
    measured_value = parse_measured_intensity(padded[20:22])
    best_value = measured_value if measured_value is not None else class_value
    if best_value is None and class_value is None:
        return None

    station = station_index.get(code, event.event_time_key())
    latitude = station.latitude if station is not None else None
    longitude = station.longitude if station is not None else None
    distance = haversine_km(event.latitude, event.longitude, latitude, longitude)
    return Observation(
        code=code,
        latitude=latitude,
        longitude=longitude,
        distance_km=distance,
        class_intensity=class_value,
        measured_intensity=measured_value,
        best_intensity=best_value,
    )


def parse_catalog_observations(
    data_dir: Path,
    years: list[int],
    station_index: StationIndex,
) -> list[EventObservations]:
    events: list[EventObservations] = []
    for year in years:
        path = data_dir / f"i{year}.zip"
        if not path.exists():
            print(f"Skipping missing catalog: {path}")
            continue
        with zipfile.ZipFile(path) as zf:
            names = [name for name in zf.namelist() if not name.endswith("/")]
            if not names:
                continue
            data = zf.read(names[0])

        current = None
        obs_by_code: dict[str, Observation] = {}
        event_counter = 0

        def finish_current() -> None:
            if current is None or not obs_by_code:
                return
            values = [obs.best_intensity for obs in obs_by_code.values() if obs.best_intensity is not None]
            max_obs = max(values) if values else None
            max_value = current.max_intensity_header_value if current.max_intensity_header_value is not None else max_obs
            events.append(
                EventObservations(
                    event_id=current.event_id,
                    year=current.year,
                    month=current.month,
                    day=current.day,
                    hour=current.hour,
                    minute=current.minute,
                    latitude=current.latitude,
                    longitude=current.longitude,
                    depth_km=current.depth_km,
                    magnitude=current.magnitude,
                    max_intensity_value=max_value,
                    hypocenter_region=current.hypocenter_region,
                    observations=list(obs_by_code.values()),
                )
            )

        for line_number, raw_line in enumerate(data.splitlines(), start=1):
            if not raw_line:
                continue
            if raw_line[:1] in SOURCE_RECORD_TYPES:
                if current is None or obs_by_code:
                    finish_current()
                    event_counter += 1
                    current = parse_source_record(raw_line, f"{path.stem}_{event_counter:06d}", line_number)
                    obs_by_code = {}
                else:
                    current.source_record_count += 1
                continue
            if current is None:
                continue
            obs = parse_observation(raw_line, current, station_index)
            if obs is None:
                continue
            previous = obs_by_code.get(obs.code)
            if previous is None:
                obs_by_code[obs.code] = obs
                continue
            prev_best = previous.best_intensity if previous.best_intensity is not None else -math.inf
            new_best = obs.best_intensity if obs.best_intensity is not None else -math.inf
            if new_best >= prev_best:
                obs_by_code[obs.code] = obs

        finish_current()
        print(f"{year}: {len(events):,} cumulative events with station observations")
    return events


def selected_events(events: list[EventObservations]) -> list[EventObservations]:
    return [
        event
        for event in events
        if event.period is not None
        and event.mag_bin is not None
        and event.is_land_region_approx
        and event.is_shallow_20km
        and event.latitude is not None
        and event.longitude is not None
    ]


def station_code_sets(station_index: StationIndex) -> tuple[dict[str, set[str]], pd.DataFrame]:
    active_by_year = {
        year: {
            rec.code
            for rec in active_station_records_at_year_end(station_index, year)
            if rec.latitude is not None and rec.longitude is not None
        }
        for year in [1980, 1995, 2003, 2022]
    }
    code_sets = {
        "legacy_1980_codes": active_by_year[1980],
        "legacy_1995_codes": active_by_year[1995],
        "common_1980_2022_codes": active_by_year[1980] & active_by_year[2022],
        "common_1995_2022_codes": active_by_year[1995] & active_by_year[2022],
    }
    rows = []
    for name, codes in code_sets.items():
        rows.append({"station_set": name, "n_station_codes": len(codes)})
    return code_sets, pd.DataFrame(rows)


def summarize_values(values: Iterable[float | None]) -> dict[str, float | int | None]:
    vals = [v for v in values if v is not None and math.isfinite(v)]
    return {
        "n": len(vals),
        "mean": mean(vals),
        "sd": stdev(vals),
        "median": qtile(vals, 0.50),
        "q25": qtile(vals, 0.25),
        "q75": qtile(vals, 0.75),
    }


def distance_band_label(distance_km: float | None) -> str | None:
    if distance_km is None:
        return None
    for low, high, label in DISTANCE_BANDS:
        if low <= distance_km < high:
            return label
    return None


def make_distance_and_fixed_event_rows(
    events: list[EventObservations],
    code_sets: dict[str, set[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    distance_rows: list[dict[str, object]] = []
    fixed_rows: list[dict[str, object]] = []

    for event in selected_events(events):
        base = {
            "event_id": event.event_id,
            "year": event.year,
            "period": event.period,
            "mag_bin": event.mag_bin,
            "magnitude": event.magnitude,
            "depth_km": event.depth_km,
            "max_intensity_value": event.max_intensity_value,
        }

        for low, high, label in DISTANCE_BANDS:
            obs = [
                item
                for item in event.observations
                if item.distance_km is not None
                and low <= item.distance_km < high
                and item.best_intensity is not None
            ]
            if not obs:
                continue
            values = [item.best_intensity for item in obs]
            distance_rows.append(
                {
                    **base,
                    "distance_band_km": label,
                    "distance_low_km": low,
                    "distance_high_km": high,
                    "n_observations": len(values),
                    "mean_intensity": mean(values),
                    "median_intensity": qtile(values, 0.50),
                    "max_intensity": max(values),
                }
            )

        all_values = [item.best_intensity for item in event.observations if item.best_intensity is not None]
        if all_values:
            fixed_rows.append(
                {
                    **base,
                    "station_set": "all_observed",
                    "n_observations": len(all_values),
                    "mean_intensity": mean(all_values),
                    "median_intensity": qtile(all_values, 0.50),
                    "max_intensity": max(all_values),
                }
            )

        for set_name, codes in code_sets.items():
            obs = [
                item
                for item in event.observations
                if item.code in codes and item.best_intensity is not None
            ]
            if not obs:
                continue
            values = [item.best_intensity for item in obs]
            fixed_rows.append(
                {
                    **base,
                    "station_set": set_name,
                    "n_observations": len(values),
                    "mean_intensity": mean(values),
                    "median_intensity": qtile(values, 0.50),
                    "max_intensity": max(values),
                }
            )

    return pd.DataFrame(distance_rows), pd.DataFrame(fixed_rows)


def grouped_summary(df: pd.DataFrame, group_cols: list[str], value_col: str) -> pd.DataFrame:
    rows = []
    for keys, group in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        stats = summarize_values(group[value_col])
        row = {col: key for col, key in zip(group_cols, keys)}
        row.update(
            {
                "n_events": int(group["event_id"].nunique()),
                "n_rows": len(group),
                f"mean_{value_col}": stats["mean"],
                f"median_{value_col}": stats["median"],
                f"q25_{value_col}": stats["q25"],
                f"q75_{value_col}": stats["q75"],
                "mean_n_observations": mean(group["n_observations"]),
                "median_n_observations": qtile(group["n_observations"], 0.50),
            }
        )
        if "max_intensity" in group.columns:
            row["mean_max_intensity"] = mean(group["max_intensity"])
        rows.append(row)
    return pd.DataFrame(rows)


def lonlat_to_xy_km(lat, lon, lat0: float, lon0: float):
    import numpy as np

    lat_arr = np.asarray(lat, dtype=float)
    lon_arr = np.asarray(lon, dtype=float)
    x = EARTH_RADIUS_KM * math.cos(math.radians(lat0)) * np.radians(lon_arr - lon0)
    y = EARTH_RADIUS_KM * np.radians(lat_arr - lat0)
    return x, y


def xy_to_lonlat_km(x, y, lat0: float, lon0: float):
    import numpy as np

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    lat = lat0 + np.degrees(y_arr / EARTH_RADIUS_KM)
    lon = lon0 + np.degrees(x_arr / (EARTH_RADIUS_KM * math.cos(math.radians(lat0))))
    return lon, lat


def idw_grid_for_event(
    event: EventObservations,
    radius_km: float,
    spacing_km: float,
    power: float,
    k_neighbors: int,
):
    import numpy as np
    from scipy.spatial import cKDTree

    obs = [
        item
        for item in event.observations
        if item.latitude is not None
        and item.longitude is not None
        and item.best_intensity is not None
        and item.distance_km is not None
        and item.distance_km <= radius_km
    ]
    if len(obs) < 4 or event.latitude is None or event.longitude is None:
        return None

    obs_lats = np.array([item.latitude for item in obs], dtype=float)
    obs_lons = np.array([item.longitude for item in obs], dtype=float)
    obs_values = np.array([item.best_intensity for item in obs], dtype=float)
    obs_x, obs_y = lonlat_to_xy_km(obs_lats, obs_lons, float(event.latitude), float(event.longitude))
    obs_xy = np.column_stack([obs_x, obs_y])

    axis = np.arange(-radius_km, radius_km + spacing_km * 0.5, spacing_km)
    xx, yy = np.meshgrid(axis, axis)
    mask = (xx**2 + yy**2) <= radius_km**2
    grid_xy = np.column_stack([xx[mask], yy[mask]])
    k = min(k_neighbors, len(obs))
    distances, indices = cKDTree(obs_xy).query(grid_xy, k=k)
    if k == 1:
        distances = distances[:, None]
        indices = indices[:, None]
    exact = distances[:, 0] < 0.05
    weights = 1.0 / np.maximum(distances, 0.25) ** power
    predicted = np.sum(weights * obs_values[indices], axis=1) / np.sum(weights, axis=1)
    predicted[exact] = obs_values[indices[exact, 0]]

    grid = np.full(xx.shape, np.nan, dtype=float)
    grid[mask] = predicted
    lon_grid, lat_grid = xy_to_lonlat_km(xx, yy, float(event.latitude), float(event.longitude))
    return {
        "grid": grid,
        "lon_grid": lon_grid,
        "lat_grid": lat_grid,
        "mask": mask,
        "observations": obs,
        "cell_area_km2": spacing_km * spacing_km,
    }


def footprint_rows(
    events: list[EventObservations],
    station_index: StationIndex,
    radius_km: float,
    spacing_km: float,
    power: float,
    k_neighbors: int,
) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    rows: list[dict[str, object]] = []
    grids: dict[str, dict[str, object]] = {}
    for event in selected_events(events):
        if event.magnitude is None or event.magnitude < 4.0:
            continue
        result = idw_grid_for_event(event, radius_km, spacing_km, power, k_neighbors)
        if result is None:
            continue
        grid = result["grid"]
        mask = result["mask"]
        obs = result["observations"]
        active_150 = None
        if event.latitude is not None and event.longitude is not None:
            active_150 = station_index.active_station_counts_within_radii(
                event.latitude, event.longitude, event.year and int(f"{event.year:04d}{event.month:02d}{event.day:02d}{event.hour:02d}{event.minute:02d}"),
                [int(radius_km)],
            ).get(int(radius_km))
        detected_150 = len(obs)
        support_ratio = detected_150 / active_150 if active_150 not in (None, 0) else None
        base = {
            "event_id": event.event_id,
            "year": event.year,
            "period": event.period,
            "mag_bin": event.mag_bin,
            "magnitude": event.magnitude,
            "depth_km": event.depth_km,
            "latitude": event.latitude,
            "longitude": event.longitude,
            "max_intensity_value": event.max_intensity_value,
            "n_observations_used": detected_150,
            "active_station_count_within_radius": active_150,
            "detected_to_active_ratio_within_radius": support_ratio,
            "radius_km": radius_km,
            "grid_spacing_km": spacing_km,
            "idw_power": power,
            "idw_k_neighbors": k_neighbors,
        }
        values = grid[mask]
        for threshold in FOOTPRINT_THRESHOLDS:
            rows.append(
                {
                    **base,
                    "threshold_intensity": threshold,
                    "footprint_area_km2": float((values >= threshold).sum() * result["cell_area_km2"]),
                    "observed_station_count_ge_threshold": sum(
                        item.best_intensity is not None and item.best_intensity >= threshold for item in obs
                    ),
                }
            )
        grids[event.event_id] = result
    return pd.DataFrame(rows), grids


def footprint_summary(footprints: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["period", "mag_bin", "threshold_intensity"]
    for keys, group in footprints.groupby(group_cols, dropna=False):
        period, mlabel, threshold = keys
        stats = summarize_values(group["footprint_area_km2"])
        rows.append(
            {
                "period": period,
                "mag_bin": mlabel,
                "threshold_intensity": threshold,
                "n_events": int(group["event_id"].nunique()),
                "mean_area_km2": stats["mean"],
                "median_area_km2": stats["median"],
                "q25_area_km2": stats["q25"],
                "q75_area_km2": stats["q75"],
                "mean_observed_station_count_ge_threshold": mean(group["observed_station_count_ge_threshold"]),
                "mean_n_observations_used": mean(group["n_observations_used"]),
                "mean_active_station_count_within_radius": mean(group["active_station_count_within_radius"]),
                "mean_detected_to_active_ratio_within_radius": mean(group["detected_to_active_ratio_within_radius"]),
            }
        )
    return pd.DataFrame(rows)


def plot_distance_band_summary(summary: pd.DataFrame, png_dir: Path) -> None:
    import matplotlib.pyplot as plt

    configure_journal_matplotlib(plt)
    data = summary[summary["mag_bin"].isin(["4.0<=M<5.0", "5.0<=M<6.0"])].copy()
    if data.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7), sharey=True)
    for ax, mlabel in zip(axes, ["4.0<=M<5.0", "5.0<=M<6.0"]):
        sub = data[data["mag_bin"] == mlabel]
        for idx, period in enumerate(PERIOD_ORDER):
            row = sub[sub["period"] == period].sort_values("distance_low_km")
            if row.empty:
                continue
            ax.plot(
                row["distance_band_km"],
                row["mean_mean_intensity"],
                marker="o",
                linewidth=1.5,
                markersize=3.8,
                color=PALETTE[idx],
                label=period,
            )
        ax.set_title(mlabel)
        ax.set_xlabel("Epicentral distance band (km)")
        style_journal_axis(ax)
        ax.tick_params(axis="x", rotation=35)
    axes[0].set_ylabel("Mean station intensity within band")
    axes[0].legend(loc="upper right")
    fig.tight_layout()
    save_journal_png(fig, png_dir / "distance_band_mean_intensity_by_period.png")
    plt.close(fig)


def plot_fixed_station_summary(summary: pd.DataFrame, png_dir: Path) -> None:
    import matplotlib.pyplot as plt

    configure_journal_matplotlib(plt)
    data = summary[
        (summary["mag_bin"] == "4.0<=M<5.0")
        & summary["station_set"].isin(["all_observed", "legacy_1980_codes", "legacy_1995_codes", "common_1995_2022_codes"])
    ].copy()
    if data.empty:
        return
    label_map = {
        "all_observed": "All observed",
        "legacy_1980_codes": "1980 active codes",
        "legacy_1995_codes": "1995 active codes",
        "common_1995_2022_codes": "Common 1995/2022",
    }
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    for idx, station_set in enumerate(label_map):
        row = data[data["station_set"] == station_set].copy()
        row["period"] = pd.Categorical(row["period"], PERIOD_ORDER, ordered=True)
        row = row.sort_values("period")
        if row.empty:
            continue
        ax.plot(
            row["period"].astype(str),
            row["mean_mean_intensity"],
            marker="o",
            linewidth=1.6,
            markersize=4.0,
            color=PALETTE[idx],
            label=label_map[station_set],
        )
    ax.set_xlabel("Period")
    ax.set_ylabel("Mean event-average intensity, M4-5 shallow inland")
    style_journal_axis(ax)
    ax.legend(loc="best")
    fig.tight_layout()
    save_journal_png(fig, png_dir / "fixed_station_mean_intensity_sensitivity.png")
    plt.close(fig)


def plot_footprint_summary(summary: pd.DataFrame, png_dir: Path) -> None:
    import matplotlib.pyplot as plt

    configure_journal_matplotlib(plt)
    data = summary[summary["mag_bin"].isin(["4.0<=M<5.0", "5.0<=M<6.0"])].copy()
    if data.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7), sharey=False)
    for ax, mlabel in zip(axes, ["4.0<=M<5.0", "5.0<=M<6.0"]):
        sub = data[data["mag_bin"] == mlabel]
        for idx, period in enumerate(PERIOD_ORDER):
            row = sub[sub["period"] == period].sort_values("threshold_intensity")
            if row.empty:
                continue
            ax.plot(
                row["threshold_intensity"],
                row["median_area_km2"],
                marker="o",
                linewidth=1.5,
                markersize=3.8,
                color=PALETTE[idx],
                label=period,
            )
        ax.set_title(mlabel)
        ax.set_xlabel("Intensity threshold")
        ax.set_yscale("symlog", linthresh=25)
        style_journal_axis(ax)
    axes[0].set_ylabel("Median IDW footprint area (km2)")
    axes[0].legend(loc="upper right")
    fig.tight_layout()
    save_journal_png(fig, png_dir / "idw_footprint_area_by_period_threshold.png")
    plt.close(fig)


def plot_footprint_support(footprints: pd.DataFrame, png_dir: Path) -> None:
    import matplotlib.pyplot as plt

    configure_journal_matplotlib(plt)
    data = footprints[
        (footprints["threshold_intensity"] == 3.0)
        & (footprints["mag_bin"] == "4.0<=M<5.0")
        & footprints["footprint_area_km2"].notna()
    ].copy()
    if data.empty:
        return
    fig, ax = plt.subplots(figsize=(6.7, 4.2))
    for idx, period in enumerate(PERIOD_ORDER):
        sub = data[data["period"] == period]
        if sub.empty:
            continue
        ax.scatter(
            sub["n_observations_used"],
            sub["footprint_area_km2"],
            s=14,
            alpha=0.42,
            color=PALETTE[idx],
            linewidths=0,
            label=period,
        )
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=25)
    ax.set_xlabel("Observed stations used for IDW, within 150 km")
    ax.set_ylabel("IDW area with intensity >= 3 (km2)")
    style_journal_axis(ax, grid_axis="both")
    ax.legend(loc="upper left")
    fig.tight_layout()
    save_journal_png(fig, png_dir / "idw_footprint_area_vs_station_support.png")
    plt.close(fig)


def plot_representative_footprints(
    events: list[EventObservations],
    footprints: pd.DataFrame,
    grids: dict[str, dict[str, object]],
    png_dir: Path,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    configure_journal_matplotlib(plt)
    candidate = footprints[
        (footprints["threshold_intensity"] == 3.0)
        & (footprints["mag_bin"] == "4.0<=M<5.0")
        & (footprints["n_observations_used"] >= 6)
    ].copy()
    if candidate.empty:
        return
    event_by_id = {event.event_id: event for event in events}
    selected_ids = []
    for period in PERIOD_ORDER:
        sub = candidate[candidate["period"] == period].copy()
        if sub.empty:
            continue
        sub["score"] = sub["n_observations_used"].rank(pct=True) + sub["footprint_area_km2"].rank(pct=True)
        selected_ids.append(str(sub.sort_values("score", ascending=False).iloc[0]["event_id"]))
    if not selected_ids:
        return

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2), constrained_layout=True)
    axes = axes.ravel()
    for ax, event_id in zip(axes, selected_ids):
        event = event_by_id[event_id]
        result = grids[event_id]
        grid = result["grid"]
        lon_grid = result["lon_grid"]
        lat_grid = result["lat_grid"]
        obs = result["observations"]
        im = ax.pcolormesh(lon_grid, lat_grid, grid, shading="auto", cmap="magma", vmin=0.5, vmax=5.5)
        levels = [2, 3, 4, 5]
        try:
            ax.contour(lon_grid, lat_grid, grid, levels=levels, colors="white", linewidths=0.65, alpha=0.8)
        except Exception:
            pass
        ax.scatter(
            [item.longitude for item in obs],
            [item.latitude for item in obs],
            c=[item.best_intensity for item in obs],
            cmap="magma",
            vmin=0.5,
            vmax=5.5,
            s=16,
            edgecolor="white",
            linewidth=0.25,
        )
        ax.scatter([event.longitude], [event.latitude], marker="*", s=90, color="#00A6D6", edgecolor="black", linewidth=0.4)
        ax.set_title(f"{event.period}  {event.event_id}  M{event.magnitude:.1f}, n={len(obs)}")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_aspect("equal", adjustable="box")
        style_journal_axis(ax, grid_axis="both")
    for ax in axes[len(selected_ids) :]:
        ax.axis("off")
    cbar = fig.colorbar(im, ax=axes[: len(selected_ids)], fraction=0.035, pad=0.02)
    cbar.set_label("IDW intensity")
    save_journal_png(fig, png_dir / "representative_idw_intensity_footprints.png")
    plt.close(fig)


def write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"CSV -> {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--code-p", type=Path, default=None)
    parser.add_argument("--years", default="1980:2022")
    parser.add_argument("--csv-dir", type=Path, default=Path("outputs/csv/remaining_closure"))
    parser.add_argument("--png-dir", type=Path, default=Path("outputs/png/remaining_closure"))
    parser.add_argument("--footprint-radius-km", type=float, default=150.0)
    parser.add_argument("--footprint-grid-km", type=float, default=5.0)
    parser.add_argument("--idw-power", type=float, default=2.0)
    parser.add_argument("--idw-k-neighbors", type=int, default=8)
    args = parser.parse_args()

    station_index = load_station_index(args.data_dir, args.code_p)
    if station_index is None:
        raise SystemExit("code_p.dat/code_p.zip was not found.")

    years = parse_years(args.years, args.data_dir)
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)

    events = parse_catalog_observations(args.data_dir, years, station_index)
    code_sets, code_set_summary = station_code_sets(station_index)
    write_csv(args.csv_dir / "fixed_station_code_sets.csv", code_set_summary)

    distance_events, fixed_events = make_distance_and_fixed_event_rows(events, code_sets)
    write_csv(args.csv_dir / "event_distance_band_intensity_metrics.csv", distance_events)
    write_csv(args.csv_dir / "event_fixed_station_intensity_metrics.csv", fixed_events)

    distance_summary = grouped_summary(
        distance_events,
        ["period", "mag_bin", "distance_band_km", "distance_low_km", "distance_high_km"],
        "mean_intensity",
    )
    fixed_summary = grouped_summary(
        fixed_events,
        ["period", "mag_bin", "station_set"],
        "mean_intensity",
    )
    write_csv(args.csv_dir / "period_distance_band_intensity_summary.csv", distance_summary)
    write_csv(args.csv_dir / "period_fixed_station_intensity_summary.csv", fixed_summary)

    footprints, grids = footprint_rows(
        events,
        station_index,
        radius_km=args.footprint_radius_km,
        spacing_km=args.footprint_grid_km,
        power=args.idw_power,
        k_neighbors=args.idw_k_neighbors,
    )
    footprint_summary_df = footprint_summary(footprints)
    write_csv(args.csv_dir / "event_idw_intensity_footprint_area.csv", footprints)
    write_csv(args.csv_dir / "period_idw_intensity_footprint_area_summary.csv", footprint_summary_df)

    plot_distance_band_summary(distance_summary, args.png_dir)
    plot_fixed_station_summary(fixed_summary, args.png_dir)
    plot_footprint_summary(footprint_summary_df, args.png_dir)
    plot_footprint_support(footprints, args.png_dir)
    plot_representative_footprints(events, footprints, grids, args.png_dir)

    print(f"Events parsed: {len(events):,}")
    print(f"Land/shallow events used for fixed and distance-band analyses: {len(selected_events(events)):,}")
    print(f"IDW footprint event-threshold rows: {len(footprints):,}")
    print(f"CSV directory: {args.csv_dir.resolve()}")
    print(f"PNG directory: {args.png_dir.resolve()}")


if __name__ == "__main__":
    main()
