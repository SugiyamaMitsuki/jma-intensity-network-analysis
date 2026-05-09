#!/usr/bin/env python3
"""Analyze nearest station spacing and paired seismic-intensity differences.

For each event month, the script builds the active station network from
code_p.dat, finds each station's nearest active neighbor, and then compares
intensity values only when both stations in a nearest-neighbor pair reported
the same earthquake.  The primary response variable is |dI|, the absolute
difference in the best available station intensity:

  measured instrumental intensity, if present; otherwise class intensity.

Outputs are separated under:
  outputs/csv/station_nearest_pair_intensity/
  outputs/png/station_nearest_pair_intensity/
"""

from __future__ import annotations

import argparse
import calendar
import math
import os
import tempfile
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from analyze_jma_intensity import (
    AXIS_COLOR,
    GRID_COLOR,
    INTENSITY_VALUE,
    PALETTE,
    SOURCE_RECORD_TYPES,
    StationIndex,
    StationRecord,
    add_network_epoch_guides,
    as_ascii,
    configure_journal_matplotlib,
    haversine_km,
    load_station_index,
    mean,
    numeric,
    parse_intensity_code,
    parse_measured_intensity,
    parse_source_record,
    period_label,
    qtile,
    save_journal_png,
    station_active_at_key,
    station_network_region,
    stdev,
    style_journal_axis,
)
from plot_station_maps_pygmt import configure_conda_gmt
from map_prefecture_boundaries import plot_prefecture_boundaries_pygmt


PAIR_DISTANCE_BINS = [
    (0.0, 1.0, "0-1"),
    (1.0, 2.0, "1-2"),
    (2.0, 5.0, "2-5"),
    (5.0, 10.0, "5-10"),
    (10.0, 20.0, "10-20"),
    (20.0, 50.0, "20-50"),
    (50.0, 1_000_000.0, "50+"),
]


@dataclass(frozen=True)
class StationObservation:
    code: str
    class_intensity: float | None
    measured_intensity: float | None
    best_intensity: float | None


@dataclass(frozen=True)
class NearestNeighbor:
    station: StationRecord
    neighbor: StationRecord
    distance_km: float


def month_end_key(year: int, month: int) -> int:
    last_day = calendar.monthrange(year, month)[1]
    return int(f"{year:04d}{month:02d}{last_day:02d}2359")


def pair_distance_bin(distance_km: float | None) -> str | None:
    if distance_km is None:
        return None
    for low, high, label in PAIR_DISTANCE_BINS:
        if low <= distance_km < high:
            return label
    return None


def distance_bin_order(label: str) -> int:
    for idx, (_, _, bin_label) in enumerate(PAIR_DISTANCE_BINS):
        if label == bin_label:
            return idx
    return 999


def active_station_records_at_month_end(
    station_index: StationIndex,
    year: int,
    month: int,
) -> list[StationRecord]:
    time_key = month_end_key(year, month)
    records: list[StationRecord] = []
    for records_for_code in station_index.by_code.values():
        active = [rec for rec in records_for_code if station_active_at_key(rec, time_key)]
        if not active:
            continue
        active.sort(key=lambda rec: rec.start_key or 0)
        rec = active[-1]
        if (
            rec.latitude is not None
            and rec.longitude is not None
            and 120.0 <= rec.longitude <= 150.0
            and 20.0 <= rec.latitude <= 48.0
        ):
            records.append(rec)
    return records


def nearest_neighbor_map(records: list[StationRecord]) -> dict[str, NearestNeighbor]:
    usable = [
        rec
        for rec in records
        if rec.latitude is not None and rec.longitude is not None
    ]
    if len(usable) < 2:
        return {}

    try:
        import numpy as np
        from scipy.spatial import cKDTree
    except Exception:
        out: dict[str, NearestNeighbor] = {}
        for idx, rec in enumerate(usable):
            nearest: tuple[float, StationRecord] | None = None
            for jdx, other in enumerate(usable):
                if idx == jdx:
                    continue
                distance = haversine_km(rec.latitude, rec.longitude, other.latitude, other.longitude)
                if distance is not None and (nearest is None or distance < nearest[0]):
                    nearest = (distance, other)
            if nearest is not None:
                out[rec.code] = NearestNeighbor(rec, nearest[1], nearest[0])
        return out

    lat = np.radians(np.array([float(rec.latitude) for rec in usable], dtype=float))
    lon = np.radians(np.array([float(rec.longitude) for rec in usable], dtype=float))
    xyz = np.column_stack(
        [
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ]
    )
    chord_distances, neighbor_idx = cKDTree(xyz).query(xyz, k=2)
    nearest_chord = np.clip(chord_distances[:, 1], 0.0, 2.0)
    distances = 2.0 * 6371.0088 * np.arcsin(nearest_chord / 2.0)

    out = {}
    for idx, rec in enumerate(usable):
        neighbor = usable[int(neighbor_idx[idx, 1])]
        out[rec.code] = NearestNeighbor(rec, neighbor, float(distances[idx]))
    return out


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

    years: set[int] = set()
    for part in value.split(","):
        item = part.strip()
        if not item:
            continue
        if ":" in item:
            start, end = [int(token) for token in item.split(":", 1)]
            years.update(range(start, end + 1))
        elif "-" in item and len(item) == 9:
            start, end = [int(token) for token in item.split("-", 1)]
            years.update(range(start, end + 1))
        else:
            years.add(int(item))
    return sorted(years)


def parse_observation(line: bytes) -> StationObservation | None:
    padded = line.ljust(96, b" ")
    code = as_ascii(padded[0:7]).strip()
    if not code:
        return None
    _, class_value = parse_intensity_code(padded[18:19])
    measured_value = parse_measured_intensity(padded[20:22])
    best_value = measured_value if measured_value is not None else class_value
    if best_value is None and class_value is None:
        return None
    return StationObservation(code, class_value, measured_value, best_value)


def finite_abs_diff(value1: float | None, value2: float | None) -> float | None:
    if value1 is None or value2 is None:
        return None
    if not math.isfinite(value1) or not math.isfinite(value2):
        return None
    return abs(value1 - value2)


def lower_code_pair(code1: str, code2: str) -> tuple[str, str]:
    return (code1, code2) if code1 <= code2 else (code2, code1)


def pair_row(
    event,
    obs_by_code: dict[str, StationObservation],
    pair: tuple[str, str],
    nn_map: dict[str, NearestNeighbor],
) -> dict[str, object] | None:
    code1, code2 = pair
    obs1 = obs_by_code.get(code1)
    obs2 = obs_by_code.get(code2)
    if obs1 is None or obs2 is None:
        return None

    nn1 = nn_map.get(code1)
    nn2 = nn_map.get(code2)
    if nn1 is not None and nn1.neighbor.code == code2:
        station1, station2, pair_distance = nn1.station, nn1.neighbor, nn1.distance_km
        source_station_code = code1
    elif nn2 is not None and nn2.neighbor.code == code1:
        station1, station2, pair_distance = nn2.neighbor, nn2.station, nn2.distance_km
        source_station_code = code2
    else:
        return None

    if station1.latitude is None or station1.longitude is None:
        return None
    if station2.latitude is None or station2.longitude is None:
        return None

    best_diff = finite_abs_diff(obs1.best_intensity, obs2.best_intensity)
    class_diff = finite_abs_diff(obs1.class_intensity, obs2.class_intensity)
    measured_diff = finite_abs_diff(obs1.measured_intensity, obs2.measured_intensity)
    if best_diff is None and class_diff is None and measured_diff is None:
        return None

    midpoint_lon = (float(station1.longitude) + float(station2.longitude)) / 2.0
    midpoint_lat = (float(station1.latitude) + float(station2.latitude)) / 2.0
    mutual = (
        nn1 is not None
        and nn2 is not None
        and nn1.neighbor.code == code2
        and nn2.neighbor.code == code1
    )
    year = int(event.year) if event.year is not None else None
    row = {
        "event_id": event.event_id,
        "year": year,
        "month": event.month,
        "day": event.day,
        "hour": event.hour,
        "minute": event.minute,
        "period": period_label(year),
        "magnitude": event.magnitude,
        "depth_km": event.depth_km,
        "event_latitude": event.latitude,
        "event_longitude": event.longitude,
        "hypocenter_region": event.hypocenter_region,
        "station_code_1": code1,
        "station_code_2": code2,
        "station_name_1": station1.name if station1.code == code1 else station2.name,
        "station_name_2": station2.name if station2.code == code2 else station1.name,
        "station_latitude_1": station1.latitude if station1.code == code1 else station2.latitude,
        "station_longitude_1": station1.longitude if station1.code == code1 else station2.longitude,
        "station_latitude_2": station2.latitude if station2.code == code2 else station1.latitude,
        "station_longitude_2": station2.longitude if station2.code == code2 else station1.longitude,
        "midpoint_latitude": midpoint_lat,
        "midpoint_longitude": midpoint_lon,
        "region_1": station_network_region(station1 if station1.code == code1 else station2),
        "region_2": station_network_region(station2 if station2.code == code2 else station1),
        "nearest_source_station_code": source_station_code,
        "mutual_nearest": int(mutual),
        "nearest_pair_distance_km": pair_distance,
        "nearest_pair_distance_bin": pair_distance_bin(pair_distance),
        "station_1_class_intensity": obs1.class_intensity,
        "station_2_class_intensity": obs2.class_intensity,
        "station_1_measured_intensity": obs1.measured_intensity,
        "station_2_measured_intensity": obs2.measured_intensity,
        "station_1_best_intensity": obs1.best_intensity,
        "station_2_best_intensity": obs2.best_intensity,
        "abs_class_intensity_diff": class_diff,
        "abs_measured_intensity_diff": measured_diff,
        "abs_best_intensity_diff": best_diff,
    }
    return row


def event_pair_rows(
    event,
    obs_by_code: dict[str, StationObservation],
    nn_map: dict[str, NearestNeighbor],
) -> list[dict[str, object]]:
    used_pairs: set[tuple[str, str]] = set()
    out: list[dict[str, object]] = []
    for code in obs_by_code:
        nn = nn_map.get(code)
        if nn is None:
            continue
        other_code = nn.neighbor.code
        if other_code not in obs_by_code:
            continue
        pair = lower_code_pair(code, other_code)
        if pair in used_pairs:
            continue
        used_pairs.add(pair)
        row = pair_row(event, obs_by_code, pair, nn_map)
        if row is not None:
            out.append(row)
    return out


def collect_comparisons(
    data_dir: Path,
    years: list[int],
    station_index: StationIndex,
    min_magnitude: float | None,
) -> tuple[list[dict[str, object]], dict[tuple[int, int], dict[str, NearestNeighbor]]]:
    comparisons: list[dict[str, object]] = []
    monthly_nn_cache: dict[tuple[int, int], dict[str, NearestNeighbor]] = {}

    def nn_for_month(year: int, month: int) -> dict[str, NearestNeighbor]:
        key = (year, month)
        if key not in monthly_nn_cache:
            records = active_station_records_at_month_end(station_index, year, month)
            monthly_nn_cache[key] = nearest_neighbor_map(records)
        return monthly_nn_cache[key]

    for year in years:
        zip_path = data_dir / f"i{year}.zip"
        if not zip_path.exists():
            print(f"Skipping missing catalog: {zip_path}")
            continue
        with zipfile.ZipFile(zip_path) as zf:
            names = [name for name in zf.namelist() if not name.endswith("/")]
            if not names:
                continue
            data = zf.read(names[0])

        current = None
        obs_by_code: dict[str, StationObservation] = {}
        event_counter = 0

        def finish_current() -> None:
            if current is None or not obs_by_code:
                return
            if current.year is None or current.month is None:
                return
            if min_magnitude is not None and (
                current.magnitude is None or current.magnitude < min_magnitude
            ):
                return
            rows = event_pair_rows(
                current,
                obs_by_code,
                nn_for_month(int(current.year), int(current.month)),
            )
            comparisons.extend(rows)

        for line_number, raw_line in enumerate(data.splitlines(), start=1):
            if not raw_line:
                continue
            first = raw_line[:1]
            if first in SOURCE_RECORD_TYPES:
                if current is None or obs_by_code:
                    finish_current()
                    event_counter += 1
                    event_id = f"{zip_path.stem}_{event_counter:06d}"
                    current = parse_source_record(raw_line, event_id, line_number)
                    obs_by_code = {}
                else:
                    current.source_record_count += 1
                continue

            if current is None:
                continue
            obs = parse_observation(raw_line)
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
        print(f"{year}: {len(comparisons):,} cumulative nearest-pair comparisons")

    return comparisons, monthly_nn_cache


def to_csv(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(path, index=False)
    print(f"CSV -> {path}")


def summarize_series(values: Iterable[float | None]) -> dict[str, float | int | None]:
    vals = [numeric(value) for value in values]
    vals = [value for value in vals if value is not None]
    return {
        "n": len(vals),
        "mean": mean(vals),
        "sd": stdev(vals),
        "median": qtile(vals, 0.50),
        "q25": qtile(vals, 0.25),
        "q75": qtile(vals, 0.75),
        "q90": qtile(vals, 0.90),
        "share_zero": mean(1.0 if value == 0.0 else 0.0 for value in vals),
        "share_le_0p5": mean(1.0 if value <= 0.5 else 0.0 for value in vals),
        "share_ge_1p0": mean(1.0 if value >= 1.0 else 0.0 for value in vals),
    }


def make_summaries(df):
    import pandas as pd

    def summarize_group(group) -> dict[str, object]:
        stats_best = summarize_series(group["abs_best_intensity_diff"])
        stats_measured = summarize_series(group["abs_measured_intensity_diff"])
        return {
            "n_comparisons": int(stats_best["n"]),
            "n_measured_comparisons": int(stats_measured["n"]),
            "mean_abs_best_diff": stats_best["mean"],
            "sd_abs_best_diff": stats_best["sd"],
            "median_abs_best_diff": stats_best["median"],
            "q25_abs_best_diff": stats_best["q25"],
            "q75_abs_best_diff": stats_best["q75"],
            "q90_abs_best_diff": stats_best["q90"],
            "share_abs_best_diff_zero": stats_best["share_zero"],
            "share_abs_best_diff_le_0p5": stats_best["share_le_0p5"],
            "share_abs_best_diff_ge_1p0": stats_best["share_ge_1p0"],
            "mean_abs_measured_diff": stats_measured["mean"],
            "median_abs_measured_diff": stats_measured["median"],
            "mean_nearest_pair_distance_km": mean(group["nearest_pair_distance_km"]),
            "median_nearest_pair_distance_km": qtile(group["nearest_pair_distance_km"], 0.50),
        }

    bin_rows: list[dict[str, object]] = []
    for label, group in df.groupby("nearest_pair_distance_bin", dropna=True):
        row = {"nearest_pair_distance_bin": label}
        row.update(summarize_group(group))
        bin_rows.append(row)
    bin_summary = pd.DataFrame(bin_rows)
    if not bin_summary.empty:
        bin_summary["_order"] = bin_summary["nearest_pair_distance_bin"].map(distance_bin_order)
        bin_summary = bin_summary.sort_values("_order").drop(columns=["_order"])

    year_rows = []
    for year, group in df.groupby("year", dropna=True):
        row = {"year": int(year)}
        row.update(summarize_group(group))
        year_rows.append(row)
    year_summary = pd.DataFrame(year_rows).sort_values("year")

    period_rows = []
    for period, group in df.groupby("period", dropna=True):
        row = {"period": period}
        row.update(summarize_group(group))
        period_rows.append(row)
    period_summary = pd.DataFrame(period_rows)

    pair_rows = []
    pair_cols = ["station_code_1", "station_code_2"]
    for (code1, code2), group in df.groupby(pair_cols, dropna=True):
        row = {
            "station_code_1": code1,
            "station_code_2": code2,
            "station_name_1": group["station_name_1"].iloc[-1],
            "station_name_2": group["station_name_2"].iloc[-1],
            "station_latitude_1": group["station_latitude_1"].iloc[-1],
            "station_longitude_1": group["station_longitude_1"].iloc[-1],
            "station_latitude_2": group["station_latitude_2"].iloc[-1],
            "station_longitude_2": group["station_longitude_2"].iloc[-1],
            "midpoint_latitude": group["midpoint_latitude"].median(),
            "midpoint_longitude": group["midpoint_longitude"].median(),
            "region_1": group["region_1"].iloc[-1],
            "region_2": group["region_2"].iloc[-1],
            "first_year": int(group["year"].min()),
            "last_year": int(group["year"].max()),
            "mutual_nearest_share": group["mutual_nearest"].mean(),
        }
        row.update(summarize_group(group))
        pair_rows.append(row)
    pair_summary = pd.DataFrame(pair_rows).sort_values(
        ["n_comparisons", "mean_abs_best_diff"], ascending=[False, False]
    )

    region_rows = []
    for region, group in df.groupby("region_1", dropna=True):
        row = {"region": region}
        row.update(summarize_group(group))
        region_rows.append(row)
    region_summary = pd.DataFrame(region_rows).sort_values("region")

    return {
        "distance_bin": bin_summary,
        "year": year_summary,
        "period": period_summary,
        "pair": pair_summary,
        "region": region_summary,
    }


def make_statistical_tests(df):
    import numpy as np
    import pandas as pd
    from scipy import stats

    rows: list[dict[str, object]] = []

    def add_stats(label: str, group) -> None:
        use = group[["nearest_pair_distance_km", "abs_best_intensity_diff"]].dropna()
        use = use[(use["nearest_pair_distance_km"] > 0.0) & (use["abs_best_intensity_diff"].notna())]
        if len(use) < 3:
            return
        x = use["nearest_pair_distance_km"].to_numpy(dtype=float)
        y = use["abs_best_intensity_diff"].to_numpy(dtype=float)
        logx = np.log10(x)
        pearson = stats.pearsonr(logx, y)
        spearman = stats.spearmanr(x, y)
        lin = stats.linregress(logx, y)
        rows.append(
            {
                "group": label,
                "n_comparisons": len(use),
                "pearson_r_log10_distance_vs_abs_diff": float(pearson.statistic),
                "pearson_p_value": float(pearson.pvalue),
                "spearman_rho_distance_vs_abs_diff": float(spearman.statistic),
                "spearman_p_value": float(spearman.pvalue),
                "ols_slope_abs_diff_per_log10_km": float(lin.slope),
                "ols_intercept": float(lin.intercept),
                "ols_r2": float(lin.rvalue**2),
                "ols_p_value": float(lin.pvalue),
            }
        )

    add_stats("all", df)
    for period, group in df.groupby("period", dropna=True):
        add_stats(f"period:{period}", group)
    for year, group in df.groupby("year", dropna=True):
        if len(group) >= 100:
            add_stats(f"year:{int(year)}", group)

    labels = [label for _, _, label in PAIR_DISTANCE_BINS]
    samples = [
        df.loc[df["nearest_pair_distance_bin"] == label, "abs_best_intensity_diff"].dropna().to_numpy(dtype=float)
        for label in labels
    ]
    samples = [sample for sample in samples if len(sample) > 0]
    if len(samples) >= 2:
        kw = stats.kruskal(*samples)
        rows.append(
            {
                "group": "distance_bin_kruskal_wallis",
                "n_comparisons": int(sum(len(sample) for sample in samples)),
                "pearson_r_log10_distance_vs_abs_diff": None,
                "pearson_p_value": None,
                "spearman_rho_distance_vs_abs_diff": None,
                "spearman_p_value": None,
                "ols_slope_abs_diff_per_log10_km": None,
                "ols_intercept": None,
                "ols_r2": float(kw.statistic),
                "ols_p_value": float(kw.pvalue),
            }
        )
    return pd.DataFrame(rows)


def plot_distance_bin_boxplot(df, bin_summary, png_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    configure_journal_matplotlib(plt)
    labels = [label for _, _, label in PAIR_DISTANCE_BINS if label in set(df["nearest_pair_distance_bin"])]
    data = [
        df.loc[df["nearest_pair_distance_bin"] == label, "abs_best_intensity_diff"].dropna().to_numpy()
        for label in labels
    ]
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    bp = ax.boxplot(
        data,
        patch_artist=True,
        widths=0.58,
        showfliers=False,
        medianprops={"color": "#1f1f1f", "linewidth": 1.25},
        boxprops={"facecolor": "#d7dee4", "edgecolor": AXIS_COLOR, "linewidth": 0.8},
        whiskerprops={"color": AXIS_COLOR, "linewidth": 0.8},
        capprops={"color": AXIS_COLOR, "linewidth": 0.8},
    )
    for patch, color in zip(bp["boxes"], ["#d8e1e8", "#c9d8e5", "#b7cce1", "#a5bfd8", "#90b0cd", "#7c9fc1", "#688eb4"]):
        patch.set_facecolor(color)
    means = [
        numeric(bin_summary.loc[bin_summary["nearest_pair_distance_bin"] == label, "mean_abs_best_diff"].iloc[0])
        for label in labels
    ]
    ax.plot(np.arange(1, len(labels) + 1), means, color="#b23a48", marker="o", markersize=4.0, linewidth=1.5, label="Mean")
    ax.set_xticks(np.arange(1, len(labels) + 1), labels)
    ax.set_xlabel("Nearest-station distance (km)")
    ax.set_ylabel("Absolute intensity difference |dI|")
    ax.set_ylim(-0.05, min(4.0, max(1.2, float(df["abs_best_intensity_diff"].quantile(0.995)) + 0.25)))
    ax.legend(loc="upper left")
    style_journal_axis(ax)
    for idx, label in enumerate(labels, start=1):
        n = int(bin_summary.loc[bin_summary["nearest_pair_distance_bin"] == label, "n_comparisons"].iloc[0])
        ax.text(idx, ax.get_ylim()[1] * 0.965, f"n={n:,}", ha="center", va="top", fontsize=7.2, color="#4a4a4a", rotation=90)
    fig.tight_layout()
    save_journal_png(fig, png_dir / "nearest_pair_intensity_difference_by_distance_bin.png")
    plt.close(fig)


def plot_yearly_summary(year_summary, png_dir: Path) -> None:
    import matplotlib.pyplot as plt

    configure_journal_matplotlib(plt)
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    years = year_summary["year"]
    ax.plot(
        years,
        year_summary["median_abs_best_diff"],
        color=PALETTE[0],
        marker="o",
        markersize=3.4,
        linewidth=1.55,
        label="Median |dI|",
    )
    ax.plot(
        years,
        year_summary["mean_abs_best_diff"],
        color=PALETTE[1],
        marker="s",
        markersize=3.2,
        linewidth=1.35,
        label="Mean |dI|",
    )
    ax.fill_between(
        years,
        year_summary["q25_abs_best_diff"],
        year_summary["q75_abs_best_diff"],
        color=PALETTE[0],
        alpha=0.14,
        linewidth=0,
        label="IQR",
    )
    ax.set_xlabel("Year")
    ax.set_ylabel("Nearest-pair absolute intensity difference |dI|")
    ax.set_ylim(0.0, max(1.2, float(year_summary["q90_abs_best_diff"].max()) + 0.2))
    add_network_epoch_guides(ax)
    style_journal_axis(ax)
    ax2 = ax.twinx()
    ax2.bar(
        years,
        year_summary["n_comparisons"],
        width=0.72,
        color="#aeb4bb",
        alpha=0.24,
        label="Paired observations",
    )
    ax2.set_ylabel("Number of paired observations")
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.grid(False)
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, loc="upper left")
    fig.tight_layout()
    save_journal_png(fig, png_dir / "nearest_pair_intensity_difference_by_year.png")
    plt.close(fig)


def plot_hexbin(df, bin_summary, png_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    configure_journal_matplotlib(plt)
    plot_df = df[["nearest_pair_distance_km", "abs_best_intensity_diff"]].dropna()
    fig, ax = plt.subplots(figsize=(6.8, 4.7))
    hb = ax.hexbin(
        plot_df["nearest_pair_distance_km"],
        plot_df["abs_best_intensity_diff"],
        gridsize=(54, 28),
        mincnt=1,
        bins="log",
        cmap="viridis",
        linewidths=0.0,
    )
    centers = []
    medians = []
    for label in bin_summary["nearest_pair_distance_bin"]:
        low, high, _ = next(item for item in PAIR_DISTANCE_BINS if item[2] == label)
        if high > 1000:
            high = max(plot_df["nearest_pair_distance_km"].max(), low + 1.0)
        centers.append(math.sqrt(max(low, 0.1) * high) if low > 0 else high / 2.0)
        medians.append(bin_summary.loc[bin_summary["nearest_pair_distance_bin"] == label, "median_abs_best_diff"].iloc[0])
    ax.plot(centers, medians, color="#d73027", marker="o", markersize=3.8, linewidth=1.55, label="Binned median")
    ax.set_xscale("log")
    ax.set_xlabel("Nearest-station distance (km)")
    ax.set_ylabel("Absolute intensity difference |dI|")
    ax.set_ylim(-0.05, min(4.0, max(1.2, float(plot_df["abs_best_intensity_diff"].quantile(0.995)) + 0.25)))
    style_journal_axis(ax, grid_axis="both")
    ax.legend(loc="upper left")
    cbar = fig.colorbar(hb, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("log10(count)")
    fig.tight_layout()
    save_journal_png(fig, png_dir / "nearest_pair_distance_vs_intensity_difference_hexbin.png")
    plt.close(fig)


def plot_period_ecdf(df, png_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    configure_journal_matplotlib(plt)
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    periods = ["1980-1995", "1996-2003", "2004-2010", "2011-2022"]
    for idx, period in enumerate(periods):
        values = np.sort(df.loc[df["period"] == period, "abs_best_intensity_diff"].dropna().to_numpy())
        if len(values) == 0:
            continue
        y = np.arange(1, len(values) + 1) / len(values)
        ax.step(values, y, where="post", color=PALETTE[idx], linewidth=1.6, label=f"{period} (n={len(values):,})")
    ax.set_xlabel("Absolute intensity difference |dI|")
    ax.set_ylabel("Cumulative probability")
    ax.set_xlim(0.0, min(4.0, max(1.2, float(df["abs_best_intensity_diff"].quantile(0.995)) + 0.25)))
    ax.set_ylim(0.0, 1.0)
    style_journal_axis(ax)
    ax.legend(loc="lower right")
    fig.tight_layout()
    save_journal_png(fig, png_dir / "nearest_pair_intensity_difference_ecdf_by_period.png")
    plt.close(fig)


def draw_pair_difference_map(pair_summary, png_dir: Path, min_pair_samples: int, region: list[float], projection: str) -> None:
    configure_conda_gmt()
    import pygmt

    plot_df = pair_summary[
        (pair_summary["n_comparisons"] >= min_pair_samples)
        & pair_summary["midpoint_longitude"].notna()
        & pair_summary["midpoint_latitude"].notna()
        & pair_summary["mean_abs_best_diff"].notna()
    ].copy()
    if plot_df.empty:
        return
    cap = max(0.75, min(2.0, float(plot_df["mean_abs_best_diff"].quantile(0.98))))
    cpt_file = tempfile.NamedTemporaryFile(suffix=".cpt", delete=False)
    cpt_file.close()
    cpt_path = Path(cpt_file.name)
    pygmt.makecpt(cmap="batlow", series=[0, cap, 0.1], continuous=True, output=str(cpt_path))

    fig = pygmt.Figure()
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_TICK_LENGTH_PRIMARY="2.5p",
        FORMAT_GEO_MAP="dddF",
        FONT_ANNOT_PRIMARY="8p,Helvetica,black",
        FONT_LABEL="9p,Helvetica,black",
        FONT_TITLE="13p,Helvetica-Bold,black",
        MAP_GRID_PEN_PRIMARY="0.15p,#d7d7d7",
        MAP_FRAME_PEN="0.7p,#222222",
    ):
        fig.basemap(
            region=region,
            projection=projection,
            frame=["xafg5+lLongitude", "yafg5+lLatitude", "+tNearest-station pair mean intensity difference"],
        )
        fig.coast(
            land="#f2f2ee",
            water="#e9f2f5",
            shorelines="0.35p,#555555",
            borders="1/0.2p,#888888",
            resolution="i",
            map_scale="jBL+w300k+o0.35c/0.35c+f+lkm",
        )
        plot_prefecture_boundaries_pygmt(fig, pen="0.16p,#666666")
        q70 = float(plot_df["n_comparisons"].quantile(0.70))
        q90 = float(plot_df["n_comparisons"].quantile(0.90))
        size_layers = [
            (plot_df["n_comparisons"] < q70, "c0.055c"),
            ((plot_df["n_comparisons"] >= q70) & (plot_df["n_comparisons"] < q90), "c0.080c"),
            (plot_df["n_comparisons"] >= q90, "c0.110c"),
        ]
        for mask, style in size_layers:
            subset = plot_df[mask]
            if subset.empty:
                continue
            fig.plot(
                x=subset["midpoint_longitude"],
                y=subset["midpoint_latitude"],
                style=style,
                fill=subset["mean_abs_best_diff"],
                cmap=str(cpt_path),
                pen="0.03p,white",
                transparency=7,
            )
        fig.colorbar(cmap=str(cpt_path), frame=["xaf", "y+lMean |dI|"])
        fig.text(
            x=region[0] + 0.65,
            y=region[3] - 0.75,
            text=f"Pairs with n >= {min_pair_samples:,}: {len(plot_df):,}",
            font="9p,Helvetica-Bold,#222222",
            justify="LM",
            fill="white@18",
            clearance="0.08c/0.05c",
        )
    out = png_dir / "nearest_pair_mean_intensity_difference_map.png"
    fig.savefig(out)
    cpt_path.unlink(missing_ok=True)
    print(f"PNG -> {out}")


def draw_annual_nearest_distance_map(
    station_index: StationIndex,
    year: int,
    png_dir: Path,
    region: list[float],
    projection: str,
) -> None:
    configure_conda_gmt()
    import pygmt

    records = active_station_records_at_month_end(station_index, year, 12)
    nn_map = nearest_neighbor_map(records)
    rows = [
        {
            "lon": nn.station.longitude,
            "lat": nn.station.latitude,
            "distance": nn.distance_km,
        }
        for nn in nn_map.values()
        if nn.station.longitude is not None and nn.station.latitude is not None
    ]
    if not rows:
        return
    import pandas as pd

    df = pd.DataFrame(rows)
    cap = max(5.0, min(50.0, float(df["distance"].quantile(0.98))))
    cpt_file = tempfile.NamedTemporaryFile(suffix=".cpt", delete=False)
    cpt_file.close()
    cpt_path = Path(cpt_file.name)
    pygmt.makecpt(cmap="batlow", series=[0, cap, 1.0], continuous=True, output=str(cpt_path))

    fig = pygmt.Figure()
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_TICK_LENGTH_PRIMARY="2.5p",
        FORMAT_GEO_MAP="dddF",
        FONT_ANNOT_PRIMARY="8p,Helvetica,black",
        FONT_LABEL="9p,Helvetica,black",
        FONT_TITLE="13p,Helvetica-Bold,black",
        MAP_GRID_PEN_PRIMARY="0.15p,#d7d7d7",
        MAP_FRAME_PEN="0.7p,#222222",
    ):
        fig.basemap(
            region=region,
            projection=projection,
            frame=["xafg5+lLongitude", "yafg5+lLatitude", f"+tNearest-station distance, active network at {year}-12"],
        )
        fig.coast(
            land="#f2f2ee",
            water="#e9f2f5",
            shorelines="0.35p,#555555",
            borders="1/0.2p,#888888",
            resolution="i",
            map_scale="jBL+w300k+o0.35c/0.35c+f+lkm",
        )
        plot_prefecture_boundaries_pygmt(fig, pen="0.16p,#666666")
        fig.plot(
            x=df["lon"],
            y=df["lat"],
            style="c0.055c",
            fill=df["distance"],
            cmap=str(cpt_path),
            pen="0.025p,white",
            transparency=8,
        )
        fig.colorbar(cmap=str(cpt_path), frame=["xaf", "y+lNearest distance (km)"])
        fig.text(
            x=region[0] + 0.65,
            y=region[3] - 0.75,
            text=f"N = {len(df):,}; median = {df['distance'].median():.2f} km",
            font="9p,Helvetica-Bold,#222222",
            justify="LM",
            fill="white@18",
            clearance="0.08c/0.05c",
        )
    out = png_dir / f"station_nearest_neighbor_distance_map_{year}.png"
    fig.savefig(out)
    cpt_path.unlink(missing_ok=True)
    print(f"PNG -> {out}")


def write_monthly_network_summary(monthly_nn_cache, csv_dir: Path) -> None:
    import pandas as pd

    rows = []
    for (year, month), nn_map in sorted(monthly_nn_cache.items()):
        distances = [nn.distance_km for nn in nn_map.values()]
        rows.append(
            {
                "year": year,
                "month": month,
                "active_station_count": len(nn_map),
                "mean_nearest_neighbor_distance_km": mean(distances),
                "median_nearest_neighbor_distance_km": qtile(distances, 0.50),
                "q25_nearest_neighbor_distance_km": qtile(distances, 0.25),
                "q75_nearest_neighbor_distance_km": qtile(distances, 0.75),
                "q90_nearest_neighbor_distance_km": qtile(distances, 0.90),
            }
        )
    to_csv(csv_dir / "monthly_active_network_nearest_neighbor_summary.csv", pd.DataFrame(rows))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--code-p", type=Path, default=None)
    parser.add_argument("--years", default="auto", help='Year list/range, e.g. "1980:2022".')
    parser.add_argument("--min-magnitude", type=float, default=None)
    parser.add_argument("--csv-dir", type=Path, default=Path("outputs/csv/station_nearest_pair_intensity"))
    parser.add_argument("--png-dir", type=Path, default=Path("outputs/png/station_nearest_pair_intensity"))
    parser.add_argument("--min-pair-samples-for-map", type=int, default=20)
    parser.add_argument("--map-year", type=int, default=2022)
    parser.add_argument("--projection", default="M15c")
    parser.add_argument("--region", nargs=4, type=float, default=[122.0, 146.5, 24.0, 46.5])
    args = parser.parse_args()

    station_index = load_station_index(args.data_dir, args.code_p)
    if station_index is None:
        raise SystemExit("code_p.dat/code_p.zip was not found.")

    years = parse_years(args.years, args.data_dir)
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)

    comparisons, monthly_nn_cache = collect_comparisons(
        data_dir=args.data_dir,
        years=years,
        station_index=station_index,
        min_magnitude=args.min_magnitude,
    )
    if not comparisons:
        raise SystemExit("No nearest-neighbor paired observations were found.")

    import pandas as pd

    df = pd.DataFrame(comparisons)
    df = df.sort_values(["year", "month", "event_id", "station_code_1", "station_code_2"])
    to_csv(args.csv_dir / "nearest_pair_intensity_comparisons.csv", df)

    summaries = make_summaries(df)
    to_csv(args.csv_dir / "nearest_pair_intensity_by_distance_bin.csv", summaries["distance_bin"])
    to_csv(args.csv_dir / "nearest_pair_intensity_by_year.csv", summaries["year"])
    to_csv(args.csv_dir / "nearest_pair_intensity_by_period.csv", summaries["period"])
    to_csv(args.csv_dir / "nearest_pair_intensity_by_station_pair.csv", summaries["pair"])
    to_csv(args.csv_dir / "nearest_pair_intensity_by_region.csv", summaries["region"])
    to_csv(args.csv_dir / "nearest_pair_intensity_statistical_tests.csv", make_statistical_tests(df))
    write_monthly_network_summary(monthly_nn_cache, args.csv_dir)

    plot_distance_bin_boxplot(df, summaries["distance_bin"], args.png_dir)
    plot_yearly_summary(summaries["year"], args.png_dir)
    plot_hexbin(df, summaries["distance_bin"], args.png_dir)
    plot_period_ecdf(df, args.png_dir)
    draw_pair_difference_map(
        summaries["pair"],
        args.png_dir,
        min_pair_samples=args.min_pair_samples_for_map,
        region=args.region,
        projection=args.projection,
    )
    draw_annual_nearest_distance_map(
        station_index,
        args.map_year,
        args.png_dir,
        region=args.region,
        projection=args.projection,
    )

    overall = summarize_series(df["abs_best_intensity_diff"])
    print("Overall nearest-pair intensity difference:")
    print(f"  comparisons: {int(overall['n']):,}")
    print(f"  mean |dI|: {overall['mean']:.3f}")
    print(f"  median |dI|: {overall['median']:.3f}")
    print(f"  share |dI| = 0: {overall['share_zero']:.3f}")
    print(f"  share |dI| <= 0.5: {overall['share_le_0p5']:.3f}")
    print(f"  share |dI| >= 1.0: {overall['share_ge_1p0']:.3f}")


if __name__ == "__main__":
    main()
