#!/usr/bin/env python3
"""Create spatial PyGMT visualizations for the JMA intensity analysis."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

from analyze_jma_intensity import (
    NETWORK_REGIONS,
    PERIODS,
    StationRecord,
    active_station_records_at_year_end,
    load_station_index,
    numeric,
    station_network_region,
)
from map_prefecture_boundaries import plot_prefecture_boundaries_pygmt


MAP_REGION = [122.0, 146.5, 24.0, 46.5]
PROJECTION = "M15c"
PERIOD_COLORS = {
    "1980-1995": "#0072B2",
    "1996-2003": "#D55E00",
    "2004-2010": "#009E73",
    "2011-2022": "#CC79A7",
}
REGION_COLORS = {
    "Hokkaido": "#0072B2",
    "Tohoku": "#56B4E9",
    "Kanto": "#D55E00",
    "Chubu": "#009E73",
    "Kinki": "#CC79A7",
    "Chugoku": "#E69F00",
    "Shikoku": "#F0E442",
    "Kyushu": "#6A3D9A",
    "Okinawa": "#999999",
    "Other": "#333333",
}


GRID_METRICS = [
    {
        "name": "event_count",
        "column": "n_events",
        "label": "Number of events",
        "cmap": "batlow",
        "series": None,
        "min_events": 1,
    },
    {
        "name": "mean_max_intensity",
        "column": "mean_max_intensity",
        "label": "Mean maximum intensity",
        "cmap": "lajolla",
        "series": [1, 6, 0.5],
        "min_events": 3,
    },
    {
        "name": "mean_station_intensity",
        "column": "mean_station_intensity",
        "label": "Mean station intensity",
        "cmap": "lajolla",
        "series": [0.8, 3.2, 0.2],
        "min_events": 3,
    },
    {
        "name": "median_shortest_distance",
        "column": "median_nearest_distance_km",
        "label": "Median shortest distance [km]",
        "cmap": "devon",
        "series": [0, 80, 10],
        "min_events": 3,
    },
    {
        "name": "mean_active_stations_50km",
        "column": "mean_active_stations_50km",
        "label": "Active stations within 50 km",
        "cmap": "batlow",
        "series": [0, 120, 20],
        "min_events": 3,
    },
    {
        "name": "mean_detected_stations_50km",
        "column": "mean_detected_stations_50km",
        "label": "Detected stations within 50 km",
        "cmap": "batlow",
        "series": [0, 90, 15],
        "min_events": 3,
    },
    {
        "name": "share_nearest_within_20km",
        "column": "share_nearest_within_20km",
        "label": "Share with nearest station <=20 km",
        "cmap": "lapaz",
        "series": [0, 1, 0.1],
        "min_events": 3,
    },
    {
        "name": "mean_station_records",
        "column": "mean_station_records",
        "label": "Mean station records",
        "cmap": "batlow",
        "series": [0, 250, 50],
        "min_events": 3,
    },
]

DIFFERENCE_METRICS = [
    ("delta_mean_max_intensity", "Delta mean maximum intensity", [-2, 2, 0.5]),
    ("delta_mean_station_intensity", "Delta mean station intensity", [-1.2, 1.2, 0.3]),
    ("delta_median_nearest_distance_km", "Delta median shortest distance [km]", [-80, 80, 20]),
    ("delta_mean_active_stations_50km", "Delta active stations within 50 km", [-120, 120, 30]),
    ("delta_share_nearest_within_20km", "Delta share with nearest station <=20 km", [-1, 1, 0.25]),
]


def configure_conda_gmt() -> None:
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        return
    prefix = Path(conda_prefix)
    bin_dir = prefix / "bin"
    libgmt = prefix / "lib" / "libgmt.dylib"
    if bin_dir.exists():
        os.environ["PATH"] = f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
    if libgmt.exists():
        os.environ.setdefault("GMT_LIBRARY_PATH", str(libgmt))


def period_label(year: int) -> str | None:
    for start, end, label in PERIODS:
        if start <= int(year) <= end:
            return label
    return None


def load_filtered_events(path: Path, min_magnitude: float) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[
        (df["is_land_region_approx"] == True)
        & (df["is_shallow_20km"] == True)
        & (df["magnitude"] >= min_magnitude)
        & df["latitude"].notna()
        & df["longitude"].notna()
    ].copy()
    df = df[
        (df["longitude"] >= MAP_REGION[0])
        & (df["longitude"] <= MAP_REGION[1])
        & (df["latitude"] >= MAP_REGION[2])
        & (df["latitude"] <= MAP_REGION[3])
    ].copy()
    df["analysis_period"] = df["year"].map(period_label)
    df = df[df["analysis_period"].notna()].copy()
    df["nearest_distance_capped_km"] = df["nearest_station_distance_km"].clip(upper=100.0)
    df["max_station_distance_capped_km"] = df["max_intensity_station_distance_km"].clip(upper=100.0)
    df["active_station_count_within_50km_capped"] = df["active_station_count_within_50km"].clip(upper=150.0)
    return df


def base_map(pygmt, title: str):
    fig = pygmt.Figure()
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_TICK_LENGTH_PRIMARY="2.5p",
        FORMAT_GEO_MAP="dddF",
        FONT_ANNOT_PRIMARY="8p,Helvetica,black",
        FONT_LABEL="9p,Helvetica,black",
        FONT_TITLE="12p,Helvetica-Bold,black",
        MAP_GRID_PEN_PRIMARY="0.15p,#d8d8d8",
        MAP_FRAME_PEN="0.7p,#222222",
    ):
        fig.basemap(
            region=MAP_REGION,
            projection=PROJECTION,
            frame=["xafg5+lLongitude", "yafg5+lLatitude", f"+t{title}"],
        )
        fig.coast(
            land="#f3f1ea",
            water="#e9f2f5",
            shorelines="0.35p,#555555",
            borders="1/0.2p,#888888",
            resolution="i",
        )
        plot_prefecture_boundaries_pygmt(fig, pen="0.16p,#666666")
    return fig


def add_map_scale(fig) -> None:
    fig.basemap(map_scale="jBL+w300k+o0.35c/0.35c+f+lkm")


def make_cpt(pygmt, cmap: str, series: list[float] | None, data: pd.Series | None = None, reverse: bool = False) -> None:
    if series is None:
        if data is None or data.dropna().empty:
            series = [0, 1, 0.1]
        else:
            low = float(data.quantile(0.02))
            high = float(data.quantile(0.98))
            if low == high:
                high = low + 1.0
            step = (high - low) / 8.0
            series = [low, high, step]
    pygmt.makecpt(cmap=cmap, series=series, reverse=reverse)


def save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=360)


def plot_period_epicenters(pygmt, df: pd.DataFrame, out_dir: Path) -> None:
    fig = base_map(pygmt, "M>=4 shallow inland epicenters by period")
    add_map_scale(fig)
    for label, color in PERIOD_COLORS.items():
        sub = df[df["analysis_period"] == label]
        if sub.empty:
            continue
        fig.plot(
            x=sub["longitude"],
            y=sub["latitude"],
            style="c0.075c",
            fill=color,
            pen="0.03p,white",
            transparency=28,
            label=f"{label} (n={len(sub):,})",
        )
    fig.legend(position="JTR+jTR+o0.15c/0.15c", box="+gwhite@18+p0.25p,#777777")
    save(fig, out_dir / "event_epicenters_by_period_m4plus.png")


def plot_point_metric(
    pygmt,
    df: pd.DataFrame,
    value_column: str,
    title: str,
    label: str,
    out_path: Path,
    cmap: str,
    series: list[float] | None,
    period: str | None = None,
    reverse: bool = False,
) -> None:
    data = df if period is None else df[df["analysis_period"] == period]
    data = data[data[value_column].notna()].copy()
    if data.empty:
        return
    fig = base_map(pygmt, title if period is None else f"{title}, {period}")
    add_map_scale(fig)
    make_cpt(pygmt, cmap=cmap, series=series, data=data[value_column], reverse=reverse)
    fig.plot(
        x=data["longitude"],
        y=data["latitude"],
        style="c0.075c",
        fill=data[value_column],
        cmap=True,
        pen="0.03p,white",
        transparency=20,
    )
    fig.colorbar(frame=[f"x+l{label}"])
    fig.text(
        x=MAP_REGION[0] + 0.65,
        y=MAP_REGION[3] - 0.75,
        text=f"N = {len(data):,}",
        font="10p,Helvetica-Bold,#222222",
        justify="LM",
        fill="white@18",
        clearance="0.08c/0.05c",
    )
    save(fig, out_path)


def bin_events(df: pd.DataFrame, grid_deg: float, min_events: int) -> pd.DataFrame:
    work = df.copy()
    work["grid_lon"] = ((work["longitude"] - MAP_REGION[0]) // grid_deg) * grid_deg + MAP_REGION[0] + grid_deg / 2.0
    work["grid_lat"] = ((work["latitude"] - MAP_REGION[2]) // grid_deg) * grid_deg + MAP_REGION[2] + grid_deg / 2.0
    rows = []
    for (period, lon, lat), group in work.groupby(["analysis_period", "grid_lon", "grid_lat"], observed=False):
        nearest = group["nearest_station_distance_km"]
        rows.append(
            {
                "period": period,
                "grid_lon": lon,
                "grid_lat": lat,
                "n_events": len(group),
                "mean_max_intensity": group["max_intensity_value"].mean(),
                "mean_station_intensity": group["mean_best_station_intensity"].mean(),
                "median_nearest_distance_km": nearest.median(),
                "mean_active_stations_50km": group["active_station_count_within_50km"].mean(),
                "mean_detected_stations_50km": group["detected_station_count_within_50km"].mean(),
                "share_nearest_within_20km": (nearest <= 20.0).mean(),
                "mean_station_records": group["n_valid_intensity"].mean(),
            }
        )
    out = pd.DataFrame(rows)
    for metric in GRID_METRICS:
        column = metric["column"]
        if column == "n_events":
            continue
        out.loc[out["n_events"] < max(min_events, metric["min_events"]), column] = pd.NA
    return out


def plot_grid_metric(pygmt, grid: pd.DataFrame, metric: dict[str, object], period: str, out_path: Path) -> None:
    column = str(metric["column"])
    data = grid[(grid["period"] == period) & grid[column].notna()].copy()
    if data.empty:
        return
    fig = base_map(pygmt, f"{metric['label']}, {period}")
    add_map_scale(fig)
    make_cpt(
        pygmt,
        cmap=str(metric["cmap"]),
        series=metric["series"],
        data=data[column],
        reverse=column == "median_nearest_distance_km",
    )
    fig.plot(
        x=data["grid_lon"],
        y=data["grid_lat"],
        style="s0.18c",
        fill=data[column],
        cmap=True,
        pen="0.03p,white",
        transparency=2,
    )
    fig.colorbar(frame=[f"x+l{metric['label']}"])
    fig.text(
        x=MAP_REGION[0] + 0.65,
        y=MAP_REGION[3] - 0.75,
        text=f"Cells = {len(data):,}",
        font="10p,Helvetica-Bold,#222222",
        justify="LM",
        fill="white@18",
        clearance="0.08c/0.05c",
    )
    save(fig, out_path)


def difference_grid(grid: pd.DataFrame, pre_period: str, post_period: str) -> pd.DataFrame:
    pre = grid[grid["period"] == pre_period].copy()
    post = grid[grid["period"] == post_period].copy()
    merged = pre.merge(post, on=["grid_lon", "grid_lat"], suffixes=("_pre", "_post"))
    if merged.empty:
        return merged
    merged["delta_mean_max_intensity"] = merged["mean_max_intensity_post"] - merged["mean_max_intensity_pre"]
    merged["delta_mean_station_intensity"] = merged["mean_station_intensity_post"] - merged["mean_station_intensity_pre"]
    merged["delta_median_nearest_distance_km"] = (
        merged["median_nearest_distance_km_post"] - merged["median_nearest_distance_km_pre"]
    )
    merged["delta_mean_active_stations_50km"] = (
        merged["mean_active_stations_50km_post"] - merged["mean_active_stations_50km_pre"]
    )
    merged["delta_share_nearest_within_20km"] = (
        merged["share_nearest_within_20km_post"] - merged["share_nearest_within_20km_pre"]
    )
    return merged


def plot_difference_metric(pygmt, diff: pd.DataFrame, column: str, label: str, series: list[float], out_path: Path) -> None:
    data = diff[diff[column].notna()].copy()
    if data.empty:
        return
    fig = base_map(pygmt, f"2011-2022 minus 1980-1995: {label}")
    add_map_scale(fig)
    make_cpt(pygmt, cmap="polar", series=series, data=data[column])
    fig.plot(
        x=data["grid_lon"],
        y=data["grid_lat"],
        style="s0.2c",
        fill=data[column],
        cmap=True,
        pen="0.03p,white",
    )
    fig.colorbar(frame=[f"x+l{label}"])
    save(fig, out_path)


def active_records_for_year(station_index, year: int) -> list[StationRecord]:
    return [
        rec
        for rec in active_station_records_at_year_end(station_index, year)
        if rec.latitude is not None
        and rec.longitude is not None
        and MAP_REGION[0] <= rec.longitude <= MAP_REGION[1]
        and MAP_REGION[2] <= rec.latitude <= MAP_REGION[3]
    ]


def bin_station_density(records: list[StationRecord], grid_deg: float) -> pd.DataFrame:
    rows = []
    for rec in records:
        lon = ((float(rec.longitude) - MAP_REGION[0]) // grid_deg) * grid_deg + MAP_REGION[0] + grid_deg / 2.0
        lat = ((float(rec.latitude) - MAP_REGION[2]) // grid_deg) * grid_deg + MAP_REGION[2] + grid_deg / 2.0
        rows.append({"grid_lon": lon, "grid_lat": lat, "region": station_network_region(rec)})
    if not rows:
        return pd.DataFrame(columns=["grid_lon", "grid_lat", "station_count"])
    return pd.DataFrame(rows).groupby(["grid_lon", "grid_lat"], as_index=False).size().rename(columns={"size": "station_count"})


def plot_station_density_grid(pygmt, stations: pd.DataFrame, year: int, out_path: Path) -> None:
    if stations.empty:
        return
    fig = base_map(pygmt, f"Active station density by grid cell, {year}")
    add_map_scale(fig)
    make_cpt(pygmt, cmap="batlow", series=[0, 40, 5], data=stations["station_count"])
    fig.plot(
        x=stations["grid_lon"],
        y=stations["grid_lat"],
        style="s0.18c",
        fill=stations["station_count"],
        cmap=True,
        pen="0.03p,white",
    )
    fig.colorbar(frame=["x+lActive stations per 0.5 deg cell"])
    save(fig, out_path)


def plot_station_start_periods(pygmt, station_index, year: int, out_path: Path) -> None:
    records = active_records_for_year(station_index, year)
    rows = []
    for rec in records:
        start = rec.start_key or 0
        start_year = int(str(start)[:4]) if start else 0
        if start_year <= 1995:
            label = "<=1995"
        elif start_year <= 2003:
            label = "1996-2003"
        elif start_year <= 2010:
            label = "2004-2010"
        else:
            label = "2011-2022"
        rows.append({"lon": rec.longitude, "lat": rec.latitude, "label": label})
    data = pd.DataFrame(rows)
    colors = {"<=1995": "#0072B2", "1996-2003": "#D55E00", "2004-2010": "#009E73", "2011-2022": "#CC79A7"}
    fig = base_map(pygmt, f"Active stations in {year} by installation period")
    add_map_scale(fig)
    for label, color in colors.items():
        sub = data[data["label"] == label]
        if sub.empty:
            continue
        fig.plot(
            x=sub["lon"],
            y=sub["lat"],
            style="c0.055c",
            fill=color,
            pen="0.02p,white",
            transparency=10,
            label=f"{label} (n={len(sub):,})",
        )
    fig.legend(position="JTR+jTR+o0.15c/0.15c", box="+gwhite@18+p0.25p,#777777")
    save(fig, out_path)


def main() -> None:
    configure_conda_gmt()
    import pygmt

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--event-csv", type=Path, default=Path("outputs/csv/event_intensity_summary.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/png/spatial_intensity_maps"))
    parser.add_argument("--csv-dir", type=Path, default=Path("outputs/csv"))
    parser.add_argument("--min-magnitude", type=float, default=4.0)
    parser.add_argument("--grid-deg", type=float, default=0.5)
    parser.add_argument("--min-cell-events", type=int, default=3)
    parser.add_argument("--station-years", default="1980,1995,2003,2022")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    df = load_filtered_events(args.event_csv, args.min_magnitude)
    print(f"Filtered events: {len(df):,}")

    plot_period_epicenters(pygmt, df, args.output_dir / "event_points")

    point_specs = [
        (
            "max_intensity",
            "max_intensity_value",
            "Observed maximum intensity",
            "Maximum intensity",
            "lajolla",
            [1, 7, 1],
            False,
        ),
        (
            "mean_station_intensity",
            "mean_best_station_intensity",
            "Event-average station intensity",
            "Mean station intensity",
            "lajolla",
            [0.5, 3.5, 0.5],
            False,
        ),
        (
            "shortest_epicentral_distance",
            "nearest_distance_capped_km",
            "Shortest detected epicentral distance",
            "Shortest distance [km], capped at 100",
            "devon",
            [0, 100, 10],
            True,
        ),
        (
            "max_intensity_station_distance",
            "max_station_distance_capped_km",
            "Distance to maximum-intensity station",
            "Distance [km], capped at 100",
            "devon",
            [0, 100, 10],
            True,
        ),
        (
            "active_station_support_50km",
            "active_station_count_within_50km_capped",
            "Station support within 50 km of epicenter",
            "Active stations within 50 km, capped at 150",
            "batlow",
            [0, 150, 25],
            False,
        ),
        (
            "detected_station_support_50km",
            "detected_station_count_within_50km",
            "Detected stations within 50 km of epicenter",
            "Detected stations within 50 km",
            "batlow",
            [0, 100, 20],
            False,
        ),
        (
            "detected_to_active_ratio_50km",
            "detected_to_active_station_ratio_within_50km",
            "Detected-to-active station ratio within 50 km",
            "Detected / active",
            "lapaz",
            [0, 1, 0.1],
            False,
        ),
    ]
    for name, column, title, label, cmap, series, reverse in point_specs:
        plot_point_metric(
            pygmt,
            df,
            column,
            f"M>=4 shallow inland events: {title}",
            label,
            args.output_dir / "event_points" / f"event_{name}_m4plus_all.png",
            cmap,
            series,
            reverse=reverse,
        )
        for period in [label_ for _, _, label_ in PERIODS]:
            plot_point_metric(
                pygmt,
                df,
                column,
                title,
                label,
                args.output_dir / "event_points" / f"event_{name}_{period.replace('-', '_')}.png",
                cmap,
                series,
                period=period,
                reverse=reverse,
            )

    grid = bin_events(df, args.grid_deg, args.min_cell_events)
    grid_csv = args.csv_dir / "spatial_intensity_grid_summary_land_approx_depth20_m4plus.csv"
    grid.to_csv(grid_csv, index=False)
    print(f"Grid summary -> {grid_csv}")
    for period in [label for _, _, label in PERIODS]:
        for metric in GRID_METRICS:
            plot_grid_metric(
                pygmt,
                grid,
                metric,
                period,
                args.output_dir / "grid_maps" / f"grid_{metric['name']}_{period.replace('-', '_')}.png",
            )

    diff = difference_grid(grid, "1980-1995", "2011-2022")
    diff_csv = args.csv_dir / "spatial_intensity_grid_difference_2011_2022_vs_1980_1995_land_approx_depth20_m4plus.csv"
    diff.to_csv(diff_csv, index=False)
    print(f"Grid difference -> {diff_csv}")
    for column, label, series in DIFFERENCE_METRICS:
        plot_difference_metric(
            pygmt,
            diff,
            column,
            label,
            series,
            args.output_dir / "difference_maps" / f"difference_{column}.png",
        )

    station_index = load_station_index(args.data_dir)
    if station_index is not None:
        station_years = [int(year.strip()) for year in args.station_years.split(",") if year.strip()]
        for year in station_years:
            records = active_records_for_year(station_index, year)
            density = bin_station_density(records, args.grid_deg)
            density.to_csv(args.csv_dir / f"station_density_grid_{year}.csv", index=False)
            plot_station_density_grid(
                pygmt,
                density,
                year,
                args.output_dir / "station_network" / f"station_density_grid_{year}.png",
            )
        plot_station_start_periods(
            pygmt,
            station_index,
            max(station_years),
            args.output_dir / "station_network" / f"station_start_periods_active_{max(station_years)}.png",
        )

    print(f"PNG maps written to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
