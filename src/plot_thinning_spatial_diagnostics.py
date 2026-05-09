#!/usr/bin/env python3
"""Plot spatial diagnostics for station-thinning interpolation tests."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from map_prefecture_boundaries import DEFAULT_PREFECTURE_BOUNDARY, plot_prefecture_boundaries_cartopy


DEFAULT_CSV_DIR = Path("outputs/csv/station_thinning_interpolation")
DEFAULT_PNG_DIR = Path("outputs/png/station_thinning_interpolation")
DEFAULT_TARGET_CSV = Path("outputs/csv/hypocenter_catalog/target_events_intensity_6upper_plus_with_hypocenter.csv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create map panels showing held-out station errors and local ground-complexity diagnostics."
    )
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--target-events-csv", type=Path, default=DEFAULT_TARGET_CSV)
    parser.add_argument("--events", default="i1995_000072,i2004_001047,i2016_000476")
    parser.add_argument("--method", default="idw")
    parser.add_argument("--keep-fractions", default="0.1,0.5")
    parser.add_argument("--random-index", type=int, default=0)
    parser.add_argument("--prefecture-boundary", type=Path, default=DEFAULT_PREFECTURE_BOUNDARY)
    return parser.parse_args()


def parse_list(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def event_metadata(path: Path) -> pd.DataFrame:
    cols = [
        "event_id",
        "intensity_origin_time",
        "hyp_region_name",
        "analysis_latitude",
        "analysis_longitude",
        "analysis_magnitude",
    ]
    df = pd.read_csv(path, usecols=lambda c: c in cols)
    df["event_label"] = (
        df["intensity_origin_time"].astype(str).str.slice(0, 10)
        + " "
        + df["hyp_region_name"].astype(str)
        + " M"
        + df["analysis_magnitude"].round(1).astype(str)
    )
    return df


def panel_extent(stations: pd.DataFrame, event_row: pd.Series) -> tuple[float, float, float, float]:
    lon = stations["longitude"].to_numpy(dtype=float)
    lat = stations["latitude"].to_numpy(dtype=float)
    lon = lon[np.isfinite(lon)]
    lat = lat[np.isfinite(lat)]
    lon_min, lon_max = np.nanquantile(lon, [0.02, 0.98])
    lat_min, lat_max = np.nanquantile(lat, [0.02, 0.98])
    lon_min = min(lon_min, float(event_row["analysis_longitude"]))
    lon_max = max(lon_max, float(event_row["analysis_longitude"]))
    lat_min = min(lat_min, float(event_row["analysis_latitude"]))
    lat_max = max(lat_max, float(event_row["analysis_latitude"]))
    pad_lon = max(0.35, (lon_max - lon_min) * 0.08)
    pad_lat = max(0.30, (lat_max - lat_min) * 0.08)
    return lon_min - pad_lon, lon_max + pad_lon, lat_min - pad_lat, lat_max + pad_lat


def setup_map(ax, extent: tuple[float, float, float, float], title: str, prefecture_boundary: Path | None) -> None:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#f3f1eb", edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#eef4f6", edgecolor="none", zorder=0)
    ax.coastlines(resolution="10m", color="#555555", linewidth=0.45, zorder=1)
    plot_prefecture_boundaries_cartopy(ax, prefecture_boundary, color="#686868", linewidth=0.28, alpha=0.82, zorder=2)
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.35,
        color="#b9b9b9",
        alpha=0.7,
        linestyle="-",
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 6, "color": "#444444"}
    gl.ylabel_style = {"size": 6, "color": "#444444"}
    ax.set_title(title, loc="left", fontsize=8.5, pad=3)


def train_subset(stations: pd.DataFrame, test: pd.DataFrame) -> pd.DataFrame:
    test_codes = set(test["station_code"].astype(str))
    return stations[~stations["station_code"].astype(str).isin(test_codes)].copy()


def plot_error_maps(args: argparse.Namespace, events: list[str], keep_fractions: list[float]) -> Path:
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    args.png_dir.mkdir(parents=True, exist_ok=True)
    pred_cols = [
        "station_code",
        "station_name",
        "latitude",
        "longitude",
        "observed_intensity",
        "event_id",
        "method",
        "keep_fraction",
        "random_index",
        "n_train",
        "train_density_per_10000km2",
        "predicted_intensity",
        "prediction_error",
        "abs_error",
    ]
    pred = pd.read_csv(args.csv_dir / "thinning_prediction_errors.csv.gz", usecols=pred_cols)
    pred = pred[
        pred["event_id"].isin(events)
        & (pred["method"] == args.method)
        & pred["keep_fraction"].isin(keep_fractions)
        & (pred["random_index"] == args.random_index)
    ].copy()
    meta = event_metadata(args.target_events_csv).set_index("event_id")
    fig, axes = plt.subplots(
        len(keep_fractions),
        len(events),
        figsize=(4.0 * len(events), 3.5 * len(keep_fractions)),
        subplot_kw={"projection": ccrs.PlateCarree()},
        squeeze=False,
    )
    scatter_for_cbar = None
    for col, event_id in enumerate(events):
        stations = pd.read_csv(args.csv_dir / "stations" / f"{event_id}_station_observations_ground_complexity.csv")
        event_row = meta.loc[event_id]
        extent = panel_extent(stations, event_row)
        for row, keep_fraction in enumerate(keep_fractions):
            ax = axes[row, col]
            test = pred[(pred["event_id"] == event_id) & (pred["keep_fraction"] == keep_fraction)].copy()
            train = train_subset(stations, test)
            title = (
                f"{event_row['event_label']} | keep={keep_fraction:g}, "
                f"n={int(test['n_train'].iloc[0]) if len(test) else 0}"
            )
            setup_map(ax, extent, title, args.prefecture_boundary)
            ax.scatter(
                train["longitude"],
                train["latitude"],
                s=4,
                c="#1f1f1f",
                alpha=0.45,
                linewidths=0,
                transform=ccrs.PlateCarree(),
                zorder=2,
            )
            scatter_for_cbar = ax.scatter(
                test["longitude"],
                test["latitude"],
                s=15,
                c=test["prediction_error"],
                cmap="RdBu_r",
                vmin=-1.5,
                vmax=1.5,
                edgecolors="white",
                linewidths=0.15,
                transform=ccrs.PlateCarree(),
                zorder=3,
            )
            ax.scatter(
                [event_row["analysis_longitude"]],
                [event_row["analysis_latitude"]],
                marker="*",
                s=70,
                c="#111111",
                edgecolors="white",
                linewidths=0.45,
                transform=ccrs.PlateCarree(),
                zorder=4,
            )
    fig.subplots_adjust(left=0.035, right=0.965, top=0.93, bottom=0.10, wspace=0.11, hspace=0.15)
    if scatter_for_cbar is not None:
        cbar = fig.colorbar(scatter_for_cbar, ax=axes.ravel().tolist(), orientation="horizontal", fraction=0.035, pad=0.045)
        cbar.set_label("Prediction error at held-out station (predicted - observed intensity)", fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    fig.suptitle("Station thinning diagnostics: spatial distribution of held-out errors", fontsize=11, y=0.985)
    out = args.png_dir / f"thinning_spatial_error_maps_{args.method}.png"
    fig.savefig(out, dpi=320)
    plt.close(fig)
    return out


def plot_complexity_maps(args: argparse.Namespace, events: list[str]) -> Path:
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    meta = event_metadata(args.target_events_csv).set_index("event_id")
    fig, axes = plt.subplots(
        1,
        len(events),
        figsize=(4.0 * len(events), 3.5),
        subplot_kw={"projection": ccrs.PlateCarree()},
        squeeze=False,
    )
    scatter_for_cbar = None
    for col, event_id in enumerate(events):
        stations = pd.read_csv(args.csv_dir / "stations" / f"{event_id}_station_observations_ground_complexity.csv")
        event_row = meta.loc[event_id]
        complexity_col = [c for c in stations.columns if c.startswith("local_site_delta_std_")][0]
        extent = panel_extent(stations, event_row)
        ax = axes[0, col]
        setup_map(ax, extent, event_row["event_label"], args.prefecture_boundary)
        scatter_for_cbar = ax.scatter(
            stations["longitude"],
            stations["latitude"],
            s=8,
            c=stations[complexity_col],
            cmap="viridis",
            vmin=0.05,
            vmax=0.45,
            edgecolors="none",
            transform=ccrs.PlateCarree(),
            zorder=3,
        )
        ax.scatter(
            [event_row["analysis_longitude"]],
            [event_row["analysis_latitude"]],
            marker="*",
            s=80,
            c="#111111",
            edgecolors="white",
            linewidths=0.45,
            transform=ccrs.PlateCarree(),
            zorder=4,
        )
    fig.subplots_adjust(left=0.035, right=0.965, top=0.88, bottom=0.17, wspace=0.11)
    if scatter_for_cbar is not None:
        cbar = fig.colorbar(scatter_for_cbar, ax=axes.ravel().tolist(), orientation="horizontal", fraction=0.06, pad=0.08)
        cbar.set_label("Local std. of site-intensity correction within 20 km", fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    fig.suptitle("Local J-SHIS ground-complexity around intensity stations", fontsize=11, y=0.98)
    out = args.png_dir / "station_local_ground_complexity_maps.png"
    fig.savefig(out, dpi=320)
    plt.close(fig)
    return out


def main() -> None:
    args = parse_args()
    events = parse_list(args.events)
    keep_fractions = [float(v) for v in parse_list(args.keep_fractions)]
    err_map = plot_error_maps(args, events, keep_fractions)
    complexity_map = plot_complexity_maps(args, events)
    print(f"Saved {err_map}")
    print(f"Saved {complexity_map}")


if __name__ == "__main__":
    main()
