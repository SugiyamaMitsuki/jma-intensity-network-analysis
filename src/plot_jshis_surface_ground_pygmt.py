#!/usr/bin/env python3
"""Plot J-SHIS V4 AVS30 and site amplification contour maps with PyGMT."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from plot_station_maps_pygmt import configure_conda_gmt
from map_prefecture_boundaries import DEFAULT_PREFECTURE_BOUNDARY, plot_prefecture_boundaries_pygmt


DEFAULT_JSHIS_CSV = Path(
    "data(j-shis)/Z-V4-JAPAN-AMP-VS400_M250/Z-V4-JAPAN-AMP-VS400_M250.csv"
)
DEFAULT_EVENT_SUMMARY = Path("outputs/csv/hypocenter_catalog/jma_intensity_events_with_hypocenter.csv")
FALLBACK_EVENT_SUMMARY = Path("outputs/csv/event_intensity_summary.csv")
DEFAULT_CSV_DIR = Path("outputs/csv/jshis_surface_ground")
DEFAULT_PNG_DIR = Path("outputs/png/jshis_surface_ground")
DEFAULT_REGION = [122.0, 146.5, 24.0, 46.5]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create PyGMT contour maps of J-SHIS surface AVS30 and "
            "site amplification factor. The 250 m mesh is aggregated to a "
            "coarser regular grid for national-scale plotting."
        )
    )
    parser.add_argument("--jshis-csv", type=Path, default=DEFAULT_JSHIS_CSV)
    parser.add_argument("--event-summary", type=Path, default=DEFAULT_EVENT_SUMMARY)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument(
        "--region",
        type=float,
        nargs=4,
        default=DEFAULT_REGION,
        metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
    )
    parser.add_argument(
        "--spacing",
        type=float,
        default=0.02,
        help="Output grid spacing in degrees for plotting and contouring.",
    )
    parser.add_argument("--chunksize", type=int, default=700_000)
    parser.add_argument(
        "--min-intensity",
        type=float,
        default=6.5,
        help="Minimum max observed intensity value for target event overlay.",
    )
    parser.add_argument(
        "--no-events",
        action="store_true",
        help="Do not overlay JMA events whose observed intensity is 6 upper or larger.",
    )
    parser.add_argument("--prefecture-boundary", type=Path, default=DEFAULT_PREFECTURE_BOUNDARY)
    return parser.parse_args()


def read_jshis_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("#"):
                break
            text = line[1:].strip()
            if "=" in text:
                key, value = [part.strip() for part in text.split("=", 1)]
                metadata[key] = value
            elif text and "," in text:
                metadata["columns"] = text
    return metadata


def decode_250m_mesh_center(codes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return center latitude/longitude from 10-digit J-SHIS 250 m mesh codes."""

    code = codes.astype(np.int64, copy=False)
    lat = (
        (code // 100_000_000) / 1.5
        + ((code // 100_000) % 10) / 12.0
        + ((code // 1_000) % 10) / 120.0
        + (((code // 10) % 10) - 1) * (7.5 / 3600.0)
        + (3.75 / 3600.0)
    )
    lon = (
        ((code // 1_000_000) % 100)
        + 100.0
        + ((code // 10_000) % 10) / 8.0
        + ((code // 100) % 10) / 80.0
        + ((code % 10) - 1) * (11.25 / 3600.0)
        + (5.625 / 3600.0)
    )
    return lat, lon


def init_accumulators(region: list[float], spacing: float) -> tuple[np.ndarray, ...]:
    lon_min, lon_max, lat_min, lat_max = region
    nx = int(math.ceil((lon_max - lon_min) / spacing))
    ny = int(math.ceil((lat_max - lat_min) / spacing))
    shape = ny * nx
    return (
        np.zeros(shape, dtype=np.float64),
        np.zeros(shape, dtype=np.float64),
        np.zeros(shape, dtype=np.int32),
        np.zeros(shape, dtype=np.int32),
        np.zeros(shape, dtype=np.int32),
        np.array([nx, ny], dtype=np.int32),
    )


def aggregate_jshis_grid(
    path: Path,
    region: list[float],
    spacing: float,
    chunksize: int,
) -> tuple[pd.DataFrame, dict[str, float | int]]:
    names = ["CODE", "JCODE", "AVS", "ARV", "AVS_EB", "AVS_REF"]
    dtypes = {
        "CODE": "int64",
        "JCODE": "int16",
        "AVS": "float32",
        "ARV": "float32",
        "AVS_EB": "string",
        "AVS_REF": "int8",
    }
    lon_min, lon_max, lat_min, lat_max = region
    sum_avs, sum_arv, counts, jcode_water, jcode_lowland, dims = init_accumulators(region, spacing)
    nx, ny = int(dims[0]), int(dims[1])

    total_rows = 0
    valid_rows = 0
    in_region_rows = 0
    avs_values: list[np.ndarray] = []
    arv_values: list[np.ndarray] = []

    for chunk in pd.read_csv(
        path,
        comment="#",
        names=names,
        skipinitialspace=True,
        dtype=dtypes,
        chunksize=chunksize,
    ):
        total_rows += len(chunk)
        valid = (chunk["JCODE"] > 0) & (chunk["AVS"] > 0) & (chunk["ARV"] > 0)
        if not valid.any():
            continue

        data = chunk.loc[valid, ["CODE", "JCODE", "AVS", "ARV"]]
        valid_rows += len(data)
        lat, lon = decode_250m_mesh_center(data["CODE"].to_numpy(np.int64))
        in_region = (
            (lon >= lon_min)
            & (lon < lon_max)
            & (lat >= lat_min)
            & (lat < lat_max)
        )
        if not np.any(in_region):
            continue

        in_region_rows += int(in_region.sum())
        lon = lon[in_region]
        lat = lat[in_region]
        avs = data["AVS"].to_numpy(np.float64)[in_region]
        arv = data["ARV"].to_numpy(np.float64)[in_region]
        jcode = data["JCODE"].to_numpy(np.int16)[in_region]

        ix = np.floor((lon - lon_min) / spacing).astype(np.int64)
        iy = np.floor((lat - lat_min) / spacing).astype(np.int64)
        inside = (0 <= ix) & (ix < nx) & (0 <= iy) & (iy < ny)
        flat = iy[inside] * nx + ix[inside]
        avs_inside = avs[inside]
        arv_inside = arv[inside]
        jcode_inside = jcode[inside]

        sum_avs += np.bincount(flat, weights=avs_inside, minlength=nx * ny)
        sum_arv += np.bincount(flat, weights=arv_inside, minlength=nx * ny)
        counts += np.bincount(flat, minlength=nx * ny).astype(np.int32)
        jcode_water += np.bincount(
            flat, weights=(jcode_inside == 24).astype(np.int8), minlength=nx * ny
        ).astype(np.int32)
        jcode_lowland += np.bincount(
            flat, weights=((10 <= jcode_inside) & (jcode_inside <= 20)).astype(np.int8),
            minlength=nx * ny,
        ).astype(np.int32)

        avs_values.append(avs_inside.astype(np.float32, copy=False))
        arv_values.append(arv_inside.astype(np.float32, copy=False))

    mean_avs = np.full(nx * ny, np.nan, dtype=np.float32)
    mean_arv = np.full(nx * ny, np.nan, dtype=np.float32)
    nonzero = counts > 0
    mean_avs[nonzero] = (sum_avs[nonzero] / counts[nonzero]).astype(np.float32)
    mean_arv[nonzero] = (sum_arv[nonzero] / counts[nonzero]).astype(np.float32)

    ys, xs = np.divmod(np.arange(nx * ny), nx)
    lon_center = lon_min + (xs + 0.5) * spacing
    lat_center = lat_min + (ys + 0.5) * spacing
    grid_df = pd.DataFrame(
        {
            "longitude": lon_center[nonzero],
            "latitude": lat_center[nonzero],
            "avs30_m_s": mean_avs[nonzero],
            "amplification_vs400": mean_arv[nonzero],
            "n_250m_mesh": counts[nonzero],
            "n_waterbody_jcode24": jcode_water[nonzero],
            "n_lowland_jcode10_20": jcode_lowland[nonzero],
        }
    )

    raw_avs = np.concatenate(avs_values) if avs_values else np.array([], dtype=np.float32)
    raw_arv = np.concatenate(arv_values) if arv_values else np.array([], dtype=np.float32)
    stats: dict[str, float | int] = {
        "input_rows": int(total_rows),
        "valid_surface_mesh_rows": int(valid_rows),
        "valid_surface_mesh_rows_in_region": int(in_region_rows),
        "output_grid_spacing_deg": float(spacing),
        "output_grid_cells_with_data": int(nonzero.sum()),
        "output_grid_total_cells": int(nx * ny),
        "region_lon_min": float(lon_min),
        "region_lon_max": float(lon_max),
        "region_lat_min": float(lat_min),
        "region_lat_max": float(lat_max),
    }
    if raw_avs.size:
        for prefix, values in [("avs30_m_s", raw_avs), ("amplification_vs400", raw_arv)]:
            stats[f"{prefix}_min"] = float(np.nanmin(values))
            stats[f"{prefix}_p01"] = float(np.nanquantile(values, 0.01))
            stats[f"{prefix}_p05"] = float(np.nanquantile(values, 0.05))
            stats[f"{prefix}_median"] = float(np.nanquantile(values, 0.50))
            stats[f"{prefix}_p95"] = float(np.nanquantile(values, 0.95))
            stats[f"{prefix}_p99"] = float(np.nanquantile(values, 0.99))
            stats[f"{prefix}_max"] = float(np.nanmax(values))
    return grid_df, stats


def grid_to_dataarray(
    grid_df: pd.DataFrame,
    region: list[float],
    spacing: float,
    column: str,
    name: str,
    attrs: dict[str, str],
) -> xr.DataArray:
    lon_min, lon_max, lat_min, lat_max = region
    nx = int(math.ceil((lon_max - lon_min) / spacing))
    ny = int(math.ceil((lat_max - lat_min) / spacing))
    data = np.full((ny, nx), np.nan, dtype=np.float32)
    ix = np.floor((grid_df["longitude"].to_numpy() - lon_min) / spacing).astype(int)
    iy = np.floor((grid_df["latitude"].to_numpy() - lat_min) / spacing).astype(int)
    data[iy, ix] = grid_df[column].to_numpy(np.float32)
    lon = lon_min + (np.arange(nx) + 0.5) * spacing
    lat = lat_min + (np.arange(ny) + 0.5) * spacing
    return xr.DataArray(data, coords={"lat": lat, "lon": lon}, dims=("lat", "lon"), name=name, attrs=attrs)


def rounded_series(values: np.ndarray, step: float, low_q: float, high_q: float) -> list[float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return [0.0, 1.0, step]
    lo = math.floor(float(np.nanquantile(finite, low_q)) / step) * step
    hi = math.ceil(float(np.nanquantile(finite, high_q)) / step) * step
    if lo >= hi:
        lo = math.floor(float(np.nanmin(finite)) / step) * step
        hi = math.ceil(float(np.nanmax(finite)) / step) * step
    return [lo, hi, step]


def load_target_events(path: Path, min_intensity: float, region: list[float]) -> pd.DataFrame:
    if not path.exists():
        path = FALLBACK_EVENT_SUMMARY
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, low_memory=False)
    lat_col = "analysis_latitude" if "analysis_latitude" in df.columns else "latitude"
    lon_col = "analysis_longitude" if "analysis_longitude" in df.columns else "longitude"
    required = {lat_col, lon_col, "max_intensity_value"}
    if not required.issubset(df.columns):
        return pd.DataFrame()
    lon_min, lon_max, lat_min, lat_max = region
    target = df[
        (df["max_intensity_value"] >= min_intensity)
        & df[lon_col].between(lon_min, lon_max)
        & df[lat_col].between(lat_min, lat_max)
    ].copy()
    target["plot_latitude"] = target[lat_col]
    target["plot_longitude"] = target[lon_col]
    cols = [
        "event_id",
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second",
        "plot_latitude",
        "plot_longitude",
        "latitude",
        "longitude",
        "hyp_latitude",
        "hyp_longitude",
        "analysis_coordinate_source",
        "depth_km",
        "analysis_depth_km",
        "magnitude",
        "analysis_magnitude",
        "max_intensity_value",
        "n_valid_intensity",
        "n_unique_stations",
        "hypocenter_region",
    ]
    return target[[col for col in cols if col in target.columns]].sort_values(
        ["year", "month", "day", "hour", "minute", "second"]
    )


def draw_contour_map(
    pygmt,
    grid: xr.DataArray,
    output_path: Path,
    region: list[float],
    title: str,
    cmap: str,
    cpt_series: list[float],
    colorbar_label: str,
    contour_interval: float,
    contour_annotation: float,
    target_events: pd.DataFrame | None,
    prefecture_boundary: Path | None = DEFAULT_PREFECTURE_BOUNDARY,
) -> None:
    projection = "B134/35/25/47/15c"
    fig = pygmt.Figure()
    pygmt.makecpt(cmap=cmap, series=cpt_series, continuous=True)

    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN="0.7p,#202020",
        MAP_TICK_LENGTH_PRIMARY="2.5p",
        MAP_GRID_PEN_PRIMARY="0.12p,#d4d4d4",
        FORMAT_GEO_MAP="dddF",
        FONT_ANNOT_PRIMARY="7.5p,Helvetica,#202020",
        FONT_LABEL="8.5p,Helvetica,#202020",
        FONT_TITLE="12p,Helvetica-Bold,#111111",
        MAP_TITLE_OFFSET="6p",
        COLOR_NAN="white",
    ):
        fig.basemap(
            region=region,
            projection=projection,
            frame=["xafg5+lLongitude", "yafg5+lLatitude", f"+t{title}"],
        )
        fig.coast(
            region=region,
            projection=projection,
            water="#eef4f6",
            land="#fbfbf8",
            shorelines="0.3p,#666666",
            resolution="i",
        )
        fig.grdimage(
            grid=grid,
            region=region,
            projection=projection,
            cmap=True,
            nan_transparent=True,
        )
        fig.grdcontour(
            grid=grid,
            region=region,
            projection=projection,
            levels=contour_interval,
            annotation=f"{contour_annotation}+f6.5p,Helvetica,#222222",
            pen="0.16p,white@35",
        )
        fig.coast(
            region=region,
            projection=projection,
            shorelines="0.35p,#333333",
            borders="1/0.18p,#777777",
            resolution="i",
            map_scale="jBL+w300k+o0.35c/0.35c+f+lkm",
        )
        plot_prefecture_boundaries_pygmt(fig, prefecture_boundary, pen="0.16p,#5f5f5f")
        if target_events is not None and not target_events.empty:
            fig.plot(
                x=target_events["plot_longitude"] if "plot_longitude" in target_events else target_events["longitude"],
                y=target_events["plot_latitude"] if "plot_latitude" in target_events else target_events["latitude"],
                style="a0.20c",
                fill="#111111",
                pen="0.20p,white",
            )
            fig.text(
                x=region[0] + 0.55,
                y=region[3] - 0.55,
                text=f"Stars: JMA events with Imax >= 6-upper ({len(target_events)})",
                font="7p,Helvetica,#111111",
                justify="TL",
                fill="white@20",
                clearance="2p/2p",
            )
        fig.colorbar(
            position="JBC+w10.5c/0.28c+o0c/0.65c+h",
            frame=[f"x+l{colorbar_label}", "y"],
        )
    fig.savefig(output_path, crop=True, dpi=300)


def write_summary(
    csv_dir: Path,
    metadata: dict[str, str],
    stats: dict[str, float | int],
    source_path: Path,
    target_events: pd.DataFrame,
    min_intensity: float,
) -> Path:
    rows: list[dict[str, str | float | int]] = [
        {"key": "source_file", "value": str(source_path)},
    ]
    rows.extend({"key": f"jshis_{key.lower().replace(' ', '_')}", "value": value} for key, value in metadata.items())
    rows.extend({"key": key, "value": value} for key, value in stats.items())
    rows.append({"key": "target_event_min_intensity_value", "value": min_intensity})
    rows.append({"key": "target_event_count", "value": int(len(target_events))})
    out = csv_dir / "jshis_surface_ground_processing_summary.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    return out


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_jshis_metadata(args.jshis_csv)

    grid_df, stats = aggregate_jshis_grid(
        args.jshis_csv,
        args.region,
        args.spacing,
        args.chunksize,
    )
    grid_csv = args.csv_dir / f"jshis_surface_ground_grid_{str(args.spacing).replace('.', 'p')}deg.csv"
    grid_df.to_csv(grid_csv, index=False)

    target_events = pd.DataFrame()
    if not args.no_events:
        target_events = load_target_events(args.event_summary, args.min_intensity, args.region)
        target_csv = args.csv_dir / "target_events_intensity_6upper_plus.csv"
        target_events.to_csv(target_csv, index=False)
    else:
        target_csv = args.csv_dir / "target_events_intensity_6upper_plus.csv"

    summary_csv = write_summary(
        args.csv_dir,
        metadata,
        stats,
        args.jshis_csv,
        target_events,
        args.min_intensity,
    )

    avs_grid = grid_to_dataarray(
        grid_df,
        args.region,
        args.spacing,
        "avs30_m_s",
        "AVS30",
        {"long_name": "Average S-wave velocity in the upper 30 m", "units": "m/s"},
    )
    arv_grid = grid_to_dataarray(
        grid_df,
        args.region,
        args.spacing,
        "amplification_vs400",
        "ARV",
        {"long_name": "Site amplification factor from Vs=400 m/s to surface", "units": "ratio"},
    )

    configure_conda_gmt()
    import pygmt

    amp_series = rounded_series(arv_grid.values, step=0.1, low_q=0.01, high_q=0.995)
    avs_series = rounded_series(avs_grid.values, step=50.0, low_q=0.005, high_q=0.99)

    amp_png = args.png_dir / "jshis_surface_amplification_vs400_contour_map.png"
    avs_png = args.png_dir / "jshis_avs30_contour_map.png"
    draw_contour_map(
        pygmt=pygmt,
        grid=arv_grid,
        output_path=amp_png,
        region=args.region,
        title="J-SHIS V4 Surface Amplification, Vs=400 m/s to Surface",
        cmap="batlow",
        cpt_series=amp_series,
        colorbar_label="Site amplification factor",
        contour_interval=0.2,
        contour_annotation=0.4,
        target_events=target_events,
        prefecture_boundary=args.prefecture_boundary,
    )
    draw_contour_map(
        pygmt=pygmt,
        grid=avs_grid,
        output_path=avs_png,
        region=args.region,
        title="J-SHIS V4 AVS30",
        cmap="viridis",
        cpt_series=avs_series,
        colorbar_label="AVS30 (m/s)",
        contour_interval=100.0,
        contour_annotation=200.0,
        target_events=target_events,
        prefecture_boundary=args.prefecture_boundary,
    )

    print(f"Saved grid CSV: {grid_csv}")
    print(f"Saved summary CSV: {summary_csv}")
    if not args.no_events:
        print(f"Saved target event CSV: {target_csv}")
    print(f"Saved PNG: {amp_png}")
    print(f"Saved PNG: {avs_png}")


if __name__ == "__main__":
    main()
