#!/usr/bin/env python3
"""Map local J-SHIS ground-amplification variability and related diagnostics."""

from __future__ import annotations

import argparse
import gzip
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy import ndimage

from map_prefecture_boundaries import DEFAULT_PREFECTURE_BOUNDARY, plot_prefecture_boundaries_pygmt
from plot_jshis_surface_ground_pygmt import rounded_series
from plot_station_maps_pygmt import configure_conda_gmt


DEFAULT_GRID_CSV = Path("outputs/csv/jshis_surface_ground/jshis_surface_ground_grid_0p02deg.csv")
DEFAULT_REGION = [122.0, 146.5, 24.0, 46.5]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-csv", type=Path, default=DEFAULT_GRID_CSV)
    parser.add_argument("--csv-dir", type=Path, default=Path("outputs/csv/jshis_ground_variability"))
    parser.add_argument("--png-dir", type=Path, default=Path("outputs/png/jshis_ground_variability"))
    parser.add_argument("--derived-csv-gz", type=Path, default=Path("data/derived/jshis_ground_variability_metrics_0p02deg.csv.gz"))
    parser.add_argument("--region", nargs=4, type=float, default=DEFAULT_REGION)
    parser.add_argument("--spacing", type=float, default=0.02)
    parser.add_argument("--radius-km", type=float, default=20.0)
    parser.add_argument("--amp-intensity-coef", type=float, default=1.72)
    parser.add_argument("--prefecture-boundary", type=Path, default=DEFAULT_PREFECTURE_BOUNDARY)
    parser.add_argument("--assets-ja", type=Path, default=Path("paper/assets_ja/fig04c_jshis_ground_variability_metrics.png"))
    parser.add_argument("--assets-en", type=Path, default=Path("paper/assets_en/fig04c_jshis_ground_variability_metrics.png"))
    return parser.parse_args()


def grid_shape(region: list[float], spacing: float) -> tuple[int, int]:
    lon_min, lon_max, lat_min, lat_max = region
    nx = int(math.ceil((lon_max - lon_min) / spacing))
    ny = int(math.ceil((lat_max - lat_min) / spacing))
    return ny, nx


def dataframe_to_array(df: pd.DataFrame, region: list[float], spacing: float, column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon_min, lon_max, lat_min, lat_max = region
    ny, nx = grid_shape(region, spacing)
    arr = np.full((ny, nx), np.nan, dtype=np.float32)
    ix = np.floor((df["longitude"].to_numpy() - lon_min) / spacing).astype(int)
    iy = np.floor((df["latitude"].to_numpy() - lat_min) / spacing).astype(int)
    inside = (0 <= ix) & (ix < nx) & (0 <= iy) & (iy < ny)
    arr[iy[inside], ix[inside]] = df[column].to_numpy(np.float32)[inside]
    lon = lon_min + (np.arange(nx) + 0.5) * spacing
    lat = lat_min + (np.arange(ny) + 0.5) * spacing
    return arr, lat, lon


def circular_kernel(radius_cells: int) -> np.ndarray:
    y, x = np.ogrid[-radius_cells : radius_cells + 1, -radius_cells : radius_cells + 1]
    kernel = (x * x + y * y) <= radius_cells * radius_cells
    return kernel.astype(np.float32)


def local_mean_std(values: np.ndarray, kernel: np.ndarray, min_count: int = 5) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    valid = np.isfinite(values).astype(np.float32)
    clean = np.nan_to_num(values, nan=0.0).astype(np.float32)
    count = ndimage.convolve(valid, kernel, mode="constant", cval=0.0)
    total = ndimage.convolve(clean, kernel, mode="constant", cval=0.0)
    total2 = ndimage.convolve(clean * clean, kernel, mode="constant", cval=0.0)
    mean = np.full(values.shape, np.nan, dtype=np.float32)
    std = np.full(values.shape, np.nan, dtype=np.float32)
    ok = count >= min_count
    mean[ok] = total[ok] / count[ok]
    variance = np.maximum(total2[ok] / count[ok] - mean[ok] * mean[ok], 0.0)
    std[ok] = np.sqrt(variance)
    return mean, std, count.astype(np.float32)


def gradient_per_10km(values: np.ndarray, lat: np.ndarray, spacing: float) -> np.ndarray:
    out = np.full(values.shape, np.nan, dtype=np.float32)
    north = values[2:, 1:-1]
    south = values[:-2, 1:-1]
    east = values[1:-1, 2:]
    west = values[1:-1, :-2]
    center_lat = lat[1:-1]
    dy = 2.0 * spacing * 111.32
    dx = 2.0 * spacing * 111.32 * np.cos(np.deg2rad(center_lat))[:, None]
    valid = np.isfinite(north) & np.isfinite(south) & np.isfinite(east) & np.isfinite(west)
    gy = (north - south) / dy
    gx = (east - west) / dx
    grad = np.sqrt(gx * gx + gy * gy) * 10.0
    out[1:-1, 1:-1][valid] = grad[valid]
    return out


def to_dataarray(values: np.ndarray, lat: np.ndarray, lon: np.ndarray, name: str, units: str) -> xr.DataArray:
    return xr.DataArray(
        values.astype(np.float32),
        coords={"lat": lat, "lon": lon},
        dims=("lat", "lon"),
        name=name,
        attrs={"units": units},
    )


def flatten_metrics(lat: np.ndarray, lon: np.ndarray, arrays: dict[str, np.ndarray]) -> pd.DataFrame:
    lon2, lat2 = np.meshgrid(lon, lat)
    mask = np.zeros_like(lat2, dtype=bool)
    for values in arrays.values():
        mask |= np.isfinite(values)
    data = {
        "longitude": lon2[mask],
        "latitude": lat2[mask],
    }
    for name, values in arrays.items():
        data[name] = values[mask]
    return pd.DataFrame(data)


def write_gzip_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
        df.to_csv(f, index=False, float_format="%.6g")


def panel(
    fig,
    pygmt,
    grid: xr.DataArray,
    *,
    region: list[float],
    title: str,
    cmap: str,
    series: list[float],
    colorbar_label: str,
    prefecture_boundary: Path | None,
) -> None:
    projection = "M7.2c"
    pygmt.makecpt(cmap=cmap, series=series, continuous=True)
    fig.basemap(region=region, projection=projection, frame=["xafg5", "yafg5", f"+t{title}"])
    fig.coast(region=region, projection=projection, water="#eef4f6", land="#fbfbf8", shorelines="0.24p,#666666", resolution="i")
    fig.grdimage(grid=grid, region=region, projection=projection, cmap=True, nan_transparent=True)
    fig.coast(region=region, projection=projection, shorelines="0.28p,#333333", borders="1/0.15p,#777777", resolution="i")
    plot_prefecture_boundaries_pygmt(fig, prefecture_boundary, pen="0.12p,#5f5f5f")
    fig.colorbar(position="JMR+w4.4c/0.16c+o0.28c/0c+v", frame=[f"x+l{colorbar_label}", "y"])


def draw_maps(
    pygmt,
    grids: dict[str, xr.DataArray],
    output_path: Path,
    region: list[float],
    prefecture_boundary: Path | None,
) -> None:
    finite_delta = grids["local_site_delta_std_20km"].values
    finite_cv = grids["local_amp_cv_20km"].values
    finite_avs = grids["local_avs30_std_20km"].values
    finite_grad = grids["site_delta_gradient_per_10km"].values

    specs = [
        (
            "local_site_delta_std_20km",
            "a  Site-term std. (20 km)",
            "batlow",
            rounded_series(finite_delta, 0.02, 0.01, 0.995),
            "std dIsite",
        ),
        (
            "local_amp_cv_20km",
            "b  Amplification CV (20 km)",
            "lajolla",
            rounded_series(finite_cv, 0.02, 0.01, 0.995),
            "CV(A)",
        ),
        (
            "local_avs30_std_20km",
            "c  AVS30 std. (20 km)",
            "devon",
            rounded_series(finite_avs, 20.0, 0.01, 0.995),
            "std AVS30 (m/s)",
        ),
        (
            "site_delta_gradient_per_10km",
            "d  Site-term gradient",
            "roma",
            rounded_series(finite_grad, 0.02, 0.01, 0.995),
            "|grad dIsite| / 10 km",
        ),
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN="0.55p,#202020",
        MAP_TICK_LENGTH_PRIMARY="2.0p",
        MAP_GRID_PEN_PRIMARY="0.10p,#d5d5d5",
        FORMAT_GEO_MAP="dddF",
        FONT_ANNOT_PRIMARY="6.0p,Helvetica,#202020",
        FONT_LABEL="6.2p,Helvetica,#202020",
        FONT_TITLE="7.2p,Helvetica-Bold,#111111",
        MAP_TITLE_OFFSET="3p",
        COLOR_NAN="white",
    ):
        fig = pygmt.Figure()
        with fig.subplot(nrows=2, ncols=2, figsize=("19.0c", "16.3c"), margins=["1.18c", "1.08c"]):
            for idx, (key, title, cmap, series, colorbar_label) in enumerate(specs):
                with fig.set_panel(panel=idx):
                    panel(
                        fig,
                        pygmt,
                        grids[key],
                        region=region,
                        title=title,
                        cmap=cmap,
                        series=series,
                        colorbar_label=colorbar_label,
                        prefecture_boundary=prefecture_boundary,
                    )
        fig.savefig(output_path, crop=True, dpi=420)


def write_summary(metrics: pd.DataFrame, path: Path, radius_km: float, spacing: float) -> None:
    rows: list[dict[str, object]] = [
        {"key": "radius_km", "value": radius_km},
        {"key": "spacing_deg", "value": spacing},
        {"key": "n_cells", "value": len(metrics)},
    ]
    for column in [
        "local_site_delta_std_20km",
        "local_amp_cv_20km",
        "local_avs30_std_20km",
        "site_delta_gradient_per_10km",
    ]:
        values = metrics[column].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        rows.extend(
            [
                {"key": f"{column}_median", "value": float(np.nanquantile(values, 0.50))},
                {"key": f"{column}_p90", "value": float(np.nanquantile(values, 0.90))},
                {"key": f"{column}_p99", "value": float(np.nanquantile(values, 0.99))},
            ]
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)

    grid = pd.read_csv(args.grid_csv)
    amp, lat, lon = dataframe_to_array(grid, args.region, args.spacing, "amplification_vs400")
    avs30, _, _ = dataframe_to_array(grid, args.region, args.spacing, "avs30_m_s")
    site_delta = args.amp_intensity_coef * np.log10(np.clip(amp, 0.05, None))

    radius_cells = max(1, int(round(args.radius_km / (111.32 * args.spacing))))
    kernel = circular_kernel(radius_cells)
    amp_mean, amp_std, cell_count = local_mean_std(amp, kernel)
    _, avs30_std, _ = local_mean_std(avs30, kernel)
    _, site_delta_std, _ = local_mean_std(site_delta, kernel)
    amp_cv = amp_std / amp_mean
    site_delta_gradient = gradient_per_10km(site_delta, lat, args.spacing)

    arrays = {
        "site_intensity_delta": site_delta,
        "local_site_delta_std_20km": site_delta_std,
        "local_amp_cv_20km": amp_cv.astype(np.float32),
        "local_avs30_std_20km": avs30_std,
        "site_delta_gradient_per_10km": site_delta_gradient,
        "local_ground_cell_count_20km": cell_count,
    }
    metrics = flatten_metrics(lat, lon, arrays)
    metrics_csv = args.csv_dir / "jshis_ground_variability_metrics_0p02deg.csv"
    summary_csv = args.csv_dir / "jshis_ground_variability_summary.csv"
    metrics.to_csv(metrics_csv, index=False, float_format="%.6g")
    write_gzip_csv(metrics, args.derived_csv_gz)
    write_summary(metrics, summary_csv, args.radius_km, args.spacing)

    dataarrays = {
        key: to_dataarray(values, lat, lon, key, "")
        for key, values in arrays.items()
        if key != "local_ground_cell_count_20km" and key != "site_intensity_delta"
    }

    configure_conda_gmt()
    import pygmt

    combined = args.png_dir / "jshis_ground_variability_metrics_map.png"
    draw_maps(pygmt, dataarrays, combined, args.region, args.prefecture_boundary)
    for asset_path in (args.assets_ja, args.assets_en):
        asset_path.parent.mkdir(parents=True, exist_ok=True)
        asset_path.write_bytes(combined.read_bytes())

    print(f"Metrics CSV -> {metrics_csv}")
    print(f"Derived CSV.gz -> {args.derived_csv_gz}")
    print(f"Summary CSV -> {summary_csv}")
    print(f"PNG -> {combined}")


if __name__ == "__main__":
    main()
