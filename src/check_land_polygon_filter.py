#!/usr/bin/env python3
"""Cross-check the region-name land approximation against land polygons."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from shapely.geometry import Point
from shapely.strtree import STRtree

from analyze_jma_intensity import configure_journal_matplotlib, save_journal_png, style_journal_axis


MAP_REGION = [122.0, 146.5, 24.0, 46.5]


def load_land_geometries(resolution: str):
    import cartopy.io.shapereader as shpreader

    path = shpreader.natural_earth(resolution=resolution, category="physical", name="land")
    return list(shpreader.Reader(path).geometries())


def classify_land_points(df: pd.DataFrame, geoms: list) -> pd.Series:
    tree = STRtree(geoms)
    values = []
    for lon, lat in zip(df["longitude"], df["latitude"]):
        if pd.isna(lon) or pd.isna(lat):
            values.append(pd.NA)
            continue
        pt = Point(float(lon), float(lat))
        on_land = False
        for item in tree.query(pt):
            geom = geoms[int(item)] if isinstance(item, (int, float)) or str(item).isdigit() else item
            if geom.contains(pt) or geom.touches(pt):
                on_land = True
                break
        values.append(on_land)
    return pd.Series(values, index=df.index, dtype="boolean")


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    filters = {
        "all_with_location": df,
        "m4plus_shallow20": df[(df["magnitude"] >= 4.0) & (df["depth_km"] <= 20.0)],
        "m4_5_shallow20": df[(df["magnitude"] >= 4.0) & (df["magnitude"] < 5.0) & (df["depth_km"] <= 20.0)],
    }
    for label, sub in filters.items():
        sub = sub[sub["land_polygon"].notna() & sub["is_land_region_approx"].notna()]
        n = len(sub)
        if n == 0:
            continue
        approx = sub["is_land_region_approx"].astype(bool)
        poly = sub["land_polygon"].astype(bool)
        rows.append(
            {
                "sample": label,
                "n_events": n,
                "both_land": int((approx & poly).sum()),
                "both_nonland": int((~approx & ~poly).sum()),
                "approx_land_polygon_nonland": int((approx & ~poly).sum()),
                "approx_nonland_polygon_land": int((~approx & poly).sum()),
                "agreement_share": float((approx == poly).mean()),
                "approx_land_share": float(approx.mean()),
                "polygon_land_share": float(poly.mean()),
            }
        )
    return pd.DataFrame(rows)


def plot_mismatch_map(df: pd.DataFrame, out_path: Path) -> None:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt

    configure_journal_matplotlib(plt)
    data = df[
        df["land_polygon"].notna()
        & df["is_land_region_approx"].notna()
        & (df["magnitude"] >= 4.0)
        & (df["depth_km"] <= 20.0)
        & (df["longitude"].between(MAP_REGION[0], MAP_REGION[1]))
        & (df["latitude"].between(MAP_REGION[2], MAP_REGION[3]))
    ].copy()
    if data.empty:
        return
    approx = data["is_land_region_approx"].astype(bool)
    poly = data["land_polygon"].astype(bool)
    data["mismatch_type"] = "agree"
    data.loc[approx & ~poly, "mismatch_type"] = "region land / polygon sea"
    data.loc[~approx & poly, "mismatch_type"] = "region sea / polygon land"

    fig = plt.figure(figsize=(7.0, 7.0))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(MAP_REGION, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#f2f0e8", edgecolor="none")
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#e9f2f5", edgecolor="none")
    ax.coastlines(resolution="10m", linewidth=0.45, color="#555555")
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.25, edgecolor="#777777")
    gl = ax.gridlines(draw_labels=True, linewidth=0.25, color="#d7d7d7", alpha=0.7)
    gl.top_labels = False
    gl.right_labels = False

    agree = data[data["mismatch_type"] == "agree"]
    ax.scatter(
        agree["longitude"],
        agree["latitude"],
        s=4,
        color="#8f969e",
        alpha=0.18,
        transform=ccrs.PlateCarree(),
        label=f"agree (n={len(agree):,})",
    )
    colors = {
        "region land / polygon sea": "#D55E00",
        "region sea / polygon land": "#0072B2",
    }
    for label, color in colors.items():
        sub = data[data["mismatch_type"] == label]
        if sub.empty:
            continue
        ax.scatter(
            sub["longitude"],
            sub["latitude"],
            s=14,
            color=color,
            alpha=0.72,
            edgecolor="white",
            linewidth=0.2,
            transform=ccrs.PlateCarree(),
            label=f"{label} (n={len(sub):,})",
        )
    ax.set_title("Land-filter cross-check, M>=4 shallow events")
    ax.legend(loc="lower left", frameon=True, fontsize=7)
    save_journal_png(fig, out_path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-csv", type=Path, default=Path("outputs/csv/event_intensity_summary.csv"))
    parser.add_argument("--csv-dir", type=Path, default=Path("outputs/csv/remaining_closure"))
    parser.add_argument("--png-dir", type=Path, default=Path("outputs/png/remaining_closure"))
    parser.add_argument("--resolution", default="10m", choices=["10m", "50m", "110m"])
    args = parser.parse_args()

    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.event_csv)
    df = df[df["latitude"].notna() & df["longitude"].notna()].copy()
    geoms = load_land_geometries(args.resolution)
    df["land_polygon"] = classify_land_points(df, geoms)
    df["land_filter_agrees"] = df["land_polygon"].astype("boolean") == df["is_land_region_approx"].astype("boolean")
    keep_cols = [
        "event_id",
        "year",
        "magnitude",
        "depth_km",
        "latitude",
        "longitude",
        "hypocenter_region",
        "is_land_region_approx",
        "land_polygon",
        "land_filter_agrees",
    ]
    event_out = args.csv_dir / "event_land_polygon_crosscheck.csv"
    summary_out = args.csv_dir / "land_polygon_crosscheck_summary.csv"
    df[keep_cols].to_csv(event_out, index=False)
    summary = summarize(df)
    summary.to_csv(summary_out, index=False)
    plot_mismatch_map(df, args.png_dir / "land_filter_mismatch_map.png")
    print(f"CSV -> {event_out}")
    print(f"CSV -> {summary_out}")
    print(f"PNG -> {args.png_dir / 'land_filter_mismatch_map.png'}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
