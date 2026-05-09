#!/usr/bin/env python3
"""Create yearly PyGMT maps of active seismic intensity stations from code_p.dat."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from analyze_jma_intensity import (
    StationRecord,
    active_station_records_at_year_end,
    load_station_index,
    station_network_region,
)
from map_prefecture_boundaries import DEFAULT_PREFECTURE_BOUNDARY, plot_prefecture_boundaries_pygmt


REGION_STYLE = {
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


def configure_conda_gmt() -> None:
    """Prefer GMT from the active conda environment before importing PyGMT."""

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


def parse_years(value: str | None, data_dir: Path) -> list[int]:
    if value is None or value.lower() == "auto":
        years = sorted(
            int(path.stem[1:])
            for path in data_dir.glob("i[0-9][0-9][0-9][0-9].zip")
            if path.stem[1:].isdigit()
        )
        if not years:
            raise SystemExit(f"No iYYYY.zip files found in {data_dir}")
        return list(range(min(years), max(years) + 1))

    out: set[int] = set()
    for part in value.split(","):
        item = part.strip()
        if not item:
            continue
        if ":" in item or "-" in item:
            sep = ":" if ":" in item else "-"
            start, end = [int(token) for token in item.split(sep, 1)]
            out.update(range(start, end + 1))
        else:
            out.add(int(item))
    return sorted(out)


def active_records_for_year(station_index, year: int) -> list[StationRecord]:
    records = active_station_records_at_year_end(station_index, year)
    return [
        rec
        for rec in records
        if rec.latitude is not None
        and rec.longitude is not None
        and 120.0 <= rec.longitude <= 150.0
        and 20.0 <= rec.latitude <= 48.0
    ]


def draw_station_map(
    pygmt,
    records: list[StationRecord],
    year: int,
    output_path: Path,
    region: list[float],
    projection: str,
    color_by_region: bool,
    prefecture_boundary: Path | None = DEFAULT_PREFECTURE_BOUNDARY,
) -> None:
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
            frame=[
                "xafg5+lLongitude",
                "yafg5+lLatitude",
                f"+tActive seismic intensity stations, {year}",
            ],
        )
        fig.coast(
            land="#f2f2ee",
            water="#e9f2f5",
            shorelines="0.35p,#555555",
            borders="1/0.2p,#888888",
            resolution="i",
            map_scale="jBL+w300k+o0.35c/0.35c+f+lkm",
        )
        plot_prefecture_boundaries_pygmt(fig, prefecture_boundary, pen="0.16p,#666666")

        if color_by_region:
            for region_name, color in REGION_STYLE.items():
                subset = [rec for rec in records if station_network_region(rec) == region_name]
                if not subset:
                    continue
                fig.plot(
                    x=[rec.longitude for rec in subset],
                    y=[rec.latitude for rec in subset],
                    style="c0.055c",
                    fill=color,
                    pen="0.03p,white",
                    transparency=8,
                    label=region_name,
                )
            fig.legend(position="JTR+jTR+o0.15c/0.15c", box="+gwhite@18+p0.25p,#777777")
        else:
            fig.plot(
                x=[rec.longitude for rec in records],
                y=[rec.latitude for rec in records],
                style="c0.052c",
                fill="#005f73",
                pen="0.03p,white",
                transparency=10,
            )

        fig.text(
            x=region[0] + 0.65,
            y=region[3] - 0.75,
            text=f"N = {len(records):,}",
            font="10p,Helvetica-Bold,#222222",
            justify="LM",
            fill="white@18",
            clearance="0.08c/0.05c",
        )
    fig.savefig(output_path)


def main() -> None:
    configure_conda_gmt()
    import pygmt

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--code-p", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/png/station_maps_pygmt"))
    parser.add_argument("--years", default="auto", help='Year list/range, e.g. "1980:2022" or "1995,2003,2022".')
    parser.add_argument("--projection", default="M15c")
    parser.add_argument("--region", nargs=4, type=float, default=[122.0, 146.5, 24.0, 46.5])
    parser.add_argument("--color-by-region", action="store_true")
    parser.add_argument("--prefecture-boundary", type=Path, default=DEFAULT_PREFECTURE_BOUNDARY)
    args = parser.parse_args()

    station_index = load_station_index(args.data_dir, args.code_p)
    if station_index is None:
        raise SystemExit("code_p.dat/code_p.zip was not found.")

    years = parse_years(args.years, args.data_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for year in years:
        records = active_records_for_year(station_index, year)
        out = args.output_dir / f"station_distribution_{year}.png"
        draw_station_map(
            pygmt=pygmt,
            records=records,
            year=year,
            output_path=out,
            region=args.region,
            projection=args.projection,
            color_by_region=args.color_by_region,
            prefecture_boundary=args.prefecture_boundary,
        )
        print(f"{year}: {len(records):,} stations -> {out}")


if __name__ == "__main__":
    main()
