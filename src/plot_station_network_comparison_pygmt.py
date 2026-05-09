#!/usr/bin/env python3
"""Plot side-by-side year-end snapshots of the JMA intensity station network."""

from __future__ import annotations

import argparse
import calendar
import csv
from pathlib import Path

from analyze_jma_intensity import StationRecord, load_station_index, station_active_at_key
from map_prefecture_boundaries import DEFAULT_PREFECTURE_BOUNDARY, plot_prefecture_boundaries_pygmt
from plot_station_maps_pygmt import configure_conda_gmt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--code-p", type=Path, default=None)
    parser.add_argument("--first-year", type=int, default=1994)
    parser.add_argument("--first-month", type=int, default=12)
    parser.add_argument("--latest-year", default="auto")
    parser.add_argument("--latest-month", type=int, default=12)
    parser.add_argument("--region", nargs=4, type=float, default=[122.0, 146.5, 24.0, 46.5])
    parser.add_argument("--projection", default="M8.5c")
    parser.add_argument("--output", type=Path, default=Path("outputs/png/station_maps_pygmt/station_network_1994_vs_latest_year_end.png"))
    parser.add_argument("--csv", type=Path, default=Path("outputs/csv/station_network_1994_vs_latest_year_end.csv"))
    parser.add_argument("--assets-ja", type=Path, default=Path("paper/assets_ja/fig03b_station_network_1994_vs_2022.png"))
    parser.add_argument("--assets-en", type=Path, default=Path("paper/assets_en/fig03b_station_network_1994_vs_2022.png"))
    parser.add_argument("--prefecture-boundary", type=Path, default=DEFAULT_PREFECTURE_BOUNDARY)
    return parser.parse_args()


def latest_intensity_year(data_dir: Path) -> int:
    years: set[int] = set()
    for pattern in ("i[0-9][0-9][0-9][0-9].zip", "i[0-9][0-9][0-9][0-9].dat"):
        for path in data_dir.glob(pattern):
            year_text = path.stem[1:5]
            if year_text.isdigit():
                years.add(int(year_text))
    if not years:
        raise SystemExit(f"No iYYYY.zip or iYYYY.dat files found in {data_dir}")
    return max(years)


def month_end_key(year: int, month: int) -> int:
    day = calendar.monthrange(year, month)[1]
    return int(f"{year:04d}{month:02d}{day:02d}2359")


def label_for_month_end(year: int, month: int) -> str:
    day = calendar.monthrange(year, month)[1]
    return f"{year:04d}-{month:02d}-{day:02d} 23:59 JST"


def active_records_for_time_key(station_index, time_key: int) -> list[StationRecord]:
    active: list[StationRecord] = []
    for records_for_code in station_index.by_code.values():
        records = [rec for rec in records_for_code if station_active_at_key(rec, time_key)]
        if not records:
            continue
        records.sort(key=lambda rec: rec.start_key or 0)
        rec = records[-1]
        if (
            rec.latitude is not None
            and rec.longitude is not None
            and 120.0 <= rec.longitude <= 150.0
            and 20.0 <= rec.latitude <= 48.0
        ):
            active.append(rec)
    return active


def draw_panel(
    fig,
    records: list[StationRecord],
    *,
    region: list[float],
    projection: str,
    title: str,
    panel_label: str,
    fill: str,
    symbol: str,
    pen: str,
    prefecture_boundary: Path | None,
) -> None:
    fig.basemap(
        region=region,
        projection=projection,
        frame=[
            "xafg5+lLongitude",
            "yafg5+lLatitude",
            f"+t{title}",
        ],
    )
    fig.coast(
        land="#f3f3f0",
        water="#e8f0f3",
        shorelines="0.35p,#5a5a5a",
        borders="1/0.18p,#8a8a8a",
        resolution="i",
    )
    plot_prefecture_boundaries_pygmt(fig, prefecture_boundary, pen="0.14p,#707070")
    fig.plot(
        x=[rec.longitude for rec in records],
        y=[rec.latitude for rec in records],
        style=symbol,
        fill=fill,
        pen=pen,
        transparency=8,
    )
    fig.text(
        x=region[0] + 0.55,
        y=region[3] - 0.55,
        text=panel_label,
        font="12p,Helvetica-Bold,#111111",
        justify="LM",
        fill="white@10",
        clearance="0.08c/0.05c",
    )
    fig.text(
        x=region[0] + 0.55,
        y=region[3] - 1.15,
        text=f"N = {len(records):,}",
        font="9.5p,Helvetica-Bold,#111111",
        justify="LM",
        fill="white@10",
        clearance="0.08c/0.05c",
    )


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["snapshot", "year", "month", "as_of", "time_key", "active_station_count"],
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    configure_conda_gmt()
    import pygmt

    args = parse_args()
    latest_year = latest_intensity_year(args.data_dir) if args.latest_year == "auto" else int(args.latest_year)

    station_index = load_station_index(args.data_dir, args.code_p)
    if station_index is None:
        raise SystemExit("code_p.dat/code_p.zip was not found.")

    snapshots = [
        ("a", "pre_expansion", args.first_year, args.first_month, "#d73027", "t0.075c", "0.05p,white"),
        ("b", "latest", latest_year, args.latest_month, "#1f78b4", "c0.040c", "0.025p,white"),
    ]

    rows: list[dict[str, object]] = []
    panels: list[tuple[str, str, list[StationRecord], str, str, str]] = []
    for panel_label, snapshot, year, month, fill, symbol, pen in snapshots:
        time_key = month_end_key(year, month)
        records = active_records_for_time_key(station_index, time_key)
        label = label_for_month_end(year, month)
        rows.append(
            {
                "snapshot": snapshot,
                "year": year,
                "month": month,
                "as_of": label,
                "time_key": time_key,
                "active_station_count": len(records),
            }
        )
        panels.append((panel_label, label, records, fill, symbol, pen))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_TICK_LENGTH_PRIMARY="2.3p",
        FORMAT_GEO_MAP="dddF",
        FONT_ANNOT_PRIMARY="7.5p,Helvetica,black",
        FONT_LABEL="8p,Helvetica,black",
        FONT_TITLE="8.8p,Helvetica-Bold,black",
        MAP_GRID_PEN_PRIMARY="0.12p,#d6d6d6",
        MAP_FRAME_PEN="0.65p,#222222",
    ):
        fig = pygmt.Figure()
        with fig.subplot(nrows=1, ncols=2, figsize=("18.0c", "10.8c"), margins=["0.22c", "0.05c"]):
            for idx, (panel_label, label, records, fill, symbol, pen) in enumerate(panels):
                with fig.set_panel(panel=idx):
                    draw_panel(
                        fig,
                        records,
                        region=args.region,
                        projection=args.projection,
                        title=label,
                        panel_label=panel_label,
                        fill=fill,
                        symbol=symbol,
                        pen=pen,
                        prefecture_boundary=args.prefecture_boundary,
                    )
        fig.savefig(args.output, dpi=450)

    write_summary(args.csv, rows)
    for asset_path in (args.assets_ja, args.assets_en):
        asset_path.parent.mkdir(parents=True, exist_ok=True)
        asset_path.write_bytes(args.output.read_bytes())
    print(f"Figure -> {args.output}")
    print(f"CSV -> {args.csv}")
    for row in rows:
        print(f"{row['snapshot']}: {row['as_of']}, N={row['active_station_count']:,}")


if __name__ == "__main__":
    main()
