#!/usr/bin/env python3
"""Create monthly PyGMT station-distribution frames and an MP4 movie."""

from __future__ import annotations

import argparse
import calendar
import csv
import subprocess
from pathlib import Path

from analyze_jma_intensity import (
    StationRecord,
    load_station_index,
    station_active_at_key,
)
from plot_station_maps_pygmt import configure_conda_gmt, draw_station_map
from map_prefecture_boundaries import DEFAULT_PREFECTURE_BOUNDARY


def parse_months(value: str, data_dir: Path) -> list[tuple[int, int]]:
    if value.lower() == "auto":
        years = sorted(
            int(path.stem[1:])
            for path in data_dir.glob("i[0-9][0-9][0-9][0-9].zip")
            if path.stem[1:].isdigit()
        )
        if not years:
            raise SystemExit(f"No iYYYY.zip files found in {data_dir}")
        return [(year, month) for year in range(min(years), max(years) + 1) for month in range(1, 13)]

    out: list[tuple[int, int]] = []
    for part in value.split(","):
        item = part.strip()
        if not item:
            continue
        if ":" in item:
            start, end = item.split(":", 1)
            sy, sm = [int(token) for token in start.split("-", 1)]
            ey, em = [int(token) for token in end.split("-", 1)]
            year, month = sy, sm
            while (year, month) <= (ey, em):
                out.append((year, month))
                month += 1
                if month == 13:
                    year += 1
                    month = 1
        else:
            year, month = [int(token) for token in item.split("-", 1)]
            out.append((year, month))
    return out


def month_end_key(year: int, month: int) -> int:
    day = calendar.monthrange(year, month)[1]
    return int(f"{year:04d}{month:02d}{day:02d}2359")


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


def write_counts_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["year", "month", "label", "active_station_count", "frame"])
        writer.writeheader()
        writer.writerows(rows)


def make_video(frame_dir: Path, output_path: Path, fps: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pattern = str(frame_dir / "station_distribution_monthly_%06d.png")
    cmd = [
        "ffmpeg",
        "-y",
        "-framerate",
        str(fps),
        "-i",
        pattern,
        "-vf",
        "pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p",
        "-c:v",
        "libx264",
        "-crf",
        "18",
        "-pix_fmt",
        "yuv420p",
        str(output_path),
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    configure_conda_gmt()
    import pygmt

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--code-p", type=Path, default=None)
    parser.add_argument("--frame-dir", type=Path, default=Path("outputs/png/station_maps_pygmt/monthly_frames"))
    parser.add_argument("--video-path", type=Path, default=Path("outputs/png/station_maps_pygmt/station_distribution_monthly_198001_202212.mp4"))
    parser.add_argument("--counts-csv", type=Path, default=Path("outputs/csv/station_distribution_monthly_counts.csv"))
    parser.add_argument("--months", default="auto", help='Month list/range, e.g. "1980-01:2022-12" or "1995-01,2022-12".')
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--projection", default="M15c")
    parser.add_argument("--region", nargs=4, type=float, default=[122.0, 146.5, 24.0, 46.5])
    parser.add_argument("--color-by-region", action="store_true")
    parser.add_argument("--prefecture-boundary", type=Path, default=DEFAULT_PREFECTURE_BOUNDARY)
    parser.add_argument("--overwrite-frames", action="store_true")
    parser.add_argument("--write-dated-copies", action="store_true")
    parser.add_argument("--skip-video", action="store_true")
    args = parser.parse_args()

    station_index = load_station_index(args.data_dir, args.code_p)
    if station_index is None:
        raise SystemExit("code_p.dat/code_p.zip was not found.")

    months = parse_months(args.months, args.data_dir)
    args.frame_dir.mkdir(parents=True, exist_ok=True)
    counts: list[dict[str, object]] = []
    for idx, (year, month) in enumerate(months, start=1):
        label = f"{year:04d}-{month:02d}"
        time_key = month_end_key(year, month)
        records = active_records_for_time_key(station_index, time_key)
        frame = args.frame_dir / f"station_distribution_monthly_{idx:06d}.png"
        if args.overwrite_frames or not frame.exists():
            draw_station_map(
                pygmt=pygmt,
                records=records,
                year=label,
                output_path=frame,
                region=args.region,
                projection=args.projection,
                color_by_region=args.color_by_region,
                prefecture_boundary=args.prefecture_boundary,
            )
        if args.write_dated_copies:
            dated_frame = args.frame_dir / f"station_distribution_{year:04d}_{month:02d}.png"
            dated_frame.write_bytes(frame.read_bytes())
        counts.append(
            {
                "year": year,
                "month": month,
                "label": label,
                "active_station_count": len(records),
                "frame": str(frame),
            }
        )
        print(f"{idx:04d}/{len(months):04d} {label}: {len(records):,} stations -> {frame}")

    write_counts_csv(args.counts_csv, counts)
    print(f"Counts CSV -> {args.counts_csv}")
    if not args.skip_video:
        make_video(args.frame_dir, args.video_path, args.fps)
        print(f"Video -> {args.video_path}")


if __name__ == "__main__":
    main()
