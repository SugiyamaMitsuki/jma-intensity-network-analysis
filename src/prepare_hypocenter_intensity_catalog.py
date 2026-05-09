#!/usr/bin/env python3
"""Prepare JMA hypocenter catalog tables and match them to intensity events."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable, Iterator

import pandas as pd


DEFAULT_DATA_DIR = Path("data")
DEFAULT_EVENT_SUMMARY = Path("outputs/csv/event_intensity_summary.csv")
DEFAULT_OUT_DIR = Path("outputs/csv/hypocenter_catalog")
PREFERRED_RECORD_TYPES = {"J"}
PREFERRED_SOURCE_FLAGS = {"K", "k", "A"}
EPOCH = datetime(1900, 1, 1)


@dataclass(frozen=True)
class HypocenterEvent:
    hypocenter_id: str
    source_archive: str
    source_member: str
    source_line_number: int
    record_type: str
    source_flag: str
    origin_time: datetime
    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: float
    origin_epoch_s: float
    latitude: float
    longitude: float
    depth_km: float | None
    magnitude: float | None
    magnitude_type: str
    magnitude2: float | None
    magnitude2_type: str
    origin_time_error_s: float | None
    latitude_error_min: float | None
    longitude_error_min: float | None
    depth_error_km: float | None
    travel_time_table: str
    hypocenter_eval: str
    hypocenter_info: str
    max_intensity_code: str
    damage_scale: str
    tsunami_scale: str
    large_area_code: str
    small_area_code: str
    region_code: str
    region_name: str
    station_count: int | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse JMA Monthly Earthquake Catalog hypocenter files and match "
            "high precision hypocenters to parsed JMA intensity event summaries."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--event-summary", type=Path, default=DEFAULT_EVENT_SUMMARY)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--time-window-s",
        type=float,
        default=10.0,
        help="Maximum origin-time difference for an automatic match.",
    )
    parser.add_argument(
        "--candidate-minute-window",
        type=int,
        default=2,
        help="Minute window used while streaming the hypocenter catalog.",
    )
    parser.add_argument(
        "--write-full-catalog",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write the full parsed hypocenter catalog as csv.gz.",
    )
    return parser.parse_args()


def _parse_int(s: str) -> int:
    return int(s.strip().replace(" ", ""))


def _parse_optional_int(s: str) -> int | None:
    if not s.strip():
        return None
    return int(s.strip().replace(" ", ""))


def _parse_scaled_int(s: str, scale: float) -> float:
    return int(s.strip().replace(" ", "")) / scale


def _parse_optional_scaled_int(s: str, scale: float) -> float | None:
    if not s.strip():
        return None
    return int(s.strip().replace(" ", "")) / scale


def deg_min_to_decimal(deg: int, minutes: float) -> float:
    sign = -1 if deg < 0 else 1
    return sign * (abs(deg) + minutes / 60.0)


def parse_depth(s: str) -> float | None:
    if not s.strip():
        return None
    if s[3:5] == "  ":
        return float(int(s[0:3]))
    return int(s.replace(" ", "0")) / 100.0


def parse_magnitude(s: str) -> float | None:
    if not s.strip():
        return None
    if s[0] == "-":
        return -int(s[1]) / 10.0
    if s[0].isalpha():
        base = -(ord(s[0].upper()) - ord("A") + 1)
        return base - int(s[1]) / 10.0
    return int(s) / 10.0


def origin_from_parts(year: int, month: int, day: int, hour: int, minute: int, second: float) -> datetime:
    return datetime(year, month, day, hour, minute) + timedelta(seconds=second)


def datetime_to_epoch_s(value: datetime) -> float:
    return (value - EPOCH).total_seconds()


def minute_key_from_epoch(epoch_s: float) -> int:
    return math.floor(epoch_s / 60.0)


def parse_hypocenter_line(
    line: bytes,
    hypocenter_id: str,
    source_archive: str,
    source_member: str,
    source_line_number: int,
) -> HypocenterEvent:
    raw = line.rstrip(b"\r\n").decode("ascii", errors="replace").ljust(96)
    year = _parse_int(raw[1:5])
    month = _parse_int(raw[5:7])
    day = _parse_int(raw[7:9])
    hour = _parse_int(raw[9:11])
    minute = _parse_int(raw[11:13])
    second = _parse_scaled_int(raw[13:17], 100)
    origin = origin_from_parts(year, month, day, hour, minute, second)

    lat_deg = _parse_int(raw[21:24])
    lat_min = _parse_scaled_int(raw[24:28], 100)
    lon_deg = _parse_int(raw[32:36])
    lon_min = _parse_scaled_int(raw[36:40], 100)
    latitude = deg_min_to_decimal(lat_deg, lat_min)
    longitude = deg_min_to_decimal(lon_deg, lon_min)

    area_large = raw[64:65].strip()
    area_small = raw[65:68].strip().zfill(3)
    return HypocenterEvent(
        hypocenter_id=hypocenter_id,
        source_archive=source_archive,
        source_member=source_member,
        source_line_number=source_line_number,
        record_type=raw[0:1],
        source_flag=raw[95:96],
        origin_time=origin,
        year=year,
        month=month,
        day=day,
        hour=hour,
        minute=minute,
        second=second,
        origin_epoch_s=datetime_to_epoch_s(origin),
        latitude=latitude,
        longitude=longitude,
        depth_km=parse_depth(raw[44:49]),
        magnitude=parse_magnitude(raw[52:54]),
        magnitude_type=raw[54:55].strip(),
        magnitude2=parse_magnitude(raw[55:57]),
        magnitude2_type=raw[57:58].strip(),
        origin_time_error_s=_parse_optional_scaled_int(raw[17:21], 100),
        latitude_error_min=_parse_optional_scaled_int(raw[28:32], 100),
        longitude_error_min=_parse_optional_scaled_int(raw[40:44], 100),
        depth_error_km=_parse_optional_scaled_int(raw[49:52], 100),
        travel_time_table=raw[58:59].strip(),
        hypocenter_eval=raw[59:60].strip(),
        hypocenter_info=raw[60:61].strip(),
        max_intensity_code=raw[61:62].strip(),
        damage_scale=raw[62:63].strip(),
        tsunami_scale=raw[63:64].strip(),
        large_area_code=area_large,
        small_area_code=area_small,
        region_code=f"{area_large}-{area_small}",
        region_name=raw[68:92].strip(),
        station_count=_parse_optional_int(raw[92:95]),
    )


def discover_hypocenter_files(data_dir: Path) -> list[Path]:
    paths = [path for path in data_dir.glob("h*.zip") if path.is_file()]
    return sorted(paths, key=lambda p: p.name)


def iter_hypocenter_events(paths: Iterable[Path]) -> Iterator[HypocenterEvent]:
    for path in paths:
        with zipfile.ZipFile(path) as zf:
            members = [name for name in zf.namelist() if not name.endswith("/")]
            if not members:
                continue
            member = members[0]
            with zf.open(member) as handle:
                for line_number, line in enumerate(handle, start=1):
                    if not line.strip():
                        continue
                    yield parse_hypocenter_line(
                        line=line,
                        hypocenter_id=f"{path.stem}_{line_number:07d}",
                        source_archive=path.name,
                        source_member=member,
                        source_line_number=line_number,
                    )


def event_origin(row: pd.Series) -> datetime:
    return origin_from_parts(
        int(row["year"]),
        int(row["month"]),
        int(row["day"]),
        int(row["hour"]),
        int(row["minute"]),
        float(row["second"]),
    )


def build_candidate_minute_keys(events: pd.DataFrame, minute_window: int) -> set[int]:
    keys: set[int] = set()
    for origin in events["intensity_origin_time_dt"]:
        center = minute_key_from_epoch(datetime_to_epoch_s(origin))
        for delta in range(-minute_window, minute_window + 1):
            keys.add(center + delta)
    return keys


def hyp_to_csv_row(event: HypocenterEvent) -> dict[str, object]:
    row = asdict(event)
    row["origin_time"] = event.origin_time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4]
    return row


def write_full_catalog_header(writer: csv.DictWriter, event: HypocenterEvent) -> None:
    writer.fieldnames = list(hyp_to_csv_row(event).keys())
    writer.writeheader()


def haversine_km(lat1: float | None, lon1: float | None, lat2: float | None, lon2: float | None) -> float | None:
    if any(value is None or not math.isfinite(float(value)) for value in [lat1, lon1, lat2, lon2]):
        return None
    radius = 6371.0088
    phi1 = math.radians(float(lat1))
    phi2 = math.radians(float(lat2))
    dphi = math.radians(float(lat2) - float(lat1))
    dlambda = math.radians(float(lon2) - float(lon1))
    a = math.sin(dphi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    return 2.0 * radius * math.asin(math.sqrt(a))


def numeric(value: object) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(out):
        return None
    return out


def candidate_score(intensity_row: pd.Series, event: HypocenterEvent) -> tuple[float, ...]:
    intensity_epoch = float(intensity_row["intensity_origin_epoch_s"])
    time_diff = abs(event.origin_epoch_s - intensity_epoch)
    distance = haversine_km(
        numeric(intensity_row.get("latitude")),
        numeric(intensity_row.get("longitude")),
        event.latitude,
        event.longitude,
    )
    depth_i = numeric(intensity_row.get("depth_km"))
    mag_i = numeric(intensity_row.get("magnitude"))
    depth_diff = abs(event.depth_km - depth_i) if event.depth_km is not None and depth_i is not None else 999.0
    mag_diff = abs(event.magnitude - mag_i) if event.magnitude is not None and mag_i is not None else 99.0
    official_penalty = 0 if event.record_type in PREFERRED_RECORD_TYPES and event.source_flag in PREFERRED_SOURCE_FLAGS else 1
    distance_score = distance if distance is not None else 999.0
    return (time_diff, official_penalty, mag_diff, distance_score, depth_diff)


def classify_match(
    time_diff_s: float | None,
    distance_km: float | None,
    magnitude_diff: float | None,
    source_flag: str | None,
) -> tuple[str, str]:
    if time_diff_s is None:
        return "unmatched", "none"
    official = source_flag in PREFERRED_SOURCE_FLAGS
    distance_ok = distance_km is None or distance_km <= 20.0
    mag_ok = magnitude_diff is None or magnitude_diff <= 0.3
    if time_diff_s <= 1.0 and distance_ok and mag_ok and official:
        return "matched", "good"
    if time_diff_s <= 5.0 and (distance_km is None or distance_km <= 50.0) and (
        magnitude_diff is None or magnitude_diff <= 0.6
    ):
        return "matched", "fair"
    return "matched", "weak"


def empty_match_fields() -> dict[str, object]:
    return {
        "hyp_match_status": "unmatched",
        "hyp_match_quality": "none",
        "hyp_match_time_diff_s": None,
        "hyp_match_distance_km": None,
        "hyp_match_depth_diff_km": None,
        "hyp_match_magnitude_diff": None,
        "hypocenter_id": None,
        "hyp_source_archive": None,
        "hyp_source_line_number": None,
        "hyp_record_type": None,
        "hyp_source_flag": None,
        "hyp_origin_time": None,
        "hyp_latitude": None,
        "hyp_longitude": None,
        "hyp_depth_km": None,
        "hyp_magnitude": None,
        "hyp_magnitude_type": None,
        "hyp_magnitude2": None,
        "hyp_magnitude2_type": None,
        "hyp_origin_time_error_s": None,
        "hyp_latitude_error_min": None,
        "hyp_longitude_error_min": None,
        "hyp_depth_error_km": None,
        "hyp_region_code": None,
        "hyp_region_name": None,
        "hyp_station_count": None,
        "hyp_max_intensity_code": None,
        "hyp_travel_time_table": None,
        "hyp_hypocenter_eval": None,
        "hyp_hypocenter_info": None,
        "analysis_latitude": None,
        "analysis_longitude": None,
        "analysis_depth_km": None,
        "analysis_magnitude": None,
        "analysis_coordinate_source": "intensity_source",
    }


def build_match_fields(intensity_row: pd.Series, event: HypocenterEvent) -> dict[str, object]:
    time_diff = abs(event.origin_epoch_s - float(intensity_row["intensity_origin_epoch_s"]))
    distance = haversine_km(
        numeric(intensity_row.get("latitude")),
        numeric(intensity_row.get("longitude")),
        event.latitude,
        event.longitude,
    )
    depth_i = numeric(intensity_row.get("depth_km"))
    mag_i = numeric(intensity_row.get("magnitude"))
    depth_diff = abs(event.depth_km - depth_i) if event.depth_km is not None and depth_i is not None else None
    mag_diff = abs(event.magnitude - mag_i) if event.magnitude is not None and mag_i is not None else None
    status, quality = classify_match(time_diff, distance, mag_diff, event.source_flag)
    return {
        "hyp_match_status": status,
        "hyp_match_quality": quality,
        "hyp_match_time_diff_s": round(time_diff, 3),
        "hyp_match_distance_km": round(distance, 3) if distance is not None else None,
        "hyp_match_depth_diff_km": round(depth_diff, 3) if depth_diff is not None else None,
        "hyp_match_magnitude_diff": round(mag_diff, 3) if mag_diff is not None else None,
        "hypocenter_id": event.hypocenter_id,
        "hyp_source_archive": event.source_archive,
        "hyp_source_line_number": event.source_line_number,
        "hyp_record_type": event.record_type,
        "hyp_source_flag": event.source_flag,
        "hyp_origin_time": event.origin_time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4],
        "hyp_latitude": event.latitude,
        "hyp_longitude": event.longitude,
        "hyp_depth_km": event.depth_km,
        "hyp_magnitude": event.magnitude,
        "hyp_magnitude_type": event.magnitude_type,
        "hyp_magnitude2": event.magnitude2,
        "hyp_magnitude2_type": event.magnitude2_type,
        "hyp_origin_time_error_s": event.origin_time_error_s,
        "hyp_latitude_error_min": event.latitude_error_min,
        "hyp_longitude_error_min": event.longitude_error_min,
        "hyp_depth_error_km": event.depth_error_km,
        "hyp_region_code": event.region_code,
        "hyp_region_name": event.region_name,
        "hyp_station_count": event.station_count,
        "hyp_max_intensity_code": event.max_intensity_code,
        "hyp_travel_time_table": event.travel_time_table,
        "hyp_hypocenter_eval": event.hypocenter_eval,
        "hyp_hypocenter_info": event.hypocenter_info,
        "analysis_latitude": event.latitude,
        "analysis_longitude": event.longitude,
        "analysis_depth_km": event.depth_km if event.depth_km is not None else depth_i,
        "analysis_magnitude": event.magnitude if event.magnitude is not None else mag_i,
        "analysis_coordinate_source": "hypocenter_catalog",
    }


def match_intensity_events(
    intensity: pd.DataFrame,
    candidates_by_minute: dict[int, list[HypocenterEvent]],
    time_window_s: float,
) -> pd.DataFrame:
    matched_rows: list[dict[str, object]] = []
    for _, row in intensity.iterrows():
        origin_epoch = float(row["intensity_origin_epoch_s"])
        center = minute_key_from_epoch(origin_epoch)
        candidates: list[HypocenterEvent] = []
        for delta in range(-1, 2):
            candidates.extend(candidates_by_minute.get(center + delta, []))
        candidates = [
            event
            for event in candidates
            if abs(event.origin_epoch_s - origin_epoch) <= time_window_s
        ]
        base = row.drop(labels=["intensity_origin_time_dt"]).to_dict()
        if not candidates:
            fields = empty_match_fields()
            fields["analysis_latitude"] = row.get("latitude")
            fields["analysis_longitude"] = row.get("longitude")
            fields["analysis_depth_km"] = row.get("depth_km")
            fields["analysis_magnitude"] = row.get("magnitude")
            matched_rows.append({**base, **fields})
            continue
        best = min(candidates, key=lambda event: candidate_score(row, event))
        matched_rows.append({**base, **build_match_fields(row, best)})
    return pd.DataFrame(matched_rows)


def summarize_matches(matched: pd.DataFrame, hyp_years: list[int], hyp_counter: Counter, out_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    rows.append({"metric": "hypocenter_year_min", "value": min(hyp_years) if hyp_years else None})
    rows.append({"metric": "hypocenter_year_max", "value": max(hyp_years) if hyp_years else None})
    rows.append({"metric": "hypocenter_records_total", "value": int(hyp_counter["total"])})
    for key, count in sorted(hyp_counter.items()):
        if key == "total":
            continue
        rows.append({"metric": f"hypocenter_records_{key}", "value": int(count)})

    rows.append({"metric": "intensity_events_total", "value": len(matched)})
    rows.append({"metric": "intensity_events_matched", "value": int((matched["hyp_match_status"] == "matched").sum())})
    rows.append({"metric": "intensity_events_unmatched", "value": int((matched["hyp_match_status"] == "unmatched").sum())})
    for quality, count in matched["hyp_match_quality"].value_counts(dropna=False).sort_index().items():
        rows.append({"metric": f"match_quality_{quality}", "value": int(count)})
    if "year" in matched.columns:
        by_year = (
            matched.groupby("year", dropna=False)["hyp_match_status"]
            .agg(total="size", matched=lambda s: int((s == "matched").sum()))
            .reset_index()
        )
        by_year["match_rate"] = by_year["matched"] / by_year["total"]
        by_year.to_csv(out_dir / "intensity_hypocenter_match_summary_by_year.csv", index=False)
    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "intensity_hypocenter_match_summary.csv", index=False)
    return out


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    intensity = pd.read_csv(args.event_summary)
    intensity["intensity_origin_time_dt"] = intensity.apply(event_origin, axis=1)
    intensity["intensity_origin_time"] = intensity["intensity_origin_time_dt"].map(
        lambda dt: dt.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4]
    )
    intensity["intensity_origin_epoch_s"] = intensity["intensity_origin_time_dt"].map(datetime_to_epoch_s)

    candidate_keys = build_candidate_minute_keys(intensity, args.candidate_minute_window)
    hypocenter_paths = discover_hypocenter_files(args.data_dir)
    if not hypocenter_paths:
        raise SystemExit(f"No h*.zip hypocenter files found under {args.data_dir}")

    full_path = args.out_dir / "jma_hypocenter_catalog.csv.gz"
    candidates_by_minute: dict[int, list[HypocenterEvent]] = defaultdict(list)
    years: set[int] = set()
    counter: Counter = Counter()
    writer: csv.DictWriter | None = None
    gzip_handle = None

    try:
        if args.write_full_catalog:
            gzip_handle = gzip.open(full_path, "wt", encoding="utf-8", newline="")
        for event in iter_hypocenter_events(hypocenter_paths):
            years.add(event.year)
            counter["total"] += 1
            counter[f"record_type_{event.record_type or 'blank'}"] += 1
            counter[f"source_flag_{event.source_flag or 'blank'}"] += 1
            key = minute_key_from_epoch(event.origin_epoch_s)
            if key in candidate_keys:
                candidates_by_minute[key].append(event)
            if gzip_handle is not None:
                if writer is None:
                    writer = csv.DictWriter(gzip_handle, fieldnames=list(hyp_to_csv_row(event).keys()))
                    writer.writeheader()
                writer.writerow(hyp_to_csv_row(event))
    finally:
        if gzip_handle is not None:
            gzip_handle.close()

    matched = match_intensity_events(intensity, candidates_by_minute, args.time_window_s)
    matched_path = args.out_dir / "jma_intensity_events_with_hypocenter.csv"
    matched.to_csv(matched_path, index=False)

    target_5upper = matched[matched["max_intensity_value"] >= 5.5].copy()
    target_5upper_path = args.out_dir / "target_events_intensity_5upper_plus_with_hypocenter.csv"
    target_5upper.to_csv(target_5upper_path, index=False)

    target_6lower = matched[matched["max_intensity_value"] >= 6.0].copy()
    target_6lower_path = args.out_dir / "target_events_intensity_6lower_plus_with_hypocenter.csv"
    target_6lower.to_csv(target_6lower_path, index=False)

    target = matched[matched["max_intensity_value"] >= 6.5].copy()
    target_path = args.out_dir / "target_events_intensity_6upper_plus_with_hypocenter.csv"
    target.to_csv(target_path, index=False)

    summary = summarize_matches(matched, sorted(years), counter, args.out_dir)
    print(f"Parsed hypocenter files: {len(hypocenter_paths)}")
    if args.write_full_catalog:
        print(f"Saved full hypocenter catalog: {full_path}")
    print(f"Saved matched intensity catalog: {matched_path}")
    print(f"Saved target matched catalog >= 5 upper: {target_5upper_path}")
    print(f"Saved target matched catalog >= 6 lower: {target_6lower_path}")
    print(f"Saved target matched catalog >= 6 upper: {target_path}")
    print(f"Saved summary: {args.out_dir / 'intensity_hypocenter_match_summary.csv'}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
