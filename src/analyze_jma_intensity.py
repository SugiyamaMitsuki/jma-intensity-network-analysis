#!/usr/bin/env python3
"""Analyze JMA seismic intensity catalog files.

The input data are JMA fixed-width intensity/acceleration files packed as
data/iYYYY.zip.  This script parses source records and station intensity
records, builds event-level summaries, and aggregates the relation between
magnitude and average observed intensity by year.
If data/code_p.zip is available, station names and epicentral distances are
also attached to the event-level summaries.

Outputs are written under the output directory:
  csv/
  - event_intensity_summary.csv
  - yearly_magbin_summary.csv
  - period_magbin_summary.csv
  - yearly_regression_summary.csv
  - yearly_distance_detection_summary.csv
  - period_distance_detection_histogram.csv
  - yearly_shortest_epicentral_distance_m4plus.csv
  - yearly_spatial_station_count_summary.csv
  - period_spatial_station_count_summary.csv
  - station_network_density_by_year_region.csv
  - station_network_density_period_summary.csv
  - station_network_density_selected_years.csv
  png/
  - yearly_mean_intensity_by_mag.png
  - yearly_max_intensity_by_mag.png
  - mean_intensity_heatmap.png
  - magnitude_relation_by_period.png
  - paper_like_yearly_max_intensity.png
  - yearly_mean_station_intensity_by_magbin.png
  - yearly_shortest_epicentral_distance.png
  - yearly_shortest_epicentral_distance_m4plus.png
  - shortest_epicentral_distance_histogram.png
  - yearly_station_counts_within_radius_m4_5.png
  - station_density_by_region.png
  - station_nearest_neighbor_distance_by_region.png
  - station_density_region_snapshot_1995_2022.png
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import statistics
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


SOURCE_RECORD_TYPES = {b"A", b"B", b"D"}


INTENSITY_VALUE = {
    "1": 1.0,
    "2": 2.0,
    "3": 3.0,
    "4": 4.0,
    # Before Oct. 1996, 5 and 6 were not split.  After the split, A-D are
    # plotted as half-step classes to keep the ordering visible.
    "5": 5.0,
    "6": 6.0,
    "7": 7.0,
    "A": 5.0,  # 5-lower
    "B": 5.5,  # 5-upper
    "C": 6.0,  # 6-lower
    "D": 6.5,  # 6-upper
}


MAG_BINS = [
    (3.0, 4.0, "3.0<=M<4.0"),
    (4.0, 5.0, "4.0<=M<5.0"),
    (5.0, 6.0, "5.0<=M<6.0"),
    (6.0, 7.0, "6.0<=M<7.0"),
    (7.0, 10.0, "7.0<=M"),
]


PERIODS = [
    (1980, 1995, "1980-1995"),
    (1996, 2003, "1996-2003"),
    (2004, 2010, "2004-2010"),
    (2011, 2022, "2011-2022"),
]

FIVE_YEAR_PERIODS = [
    (1980, 1985, "1980-1985"),
    (1986, 1990, "1986-1990"),
    (1991, 1995, "1991-1995"),
    (1996, 2000, "1996-2000"),
    (2001, 2005, "2001-2005"),
    (2006, 2010, "2006-2010"),
    (2011, 2015, "2011-2015"),
    (2016, 2020, "2016-2020"),
    (2021, 2022, "2021-2022"),
]

DISTANCE_BINS = [
    (0.0, 10.0, "0-10"),
    (10.0, 20.0, "10-20"),
    (20.0, 30.0, "20-30"),
    (30.0, 40.0, "30-40"),
    (40.0, 50.0, "40-50"),
    (50.0, 100.0, "50-100"),
    (100.0, 1_000_000.0, "100+"),
]

STATION_COUNT_RADII_KM = [10, 20, 50, 100]


@dataclass(frozen=True)
class NetworkRegion:
    name: str
    prefixes: tuple[int, ...]
    area_km2: float


# Approximate regional land areas.  They are used only to convert active station
# counts to a comparable density scale; exact boundary work should use GIS
# polygons outside this catalog parser.
NETWORK_REGIONS = [
    NetworkRegion("Hokkaido", tuple(range(10, 17)), 83424.0),
    NetworkRegion("Tohoku", tuple(range(20, 26)), 66948.0),
    NetworkRegion("Kanto", tuple(range(30, 37)), 32423.0),
    NetworkRegion("Chubu", tuple(range(37, 47)), 72572.0),
    NetworkRegion("Kinki", tuple(range(50, 56)), 33112.0),
    NetworkRegion("Chugoku", (56, 57, 58, 59, 70), 31922.0),
    NetworkRegion("Shikoku", tuple(range(60, 64)), 18804.0),
    NetworkRegion("Kyushu", tuple(range(71, 78)), 42231.0),
    NetworkRegion("Okinawa", (80,), 2282.0),
]
JAPAN_AREA_KM2 = sum(region.area_km2 for region in NETWORK_REGIONS)
NETWORK_REGION_BY_PREFIX = {
    prefix: region.name for region in NETWORK_REGIONS for prefix in region.prefixes
}
NETWORK_REGION_AREA = {region.name: region.area_km2 for region in NETWORK_REGIONS}
NETWORK_REGION_AREA["Japan"] = JAPAN_AREA_KM2


PALETTE = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#F0E442",
    "#000000",
    "#999999",
]
JOURNAL_DPI = 360
AXIS_COLOR = "#222222"
GRID_COLOR = "#d8d8d8"
BAR_COLOR = "#a9b0b6"
EXPANSION_LINE_COLOR = "#6b6b6b"


@dataclass
class StationRecord:
    code: str
    name: str
    latitude: float | None
    longitude: float | None
    start_key: int | None
    end_key: int | None


class StationIndex:
    def __init__(self, records: list[StationRecord]):
        by_code: dict[str, list[StationRecord]] = defaultdict(list)
        for rec in records:
            by_code[rec.code].append(rec)
        for records_for_code in by_code.values():
            records_for_code.sort(key=lambda rec: rec.start_key or 0)
        self.by_code = by_code
        self.records = records
        try:
            import numpy as np

            usable = [
                rec
                for rec in records
                if rec.latitude is not None and rec.longitude is not None
            ]
            self._np = np
            self._lats = np.array([rec.latitude for rec in usable], dtype=float)
            self._lons = np.array([rec.longitude for rec in usable], dtype=float)
            self._starts = np.array(
                [rec.start_key if rec.start_key is not None else -10**12 for rec in usable],
                dtype=np.int64,
            )
            self._ends = np.array(
                [rec.end_key if rec.end_key is not None else 10**12 for rec in usable],
                dtype=np.int64,
            )
        except Exception:
            self._np = None
            self._lats = None
            self._lons = None
            self._starts = None
            self._ends = None

    def get(self, code: str, event_time_key: int | None) -> StationRecord | None:
        records = self.by_code.get(code)
        if not records:
            return None
        if event_time_key is not None:
            for rec in records:
                if (rec.start_key is None or rec.start_key <= event_time_key) and (
                    rec.end_key is None or event_time_key < rec.end_key
                ):
                    return rec
        for rec in reversed(records):
            if rec.latitude is not None and rec.longitude is not None:
                return rec
        return records[-1]

    def active_station_counts_within_radii(
        self,
        latitude: float | None,
        longitude: float | None,
        event_time_key: int | None,
        radii_km: list[int],
    ) -> dict[int, int | None]:
        if latitude is None or longitude is None or event_time_key is None:
            return {radius: None for radius in radii_km}
        if self._np is not None and self._lats is not None and self._lons is not None:
            np = self._np
            active = (self._starts <= event_time_key) & (event_time_key < self._ends)
            if not bool(active.any()):
                return {radius: 0 for radius in radii_km}
            lat_arr = self._lats[active]
            lon_arr = self._lons[active]
            phi1 = math.radians(latitude)
            lambda1 = math.radians(longitude)
            phi2 = np.radians(lat_arr)
            lambda2 = np.radians(lon_arr)
            dphi = phi2 - phi1
            dlambda = lambda2 - lambda1
            a = np.sin(dphi / 2.0) ** 2 + math.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2.0) ** 2
            distances = 2.0 * 6371.0088 * np.arcsin(np.sqrt(a))
            return {radius: int(np.count_nonzero(distances <= radius)) for radius in radii_km}

        counts = {radius: 0 for radius in radii_km}
        for code in self.by_code:
            rec = self.get(code, event_time_key)
            if rec is None or rec.latitude is None or rec.longitude is None:
                continue
            distance = haversine_km(latitude, longitude, rec.latitude, rec.longitude)
            if distance is None:
                continue
            for radius in radii_km:
                if distance <= radius:
                    counts[radius] += 1
        return counts


@dataclass
class EventBuilder:
    event_id: str
    source_record_count: int
    source_line_number: int
    record_type: str
    year: int | None
    month: int | None
    day: int | None
    hour: int | None
    minute: int | None
    second: float | None
    latitude: float | None
    longitude: float | None
    depth_km: float | None
    magnitude: float | None
    magnitude_type: str
    magnitude2: float | None
    magnitude2_type: str
    max_intensity_header_code: str
    max_intensity_header_value: float | None
    header_station_count: int | None
    large_area_code: int | None
    small_area_code: int | None
    hypocenter_region: str
    source_flag: str
    obs_intensity_values: list[float] = field(default_factory=list)
    obs_measured_values: list[float] = field(default_factory=list)
    obs_best_values: list[float] = field(default_factory=list)
    obs_pga_total_values: list[float] = field(default_factory=list)
    obs_station_codes: list[str] = field(default_factory=list)
    obs_station_distance_by_code: dict[str, float] = field(default_factory=dict)
    unknown_intensity_count: int = 0
    nearest_station_distance_km: float | None = None
    max_intensity_station_code: str | None = None
    max_intensity_station_name: str | None = None
    max_intensity_station_distance_km: float | None = None
    max_intensity_station_score: float | None = None

    def finish(self, station_index: StationIndex | None = None) -> dict[str, object]:
        obs_count = len(self.obs_intensity_values) + self.unknown_intensity_count
        mean_class = mean(self.obs_intensity_values)
        max_from_obs = max(self.obs_intensity_values) if self.obs_intensity_values else None
        max_value = self.max_intensity_header_value
        if max_value is None:
            max_value = max_from_obs
        if station_index is None:
            active_station_counts = {radius: None for radius in STATION_COUNT_RADII_KM}
        else:
            active_station_counts = station_index.active_station_counts_within_radii(
                self.latitude,
                self.longitude,
                self.event_time_key(),
                STATION_COUNT_RADII_KM,
            )
        detected_station_counts = {
            radius: sum(1 for distance in self.obs_station_distance_by_code.values() if distance <= radius)
            for radius in STATION_COUNT_RADII_KM
        }

        row = {
            "event_id": self.event_id,
            "source_record_count": self.source_record_count,
            "source_line_number": self.source_line_number,
            "record_type": self.record_type,
            "year": self.year,
            "month": self.month,
            "day": self.day,
            "hour": self.hour,
            "minute": self.minute,
            "second": self.second,
            "latitude": self.latitude,
            "longitude": self.longitude,
            "depth_km": self.depth_km,
            "magnitude": self.magnitude,
            "magnitude_type": self.magnitude_type,
            "magnitude2": self.magnitude2,
            "magnitude2_type": self.magnitude2_type,
            "max_intensity_header_code": self.max_intensity_header_code,
            "max_intensity_header_value": self.max_intensity_header_value,
            "max_intensity_from_obs": max_from_obs,
            "max_intensity_value": max_value,
            "mean_intensity_class": mean_class,
            "mean_measured_intensity": mean(self.obs_measured_values),
            "mean_best_station_intensity": mean(self.obs_best_values),
            "max_best_station_intensity": max(self.obs_best_values) if self.obs_best_values else None,
            "n_measured_intensity": len(self.obs_measured_values),
            "n_station_records": obs_count,
            "n_valid_intensity": len(self.obs_intensity_values),
            "n_unknown_intensity": self.unknown_intensity_count,
            "n_unique_stations": len(set(self.obs_station_codes)),
            "mean_pga_total_gal": mean(self.obs_pga_total_values),
            "max_pga_total_gal": max(self.obs_pga_total_values) if self.obs_pga_total_values else None,
            "max_intensity_station_code": self.max_intensity_station_code,
            "max_intensity_station_name": self.max_intensity_station_name,
            "max_intensity_station_distance_km": self.max_intensity_station_distance_km,
            "nearest_station_distance_km": self.nearest_station_distance_km,
            "header_station_count": self.header_station_count,
            "large_area_code": self.large_area_code,
            "small_area_code": self.small_area_code,
            "hypocenter_region": self.hypocenter_region,
            "source_flag": self.source_flag,
        }
        for radius in STATION_COUNT_RADII_KM:
            active_count = active_station_counts.get(radius)
            detected_count = detected_station_counts.get(radius)
            ratio = None
            if active_count not in (None, 0):
                ratio = detected_count / active_count
            row[f"active_station_count_within_{radius}km"] = active_count
            row[f"detected_station_count_within_{radius}km"] = detected_count
            row[f"detected_to_active_station_ratio_within_{radius}km"] = ratio
        return row

    def event_time_key(self) -> int | None:
        if None in (self.year, self.month, self.day, self.hour, self.minute):
            return None
        return int(f"{self.year:04d}{self.month:02d}{self.day:02d}{self.hour:02d}{self.minute:02d}")


def as_ascii(raw: bytes) -> str:
    return raw.decode("ascii", errors="ignore")


def parse_int(raw: bytes) -> int | None:
    text = as_ascii(raw).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def parse_fixed_2(raw: bytes) -> float | None:
    text = as_ascii(raw).strip()
    if not text:
        return None
    try:
        return int(text) / 100.0
    except ValueError:
        return None


def parse_second(raw: bytes) -> float | None:
    text = as_ascii(raw).strip()
    if not text:
        return None
    try:
        return int(text) / 100.0
    except ValueError:
        return None


def parse_magnitude(raw: bytes) -> float | None:
    text = as_ascii(raw)
    if not text.strip():
        return None
    text = text.strip().upper()
    if len(text) == 1 and text.isdigit():
        return int(text) / 10.0
    if len(text) != 2:
        return None
    if text[0] == "-" and text[1].isdigit():
        return -int(text[1]) / 10.0
    if text[0].isalpha() and text[1].isdigit():
        return -(ord(text[0]) - ord("A") + 1) - int(text[1]) / 10.0
    if text.isdigit():
        return int(text) / 10.0
    return None


def parse_depth(raw: bytes) -> float | None:
    text = as_ascii(raw)
    if not text.strip():
        return None
    if text.endswith("  "):
        return parse_int(raw[:3])
    try:
        return int(text.strip()) / 100.0
    except ValueError:
        try:
            return float(text.strip())
        except ValueError:
            return None


def parse_lat_lon(deg_raw: bytes, min_raw: bytes) -> float | None:
    deg = parse_int(deg_raw)
    minute = parse_fixed_2(min_raw)
    if deg is None or minute is None:
        return None
    return deg + minute / 60.0


def parse_intensity_code(raw: bytes) -> tuple[str, float | None]:
    code = as_ascii(raw).strip()
    if not code:
        return "", None
    return code, INTENSITY_VALUE.get(code)


def parse_measured_intensity(raw: bytes) -> float | None:
    text = as_ascii(raw).strip()
    if not text or "/" in text:
        return None
    try:
        return int(text) / 10.0
    except ValueError:
        return None


def parse_scaled_int(raw: bytes, scale: float) -> float | None:
    text = as_ascii(raw).strip()
    if not text or "/" in text:
        return None
    try:
        return int(text) / scale
    except ValueError:
        return None


def decode_region(raw: bytes) -> str:
    return raw.decode("cp932", errors="replace").strip()


def parse_code_p_lat_lon(value: str, deg_digits: int) -> float | None:
    value = value.strip()
    if not value or not value.isdigit():
        return None
    try:
        deg = int(value[:deg_digits])
        minute = int(value[deg_digits : deg_digits + 2])
        return deg + minute / 60.0
    except ValueError:
        return None


def parse_code_p_time_key(value: str) -> int | None:
    value = value.strip()
    if not value or "9999" in value or len(value) < 12:
        return None
    try:
        return int(value[:12])
    except ValueError:
        return None


def parse_code_p_bytes(data: bytes) -> StationIndex:
    text = data.decode("euc_jp", errors="replace")
    records: list[StationRecord] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        code = parts[0].strip()
        if not code or code.startswith("._"):
            continue
        records.append(
            StationRecord(
                code=code,
                name=parts[1].strip(),
                latitude=parse_code_p_lat_lon(parts[2], 2),
                longitude=parse_code_p_lat_lon(parts[3], 3),
                start_key=parse_code_p_time_key(parts[4]),
                end_key=parse_code_p_time_key(parts[5]) if len(parts) > 5 else None,
            )
        )
    return StationIndex(records)


def load_station_index(data_dir: Path, code_p_path: Path | None = None) -> StationIndex | None:
    candidates = [code_p_path] if code_p_path is not None else [data_dir / "code_p.zip", data_dir / "code_p.dat"]
    for path in candidates:
        if path is None or not path.exists():
            continue
        if path.suffix.lower() == ".zip":
            with zipfile.ZipFile(path) as zf:
                for name in zf.namelist():
                    if Path(name).name.lower() == "code_p.dat":
                        return parse_code_p_bytes(zf.read(name))
        else:
            return parse_code_p_bytes(path.read_bytes())
    return None


def haversine_km(lat1: float | None, lon1: float | None, lat2: float | None, lon2: float | None) -> float | None:
    if None in (lat1, lon1, lat2, lon2):
        return None
    radius_km = 6371.0088
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    return 2.0 * radius_km * math.asin(math.sqrt(a))


def station_active_at_key(rec: StationRecord, time_key: int) -> bool:
    return (rec.start_key is None or rec.start_key <= time_key) and (
        rec.end_key is None or time_key < rec.end_key
    )


def active_station_records_at_year_end(station_index: StationIndex, year: int) -> list[StationRecord]:
    time_key = int(f"{year:04d}12312359")
    active: list[StationRecord] = []
    for records_for_code in station_index.by_code.values():
        records = [rec for rec in records_for_code if station_active_at_key(rec, time_key)]
        if not records:
            continue
        records.sort(key=lambda rec: rec.start_key or 0)
        active.append(records[-1])
    return active


def station_network_region(rec: StationRecord) -> str:
    try:
        prefix = int(rec.code[:2])
    except ValueError:
        return "Other"
    return NETWORK_REGION_BY_PREFIX.get(prefix, "Other")


def nearest_neighbor_distances_km(records: list[StationRecord]) -> list[float]:
    coords = [
        (float(rec.latitude), float(rec.longitude))
        for rec in records
        if rec.latitude is not None and rec.longitude is not None
    ]
    if len(coords) < 2:
        return []
    try:
        import numpy as np
        from scipy.spatial import cKDTree
    except Exception:
        distances: list[float] = []
        for idx, (lat1, lon1) in enumerate(coords):
            nearest = None
            for jdx, (lat2, lon2) in enumerate(coords):
                if idx == jdx:
                    continue
                distance = haversine_km(lat1, lon1, lat2, lon2)
                if distance is not None and (nearest is None or distance < nearest):
                    nearest = distance
            if nearest is not None:
                distances.append(nearest)
        return distances

    lat = np.radians(np.array([lat for lat, _ in coords], dtype=float))
    lon = np.radians(np.array([lon for _, lon in coords], dtype=float))
    xyz = np.column_stack(
        [
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ]
    )
    chord_distances, _ = cKDTree(xyz).query(xyz, k=2)
    nearest_chord = np.clip(chord_distances[:, 1], 0.0, 2.0)
    return list(2.0 * 6371.0088 * np.arcsin(nearest_chord / 2.0))


def station_network_density_summaries(
    station_index: StationIndex,
    years: list[int],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    yearly_rows: list[dict[str, object]] = []
    region_order = ["Japan", *[region.name for region in NETWORK_REGIONS]]

    for year in years:
        active = [
            rec
            for rec in active_station_records_at_year_end(station_index, year)
            if rec.latitude is not None and rec.longitude is not None
        ]
        by_region: dict[str, list[StationRecord]] = {region_name: [] for region_name in region_order}
        by_region["Japan"] = list(active)
        for rec in active:
            region_name = station_network_region(rec)
            if region_name in by_region:
                by_region[region_name].append(rec)

        for region_name in region_order:
            rows = by_region[region_name]
            area = NETWORK_REGION_AREA[region_name]
            nearest_distances = nearest_neighbor_distances_km(rows)
            yearly_rows.append(
                {
                    "year": year,
                    "region": region_name,
                    "active_station_count": len(rows),
                    "area_km2": area,
                    "density_per_10000km2": len(rows) / area * 10000.0 if area else None,
                    "mean_nearest_neighbor_distance_km": mean(nearest_distances),
                    "median_nearest_neighbor_distance_km": qtile(nearest_distances, 0.50),
                    "q25_nearest_neighbor_distance_km": qtile(nearest_distances, 0.25),
                    "q75_nearest_neighbor_distance_km": qtile(nearest_distances, 0.75),
                    "q90_nearest_neighbor_distance_km": qtile(nearest_distances, 0.90),
                }
            )

    period_rows: list[dict[str, object]] = []
    for _, _, plabel in PERIODS:
        start_year, end_year = [int(value) for value in plabel.split("-")]
        for region_name in region_order:
            rows = [
                row
                for row in yearly_rows
                if row["region"] == region_name and start_year <= int(row["year"]) <= end_year
            ]
            if not rows:
                continue
            first = min(rows, key=lambda row: int(row["year"]))
            last = max(rows, key=lambda row: int(row["year"]))
            period_rows.append(
                {
                    "period": plabel,
                    "region": region_name,
                    "first_year": first["year"],
                    "last_year": last["year"],
                    "active_station_count_first_year": first["active_station_count"],
                    "active_station_count_last_year": last["active_station_count"],
                    "active_station_count_change": (
                        int(last["active_station_count"]) - int(first["active_station_count"])
                    ),
                    "density_first_year_per_10000km2": first["density_per_10000km2"],
                    "density_last_year_per_10000km2": last["density_per_10000km2"],
                    "density_change_per_10000km2": (
                        numeric(last["density_per_10000km2"]) - numeric(first["density_per_10000km2"])
                        if numeric(last["density_per_10000km2"]) is not None
                        and numeric(first["density_per_10000km2"]) is not None
                        else None
                    ),
                    "median_nn_distance_first_year_km": first["median_nearest_neighbor_distance_km"],
                    "median_nn_distance_last_year_km": last["median_nearest_neighbor_distance_km"],
                    "mean_active_station_count": mean(numeric(row["active_station_count"]) for row in rows),
                    "mean_density_per_10000km2": mean(numeric(row["density_per_10000km2"]) for row in rows),
                    "mean_median_nearest_neighbor_distance_km": mean(
                        numeric(row["median_nearest_neighbor_distance_km"]) for row in rows
                    ),
                }
            )

    selected_years = [year for year in [1980, 1995, 1996, 2003, 2010, 2022] if year in years]
    selected_rows = [
        row
        for row in yearly_rows
        if int(row["year"]) in selected_years and row["region"] in region_order
    ]
    return yearly_rows, period_rows, selected_rows


def parse_source_record(line: bytes, event_id: str, line_number: int) -> EventBuilder:
    padded = line.ljust(96, b" ")
    max_code, max_value = parse_intensity_code(padded[61:62])

    return EventBuilder(
        event_id=event_id,
        source_record_count=1,
        source_line_number=line_number,
        record_type=as_ascii(padded[0:1]),
        year=parse_int(padded[1:5]),
        month=parse_int(padded[5:7]),
        day=parse_int(padded[7:9]),
        hour=parse_int(padded[9:11]),
        minute=parse_int(padded[11:13]),
        second=parse_second(padded[13:17]),
        latitude=parse_lat_lon(padded[21:24], padded[24:28]),
        longitude=parse_lat_lon(padded[32:36], padded[36:40]),
        depth_km=parse_depth(padded[44:49]),
        magnitude=parse_magnitude(padded[52:54]),
        magnitude_type=as_ascii(padded[54:55]).strip(),
        magnitude2=parse_magnitude(padded[55:57]),
        magnitude2_type=as_ascii(padded[57:58]).strip(),
        max_intensity_header_code=max_code,
        max_intensity_header_value=max_value,
        header_station_count=parse_int(padded[90:95]),
        large_area_code=parse_int(padded[64:65]),
        small_area_code=parse_int(padded[65:68]),
        hypocenter_region=decode_region(padded[68:90]),
        source_flag=as_ascii(padded[95:96]).strip(),
    )


def add_observation_record(event: EventBuilder, line: bytes, station_index: StationIndex | None = None) -> None:
    padded = line.ljust(96, b" ")
    station_code = as_ascii(padded[0:7]).strip()
    if station_code:
        event.obs_station_codes.append(station_code)

    _, intensity_value = parse_intensity_code(padded[18:19])
    if intensity_value is None:
        event.unknown_intensity_count += 1
    else:
        event.obs_intensity_values.append(intensity_value)

    measured = parse_measured_intensity(padded[20:22])
    if measured is not None:
        event.obs_measured_values.append(measured)

    best_value = measured if measured is not None else intensity_value
    if best_value is not None:
        event.obs_best_values.append(best_value)

    pga_total = parse_scaled_int(padded[29:34], 10.0)
    if pga_total is not None:
        event.obs_pga_total_values.append(pga_total)

    if station_index is not None and station_code:
        station = station_index.get(station_code, event.event_time_key())
        if station is not None:
            distance = haversine_km(event.latitude, event.longitude, station.latitude, station.longitude)
            if distance is not None:
                previous_distance = event.obs_station_distance_by_code.get(station.code)
                if previous_distance is None or distance < previous_distance:
                    event.obs_station_distance_by_code[station.code] = distance
                if event.nearest_station_distance_km is None or distance < event.nearest_station_distance_km:
                    event.nearest_station_distance_km = distance
                if best_value is not None and (
                    event.max_intensity_station_score is None
                    or best_value > event.max_intensity_station_score
                    or (
                        abs(best_value - event.max_intensity_station_score) < 1e-9
                        and (
                            event.max_intensity_station_distance_km is None
                            or distance < event.max_intensity_station_distance_km
                        )
                    )
                ):
                    event.max_intensity_station_score = best_value
                    event.max_intensity_station_distance_km = distance
                    event.max_intensity_station_code = station.code
                    event.max_intensity_station_name = station.name


def parse_catalog_zip(path: Path, station_index: StationIndex | None = None) -> list[dict[str, object]]:
    events: list[dict[str, object]] = []
    with zipfile.ZipFile(path) as zf:
        names = [name for name in zf.namelist() if not name.endswith("/")]
        if not names:
            return events
        data = zf.read(names[0])

    current: EventBuilder | None = None
    event_counter = 0
    for line_number, raw_line in enumerate(data.splitlines(), start=1):
        if not raw_line:
            continue
        first = raw_line[:1]
        if first in SOURCE_RECORD_TYPES:
            if current is None or current.obs_intensity_values or current.unknown_intensity_count:
                if current is not None:
                    events.append(current.finish(station_index))
                event_counter += 1
                event_id = f"{path.stem}_{event_counter:06d}"
                current = parse_source_record(raw_line, event_id, line_number)
            else:
                current.source_record_count += 1
            continue

        if current is None:
            continue
        add_observation_record(current, raw_line, station_index)

    if current is not None:
        events.append(current.finish(station_index))
    return events


def mean(values: Iterable[float | None]) -> float | None:
    vals = [v for v in values if v is not None and math.isfinite(v)]
    if not vals:
        return None
    return sum(vals) / len(vals)


def stdev(values: Iterable[float | None]) -> float | None:
    vals = [v for v in values if v is not None and math.isfinite(v)]
    if len(vals) < 2:
        return None
    return statistics.stdev(vals)


def qtile(values: Iterable[float | None], q: float) -> float | None:
    vals = sorted(v for v in values if v is not None and math.isfinite(v))
    if not vals:
        return None
    idx = (len(vals) - 1) * q
    low = int(math.floor(idx))
    high = int(math.ceil(idx))
    if low == high:
        return vals[low]
    return vals[low] * (high - idx) + vals[high] * (idx - low)


def ratio_le(values: Iterable[float | None], threshold: float) -> float | None:
    vals = [v for v in values if v is not None and math.isfinite(v)]
    if not vals:
        return None
    return sum(v <= threshold for v in vals) / len(vals)


def mag_bin(mag: float | None) -> str | None:
    if mag is None:
        return None
    for low, high, label in MAG_BINS:
        if low <= mag < high:
            return label
    return None


def period_label(year: int | None) -> str | None:
    if year is None:
        return None
    for start, end, label in PERIODS:
        if start <= year <= end:
            return label
    return None


def five_year_period_label(year: int | None) -> str | None:
    if year is None:
        return None
    for start, end, label in FIVE_YEAR_PERIODS:
        if start <= year <= end:
            return label
    return None


def distance_bin(distance_km: float | None) -> str | None:
    if distance_km is None:
        return None
    for low, high, label in DISTANCE_BINS:
        if low <= distance_km < high:
            return label
    return None


def is_land_region_approx(region: object) -> bool:
    """Approximate inland epicentral regions from the JMA region name.

    Sugiyama et al. (2020) separate inland and offshore events.  The local
    catalog files include the Japanese hypocentral region name, but no polygon
    layer.  For a reproducible lightweight filter, terms that denote seas,
    offshore areas, bays, straits, and near-island offshore regions are treated
    as non-land.  Use GIS polygons if a strict inland/offshore split is needed.
    """

    text = "" if region is None else str(region)
    sea_terms = (
        "沖",
        "灘",
        "湾",
        "海峡",
        "海域",
        "近海",
        "はるか沖",
        "遠地",
        "列島",
        "島近海",
        "島沖",
        "南方沖",
        "東方沖",
        "西方沖",
        "北方沖",
    )
    return not any(term in text for term in sea_terms)


def numeric(value: object) -> float | None:
    if value is None or value == "":
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def summarize_events(events: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_year_bin: dict[tuple[int, str], list[dict[str, object]]] = defaultdict(list)
    by_period_bin: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)

    for row in events:
        year = row["year"]
        mlabel = mag_bin(numeric(row["magnitude"]))
        plabel = period_label(int(year)) if isinstance(year, int) else None
        if isinstance(year, int) and mlabel is not None:
            by_year_bin[(year, mlabel)].append(row)
        if plabel is not None and mlabel is not None:
            by_period_bin[(plabel, mlabel)].append(row)

    yearly = [
        make_group_summary({"year": year, "mag_bin": mlabel}, rows)
        for (year, mlabel), rows in sorted(by_year_bin.items())
    ]
    period = [
        make_group_summary({"period": plabel, "mag_bin": mlabel}, rows)
        for (plabel, mlabel), rows in sorted(
            by_period_bin.items(), key=lambda item: (period_sort_key(item[0][0]), item[0][1])
        )
    ]
    return yearly, period


def make_group_summary(keys: dict[str, object], rows: list[dict[str, object]]) -> dict[str, object]:
    mean_values = [numeric(r["mean_intensity_class"]) for r in rows]
    measured_values = [numeric(r["mean_measured_intensity"]) for r in rows]
    best_values = [numeric(r["mean_best_station_intensity"]) for r in rows]
    max_values = [numeric(r["max_intensity_value"]) for r in rows]
    max_best_values = [numeric(r["max_best_station_intensity"]) for r in rows]
    mags = [numeric(r["magnitude"]) for r in rows]
    depths = [numeric(r["depth_km"]) for r in rows]
    station_counts = [numeric(r["n_valid_intensity"]) for r in rows]
    nearest_distances = [numeric(r["nearest_station_distance_km"]) for r in rows]
    max_station_distances = [numeric(r["max_intensity_station_distance_km"]) for r in rows]

    out = dict(keys)
    out.update(
        {
            "n_events": len(rows),
            "mean_magnitude": mean(mags),
            "mean_depth_km": mean(depths),
            "mean_event_mean_intensity": mean(mean_values),
            "sd_event_mean_intensity": stdev(mean_values),
            "q25_event_mean_intensity": qtile(mean_values, 0.25),
            "q75_event_mean_intensity": qtile(mean_values, 0.75),
            "mean_event_mean_measured_intensity": mean(measured_values),
            "n_events_with_measured_intensity": sum(v is not None for v in measured_values),
            "mean_event_mean_best_station_intensity": mean(best_values),
            "sd_event_mean_best_station_intensity": stdev(best_values),
            "mean_event_max_intensity": mean(max_values),
            "sd_event_max_intensity": stdev(max_values),
            "mean_event_max_best_station_intensity": mean(max_best_values),
            "sd_event_max_best_station_intensity": stdev(max_best_values),
            "mean_valid_station_records": mean(station_counts),
            "mean_nearest_station_distance_km": mean(nearest_distances),
            "median_nearest_station_distance_km": qtile(nearest_distances, 0.50),
            "q25_nearest_station_distance_km": qtile(nearest_distances, 0.25),
            "q75_nearest_station_distance_km": qtile(nearest_distances, 0.75),
            "share_nearest_station_within_10km": ratio_le(nearest_distances, 10.0),
            "share_nearest_station_within_20km": ratio_le(nearest_distances, 20.0),
            "share_nearest_station_within_50km": ratio_le(nearest_distances, 50.0),
            "mean_max_intensity_station_distance_km": mean(max_station_distances),
            "median_max_intensity_station_distance_km": qtile(max_station_distances, 0.50),
        }
    )
    for radius in STATION_COUNT_RADII_KM:
        active_counts = [numeric(r.get(f"active_station_count_within_{radius}km")) for r in rows]
        detected_counts = [numeric(r.get(f"detected_station_count_within_{radius}km")) for r in rows]
        detection_ratios = [numeric(r.get(f"detected_to_active_station_ratio_within_{radius}km")) for r in rows]
        out.update(
            {
                f"mean_active_station_count_within_{radius}km": mean(active_counts),
                f"median_active_station_count_within_{radius}km": qtile(active_counts, 0.50),
                f"mean_detected_station_count_within_{radius}km": mean(detected_counts),
                f"median_detected_station_count_within_{radius}km": qtile(detected_counts, 0.50),
                f"mean_detected_to_active_station_ratio_within_{radius}km": mean(detection_ratios),
            }
        )
    return out


def select_spatial_station_count_columns(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        selected = {
            key: row.get(key)
            for key in ("year", "period", "mag_bin", "n_events", "mean_magnitude", "mean_depth_km")
            if key in row
        }
        for radius in STATION_COUNT_RADII_KM:
            for name in (
                "mean_active_station_count",
                "median_active_station_count",
                "mean_detected_station_count",
                "median_detected_station_count",
                "mean_detected_to_active_station_ratio",
            ):
                full_name = f"{name}_within_{radius}km"
                selected[full_name] = row.get(full_name)
        out.append(selected)
    return out


def period_sort_key(label: str) -> int:
    try:
        return int(label.split("-")[0])
    except ValueError:
        return 9999


def fit_line(rows: list[dict[str, object]], x_name: str, y_name: str) -> dict[str, object]:
    pairs = []
    for row in rows:
        x = numeric(row.get(x_name))
        y = numeric(row.get(y_name))
        if x is not None and y is not None:
            pairs.append((x, y))
    n = len(pairs)
    if n < 2:
        return {"n": n, "intercept": None, "slope": None, "r2": None}
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    syy = sum((y - ybar) ** 2 for y in ys)
    sxy = sum((x - xbar) * (y - ybar) for x, y in pairs)
    if sxx == 0:
        return {"n": n, "intercept": None, "slope": None, "r2": None}
    slope = sxy / sxx
    intercept = ybar - slope * xbar
    r2 = (sxy * sxy) / (sxx * syy) if syy else None
    return {"n": n, "intercept": intercept, "slope": slope, "r2": r2}


def regression_summaries(events: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for year in sorted({r["year"] for r in events if isinstance(r["year"], int)}):
        rows = [
            r
            for r in events
            if r["year"] == year
            and numeric(r["magnitude"]) is not None
            and numeric(r["mean_intensity_class"]) is not None
        ]
        mean_intensity = mean(numeric(r["mean_intensity_class"]) for r in rows)
        mean_best = mean(numeric(r["mean_best_station_intensity"]) for r in rows)
        mean_max = mean(numeric(r["max_intensity_value"]) for r in rows)
        fit_mean = fit_line(rows, "magnitude", "mean_intensity_class")
        fit_best = fit_line(rows, "magnitude", "mean_best_station_intensity")
        fit_max = fit_line(rows, "magnitude", "max_intensity_value")
        out.append(
            {
                "year": year,
                "n_events": len(rows),
                "mean_event_mean_intensity": mean_intensity,
                "mean_event_mean_best_station_intensity": mean_best,
                "mean_event_max_intensity": mean_max,
                "mean_intensity_intercept": fit_mean["intercept"],
                "mean_intensity_slope_per_magnitude": fit_mean["slope"],
                "mean_intensity_r2": fit_mean["r2"],
                "mean_best_station_intensity_intercept": fit_best["intercept"],
                "mean_best_station_intensity_slope_per_magnitude": fit_best["slope"],
                "mean_best_station_intensity_r2": fit_best["r2"],
                "max_intensity_intercept": fit_max["intercept"],
                "max_intensity_slope_per_magnitude": fit_max["slope"],
                "max_intensity_r2": fit_max["r2"],
            }
        )
    return out


def magbin_annual_trend_summaries(yearly: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    metrics = [
        ("mean_event_mean_best_station_intensity", "mean_station_intensity"),
        ("mean_event_mean_intensity", "class_mean_intensity"),
        ("mean_event_max_intensity", "observed_max_intensity"),
        ("median_nearest_station_distance_km", "median_shortest_epicentral_distance"),
        ("mean_valid_station_records", "mean_station_records"),
    ]
    for _, _, mlabel in MAG_BINS:
        group = [r for r in yearly if r["mag_bin"] == mlabel]
        for metric, metric_label in metrics:
            points = [
                (float(r["year"]), numeric(r.get(metric)))
                for r in group
                if numeric(r.get(metric)) is not None
            ]
            if len(points) < 3:
                continue
            xs = [p[0] for p in points]
            ys = [p[1] for p in points if p[1] is not None]
            xbar = sum(xs) / len(xs)
            ybar = sum(ys) / len(ys)
            sxx = sum((x - xbar) ** 2 for x in xs)
            syy = sum((y - ybar) ** 2 for y in ys)
            sxy = sum((x - xbar) * (y - ybar) for x, y in points if y is not None)
            slope = sxy / sxx if sxx else None
            intercept = ybar - slope * xbar if slope is not None else None
            r = sxy / math.sqrt(sxx * syy) if sxx and syy else None
            pre_values = [numeric(r.get(metric)) for r in group if int(r["year"]) <= 1995]
            post_values = [numeric(r.get(metric)) for r in group if int(r["year"]) >= 1996]
            rows.append(
                {
                    "mag_bin": mlabel,
                    "metric": metric_label,
                    "n_year_bins": len(points),
                    "first_year": int(min(xs)),
                    "last_year": int(max(xs)),
                    "slope_per_year": slope,
                    "slope_per_decade": slope * 10.0 if slope is not None else None,
                    "intercept": intercept,
                    "r": r,
                    "mean_1980_1995": mean(pre_values),
                    "mean_1996_2022": mean(post_values),
                    "post_minus_pre": (
                        mean(post_values) - mean(pre_values)
                        if mean(post_values) is not None and mean(pre_values) is not None
                        else None
                    ),
                }
            )
    return rows


def yearly_slope_period_summaries(yearly_regression: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = [
        ("mean_best_station_intensity_slope_per_magnitude", "mean_station_intensity"),
        ("mean_intensity_slope_per_magnitude", "class_mean_intensity"),
        ("max_intensity_slope_per_magnitude", "observed_max_intensity"),
    ]
    rows: list[dict[str, object]] = []
    for _, _, plabel in PERIODS:
        group = [
            r
            for r in yearly_regression
            if period_label(int(r["year"]) if isinstance(r["year"], int) else None) == plabel
        ]
        for metric, metric_label in metrics:
            values = [numeric(r.get(metric)) for r in group]
            values_clean = [v for v in values if v is not None]
            rows.append(
                {
                    "period": plabel,
                    "metric": metric_label,
                    "n_years": len(values_clean),
                    "mean_slope_per_magnitude": mean(values_clean),
                    "median_slope_per_magnitude": qtile(values_clean, 0.50),
                    "sd_slope_per_magnitude": stdev(values_clean),
                    "min_slope_per_magnitude": min(values_clean) if values_clean else None,
                    "max_slope_per_magnitude": max(values_clean) if values_clean else None,
                }
            )
    return rows


def period_event_regression_summaries(events: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = [
        ("mean_best_station_intensity", "mean_station_intensity"),
        ("mean_intensity_class", "class_mean_intensity"),
        ("max_intensity_value", "observed_max_intensity"),
    ]
    rows: list[dict[str, object]] = []
    for min_magnitude in [3.0, 4.0]:
        for _, _, plabel in PERIODS:
            group = [
                r
                for r in events
                if period_label(r["year"] if isinstance(r["year"], int) else None) == plabel
                and numeric(r["magnitude"]) is not None
                and numeric(r["magnitude"]) >= min_magnitude
            ]
            for metric, metric_label in metrics:
                rows_for_fit = [r for r in group if numeric(r.get(metric)) is not None]
                fit = fit_line(rows_for_fit, "magnitude", metric)
                rows.append(
                    {
                        "period": plabel,
                        "min_magnitude": min_magnitude,
                        "metric": metric_label,
                        "n_events": fit["n"],
                        "intercept": fit["intercept"],
                        "slope_per_magnitude": fit["slope"],
                        "r2": fit["r2"],
                        "mean_metric": mean(numeric(r.get(metric)) for r in rows_for_fit),
                    }
                )
    return rows


def adjusted_event_model_summaries(events: list[dict[str, object]]) -> list[dict[str, object]]:
    try:
        import numpy as np
    except Exception:
        return []

    period_names = [label for _, _, label in PERIODS]
    dummy_periods = period_names[1:]
    metrics = [
        ("mean_best_station_intensity", "mean_station_intensity"),
        ("mean_intensity_class", "class_mean_intensity"),
        ("max_intensity_value", "observed_max_intensity"),
    ]
    rows: list[dict[str, object]] = []
    base_rows = [
        r
        for r in events
        if numeric(r["magnitude"]) is not None
        and numeric(r["magnitude"]) >= 4.0
        and numeric(r.get("nearest_station_distance_km")) is not None
        and numeric(r.get("n_valid_intensity")) is not None
        and period_label(r["year"] if isinstance(r["year"], int) else None) in period_names
    ]
    for metric, metric_label in metrics:
        model_rows = [r for r in base_rows if numeric(r.get(metric)) is not None]
        if len(model_rows) < 10:
            continue
        x_rows = []
        y_vals = []
        for row in model_rows:
            plabel = period_label(row["year"] if isinstance(row["year"], int) else None)
            x_rows.append(
                [
                    1.0,
                    float(numeric(row["magnitude"])),
                    math.log1p(float(numeric(row["nearest_station_distance_km"]))),
                    math.log1p(float(numeric(row["n_valid_intensity"]))),
                    *[1.0 if plabel == dummy else 0.0 for dummy in dummy_periods],
                ]
            )
            y_vals.append(float(numeric(row[metric])))
        x = np.asarray(x_rows, dtype=float)
        y = np.asarray(y_vals, dtype=float)
        coef, *_ = np.linalg.lstsq(x, y, rcond=None)
        pred = x @ coef
        ss_res = float(((y - pred) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = None if ss_tot == 0 else 1.0 - ss_res / ss_tot
        out = {
            "metric": metric_label,
            "n_events": len(model_rows),
            "r2": r2,
            "intercept": float(coef[0]),
            "magnitude_coef": float(coef[1]),
            "log_nearest_distance_coef": float(coef[2]),
            "log_station_records_coef": float(coef[3]),
        }
        for idx, dummy in enumerate(dummy_periods, start=4):
            out[f"period_{dummy}_coef_vs_1980_1995"] = float(coef[idx])
        rows.append(out)
    return rows


def distance_detection_summaries(
    events: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_year_bin: dict[tuple[int, str], list[dict[str, object]]] = defaultdict(list)
    by_period_distance: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)

    for row in events:
        year = row["year"]
        year_int = year if isinstance(year, int) else None
        mlabel = mag_bin(numeric(row["magnitude"]))
        if year_int is not None and mlabel is not None:
            by_year_bin[(year_int, mlabel)].append(row)

        period = five_year_period_label(year_int)
        mag = numeric(row["magnitude"])
        nearest = numeric(row.get("nearest_station_distance_km"))
        dlabel = distance_bin(nearest)
        if period is not None and mag is not None and mag >= 4.0 and dlabel is not None:
            by_period_distance[(period, dlabel)].append(row)

    yearly_rows: list[dict[str, object]] = []
    for (year, mlabel), rows in sorted(by_year_bin.items()):
        nearest_distances = [numeric(r.get("nearest_station_distance_km")) for r in rows]
        max_station_distances = [numeric(r.get("max_intensity_station_distance_km")) for r in rows]
        yearly_rows.append(
            {
                "year": year,
                "mag_bin": mlabel,
                "n_events": len(rows),
                "n_events_with_distance": sum(v is not None for v in nearest_distances),
                "mean_nearest_station_distance_km": mean(nearest_distances),
                "median_nearest_station_distance_km": qtile(nearest_distances, 0.50),
                "q25_nearest_station_distance_km": qtile(nearest_distances, 0.25),
                "q75_nearest_station_distance_km": qtile(nearest_distances, 0.75),
                "share_nearest_station_within_10km": ratio_le(nearest_distances, 10.0),
                "share_nearest_station_within_20km": ratio_le(nearest_distances, 20.0),
                "share_nearest_station_within_50km": ratio_le(nearest_distances, 50.0),
                "mean_max_intensity_station_distance_km": mean(max_station_distances),
                "median_max_intensity_station_distance_km": qtile(max_station_distances, 0.50),
            }
        )

    period_rows: list[dict[str, object]] = []
    for period in [label for _, _, label in FIVE_YEAR_PERIODS]:
        period_events = [
            r
            for r in events
            if five_year_period_label(r["year"] if isinstance(r["year"], int) else None) == period
            and numeric(r["magnitude"]) is not None
            and numeric(r["magnitude"]) >= 4.0
            and numeric(r.get("nearest_station_distance_km")) is not None
        ]
        total = len(period_events)
        for _, _, dlabel in DISTANCE_BINS:
            rows = by_period_distance.get((period, dlabel), [])
            period_rows.append(
                {
                    "period": period,
                    "distance_bin_km": dlabel,
                    "n_events": len(rows),
                    "total_period_events_m4plus": total,
                    "share": len(rows) / total if total else None,
                }
            )
    return yearly_rows, period_rows


def yearly_shortest_distance_threshold_summary(
    events: list[dict[str, object]],
    min_magnitude: float = 4.0,
) -> list[dict[str, object]]:
    by_year: dict[int, list[dict[str, object]]] = defaultdict(list)
    for row in events:
        year = row["year"]
        magnitude = numeric(row.get("magnitude"))
        nearest = numeric(row.get("nearest_station_distance_km"))
        if isinstance(year, int) and magnitude is not None and magnitude >= min_magnitude and nearest is not None:
            by_year[year].append(row)

    rows: list[dict[str, object]] = []
    for year, group in sorted(by_year.items()):
        nearest_distances = [numeric(r.get("nearest_station_distance_km")) for r in group]
        rows.append(
            {
                "year": year,
                "min_magnitude": min_magnitude,
                "n_events": len(group),
                "mean_nearest_station_distance_km": mean(nearest_distances),
                "median_nearest_station_distance_km": qtile(nearest_distances, 0.50),
                "q25_nearest_station_distance_km": qtile(nearest_distances, 0.25),
                "q75_nearest_station_distance_km": qtile(nearest_distances, 0.75),
                "share_nearest_station_within_10km": ratio_le(nearest_distances, 10.0),
                "share_nearest_station_within_20km": ratio_le(nearest_distances, 20.0),
                "share_nearest_station_within_50km": ratio_le(nearest_distances, 50.0),
                "mean_event_max_intensity": mean(numeric(r.get("max_intensity_value")) for r in group),
                "mean_valid_station_records": mean(numeric(r.get("n_valid_intensity")) for r in group),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: format_cell(row.get(k)) for k in fieldnames})


def format_cell(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return value


def configure_journal_matplotlib(plt) -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": AXIS_COLOR,
            "axes.labelcolor": AXIS_COLOR,
            "axes.linewidth": 0.8,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "xtick.color": AXIS_COLOR,
            "ytick.color": AXIS_COLOR,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "legend.frameon": False,
            "legend.fontsize": 8.0,
            "legend.title_fontsize": 8.2,
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": JOURNAL_DPI,
        }
    )


def style_journal_axis(ax, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(AXIS_COLOR)
    ax.spines["bottom"].set_color(AXIS_COLOR)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.55, alpha=0.85)
    else:
        ax.grid(False)
    ax.set_axisbelow(True)


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.08,
        1.03,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color=AXIS_COLOR,
    )


def add_network_epoch_guides(ax, annotate: bool = True) -> None:
    for year, text in [(1995.5, "1995-1996"), (2003.5, "2004")]:
        ax.axvline(year, color=EXPANSION_LINE_COLOR, linewidth=0.8, linestyle=(0, (3, 3)), alpha=0.65)
        if not annotate:
            continue
        ymax = ax.get_ylim()[1]
        ymin = ax.get_ylim()[0]
        ax.text(
            year + 0.25,
            ymin + (ymax - ymin) * 0.92,
            text,
            rotation=90,
            ha="left",
            va="top",
            fontsize=7.5,
            color=EXPANSION_LINE_COLOR,
        )


def save_journal_png(fig, path: Path) -> None:
    fig.savefig(path, dpi=JOURNAL_DPI, bbox_inches="tight", facecolor="white")


def svg_escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def render_line_chart(
    path: Path,
    yearly: list[dict[str, object]],
    value_name: str,
    title: str,
    y_label: str,
    min_events: int,
) -> None:
    years = sorted({int(r["year"]) for r in yearly if isinstance(r["year"], int)})
    labels = [label for _, _, label in MAG_BINS]
    series: dict[str, list[tuple[int, float]]] = {label: [] for label in labels}
    for row in yearly:
        if int(row["n_events"]) < min_events:
            continue
        value = numeric(row.get(value_name))
        if value is None:
            continue
        series[str(row["mag_bin"])].append((int(row["year"]), value))

    width, height = 1120, 620
    margin_left, margin_right, margin_top, margin_bottom = 76, 230, 58, 76
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    min_year = min(years)
    max_year = max(years)
    y_min, y_max = 0.8, 7.1

    def xmap(year: int) -> float:
        return margin_left + (year - min_year) / (max_year - min_year) * plot_w

    def ymap(value: float) -> float:
        return margin_top + (y_max - value) / (y_max - y_min) * plot_h

    parts = svg_header(width, height, title)
    parts.append(f'<text x="{width/2}" y="30" text-anchor="middle" class="title">{svg_escape(title)}</text>')
    parts.append(f'<text x="22" y="{height/2}" transform="rotate(-90 22 {height/2})" text-anchor="middle" class="axis-label">{svg_escape(y_label)}</text>')
    parts.append(f'<text x="{margin_left + plot_w/2}" y="{height-20}" text-anchor="middle" class="axis-label">Year</text>')

    for y in range(1, 8):
        yp = ymap(float(y))
        parts.append(f'<line x1="{margin_left}" y1="{yp:.2f}" x2="{margin_left + plot_w}" y2="{yp:.2f}" class="grid"/>')
        parts.append(f'<text x="{margin_left - 10}" y="{yp + 4:.2f}" text-anchor="end" class="tick">{y}</text>')
    for year in range((min_year // 5) * 5, max_year + 1, 5):
        if year < min_year:
            continue
        xp = xmap(year)
        parts.append(f'<line x1="{xp:.2f}" y1="{margin_top}" x2="{xp:.2f}" y2="{margin_top + plot_h}" class="grid xgrid"/>')
        parts.append(f'<text x="{xp:.2f}" y="{margin_top + plot_h + 22}" text-anchor="middle" class="tick">{year}</text>')

    parts.append(f'<rect x="{margin_left}" y="{margin_top}" width="{plot_w}" height="{plot_h}" class="frame"/>')

    for idx, label in enumerate(labels):
        points = series[label]
        if not points:
            continue
        color = PALETTE[idx % len(PALETTE)]
        point_str = " ".join(f"{xmap(year):.2f},{ymap(value):.2f}" for year, value in points)
        parts.append(f'<polyline points="{point_str}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        for year, value in points:
            parts.append(f'<circle cx="{xmap(year):.2f}" cy="{ymap(value):.2f}" r="3.2" fill="{color}"><title>{year} {label}: {value:.2f}</title></circle>')
        lx, ly = margin_left + plot_w + 24, margin_top + 24 + idx * 26
        parts.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 24}" y2="{ly}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="{lx + 32}" y="{ly + 4}" class="legend">{svg_escape(label)}</text>')

    parts.append(f'<text x="{margin_left}" y="{height - 48}" class="note">Points with fewer than {min_events} events in a year-bin are omitted.</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def render_heatmap(path: Path, yearly: list[dict[str, object]], min_events: int) -> None:
    years = sorted({int(r["year"]) for r in yearly if isinstance(r["year"], int)})
    labels = [label for _, _, label in MAG_BINS]
    values = {}
    for row in yearly:
        if int(row["n_events"]) < min_events:
            continue
        value = numeric(row["mean_event_mean_best_station_intensity"])
        if value is not None:
            values[(int(row["year"]), str(row["mag_bin"]))] = value

    width, height = 1180, 430
    margin_left, margin_right, margin_top, margin_bottom = 110, 190, 58, 70
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    cell_w = plot_w / len(years)
    cell_h = plot_h / len(labels)
    vmin, vmax = 1.0, 5.0

    parts = svg_header(width, height, "Mean event-average intensity by year and magnitude bin")
    parts.append('<text x="590" y="30" text-anchor="middle" class="title">Mean event-average intensity by year and magnitude bin</text>')
    for i, label in enumerate(labels):
        y = margin_top + i * cell_h + cell_h / 2
        parts.append(f'<text x="{margin_left - 10}" y="{y + 4:.2f}" text-anchor="end" class="tick">{svg_escape(label)}</text>')
    for j, year in enumerate(years):
        x = margin_left + j * cell_w + cell_w / 2
        if year % 5 == 0 or year == years[0] or year == years[-1]:
            parts.append(f'<text x="{x:.2f}" y="{margin_top + plot_h + 22}" text-anchor="middle" class="tick">{year}</text>')

    for i, label in enumerate(labels):
        for j, year in enumerate(years):
            value = values.get((year, label))
            color = "#f4f4f4" if value is None else heat_color(value, vmin, vmax)
            x = margin_left + j * cell_w
            y = margin_top + i * cell_h
            title = "no data" if value is None else f"{year} {label}: {value:.2f}"
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{cell_w + 0.2:.2f}" height="{cell_h + 0.2:.2f}" fill="{color}" stroke="#ffffff" stroke-width="0.7"><title>{svg_escape(title)}</title></rect>')

    legend_x = margin_left + plot_w + 36
    legend_y = margin_top + 10
    parts.append(f'<text x="{legend_x}" y="{legend_y - 8}" class="legend">Intensity</text>')
    for k in range(80):
        frac = k / 79
        color = heat_color(vmin + frac * (vmax - vmin), vmin, vmax)
        y = legend_y + (1 - frac) * 220
        parts.append(f'<rect x="{legend_x}" y="{y:.2f}" width="22" height="3.2" fill="{color}" stroke="none"/>')
    for val in [1, 2, 3, 4, 5]:
        y = legend_y + (1 - (val - vmin) / (vmax - vmin)) * 220
        parts.append(f'<text x="{legend_x + 30}" y="{y + 4:.2f}" class="tick">{val}</text>')

    parts.append(f'<text x="{margin_left}" y="{height - 24}" class="note">Cells with fewer than {min_events} events are shown as gray.</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def render_period_relation(path: Path, period: list[dict[str, object]], min_events: int) -> None:
    rows = [r for r in period if int(r["n_events"]) >= min_events]
    labels = [label for _, _, label in MAG_BINS]
    centers = {label: (low + min(high, 8.0)) / 2 for low, high, label in MAG_BINS}
    periods = [label for _, _, label in PERIODS]
    series = {p: [] for p in periods}
    for row in rows:
        value = numeric(row["mean_event_mean_best_station_intensity"])
        if value is None:
            continue
        series[str(row["period"])].append((centers[str(row["mag_bin"])], value, str(row["mag_bin"]), int(row["n_events"])))

    width, height = 820, 560
    margin_left, margin_right, margin_top, margin_bottom = 78, 190, 58, 72
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    x_min, x_max = 2.9, 7.7
    y_min, y_max = 0.8, 5.2

    def xmap(x: float) -> float:
        return margin_left + (x - x_min) / (x_max - x_min) * plot_w

    def ymap(y: float) -> float:
        return margin_top + (y_max - y) / (y_max - y_min) * plot_h

    parts = svg_header(width, height, "Magnitude relation by period")
    parts.append('<text x="410" y="30" text-anchor="middle" class="title">Magnitude relation by period</text>')
    parts.append(f'<text x="{margin_left + plot_w/2}" y="{height-20}" text-anchor="middle" class="axis-label">Magnitude bin center</text>')
    parts.append(f'<text x="22" y="{height/2}" transform="rotate(-90 22 {height/2})" text-anchor="middle" class="axis-label">Mean event-average intensity</text>')

    for y in [1, 2, 3, 4, 5]:
        yp = ymap(y)
        parts.append(f'<line x1="{margin_left}" y1="{yp:.2f}" x2="{margin_left + plot_w}" y2="{yp:.2f}" class="grid"/>')
        parts.append(f'<text x="{margin_left - 10}" y="{yp + 4:.2f}" text-anchor="end" class="tick">{y}</text>')
    for x in [3, 4, 5, 6, 7]:
        xp = xmap(x)
        parts.append(f'<line x1="{xp:.2f}" y1="{margin_top}" x2="{xp:.2f}" y2="{margin_top + plot_h}" class="grid xgrid"/>')
        parts.append(f'<text x="{xp:.2f}" y="{margin_top + plot_h + 22}" text-anchor="middle" class="tick">{x}</text>')
    parts.append(f'<rect x="{margin_left}" y="{margin_top}" width="{plot_w}" height="{plot_h}" class="frame"/>')

    for idx, plabel in enumerate(periods):
        points = sorted(series[plabel])
        if not points:
            continue
        color = PALETTE[idx % len(PALETTE)]
        point_str = " ".join(f"{xmap(x):.2f},{ymap(y):.2f}" for x, y, _, _ in points)
        parts.append(f'<polyline points="{point_str}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        for x, y, label, n_events in points:
            parts.append(f'<circle cx="{xmap(x):.2f}" cy="{ymap(y):.2f}" r="4" fill="{color}"><title>{svg_escape(plabel)} {svg_escape(label)} n={n_events}: {y:.2f}</title></circle>')
        lx, ly = margin_left + plot_w + 24, margin_top + 24 + idx * 26
        parts.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 24}" y2="{ly}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="{lx + 32}" y="{ly + 4}" class="legend">{svg_escape(plabel)}</text>')

    parts.append(f'<text x="{margin_left}" y="{height - 48}" class="note">Period-bin points with fewer than {min_events} events are omitted.</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def render_reference_style_pngs(
    out_dir: Path,
    suffix: str,
    yearly: list[dict[str, object]],
    yearly_regression: list[dict[str, object]],
    min_events: int,
) -> None:
    """Create matplotlib PNGs modeled on the supplied reference package."""

    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:
        print(f"Skipping PNG plots because matplotlib/numpy is unavailable: {exc}")
        return

    configure_journal_matplotlib(plt)
    suffix_part = f"_{suffix}" if suffix else ""
    bin_labels = ["4.0<=M<5.0", "5.0<=M<6.0"]

    def rows_for(label: str) -> list[dict[str, object]]:
        return sorted(
            [r for r in yearly if r["mag_bin"] == label and int(r["n_events"]) >= min_events],
            key=lambda r: int(r["year"]),
        )

    fig, ax1 = plt.subplots(figsize=(7.2, 4.4))
    for idx, label in enumerate(bin_labels):
        rows = rows_for(label)
        if not rows:
            continue
        ax1.plot(
            [int(r["year"]) for r in rows],
            [numeric(r["mean_event_max_intensity"]) for r in rows],
            marker="o",
            markersize=3.8,
            linewidth=1.8,
            color=PALETTE[idx],
            label=f"Mean max intensity {label}",
        )
    ax1.set_xlabel("Year")
    ax1.set_ylabel("Mean observed maximum intensity")
    ax1.set_ylim(0.5, 7.2)
    style_journal_axis(ax1)
    ax2 = ax1.twinx()
    width = 0.36
    offsets = [-width / 2, width / 2]
    for idx, (offset, label) in enumerate(zip(offsets, bin_labels)):
        rows = rows_for(label)
        if not rows:
            continue
        ax2.bar(
            [int(r["year"]) + offset for r in rows],
            [int(r["n_events"]) for r in rows],
            width=width,
            alpha=0.22,
            color=PALETTE[idx],
            label=f"Events {label}",
        )
    ax2.set_ylabel("Number of events")
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.grid(False)
    add_network_epoch_guides(ax1)
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="upper left", bbox_to_anchor=(0.01, 0.99))
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"paper_like_yearly_max_intensity{suffix_part}.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for idx, (_, _, label) in enumerate(MAG_BINS):
        rows = rows_for(label) if label in bin_labels else sorted(
            [r for r in yearly if r["mag_bin"] == label and int(r["n_events"]) >= min_events],
            key=lambda r: int(r["year"]),
        )
        values = [numeric(r["mean_event_mean_best_station_intensity"]) for r in rows]
        pairs = [(int(r["year"]), v) for r, v in zip(rows, values) if v is not None]
        if not pairs:
            continue
        ax.plot(
            [p[0] for p in pairs],
            [p[1] for p in pairs],
            marker="o",
            markersize=3.6,
            linewidth=1.7,
            color=PALETTE[idx % len(PALETTE)],
            label=label,
        )
    ax.set_xlabel("Year")
    ax.set_ylabel("Mean of event-mean station intensity")
    style_journal_axis(ax)
    add_network_epoch_guides(ax)
    ax.legend(title="Magnitude bin", ncols=2, loc="upper right")
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"yearly_mean_station_intensity_by_magbin{suffix_part}.png")
    plt.close(fig)

    bin_centers = {label: (low + min(high, 8.0)) / 2.0 for low, high, label in MAG_BINS}
    years = sorted({int(r["year"]) for r in yearly})
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    for year in years:
        rows = [r for r in yearly if int(r["year"]) == year and int(r["n_events"]) >= min_events]
        points = [
            (bin_centers[str(r["mag_bin"])], numeric(r["mean_event_max_intensity"]))
            for r in rows
            if str(r["mag_bin"]) in bin_centers and numeric(r["mean_event_max_intensity"]) is not None
        ]
        if len(points) >= 2:
            points.sort()
            ax.plot([p[0] for p in points], [p[1] for p in points], alpha=0.22, color="#4c566a", linewidth=0.9)
    ax.set_xlabel("Magnitude bin center")
    ax.set_ylabel("Annual mean observed maximum intensity")
    style_journal_axis(ax)
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"magnitude_vs_annual_mean_max_intensity{suffix_part}.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    slope_specs = [
        ("max_intensity_slope_per_magnitude", "Max intensity"),
        ("mean_best_station_intensity_slope_per_magnitude", "Mean station intensity"),
        ("mean_intensity_slope_per_magnitude", "Mean class intensity"),
    ]
    for idx, (column, label) in enumerate(slope_specs):
        pairs = [
            (int(row["year"]), numeric(row.get(column)))
            for row in yearly_regression
            if numeric(row.get(column)) is not None
        ]
        if not pairs:
            continue
        pairs.sort()
        ax.plot(
            [p[0] for p in pairs],
            [p[1] for p in pairs],
            marker="o",
            markersize=3.6,
            linewidth=1.7,
            color=PALETTE[idx],
            label=label,
        )
    ax.set_xlabel("Year")
    ax.set_ylabel("Regression slope: intensity per M unit")
    style_journal_axis(ax)
    add_network_epoch_guides(ax)
    ax.legend(ncols=3, loc="upper center", bbox_to_anchor=(0.5, 1.14))
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"yearly_regression_slope_intensity_vs_magnitude{suffix_part}.png")
    plt.close(fig)


def render_standard_png_charts(
    out_dir: Path,
    suffix: str,
    yearly: list[dict[str, object]],
    period: list[dict[str, object]],
    distance_yearly: list[dict[str, object]],
    distance_period: list[dict[str, object]],
    min_events: int,
) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:
        print(f"Skipping PNG plots because matplotlib/numpy is unavailable: {exc}")
        return

    configure_journal_matplotlib(plt)
    suffix_part = f"_{suffix}" if suffix else ""
    labels = [label for _, _, label in MAG_BINS]

    def sorted_yearly(label: str, value_name: str) -> tuple[list[int], list[float]]:
        rows = [
            r
            for r in yearly
            if r["mag_bin"] == label and int(r["n_events"]) >= min_events and numeric(r.get(value_name)) is not None
        ]
        rows.sort(key=lambda r: int(r["year"]))
        return [int(r["year"]) for r in rows], [float(numeric(r[value_name])) for r in rows]

    for filename, value_name, ylabel, title, ylim in [
        (
            "yearly_mean_intensity_by_mag",
            "mean_event_mean_best_station_intensity",
            "Mean station intensity",
            "Yearly mean station intensity by magnitude bin",
            None,
        ),
        (
            "yearly_max_intensity_by_mag",
            "mean_event_max_intensity",
            "Mean observed maximum intensity",
            "Yearly mean observed maximum intensity by magnitude bin",
            (0.5, 7.2),
        ),
    ]:
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        for idx, label in enumerate(labels):
            years, values = sorted_yearly(label, value_name)
            if years:
                ax.plot(
                    years,
                    values,
                    marker="o",
                    markersize=3.6,
                    linewidth=1.7,
                    color=PALETTE[idx % len(PALETTE)],
                    label=label,
                )
        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel)
        ax.set_title(title + ("" if not suffix else f" ({suffix.replace('_', ' ')})"))
        if ylim:
            ax.set_ylim(*ylim)
        style_journal_axis(ax)
        add_network_epoch_guides(ax)
        ax.legend(title="Magnitude bin", ncols=2, loc="upper right")
        fig.tight_layout()
        save_journal_png(fig, out_dir / f"{filename}{suffix_part}.png")
        plt.close(fig)

    heat_rows = [r for r in yearly if int(r["n_events"]) >= min_events]
    years = sorted({int(r["year"]) for r in heat_rows})
    heat = np.full((len(labels), len(years)), np.nan)
    year_index = {year: idx for idx, year in enumerate(years)}
    label_index = {label: idx for idx, label in enumerate(labels)}
    for row in heat_rows:
        value = numeric(row.get("mean_event_mean_best_station_intensity"))
        if value is None:
            continue
        heat[label_index[str(row["mag_bin"])], year_index[int(row["year"])]] = value
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    cmap = plt.get_cmap("cividis").copy()
    cmap.set_bad("#f0f0f0")
    im = ax.imshow(heat, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0.8, vmax=5.0)
    ax.set_yticks(range(len(labels)), labels=labels)
    tick_positions = [i for i, year in enumerate(years) if year % 5 == 0 or i in (0, len(years) - 1)]
    ax.set_xticks(tick_positions, [str(years[i]) for i in tick_positions], rotation=45, ha="right")
    ax.set_title("Mean station intensity heatmap" + ("" if not suffix else f" ({suffix.replace('_', ' ')})"))
    ax.set_xlabel("Year")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    colorbar = fig.colorbar(im, ax=ax, label="Mean station intensity", fraction=0.034, pad=0.018)
    colorbar.outline.set_visible(False)
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"mean_intensity_heatmap{suffix_part}.png")
    plt.close(fig)

    period_labels = [label for _, _, label in PERIODS]
    bin_centers = {label: (low + min(high, 8.0)) / 2.0 for low, high, label in MAG_BINS}
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for pidx, plabel in enumerate(period_labels):
        rows = [
            r
            for r in period
            if r["period"] == plabel
            and int(r["n_events"]) >= min_events
            and numeric(r.get("mean_event_mean_best_station_intensity")) is not None
        ]
        points = sorted((bin_centers[str(r["mag_bin"])], numeric(r["mean_event_mean_best_station_intensity"])) for r in rows)
        if points:
            ax.plot(
                [p[0] for p in points],
                [p[1] for p in points],
                marker="o",
                markersize=4.0,
                linewidth=1.8,
                color=PALETTE[pidx % len(PALETTE)],
                label=plabel,
            )
    ax.set_xlabel("Magnitude bin center")
    ax.set_ylabel("Mean station intensity")
    ax.set_title("Magnitude relation by period" + ("" if not suffix else f" ({suffix.replace('_', ' ')})"))
    style_journal_axis(ax)
    ax.legend(title="Period", loc="upper left")
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"magnitude_relation_by_period{suffix_part}.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for idx, label in enumerate(labels):
        rows = [
            r
            for r in distance_yearly
            if r["mag_bin"] == label
            and int(r["n_events_with_distance"]) >= min_events
            and numeric(r.get("median_nearest_station_distance_km")) is not None
        ]
        rows.sort(key=lambda r: int(r["year"]))
        if rows:
            ax.plot(
                [int(r["year"]) for r in rows],
                [numeric(r["median_nearest_station_distance_km"]) for r in rows],
                marker="o",
                markersize=3.6,
                linewidth=1.7,
                color=PALETTE[idx % len(PALETTE)],
                label=label,
            )
    ax.set_xlabel("Year")
    ax.set_ylabel("Median shortest epicentral distance [km]")
    ax.set_title("Yearly shortest detected epicentral distance by magnitude bin" + ("" if not suffix else f" ({suffix.replace('_', ' ')})"))
    style_journal_axis(ax)
    add_network_epoch_guides(ax)
    ax.legend(title="Magnitude bin", ncols=2, loc="upper right")
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"yearly_shortest_epicentral_distance{suffix_part}.png")
    plt.close(fig)

    spatial_rows = [
        r
        for r in yearly
        if r["mag_bin"] == "4.0<=M<5.0" and int(r["n_events"]) >= min_events
    ]
    spatial_rows.sort(key=lambda r: int(r["year"]))
    if spatial_rows:
        fig, axes = plt.subplots(3, 1, figsize=(7.4, 7.0), sharex=True)
        for idx, radius in enumerate(STATION_COUNT_RADII_KM):
            color = PALETTE[idx % len(PALETTE)]
            years = [int(r["year"]) for r in spatial_rows]
            active = [numeric(r.get(f"mean_active_station_count_within_{radius}km")) for r in spatial_rows]
            detected = [numeric(r.get(f"mean_detected_station_count_within_{radius}km")) for r in spatial_rows]
            ratios = [numeric(r.get(f"mean_detected_to_active_station_ratio_within_{radius}km")) for r in spatial_rows]
            if any(v is not None for v in active):
                axes[0].plot(years, active, marker="o", markersize=3.2, linewidth=1.5, color=color, label=f"{radius} km")
            if any(v is not None for v in detected):
                axes[1].plot(years, detected, marker="o", markersize=3.2, linewidth=1.5, color=color, label=f"{radius} km")
            if any(v is not None for v in ratios):
                axes[2].plot(years, ratios, marker="o", markersize=3.2, linewidth=1.5, color=color, label=f"{radius} km")
        axes[0].set_ylabel("Active stations")
        axes[1].set_ylabel("Detected stations")
        axes[2].set_ylabel("Detected / active")
        axes[2].set_xlabel("Year")
        axes[0].set_title(
            "Stations within epicentral radius for M4.0-5.0 events"
            + ("" if not suffix else f" ({suffix.replace('_', ' ')})")
        )
        for label, ax in zip(["a", "b", "c"], axes):
            style_journal_axis(ax)
            add_panel_label(ax, label)
            add_network_epoch_guides(ax, annotate=label == "a")
        axes[0].legend(title="Radius", ncols=4, loc="upper center", bbox_to_anchor=(0.5, 1.22))
        fig.tight_layout()
        save_journal_png(fig, out_dir / f"yearly_station_counts_within_radius_m4_5{suffix_part}.png")
        plt.close(fig)

    periods = [label for _, _, label in FIVE_YEAR_PERIODS]
    distance_labels = [label for _, _, label in DISTANCE_BINS]
    x = np.arange(len(distance_labels))
    width = min(0.82 / max(len(periods), 1), 0.11)
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    for idx, plabel in enumerate(periods):
        rows = {r["distance_bin_km"]: r for r in distance_period if r["period"] == plabel}
        values = [int(rows.get(dlabel, {}).get("n_events", 0)) for dlabel in distance_labels]
        if not any(values):
            continue
        offset = (idx - (len(periods) - 1) / 2.0) * width
        ax.bar(x + offset, values, width=width, alpha=0.82, color=PALETTE[idx % len(PALETTE)], label=plabel)
    ax.set_xticks(x, distance_labels)
    ax.set_xlabel("Shortest epicentral distance bin [km]")
    ax.set_ylabel("Number of M>=4 events")
    ax.set_title("Distribution of shortest detected epicentral distance" + ("" if not suffix else f" ({suffix.replace('_', ' ')})"))
    style_journal_axis(ax)
    ax.legend(title="Period", ncols=3, loc="upper center", bbox_to_anchor=(0.5, 1.22))
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"shortest_epicentral_distance_histogram{suffix_part}.png")
    plt.close(fig)


def render_station_network_density_pngs(out_dir: Path, yearly_rows: list[dict[str, object]]) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:
        print(f"Skipping station network PNG plots because matplotlib/numpy is unavailable: {exc}")
        return

    configure_journal_matplotlib(plt)
    regions = [region.name for region in NETWORK_REGIONS]

    def plot_metric(filename: str, metric: str, ylabel: str, title: str) -> None:
        fig, ax = plt.subplots(figsize=(7.6, 4.4))
        for idx, region_name in enumerate(regions):
            rows = [
                r
                for r in yearly_rows
                if r["region"] == region_name and numeric(r.get(metric)) is not None
            ]
            rows.sort(key=lambda r: int(r["year"]))
            if not rows:
                continue
            ax.plot(
                [int(r["year"]) for r in rows],
                [numeric(r[metric]) for r in rows],
                marker="o",
                markersize=3.2,
                linewidth=1.55,
                label=region_name,
                color=PALETTE[idx % len(PALETTE)],
            )
        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        style_journal_axis(ax)
        add_network_epoch_guides(ax)
        ax.legend(title="Region", ncols=3, loc="upper center", bbox_to_anchor=(0.5, 1.24))
        fig.tight_layout()
        save_journal_png(fig, out_dir / filename)
        plt.close(fig)

    plot_metric(
        "station_density_by_region.png",
        "density_per_10000km2",
        "Active stations per 10,000 km2",
        "Active seismic intensity station density by region",
    )
    plot_metric(
        "station_nearest_neighbor_distance_by_region.png",
        "median_nearest_neighbor_distance_km",
        "Median nearest station distance [km]",
        "Median nearest-neighbor distance among active stations",
    )

    comparison_years = [1995, 2022]
    rows_by_year_region = {
        (int(row["year"]), str(row["region"])): row
        for row in yearly_rows
        if int(row["year"]) in comparison_years
    }
    x = np.arange(len(regions))
    width = 0.36
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    for idx, year in enumerate(comparison_years):
        values = [
            numeric(rows_by_year_region.get((year, region_name), {}).get("density_per_10000km2"))
            or 0.0
            for region_name in regions
        ]
        ax.bar(x + (idx - 0.5) * width, values, width=width, color=PALETTE[idx], alpha=0.88, label=str(year))
    ax.set_xticks(x, regions, rotation=30, ha="right")
    ax.set_ylabel("Active stations per 10,000 km2")
    ax.set_title("Station density change by region: 1995 vs 2022")
    style_journal_axis(ax)
    ax.legend(title="Year", loc="upper left")
    fig.tight_layout()
    save_journal_png(fig, out_dir / "station_density_region_snapshot_1995_2022.png")
    plt.close(fig)


def render_m4plus_shortest_distance_png(
    out_dir: Path,
    suffix: str,
    yearly_rows: list[dict[str, object]],
    min_events: int,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Skipping M4+ shortest-distance PNG plot because matplotlib is unavailable: {exc}")
        return

    configure_journal_matplotlib(plt)
    rows = [r for r in yearly_rows if int(r.get("n_events", 0)) >= min_events]
    rows.sort(key=lambda r: int(r["year"]))
    if not rows:
        return
    suffix_part = f"_{suffix}" if suffix else ""
    years = [int(r["year"]) for r in rows]

    fig, axes = plt.subplots(3, 1, figsize=(7.4, 7.0), sharex=True)
    axes[0].plot(
        years,
        [numeric(r.get("median_nearest_station_distance_km")) for r in rows],
        marker="o",
        markersize=3.4,
        linewidth=1.7,
        color=PALETTE[0],
        label="Median",
    )
    axes[0].fill_between(
        years,
        [numeric(r.get("q25_nearest_station_distance_km")) for r in rows],
        [numeric(r.get("q75_nearest_station_distance_km")) for r in rows],
        color=PALETTE[0],
        alpha=0.16,
        label="25-75%",
    )
    axes[0].set_ylabel("Shortest distance [km]")
    axes[0].set_title(
        "Yearly shortest detected epicentral distance for M>=4"
        + ("" if not suffix else f" ({suffix.replace('_', ' ')})")
    )
    axes[0].legend()

    for idx, (column, label) in enumerate(
        [
            ("share_nearest_station_within_10km", "<=10 km"),
            ("share_nearest_station_within_20km", "<=20 km"),
            ("share_nearest_station_within_50km", "<=50 km"),
        ]
    ):
        axes[1].plot(
            years,
            [numeric(r.get(column)) for r in rows],
            marker="o",
            markersize=3.4,
            linewidth=1.7,
            color=PALETTE[idx + 1],
            label=label,
        )
    axes[1].set_ylabel("Share of events")
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].legend(title="Nearest station")

    axes[2].bar(years, [int(r["n_events"]) for r in rows], color=BAR_COLOR, alpha=0.9)
    axes[2].set_ylabel("M>=4 events")
    axes[2].set_xlabel("Year")

    for label, ax in zip(["a", "b", "c"], axes):
        style_journal_axis(ax)
        add_panel_label(ax, label)
        add_network_epoch_guides(ax, annotate=label == "a")
    fig.tight_layout()
    save_journal_png(fig, out_dir / f"yearly_shortest_epicentral_distance_m4plus{suffix_part}.png")
    plt.close(fig)


def svg_header(width: int, height: int, title: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{svg_escape(title)}">',
        "<style>",
        "text { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; fill: #1f2933; }",
        ".title { font-size: 20px; font-weight: 700; }",
        ".axis-label { font-size: 13px; font-weight: 600; fill: #334e68; }",
        ".tick { font-size: 11px; fill: #52606d; }",
        ".legend { font-size: 12px; fill: #334e68; }",
        ".note { font-size: 11px; fill: #627d98; }",
        ".grid { stroke: #d9e2ec; stroke-width: 1; }",
        ".xgrid { stroke-dasharray: 2 4; }",
        ".frame { fill: none; stroke: #9fb3c8; stroke-width: 1; }",
        "</style>",
    ]


def heat_color(value: float, vmin: float, vmax: float) -> str:
    value = max(vmin, min(vmax, value))
    t = (value - vmin) / (vmax - vmin)
    stops = [
        (0.00, (247, 251, 255)),
        (0.25, (198, 219, 239)),
        (0.50, (107, 174, 214)),
        (0.75, (33, 113, 181)),
        (1.00, (8, 48, 107)),
    ]
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t0 <= t <= t1:
            local = (t - t0) / (t1 - t0)
            rgb = tuple(round(c0[i] * (1 - local) + c1[i] * local) for i in range(3))
            return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
    rgb = stops[-1][1]
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def write_analysis_outputs(
    csv_dir: Path,
    png_dir: Path,
    suffix: str,
    events: list[dict[str, object]],
    min_events_for_plots: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    yearly, period = summarize_events(events)
    yearly_regression = regression_summaries(events)
    distance_yearly, distance_period = distance_detection_summaries(events)
    yearly_shortest_m4plus = yearly_shortest_distance_threshold_summary(events, min_magnitude=4.0)
    magbin_trend = magbin_annual_trend_summaries(yearly)
    slope_period = yearly_slope_period_summaries(yearly_regression)
    period_event_regression = period_event_regression_summaries(events)
    adjusted_models = adjusted_event_model_summaries(events)
    yearly_spatial_counts = select_spatial_station_count_columns(yearly)
    period_spatial_counts = select_spatial_station_count_columns(period)

    suffix_part = f"_{suffix}" if suffix else ""
    write_csv(csv_dir / f"yearly_magbin_summary{suffix_part}.csv", yearly)
    write_csv(csv_dir / f"period_magbin_summary{suffix_part}.csv", period)
    write_csv(csv_dir / f"yearly_regression_summary{suffix_part}.csv", yearly_regression)
    write_csv(csv_dir / f"yearly_distance_detection_summary{suffix_part}.csv", distance_yearly)
    write_csv(csv_dir / f"period_distance_detection_histogram{suffix_part}.csv", distance_period)
    write_csv(csv_dir / f"yearly_shortest_epicentral_distance_m4plus{suffix_part}.csv", yearly_shortest_m4plus)
    write_csv(csv_dir / f"magbin_annual_trend_summary{suffix_part}.csv", magbin_trend)
    write_csv(csv_dir / f"yearly_slope_summary_by_period{suffix_part}.csv", slope_period)
    write_csv(csv_dir / f"period_event_regression{suffix_part}.csv", period_event_regression)
    write_csv(csv_dir / f"adjusted_event_models{suffix_part}.csv", adjusted_models)
    write_csv(csv_dir / f"yearly_spatial_station_count_summary{suffix_part}.csv", yearly_spatial_counts)
    write_csv(csv_dir / f"period_spatial_station_count_summary{suffix_part}.csv", period_spatial_counts)

    title_suffix = "" if not suffix else f" ({suffix.replace('_', ' ')})"
    _ = title_suffix
    render_standard_png_charts(png_dir, suffix, yearly, period, distance_yearly, distance_period, min_events_for_plots)
    render_m4plus_shortest_distance_png(png_dir, suffix, yearly_shortest_m4plus, min_events_for_plots)
    render_reference_style_pngs(png_dir, suffix, yearly, yearly_regression, min_events_for_plots)
    return yearly, period, yearly_regression


def print_key_findings(label: str, events: list[dict[str, object]], yearly: list[dict[str, object]]) -> None:
    target_rows = [
        r
        for r in yearly
        if r["mag_bin"] == "4.0<=M<5.0"
        and numeric(r["mean_event_max_intensity"]) is not None
        and int(r["n_events"]) >= 3
    ]
    pre = [r for r in target_rows if int(r["year"]) <= 1995]
    post = [r for r in target_rows if int(r["year"]) >= 1996]
    pre_max = mean(numeric(r["mean_event_max_intensity"]) for r in pre)
    post_max = mean(numeric(r["mean_event_max_intensity"]) for r in post)
    pre_mean = mean(numeric(r["mean_event_mean_intensity"]) for r in pre)
    post_mean = mean(numeric(r["mean_event_mean_intensity"]) for r in post)
    pre_best = mean(numeric(r["mean_event_mean_best_station_intensity"]) for r in pre)
    post_best = mean(numeric(r["mean_event_mean_best_station_intensity"]) for r in post)

    all_years = sorted({r["year"] for r in events if isinstance(r["year"], int)})
    print(f"\n[{label}]")
    print(f"Events: {len(events):,}")
    print(f"Years: {all_years[0]}-{all_years[-1]} ({len(all_years)} files; missing years may reflect local data availability)")
    print("M4.0-5.0 yearly bins with n>=3:")
    print(f"  mean observed maximum intensity, <=1995: {pre_max:.2f}" if pre_max is not None else "  pre max: n/a")
    print(f"  mean observed maximum intensity, >=1996: {post_max:.2f}" if post_max is not None else "  post max: n/a")
    print(f"  mean station intensity, <=1995: {pre_best:.2f}" if pre_best is not None else "  pre station mean: n/a")
    print(f"  mean station intensity, >=1996: {post_best:.2f}" if post_best is not None else "  post station mean: n/a")
    print(f"  class-only mean intensity, <=1995: {pre_mean:.2f}" if pre_mean is not None else "  pre class mean: n/a")
    print(f"  class-only mean intensity, >=1996: {post_mean:.2f}" if post_mean is not None else "  post class mean: n/a")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--code-p", type=Path, default=None, help="Optional path to code_p.dat or code_p.zip.")
    parser.add_argument("--no-station-metadata", action="store_true", help="Skip station names and distance calculation.")
    parser.add_argument("--min-events-for-plots", type=int, default=3)
    args = parser.parse_args()

    zip_paths = sorted(args.data_dir.glob("i[0-9][0-9][0-9][0-9].zip"))
    if not zip_paths:
        raise SystemExit(f"No iYYYY.zip files found in {args.data_dir}")

    station_index = None if args.no_station_metadata else load_station_index(args.data_dir, args.code_p)
    if station_index is None and not args.no_station_metadata:
        print("code_p.dat was not found; station names and distances will be blank.")

    events: list[dict[str, object]] = []
    for path in zip_paths:
        events.extend(parse_catalog_zip(path, station_index))

    events = [
        row
        for row in events
        if numeric(row["magnitude"]) is not None
        and numeric(row["mean_intensity_class"]) is not None
        and numeric(row["max_intensity_value"]) is not None
    ]

    for row in events:
        row["mag_bin"] = mag_bin(numeric(row["magnitude"]))
        row["period"] = period_label(row["year"] if isinstance(row["year"], int) else None)
        row["is_land_region_approx"] = is_land_region_approx(row["hypocenter_region"])
        depth = numeric(row["depth_km"])
        row["is_shallow_20km"] = depth is not None and depth <= 20.0

    out_dir = args.output_dir
    csv_dir = out_dir / "csv"
    png_dir = out_dir / "png"
    csv_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)
    write_csv(csv_dir / "event_intensity_summary.csv", events)

    if station_index is not None:
        zip_years = [int(path.stem[1:]) for path in zip_paths if path.stem[1:].isdigit()]
        station_years = list(range(min(zip_years), max(zip_years) + 1))
        station_yearly, station_period, station_selected = station_network_density_summaries(
            station_index,
            station_years,
        )
        write_csv(csv_dir / "station_network_density_by_year_region.csv", station_yearly)
        write_csv(csv_dir / "station_network_density_period_summary.csv", station_period)
        write_csv(csv_dir / "station_network_density_selected_years.csv", station_selected)
        render_station_network_density_pngs(png_dir, station_yearly)

    datasets = [
        ("all", "", events),
        ("land_approx", "land_approx", [r for r in events if r["is_land_region_approx"]]),
        (
            "land_approx_depth20",
            "land_approx_depth20",
            [r for r in events if r["is_land_region_approx"] and r["is_shallow_20km"]],
        ),
    ]
    for label, suffix, rows in datasets:
        yearly, _, _ = write_analysis_outputs(csv_dir, png_dir, suffix, rows, args.min_events_for_plots)
        print_key_findings(label, rows, yearly)

    print(f"CSV outputs written to: {csv_dir.resolve()}")
    print(f"PNG outputs written to: {png_dir.resolve()}")


if __name__ == "__main__":
    main()
