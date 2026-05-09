#!/usr/bin/env python3
"""Estimate how recent maximum intensities change under the 1994 station geometry."""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from analyze_jma_intensity import active_station_records_at_year_end, haversine_km, load_station_index
from analyze_station_thinning_interpolation import intensity_class_index, intensity_class_labels
from estimate_jma_intensity_distribution import event_context, load_event_station_observations, project_km


DEFAULT_TARGET_CSV = Path("outputs/csv/hypocenter_catalog/target_events_intensity_5upper_plus_with_hypocenter.csv")
DEFAULT_CSV_DIR = Path("outputs/csv/recent_1994_network_counterfactual")
DEFAULT_PNG_DIR = Path("outputs/png/recent_1994_network_counterfactual")
DEFAULT_DERIVED_DIR = Path("data/derived/recent_1994_network_counterfactual")
DEFAULT_PAPER_ASSETS_JA = Path("paper/assets_ja")
DEFAULT_PAPER_ASSETS_EN = Path("paper/assets_en")
JAPAN_LAND_AREA_KM2 = 377_975.26
SUMMARY_RADII_KM = (20.0, 50.0, 100.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resample recent strong-intensity events to the station geometry active at the end of 1994 "
            "and compare the sampled maximum intensity with the current observed maximum."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--target-events-csv", type=Path, default=DEFAULT_TARGET_CSV)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--derived-dir", type=Path, default=DEFAULT_DERIVED_DIR)
    parser.add_argument("--paper-assets-ja", type=Path, default=DEFAULT_PAPER_ASSETS_JA)
    parser.add_argument("--paper-assets-en", type=Path, default=DEFAULT_PAPER_ASSETS_EN)
    parser.add_argument("--recent-start-year", type=int, default=2011)
    parser.add_argument("--recent-end-year", type=int, default=2022)
    parser.add_argument("--counterfactual-year", type=int, default=1994)
    parser.add_argument("--match-radius-km", type=float, default=10.0)
    parser.add_argument("--min-station-intensity", type=float, default=1.0)
    parser.add_argument("--max-events", type=int, default=None)
    return parser.parse_args()


def configure_matplotlib() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 330,
            "font.family": "DejaVu Sans",
            "font.size": 7.4,
            "axes.labelsize": 7.8,
            "axes.titlesize": 8.2,
            "axes.titleweight": "bold",
            "axes.linewidth": 0.8,
            "axes.edgecolor": "#222222",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 6.8,
            "legend.frameon": False,
            "xtick.color": "#222222",
            "ytick.color": "#222222",
            "axes.labelcolor": "#222222",
            "grid.color": "#d8d8d8",
            "grid.linewidth": 0.45,
        }
    )


def class_index(value: float) -> int:
    return int(intensity_class_index(np.array([value], dtype=float))[0])


def class_label(value: float) -> str:
    idx = np.array([class_index(value)], dtype=int)
    return str(intensity_class_labels(idx)[0])


def class_label_from_index(index: int) -> str:
    return str(intensity_class_labels(np.array([index], dtype=int))[0])


def load_targets(args: argparse.Namespace) -> pd.DataFrame:
    if not args.target_events_csv.exists():
        raise SystemExit(f"Target catalog not found: {args.target_events_csv}")
    df = pd.read_csv(args.target_events_csv, low_memory=False)
    required = {"event_id", "year", "month", "day", "hour", "minute"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Target catalog is missing columns: {sorted(missing)}")
    out = df[(df["year"] >= args.recent_start_year) & (df["year"] <= args.recent_end_year)].copy()
    out = out.sort_values(["year", "month", "day", "hour", "minute", "event_id"]).reset_index(drop=True)
    if args.max_events is not None:
        out = out.head(args.max_events).copy()
    if out.empty:
        raise SystemExit("No target events after filtering.")
    return out


def observed_station_frame(event, data_dir: Path, min_intensity: float) -> pd.DataFrame:
    stations = load_event_station_observations(event, data_dir)
    stations = stations.dropna(subset=["latitude", "longitude", "observed_intensity"]).copy()
    stations["observed_intensity"] = stations["observed_intensity"].astype(float)
    stations = stations[stations["observed_intensity"] >= min_intensity].copy()
    return stations.reset_index(drop=True)


def finite_min(values: pd.Series | np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.min(arr)) if len(arr) else float("nan")


def finite_median(values: pd.Series | np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if len(arr) else float("nan")


def count_within_distance(stations: pd.DataFrame, radius_km: float) -> int:
    if stations.empty or "epicentral_distance_km" not in stations:
        return 0
    distances = stations["epicentral_distance_km"].astype(float).to_numpy()
    return int(np.sum(np.isfinite(distances) & (distances <= radius_km)))


def build_counterfactual_matches(
    event,
    current: pd.DataFrame,
    active_records,
    match_radius_km: float,
) -> pd.DataFrame:
    if current.empty:
        return pd.DataFrame()
    obs_xy = project_km(
        current["longitude"].to_numpy(dtype=float),
        current["latitude"].to_numpy(dtype=float),
        event.longitude,
        event.latitude,
    )
    tree = cKDTree(obs_xy)
    obs = current.reset_index(drop=True)
    rows: list[dict[str, object]] = []
    for rec in active_records:
        if rec.latitude is None or rec.longitude is None:
            continue
        query_xy = project_km(np.array([rec.longitude]), np.array([rec.latitude]), event.longitude, event.latitude)
        nearest_distance, nearest_idx = tree.query(query_xy, k=1)
        nearest_distance = float(np.ravel(nearest_distance)[0])
        nearest_idx = int(np.ravel(nearest_idx)[0])
        matched = nearest_distance <= match_radius_km
        src = obs.iloc[nearest_idx]
        epicentral_distance = haversine_km(event.latitude, event.longitude, rec.latitude, rec.longitude)
        rows.append(
            {
                "event_id": event.event_id,
                "counterfactual_station_code": rec.code,
                "counterfactual_station_name": rec.name,
                "counterfactual_latitude": rec.latitude,
                "counterfactual_longitude": rec.longitude,
                "counterfactual_epicentral_distance_km": epicentral_distance,
                "nearest_current_station_code": src["station_code"],
                "nearest_current_station_name": src["station_name"],
                "nearest_current_latitude": float(src["latitude"]),
                "nearest_current_longitude": float(src["longitude"]),
                "nearest_current_epicentral_distance_km": float(src["epicentral_distance_km"]),
                "match_distance_km": nearest_distance,
                "matched": bool(matched),
                "pseudo_observed_intensity": float(src["observed_intensity"]) if matched else np.nan,
            }
        )
    return pd.DataFrame(rows)


def summarize_event(row: pd.Series, args: argparse.Namespace, active_records) -> tuple[dict[str, object], pd.DataFrame]:
    event = event_context(row)
    current = observed_station_frame(event, args.data_dir, args.min_station_intensity)
    if current.empty:
        return (
            {
                "event_id": event.event_id,
                "origin_time": event.origin_time,
                "analysis_status": "no_current_station",
            },
            pd.DataFrame(),
        )

    matches = build_counterfactual_matches(event, current, active_records, args.match_radius_km)
    matched = matches[matches["matched"]].copy()
    current_max = float(current["observed_intensity"].max())
    current_max_station = current.loc[current["observed_intensity"].idxmax()]
    cf_max = float(matched["pseudo_observed_intensity"].max()) if not matched.empty else float("nan")
    if matched.empty:
        cf_max_station = None
    else:
        cf_max_station = matched.loc[matched["pseudo_observed_intensity"].idxmax()]

    current_class = class_index(current_max)
    cf_class = class_index(cf_max) if np.isfinite(cf_max) else -1
    row_out: dict[str, object] = {
        "event_id": event.event_id,
        "origin_time": event.origin_time,
        "year": event.year,
        "region_name": row.get("hyp_region_name", row.get("hypocenter_region", "")),
        "japanese_region_name": row.get("hypocenter_region", ""),
        "latitude": event.latitude,
        "longitude": event.longitude,
        "depth_km": event.depth_km,
        "magnitude": event.magnitude,
        "current_station_count": int(len(current)),
        "current_max_intensity": current_max,
        "current_max_class_index": current_class,
        "current_max_class": class_label_from_index(current_class),
        "current_max_station_code": current_max_station["station_code"],
        "current_max_station_name": current_max_station["station_name"],
        "current_max_station_epicentral_distance_km": float(current_max_station["epicentral_distance_km"]),
        "counterfactual_year_end": args.counterfactual_year,
        "counterfactual_active_station_count": int(len(active_records)),
        "counterfactual_matched_station_count": int(len(matched)),
        "counterfactual_matched_fraction": float(len(matched) / len(active_records)) if active_records else float("nan"),
        "counterfactual_max_intensity": cf_max,
        "counterfactual_max_class_index": cf_class,
        "counterfactual_max_class": class_label_from_index(cf_class) if cf_class >= 0 else "",
        "max_intensity_drop": current_max - cf_max if np.isfinite(cf_max) else float("nan"),
        "max_class_drop": current_class - cf_class if cf_class >= 0 else float("nan"),
        "match_radius_km": args.match_radius_km,
        "current_nearest_epicentral_distance_km": finite_min(current["epicentral_distance_km"]),
        "current_median_epicentral_distance_km": finite_median(current["epicentral_distance_km"]),
        "counterfactual_nearest_epicentral_distance_km": finite_min(matched["counterfactual_epicentral_distance_km"])
        if not matched.empty
        else float("nan"),
        "counterfactual_median_epicentral_distance_km": finite_median(matched["counterfactual_epicentral_distance_km"])
        if not matched.empty
        else float("nan"),
        "counterfactual_match_distance_median_km": finite_median(matched["match_distance_km"])
        if not matched.empty
        else float("nan"),
        "counterfactual_match_distance_max_km": float(matched["match_distance_km"].max()) if not matched.empty else float("nan"),
        "analysis_status": "ok" if not matched.empty else "no_counterfactual_match",
    }
    if cf_max_station is not None:
        row_out.update(
            {
                "counterfactual_max_station_code": cf_max_station["counterfactual_station_code"],
                "counterfactual_max_station_name": cf_max_station["counterfactual_station_name"],
                "counterfactual_max_station_epicentral_distance_km": float(
                    cf_max_station["counterfactual_epicentral_distance_km"]
                ),
                "counterfactual_max_source_station_code": cf_max_station["nearest_current_station_code"],
                "counterfactual_max_source_station_name": cf_max_station["nearest_current_station_name"],
                "counterfactual_max_source_match_distance_km": float(cf_max_station["match_distance_km"]),
            }
        )

    for radius in SUMMARY_RADII_KM:
        key = int(radius)
        row_out[f"current_station_count_within_{key}km"] = count_within_distance(current, radius)
        row_out[f"counterfactual_matched_count_within_{key}km"] = count_within_distance(
            matched.rename(columns={"counterfactual_epicentral_distance_km": "epicentral_distance_km"}),
            radius,
        )
    row_out["current_density_national_per_10000km2"] = len(current) / JAPAN_LAND_AREA_KM2 * 10_000.0
    row_out["counterfactual_active_density_national_per_10000km2"] = len(active_records) / JAPAN_LAND_AREA_KM2 * 10_000.0
    row_out["counterfactual_matched_density_national_per_10000km2"] = len(matched) / JAPAN_LAND_AREA_KM2 * 10_000.0

    matches["origin_time"] = event.origin_time
    matches["year"] = event.year
    return row_out, matches


def summarize_groups(event_df: pd.DataFrame) -> pd.DataFrame:
    valid = event_df[event_df["analysis_status"].eq("ok") & event_df["counterfactual_max_intensity"].notna()].copy()
    rows: list[dict[str, object]] = []
    groups: list[tuple[str, pd.DataFrame]] = [("all_recent_I5upper_plus", valid)]
    groups.append(
        (
            "support_current_ge100_and_cf_matched_ge20",
            valid[(valid["current_station_count"] >= 100) & (valid["counterfactual_matched_station_count"] >= 20)],
        )
    )
    groups.append(("support_cf_matched_ge50", valid[valid["counterfactual_matched_station_count"] >= 50]))
    groups.append(("current_max_6upper_plus", valid[valid["current_max_class_index"].astype(int) >= 8]))
    for idx in sorted(valid["current_max_class_index"].dropna().astype(int).unique()):
        groups.append((f"current_max_{class_label_from_index(idx)}", valid[valid["current_max_class_index"].astype(int) == idx]))
    for label, group in groups:
        if group.empty:
            continue
        drop = group["max_intensity_drop"].astype(float)
        class_drop = group["max_class_drop"].astype(float)
        rows.append(
            {
                "group": label,
                "n_events": int(len(group)),
                "current_max_mean": float(group["current_max_intensity"].mean()),
                "current_max_median": float(group["current_max_intensity"].median()),
                "counterfactual_max_mean": float(group["counterfactual_max_intensity"].mean()),
                "counterfactual_max_median": float(group["counterfactual_max_intensity"].median()),
                "max_intensity_drop_mean": float(drop.mean()),
                "max_intensity_drop_median": float(drop.median()),
                "max_intensity_drop_q75": float(drop.quantile(0.75)),
                "share_drop_gt_0": float((drop > 0).mean()),
                "share_drop_ge_0p5": float((drop >= 0.5).mean()),
                "share_drop_ge_1p0": float((drop >= 1.0).mean()),
                "share_class_drop_ge_1": float((class_drop >= 1).mean()),
                "current_station_count_median": float(group["current_station_count"].median()),
                "counterfactual_matched_station_count_median": float(
                    group["counterfactual_matched_station_count"].median()
                ),
                "current_nearest_epicentral_distance_median_km": float(
                    group["current_nearest_epicentral_distance_km"].median()
                ),
                "counterfactual_nearest_epicentral_distance_median_km": float(
                    group["counterfactual_nearest_epicentral_distance_km"].median()
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_summary(event_df: pd.DataFrame, png_dir: Path) -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    configure_matplotlib()
    valid = event_df[event_df["analysis_status"].eq("ok") & event_df["counterfactual_max_intensity"].notna()].copy()
    valid = valid.sort_values("max_intensity_drop", ascending=False)
    if valid.empty:
        raise SystemExit("No valid events for plotting.")

    fig = plt.figure(figsize=(10.8, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.0, 1.05], hspace=0.32, wspace=0.30)
    ax_scatter = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1])
    ax_counts = fig.add_subplot(gs[1, 0])
    ax_rank = fig.add_subplot(gs[1, 1])

    drop = valid["max_intensity_drop"].astype(float).to_numpy()
    current = valid["current_max_intensity"].astype(float).to_numpy()
    cf = valid["counterfactual_max_intensity"].astype(float).to_numpy()
    norm = Normalize(vmin=0.0, vmax=max(1.5, float(np.nanquantile(drop, 0.95))))
    sc = ax_scatter.scatter(
        current,
        cf,
        c=drop,
        cmap="magma_r",
        norm=norm,
        s=np.clip(valid["counterfactual_matched_station_count"].to_numpy(dtype=float) / 4.0, 18, 92),
        edgecolor="#242424",
        linewidth=0.35,
        alpha=0.88,
    )
    lims = [4.85, 7.05]
    ax_scatter.plot(lims, lims, color="#222222", linewidth=0.9, linestyle="--", label="1:1")
    for bound in [5.0, 5.5, 6.0, 6.5]:
        ax_scatter.axvline(bound, color="#d0d0d0", linewidth=0.55, zorder=0)
        ax_scatter.axhline(bound, color="#d0d0d0", linewidth=0.55, zorder=0)
    ax_scatter.set_xlim(lims)
    ax_scatter.set_ylim(lims)
    ax_scatter.set_xlabel("Current observed maximum intensity")
    ax_scatter.set_ylabel("1994-geometry sampled maximum")
    ax_scatter.set_title("a. Maximum intensity sampled by network geometry")
    ax_scatter.legend(loc="lower right")
    cbar = fig.colorbar(sc, ax=ax_scatter, fraction=0.046, pad=0.02)
    cbar.set_label("Current - 1994 geometry")

    bins = np.arange(-0.05, max(2.05, math.ceil(float(np.nanmax(drop)) * 2) / 2 + 0.55), 0.25)
    ax_hist.hist(drop, bins=bins, color="#4c6f8f", edgecolor="white", linewidth=0.7)
    ax_hist.axvline(float(np.nanmedian(drop)), color="#9c2f2f", linewidth=1.2, label="median")
    ax_hist.axvline(0.5, color="#222222", linewidth=0.9, linestyle=":")
    ax_hist.axvline(1.0, color="#222222", linewidth=0.9, linestyle=":")
    ax_hist.set_xlabel("Maximum-intensity reduction")
    ax_hist.set_ylabel("Number of events")
    ax_hist.set_title("b. Distribution of maximum-intensity reduction")
    ax_hist.legend(loc="upper right")

    ordered = valid.sort_values("current_station_count")
    x = np.arange(len(ordered))
    ax_counts.plot(x, ordered["current_station_count"], color="#1b567d", linewidth=1.5, label="Current reporting")
    ax_counts.plot(
        x,
        ordered["counterfactual_matched_station_count"],
        color="#a33a2b",
        linewidth=1.5,
        label="1994 geometry matched",
    )
    ax_counts.set_yscale("log")
    ax_counts.set_xlabel("Events sorted by current reporting-station count")
    ax_counts.set_ylabel("Station count")
    ax_counts.set_title("c. Station support for each event")
    ax_counts.legend(loc="upper left")
    ax_counts.grid(True, axis="y")

    top = valid.head(14).iloc[::-1].copy()
    labels = [
        f"{str(r.event_id).replace('i', '')[:4]} M{float(r.magnitude):.1f} {str(r.region_name)[:18]}"
        for r in top.itertuples()
    ]
    y = np.arange(len(top))
    ax_rank.barh(y, top["max_intensity_drop"], color="#6f5b8c", edgecolor="#2a2435", linewidth=0.35)
    ax_rank.set_yticks(y)
    ax_rank.set_yticklabels(labels)
    ax_rank.set_xlabel("Current - 1994 geometry")
    ax_rank.set_title("d. Events with the largest maximum-intensity reduction")
    ax_rank.grid(True, axis="x")

    fig.suptitle(
        "1994 year-end station geometry counterfactual",
        x=0.03,
        y=0.995,
        ha="left",
        fontsize=9.0,
        fontweight="bold",
    )
    fig.text(
        0.03,
        0.965,
        "Recent target events: 2011-2022, maximum JMA intensity 5 upper or larger; "
        "pseudo-observations use the nearest current report within 10 km.",
        ha="left",
        va="top",
        fontsize=7.4,
        color="#333333",
    )
    png_dir.mkdir(parents=True, exist_ok=True)
    out = png_dir / "recent_1994_counterfactual_max_intensity.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)
    args.derived_dir.mkdir(parents=True, exist_ok=True)

    station_index = load_station_index(args.data_dir)
    if station_index is None:
        raise SystemExit(f"Station index not found under {args.data_dir}")
    active_records = active_station_records_at_year_end(station_index, args.counterfactual_year)
    if not active_records:
        raise SystemExit(f"No active stations at the end of {args.counterfactual_year}.")

    targets = load_targets(args)
    rows: list[dict[str, object]] = []
    match_frames: list[pd.DataFrame] = []
    for _, row in targets.iterrows():
        summary, matches = summarize_event(row, args, active_records)
        rows.append(summary)
        if not matches.empty:
            match_frames.append(matches)

    event_df = pd.DataFrame(rows)
    group_df = summarize_groups(event_df)
    matches_df = pd.concat(match_frames, ignore_index=True) if match_frames else pd.DataFrame()
    top_df = (
        event_df[event_df["analysis_status"].eq("ok")]
        .sort_values("max_intensity_drop", ascending=False)
        .head(20)
        .reset_index(drop=True)
    )

    event_csv = args.csv_dir / "recent_1994_counterfactual_event_summary.csv"
    group_csv = args.csv_dir / "recent_1994_counterfactual_group_summary.csv"
    top_csv = args.csv_dir / "recent_1994_counterfactual_largest_drops.csv"
    match_csv = args.csv_dir / "recent_1994_counterfactual_station_matches.csv.gz"
    event_df.to_csv(event_csv, index=False)
    group_df.to_csv(group_csv, index=False)
    top_df.to_csv(top_csv, index=False)
    matches_df.to_csv(match_csv, index=False)

    derived_event_csv = args.derived_dir / event_csv.name
    derived_group_csv = args.derived_dir / group_csv.name
    derived_top_csv = args.derived_dir / top_csv.name
    shutil.copy2(event_csv, derived_event_csv)
    shutil.copy2(group_csv, derived_group_csv)
    shutil.copy2(top_csv, derived_top_csv)

    fig_path = plot_summary(event_df, args.png_dir)
    for assets_dir in [args.paper_assets_ja, args.paper_assets_en]:
        assets_dir.mkdir(parents=True, exist_ok=True)
        prefix = "fig13"
        shutil.copy2(fig_path, assets_dir / f"{prefix}_recent_1994_counterfactual_max_intensity.png")

    valid = event_df[event_df["analysis_status"].eq("ok")]
    print(f"Targets: {len(targets)}")
    print(f"Valid counterfactual events: {len(valid)}")
    print(f"1994 year-end active stations: {len(active_records)}")
    print(group_df.to_string(index=False))
    print(f"Wrote: {event_csv}")
    print(f"Wrote: {group_csv}")
    print(f"Wrote: {top_csv}")
    print(f"Wrote: {fig_path}")


if __name__ == "__main__":
    main()
