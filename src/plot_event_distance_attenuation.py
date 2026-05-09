#!/usr/bin/env python3
"""Plot measured JMA intensity against source distance for selected events."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from estimate_jma_intensity_distribution import event_context, load_event_station_observations


DEFAULT_EVENTS_CSV = Path("outputs/csv/hypocenter_catalog/jma_intensity_events_with_hypocenter.csv")
DEFAULT_GROUND_GRID = Path("outputs/csv/jshis_surface_ground/jshis_surface_ground_grid_0p02deg.csv")
DEFAULT_EVENT_IDS = "i2018_000842,i2018_000435,i2016_000828,i2018_001267"
DEFAULT_ASSET_JA = Path("paper/assets_ja/fig05c_event_distance_attenuation.png")
DEFAULT_CSV_DIR = Path("outputs/csv/event_distance_attenuation")
DEFAULT_PNG_DIR = Path("outputs/png/event_distance_attenuation")


@dataclass(frozen=True)
class TrendFit:
    intercept: float
    slope: float
    sigma: float

    def predict(self, distance_km: np.ndarray) -> np.ndarray:
        return self.intercept + self.slope * np.log10(np.maximum(distance_km, 1.0))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--events-csv", type=Path, default=DEFAULT_EVENTS_CSV)
    parser.add_argument("--ground-grid", type=Path, default=DEFAULT_GROUND_GRID)
    parser.add_argument("--event-ids", default=DEFAULT_EVENT_IDS)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--asset-ja", type=Path, default=DEFAULT_ASSET_JA)
    parser.add_argument("--amp-intensity-coef", type=float, default=1.72)
    parser.add_argument("--max-ground-match-km", type=float, default=5.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)
    args.asset_ja.parent.mkdir(parents=True, exist_ok=True)

    events = pd.read_csv(args.events_csv, low_memory=False)
    ground = read_ground_grid(args.ground_grid) if args.ground_grid.exists() else None
    event_ids = [token.strip() for token in args.event_ids.split(",") if token.strip()]

    frames: list[pd.DataFrame] = []
    contexts = []
    summaries: list[dict[str, object]] = []
    for event_id in event_ids:
        rows = events[events["event_id"].astype(str) == event_id]
        if rows.empty:
            raise SystemExit(f"Event id not found: {event_id}")
        ctx = event_context(rows.iloc[0])
        obs = load_event_station_observations(ctx, args.data_dir)
        if obs.empty:
            raise SystemExit(f"No observations found for {event_id}")
        obs["source_distance_km"] = np.hypot(obs["epicentral_distance_km"].to_numpy(float), ctx.depth_km)
        if ground is not None:
            obs = attach_site_term(obs, ground, args.amp_intensity_coef, args.max_ground_match_km)
        fit = fit_log_distance(obs)
        obs["log_distance_fit"] = fit.predict(obs["source_distance_km"].to_numpy(float))
        obs["log_distance_residual"] = obs["measured_intensity"] - obs["log_distance_fit"]
        obs.to_csv(args.csv_dir / f"{event_id}_distance_attenuation.csv", index=False)
        frames.append(obs)
        contexts.append(ctx)
        summaries.append(summarize_event(ctx, obs, fit))

    pd.DataFrame(summaries).to_csv(args.csv_dir / "event_distance_attenuation_summary.csv", index=False)
    draw_panel_figure(frames, contexts, args.asset_ja)
    draw_panel_figure(frames, contexts, args.png_dir / "event_distance_attenuation_selected_events.png")
    print(f"Wrote {args.asset_ja}")


def read_ground_grid(path: Path) -> pd.DataFrame:
    cols = ["longitude", "latitude", "avs30_m_s", "amplification_vs400"]
    df = pd.read_csv(path, usecols=lambda c: c in cols)
    return df.dropna(subset=["longitude", "latitude", "amplification_vs400"]).reset_index(drop=True)


def attach_site_term(
    obs: pd.DataFrame,
    ground: pd.DataFrame,
    amp_intensity_coef: float,
    max_ground_match_km: float,
) -> pd.DataFrame:
    from scipy.spatial import cKDTree

    out = obs.copy()
    lat0 = float(out["latitude"].median())
    grid_xy = local_xy_km(ground["latitude"], ground["longitude"], lat0)
    obs_xy = local_xy_km(out["latitude"], out["longitude"], lat0)
    distances, indices = cKDTree(grid_xy).query(obs_xy, k=1)
    matched = distances <= max_ground_match_km

    amp = np.full(len(out), np.nan)
    avs30 = np.full(len(out), np.nan)
    amp[matched] = ground["amplification_vs400"].to_numpy(float)[indices[matched]]
    if "avs30_m_s" in ground.columns:
        avs30[matched] = ground["avs30_m_s"].to_numpy(float)[indices[matched]]
    out["ground_nearest_distance_km"] = distances
    out["amplification_vs400"] = amp
    out["avs30_m_s"] = avs30
    out["site_intensity_delta"] = amp_intensity_coef * np.log10(amp)
    return out


def local_xy_km(lat: Iterable[float], lon: Iterable[float], lat0: float) -> np.ndarray:
    lat_arr = np.asarray(lat, dtype=float)
    lon_arr = np.asarray(lon, dtype=float)
    x = lon_arr * 111.32 * np.cos(np.radians(lat0))
    y = lat_arr * 110.57
    return np.column_stack([x, y])


def fit_log_distance(obs: pd.DataFrame) -> TrendFit:
    data = obs.dropna(subset=["measured_intensity", "source_distance_km"])
    x = np.log10(np.maximum(data["source_distance_km"].to_numpy(float), 1.0))
    y = data["measured_intensity"].to_numpy(float)
    if len(y) < 3:
        return TrendFit(float(np.nanmean(y)), 0.0, float("nan"))
    a = np.column_stack([np.ones(len(x)), x])
    coef, *_ = np.linalg.lstsq(a, y, rcond=None)
    residual = y - a @ coef
    return TrendFit(float(coef[0]), float(coef[1]), float(np.sqrt(np.mean(residual**2))))


def summarize_event(ctx, obs: pd.DataFrame, fit: TrendFit) -> dict[str, object]:
    return {
        "event_id": ctx.event_id,
        "origin_time": ctx.origin_time,
        "region_name": event_label(ctx),
        "magnitude": ctx.magnitude,
        "depth_km": ctx.depth_km,
        "n_observations": int(len(obs)),
        "max_measured_intensity": float(obs["measured_intensity"].max()),
        "median_source_distance_km": float(obs["source_distance_km"].median()),
        "min_source_distance_km": float(obs["source_distance_km"].min()),
        "log_distance_intercept": fit.intercept,
        "log_distance_slope": fit.slope,
        "log_distance_rmse": fit.sigma,
    }


def draw_panel_figure(frames: list[pd.DataFrame], contexts: list, output_path: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    configure_plot_style(font_manager, plt)
    fig, axes = plt.subplots(2, 2, figsize=(9.8, 7.4), constrained_layout=True)
    last_scatter = None
    for ax, obs, ctx in zip(axes.ravel(), frames, contexts):
        fit = fit_log_distance(obs)
        last_scatter = draw_event_panel(ax, obs, ctx, fit)
    if last_scatter is not None and "site_intensity_delta" in frames[0].columns:
        cbar = fig.colorbar(last_scatter, ax=axes.ravel().tolist(), shrink=0.86, pad=0.015)
        cbar.set_label(r"地盤項 $1.72\log_{10}(A)$")
    fig.savefig(output_path, dpi=360)
    plt.close(fig)


def draw_event_panel(ax, obs: pd.DataFrame, ctx, fit: TrendFit):
    distance = obs["source_distance_km"].to_numpy(float)
    measured = obs["measured_intensity"].to_numpy(float)
    color = obs.get("site_intensity_delta")
    if color is None or color.isna().all():
        scatter = ax.scatter(distance, measured, s=15, color="#3b6ea8", alpha=0.62, linewidths=0)
    else:
        scatter = ax.scatter(
            distance,
            measured,
            c=color.to_numpy(float),
            cmap="viridis",
            vmin=-0.2,
            vmax=0.8,
            s=16,
            alpha=0.72,
            edgecolors="#202020",
            linewidths=0.18,
        )
    draw_distance_medians(ax, distance, measured)
    xline = np.geomspace(max(1.0, np.nanmin(distance) * 0.8), min(1200.0, np.nanmax(distance) * 1.12), 240)
    yline = fit.predict(xline)
    ax.plot(xline, yline, color="#0b3c8c", lw=2.0, label="log距離回帰")
    if np.isfinite(fit.sigma):
        ax.plot(xline, yline + fit.sigma, color="#0b3c8c", lw=0.9, ls="--", alpha=0.55)
        ax.plot(xline, yline - fit.sigma, color="#0b3c8c", lw=0.9, ls="--", alpha=0.55)
    ax.set_xscale("log")
    ax.set_ylim(0.0, 7.15)
    ax.set_xlabel("震源距離 (km)")
    ax.set_ylabel("計測震度")
    ax.set_title(panel_title(ctx, obs), loc="left", fontsize=9.4)
    ax.grid(True, which="major", color="#d7dce2", lw=0.55)
    ax.grid(True, which="minor", axis="x", color="#eceff3", lw=0.35)
    for spine in ax.spines.values():
        spine.set_color("#222222")
        spine.set_linewidth(0.8)
    return scatter


def draw_distance_medians(ax, distance: np.ndarray, measured: np.ndarray) -> None:
    valid = np.isfinite(distance) & np.isfinite(measured) & (distance > 0)
    if valid.sum() < 10:
        return
    bins = np.geomspace(max(1.0, np.nanmin(distance[valid])), np.nanmax(distance[valid]), 10)
    centers: list[float] = []
    medians: list[float] = []
    for left, right in zip(bins[:-1], bins[1:]):
        mask = valid & (distance >= left) & (distance < right)
        if mask.sum() < 3:
            continue
        centers.append(float(np.sqrt(left * right)))
        medians.append(float(np.median(measured[mask])))
    if centers:
        ax.plot(centers, medians, color="#c43c39", marker="o", ms=3.4, lw=1.25, label="距離帯中央値")


def panel_title(ctx, obs: pd.DataFrame) -> str:
    origin = str(ctx.origin_time)[:16]
    mag = "" if ctx.magnitude is None else f"M{ctx.magnitude:.1f}, "
    return f"{ctx.event_id}  {origin}\n{event_label(ctx)} ({mag}h={ctx.depth_km:.1f} km, n={len(obs)})"


def event_label(ctx) -> str:
    labels = {
        "i2018_000842": "大阪府北部",
        "i2018_000435": "島根県西部",
        "i2016_000828": "熊本県熊本地方",
        "i2018_001267": "胆振地方中東部",
    }
    return labels.get(ctx.event_id, getattr(ctx, "origin_time", ""))


def configure_plot_style(font_manager, plt) -> None:
    candidates = [
        "Hiragino Sans",
        "Hiragino Kaku Gothic ProN",
        "Yu Gothic",
        "Noto Sans CJK JP",
        "IPAexGothic",
        "DejaVu Sans",
    ]
    available = {font.name for font in font_manager.fontManager.ttflist}
    for candidate in candidates:
        if candidate in available:
            plt.rcParams["font.family"] = candidate
            break
    plt.rcParams.update(
        {
            "axes.unicode_minus": False,
            "axes.linewidth": 0.8,
            "font.size": 8.7,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


if __name__ == "__main__":
    main()
