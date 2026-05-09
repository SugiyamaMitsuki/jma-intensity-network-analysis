#!/usr/bin/env python3
"""Resolve final methodological concerns for the 6-lower-plus analysis.

This script is intentionally conservative.  It does not try to rescue strong
claims; it checks which claims survive after excluding pre-1996 unsplit
intensity classes, after restricting validation to high-intensity held-out
stations, and after adding uncertainty-aware summaries.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_TARGET_EVENTS = Path("outputs/csv/hypocenter_catalog/target_events_intensity_6lower_plus_with_hypocenter.csv")
DEFAULT_FOOTPRINTS = Path(
    "outputs/csv/estimated_intensity_distribution_6lower_plus/synthesis/intensity_footprint_area_summary.csv"
)
DEFAULT_PREDICTIONS = Path(
    "outputs/csv/station_thinning_interpolation_6lower_plus_class/thinning_prediction_errors.csv.gz"
)
DEFAULT_PROB_SUMMARY = Path(
    "outputs/csv/probabilistic_uncertainty_6lower_plus/probabilistic_validation_summary.csv"
)
DEFAULT_UNCERTAINTY_COMPONENTS = Path(
    "outputs/csv/probabilistic_uncertainty_6lower_plus/uncertainty_components_by_method_density.csv"
)
DEFAULT_CSV_DIR = Path("outputs/csv/final_concern_resolution")
DEFAULT_PNG_DIR = Path("outputs/png/final_concern_resolution")

METHOD_ORDER = ["gmpe_raw", "gmpe_calibrated", "idw", "kriging", "gmpe_kriging"]
CURRENT_SCALE_START = pd.Timestamp("1996-10-01")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create final concern-resolution tables and figures.")
    parser.add_argument("--target-events", type=Path, default=DEFAULT_TARGET_EVENTS)
    parser.add_argument("--footprints", type=Path, default=DEFAULT_FOOTPRINTS)
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--probabilistic-summary", type=Path, default=DEFAULT_PROB_SUMMARY)
    parser.add_argument("--uncertainty-components", type=Path, default=DEFAULT_UNCERTAINTY_COMPONENTS)
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--png-dir", type=Path, default=DEFAULT_PNG_DIR)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260509)
    parser.add_argument("--no-png", action="store_true")
    return parser.parse_args()


def density_bins(series: pd.Series) -> pd.Categorical:
    bins = [0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000, np.inf]
    return pd.cut(series, bins=bins, right=False)


def add_event_flags(events: pd.DataFrame) -> pd.DataFrame:
    out = events.copy()
    out["origin_timestamp"] = pd.to_datetime(out["intensity_origin_time"], errors="coerce")
    out["current_jma_intensity_scale"] = out["origin_timestamp"] >= CURRENT_SCALE_START
    out["pre_1996_unsplit_intensity_scale"] = ~out["current_jma_intensity_scale"]
    out["max_intensity_class_resolved"] = np.select(
        [
            out["max_intensity_value"] >= 7.0,
            out["max_intensity_value"] >= 6.5,
            out["max_intensity_value"] >= 6.0,
        ],
        ["7", "6_upper", "6_lower_or_old_6"],
        default="below_6_lower",
    )
    return out


def event_selection_sensitivity(events: pd.DataFrame) -> pd.DataFrame:
    rows = []
    scopes = {
        "all_catalog_6lower_plus": events,
        "current_scale_only": events[events["current_jma_intensity_scale"]],
        "pre_1996_unsplit_only": events[events["pre_1996_unsplit_intensity_scale"]],
    }
    for scope, sub in scopes.items():
        rows.append(
            {
                "scope": scope,
                "n_events": len(sub),
                "n_current_scale": int(sub["current_jma_intensity_scale"].sum()),
                "n_pre_1996_unsplit": int(sub["pre_1996_unsplit_intensity_scale"].sum()),
                "n_max_6lower_or_old_6": int((sub["max_intensity_value"] == 6.0).sum()),
                "n_max_6upper": int((sub["max_intensity_value"] == 6.5).sum()),
                "n_max_7": int((sub["max_intensity_value"] >= 7.0).sum()),
                "year_min": int(sub["year"].min()) if len(sub) else np.nan,
                "year_max": int(sub["year"].max()) if len(sub) else np.nan,
                "median_n_unique_stations": float(sub["n_unique_stations"].median()) if len(sub) else np.nan,
                "median_nearest_station_distance_km": float(sub["nearest_station_distance_km"].median())
                if len(sub)
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def footprint_sensitivity(footprints: pd.DataFrame, events: pd.DataFrame) -> pd.DataFrame:
    flags = events[["event_id", "current_jma_intensity_scale"]]
    fp = footprints.merge(flags, on="event_id", how="left")
    rows = []
    for scope, sub in [
        ("all_catalog_6lower_plus", fp),
        ("current_scale_only", fp[fp["current_jma_intensity_scale"]]),
    ]:
        for method, group in sub.groupby("method", observed=True):
            row = {
                "scope": scope,
                "method": method,
                "n_events": int(group["event_id"].nunique()),
                "median_max_estimated_surface_intensity": float(group["max_estimated_surface_intensity"].median()),
            }
            for threshold in [5.5, 6.0, 6.5, 7.0]:
                key = str(threshold).replace(".", "p")
                values = group[f"area_ge_{key}_km2"].to_numpy(dtype=float)
                row[f"area_ge_{key}_positive_event_rate"] = float(np.mean(values > 0.0))
                row[f"area_ge_{key}_median_km2"] = float(np.median(values))
                row[f"area_ge_{key}_p90_km2"] = float(np.quantile(values, 0.90))
            rows.append(row)
    out = pd.DataFrame(rows)
    out["method"] = pd.Categorical(out["method"], categories=METHOD_ORDER, ordered=True)
    return out.sort_values(["scope", "method"]).reset_index(drop=True)


def bootstrap_ci(values: np.ndarray, rng: np.random.Generator, n_boot: int) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan"), float("nan")
    if len(arr) == 1:
        return float(arr[0]), float(arr[0])
    idx = rng.integers(0, len(arr), size=(n_boot, len(arr)))
    med = np.median(arr[idx], axis=1)
    return float(np.quantile(med, 0.025)), float(np.quantile(med, 0.975))


def prediction_trial_metrics(predictions_path: Path, events: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "event_id",
        "method",
        "keep_fraction",
        "random_index",
        "effective_train_density_per_10000km2",
        "effective_train_median_nn_distance_km",
        "observed_intensity_class_index",
        "class_match",
        "abs_error",
        "prediction_error",
    ]
    pred = pd.read_csv(predictions_path, usecols=cols)
    pred = pred.merge(events[["event_id", "current_jma_intensity_scale"]], on="event_id", how="left")
    pred["density_bin"] = density_bins(pred["effective_train_density_per_10000km2"]).astype(str)

    frames = [
        pred.assign(intensity_group="all"),
        pred[pred["observed_intensity_class_index"] >= 5].assign(intensity_group="I>=5-"),
        pred[pred["observed_intensity_class_index"] >= 7].assign(intensity_group="I>=6-"),
    ]
    long = pd.concat(frames, ignore_index=True)
    long["scope"] = np.where(long["current_jma_intensity_scale"], "current_scale_only", "pre_1996_unsplit_only")
    all_scope = long.assign(scope="all_catalog_6lower_plus")
    long = pd.concat([long, all_scope], ignore_index=True)

    grouped = long.groupby(
        ["scope", "event_id", "method", "keep_fraction", "random_index", "density_bin", "intensity_group"],
        observed=True,
    )
    trial = grouped.agg(
        n_predictions=("class_match", "size"),
        median_density_per_10000km2=("effective_train_density_per_10000km2", "median"),
        median_train_nn_km=("effective_train_median_nn_distance_km", "median"),
        class_accuracy=("class_match", "mean"),
        mae=("abs_error", "mean"),
        rmse=("prediction_error", lambda x: float(math.sqrt(np.mean(np.asarray(x, dtype=float) ** 2)))),
        under_rate=("prediction_error", lambda x: float(np.mean(np.asarray(x, dtype=float) < 0.0))),
        over_rate=("prediction_error", lambda x: float(np.mean(np.asarray(x, dtype=float) > 0.0))),
    ).reset_index()
    return trial


def summarize_trial_metrics(trial: pd.DataFrame, rng: np.random.Generator, n_boot: int) -> pd.DataFrame:
    rows = []
    group_cols = ["scope", "method", "density_bin", "intensity_group"]
    for keys, group in trial.groupby(group_cols, observed=True):
        scope, method, density_bin, intensity_group = keys
        lo, hi = bootstrap_ci(group["class_accuracy"].to_numpy(dtype=float), rng, n_boot)
        rows.append(
            {
                "scope": scope,
                "method": method,
                "density_bin": density_bin,
                "intensity_group": intensity_group,
                "n_trials": len(group),
                "n_events": int(group["event_id"].nunique()),
                "median_n_predictions_per_trial": float(group["n_predictions"].median()),
                "median_density_per_10000km2": float(group["median_density_per_10000km2"].median()),
                "median_train_nn_km": float(group["median_train_nn_km"].median())
                if group["median_train_nn_km"].notna().any()
                else np.nan,
                "median_class_accuracy": float(group["class_accuracy"].median()),
                "class_accuracy_ci95_low": lo,
                "class_accuracy_ci95_high": hi,
                "median_mae": float(group["mae"].median()),
                "median_rmse": float(group["rmse"].median()),
                "median_under_rate": float(group["under_rate"].median()),
                "median_over_rate": float(group["over_rate"].median()),
            }
        )
    out = pd.DataFrame(rows)
    out["method"] = pd.Categorical(out["method"], categories=METHOD_ORDER, ordered=True)
    intensity_order = pd.CategoricalDtype(categories=["all", "I>=5-", "I>=6-"], ordered=True)
    out["intensity_group"] = out["intensity_group"].astype(intensity_order)
    return out.sort_values(["scope", "intensity_group", "method", "median_density_per_10000km2"]).reset_index(drop=True)


def density_threshold_resolution(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (scope, intensity_group, method), group in summary.groupby(
        ["scope", "intensity_group", "method"], observed=True
    ):
        group = group.sort_values("median_density_per_10000km2")
        for target in [0.70, 0.80, 0.90]:
            median_ok = group[group["median_class_accuracy"] >= target]
            robust_ok = group[group["class_accuracy_ci95_low"] >= target]
            first = median_ok.iloc[0] if not median_ok.empty else None
            robust = robust_ok.iloc[0] if not robust_ok.empty else None
            rows.append(
                {
                    "scope": scope,
                    "intensity_group": intensity_group,
                    "method": method,
                    "target_class_accuracy": target,
                    "median_reached": first is not None,
                    "median_reached_density_per_10000km2": first["median_density_per_10000km2"]
                    if first is not None
                    else np.nan,
                    "median_reached_train_nn_km": first["median_train_nn_km"] if first is not None else np.nan,
                    "median_reached_accuracy": first["median_class_accuracy"] if first is not None else np.nan,
                    "robust_ci_reached": robust is not None,
                    "robust_ci_reached_density_per_10000km2": robust["median_density_per_10000km2"]
                    if robust is not None
                    else np.nan,
                    "robust_ci_low_at_reached": robust["class_accuracy_ci95_low"] if robust is not None else np.nan,
                    "final_interpretation": (
                        "supported_by_ci"
                        if robust is not None
                        else ("median_only" if first is not None else "not_reached")
                    ),
                }
            )
    out = pd.DataFrame(rows)
    out["method"] = pd.Categorical(out["method"], categories=METHOD_ORDER, ordered=True)
    return out.sort_values(["scope", "intensity_group", "target_class_accuracy", "method"]).reset_index(drop=True)


def uncertainty_selected_summary(prob_summary_path: Path, components_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    prob = pd.read_csv(prob_summary_path)
    components = pd.read_csv(components_path)
    selected_prob = prob[
        prob["observed_intensity_group"].isin(["all", "I>=5-"])
        & prob["method"].isin(["gmpe_raw", "idw", "kriging", "gmpe_kriging"])
    ].copy()
    highest_density = selected_prob.sort_values("median_density_per_10000km2").groupby(
        ["method", "observed_intensity_group"], observed=True
    ).tail(1)
    highest_components = components.sort_values("median_density_per_10000km2").groupby("method", observed=True).tail(1)
    return highest_density.reset_index(drop=True), highest_components.reset_index(drop=True)


def concern_status_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "concern": "旧震度6を現行震度6弱として扱う問題",
                "resolution": "1996-10-01以降の現行震度体系だけを別スコープ current_scale_only として集計した。",
                "claim_status": "対象定義は「旧震度6または現行震度6弱以上」と明記し、現行震度のみの感度分析を併記する。",
            },
            {
                "concern": "全観測点の70%一致が強震域精度を代表しない問題",
                "resolution": "保留観測点を all, I>=5-, I>=6- に分け、trial単位中央値とブートストラップCIを算出した。",
                "claim_status": "70%一致は全観測点限定の結果。強震点では今回の密度範囲で70%未到達。",
            },
            {
                "concern": "密度閾値が一点推定で不確実性を持たない問題",
                "resolution": "密度ビンごとのtrial単位階級一致率に95% bootstrap CIを付け、中央値到達とCI下限到達を分けた。",
                "claim_status": "密度閾値は決定的な十分条件ではなく、中央値ベースの参考値として扱う。",
            },
            {
                "concern": "距離減衰式と地盤増幅換算のばらつき",
                "resolution": "経験的予測sigmaと地盤増幅換算sigmaを導入した確率論的評価を追加した。",
                "claim_status": "点推定の階級一致だけでなく、真の階級への確率、予測区間被覆率、分散成分で報告する。",
            },
            {
                "concern": "GMPE基準場がイベントタイプ・断層面を十分に反映しない問題",
                "resolution": "現時点では物理的に完全解決できないため、GMPEを厳密な予測モデルではなく距離減衰ベースラインとして位置づける。",
                "claim_status": "距離減衰との改善量は参考比較。最終稿ではイベントタイプ別GMPE・断層面距離を今後課題として明記する。",
            },
        ]
    )


def setup_plot_style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 320,
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#222222",
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "grid.color": "#d8d8d8",
            "grid.linewidth": 0.55,
            "grid.alpha": 0.85,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
        }
    )


def savefig(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(fig)


def format_density_axis(ax, values: pd.Series) -> None:
    max_density = float(np.nanmax(values)) if len(values) else 1.0
    upper = max(1.0, max_density * 1.15)
    ticks = [0, 1, 3, 10, 30, 100, 300, 1000]
    ticks = [tick for tick in ticks if tick <= upper]
    ax.set_xscale("symlog", linthresh=1.0, linscale=0.75)
    ax.set_xlim(0, upper)
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0" if tick == 0 else f"{tick:g}" for tick in ticks])


def make_plots(
    trial_summary: pd.DataFrame,
    footprint: pd.DataFrame,
    uncertainty_components: pd.DataFrame,
    png_dir: Path,
) -> None:
    import matplotlib.pyplot as plt

    setup_plot_style()
    png_dir.mkdir(parents=True, exist_ok=True)
    colors = {
        "gmpe_raw": "#666666",
        "gmpe_calibrated": "#A6761D",
        "idw": "#D55E00",
        "kriging": "#CC79A7",
        "gmpe_kriging": "#7A3B00",
    }

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.2), sharex=True, sharey=True)
    panels = [
        ("all_catalog_6lower_plus", "all", "All catalog, all held-out stations"),
        ("all_catalog_6lower_plus", "I>=5-", "All catalog, held-out I >= 5-"),
        ("current_scale_only", "all", "Current scale only, all held-out stations"),
        ("current_scale_only", "I>=5-", "Current scale only, held-out I >= 5-"),
    ]
    for ax, (scope, intensity_group, title) in zip(axes.flat, panels):
        sub = trial_summary[(trial_summary["scope"] == scope) & (trial_summary["intensity_group"] == intensity_group)]
        for method, group in sub.groupby("method", observed=True):
            group = group.sort_values("median_density_per_10000km2")
            ax.plot(
                group["median_density_per_10000km2"],
                group["median_class_accuracy"],
                marker="o",
                linewidth=1.6,
                color=colors.get(str(method)),
                label=str(method),
            )
            ax.fill_between(
                group["median_density_per_10000km2"],
                group["class_accuracy_ci95_low"],
                group["class_accuracy_ci95_high"],
                color=colors.get(str(method)),
                alpha=0.12,
                linewidth=0,
            )
        ax.axhline(0.70, color="#222222", linestyle="--", linewidth=0.9)
        ax.set_title(title)
        format_density_axis(ax, sub["median_density_per_10000km2"])
    axes[0, 0].set_ylabel("Exact class hit rate")
    axes[1, 0].set_ylabel("Exact class hit rate")
    axes[1, 0].set_xlabel("Station density (/10,000 km$^2$)")
    axes[1, 1].set_xlabel("Station density (/10,000 km$^2$)")
    axes[0, 1].legend(loc="lower right", ncol=1)
    savefig(fig, png_dir / "class_accuracy_concern_resolution.png")

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    fp = footprint[footprint["method"].isin(["gmpe_raw", "idw", "kriging", "gmpe_kriging"])].copy()
    x_pos = np.arange(len(fp["method"].cat.categories if hasattr(fp["method"], "cat") else METHOD_ORDER))
    methods = [m for m in METHOD_ORDER if m in set(fp["method"].astype(str))]
    width = 0.36
    for i, scope in enumerate(["all_catalog_6lower_plus", "current_scale_only"]):
        sub = fp[fp["scope"] == scope].set_index(fp[fp["scope"] == scope]["method"].astype(str)).reindex(methods)
        ax.bar(
            np.arange(len(methods)) + (i - 0.5) * width,
            sub["area_ge_6p0_positive_event_rate"],
            width=width,
            label=scope,
            color="#666666" if i == 0 else "#D55E00",
            alpha=0.75,
        )
    ax.set_xticks(np.arange(len(methods)))
    ax.set_xticklabels(methods, rotation=20, ha="right")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Fraction of events with nonzero I >= 6.0 area")
    ax.set_title("Footprint sensitivity to pre-1996 unsplit intensity classes")
    ax.legend()
    savefig(fig, png_dir / "footprint_current_scale_sensitivity.png")

    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    comp = uncertainty_components[uncertainty_components["method"].isin(["gmpe_raw", "idw", "kriging", "gmpe_kriging"])]
    comp = comp.sort_values("median_density_per_10000km2").groupby("method", observed=True).tail(1)
    comp = comp.set_index("method").reindex([m for m in METHOD_ORDER if m in set(comp.index)]).reset_index()
    ax.bar(comp["method"], comp["sigma_empirical_median"], color="#555555", alpha=0.75, label="Empirical sigma")
    ax.bar(
        comp["method"],
        comp["sigma_site_conversion_median"],
        color="#0072B2",
        alpha=0.85,
        label="Site-conversion sigma",
    )
    ax.set_ylabel("Intensity sigma")
    ax.set_title("Uncertainty components at the highest density bin")
    ax.tick_params(axis="x", rotation=20)
    ax.legend()
    savefig(fig, png_dir / "highest_density_uncertainty_components.png")


def main() -> None:
    args = parse_args()
    args.csv_dir.mkdir(parents=True, exist_ok=True)
    args.png_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    events = add_event_flags(pd.read_csv(args.target_events, low_memory=False))
    footprints = pd.read_csv(args.footprints)

    event_sensitivity = event_selection_sensitivity(events)
    footprint = footprint_sensitivity(footprints, events)
    trial = prediction_trial_metrics(args.predictions, events)
    trial_summary = summarize_trial_metrics(trial, rng, args.bootstrap)
    threshold = density_threshold_resolution(trial_summary)
    selected_prob, selected_components = uncertainty_selected_summary(
        args.probabilistic_summary, args.uncertainty_components
    )
    status = concern_status_table()

    event_sensitivity.to_csv(args.csv_dir / "event_selection_sensitivity.csv", index=False)
    footprint.to_csv(args.csv_dir / "footprint_current_scale_sensitivity.csv", index=False)
    trial.to_csv(args.csv_dir / "trial_metrics_by_scope_intensity_group.csv.gz", index=False, compression="gzip")
    trial_summary.to_csv(args.csv_dir / "trial_metric_bootstrap_summary.csv", index=False)
    threshold.to_csv(args.csv_dir / "density_threshold_resolution.csv", index=False)
    selected_prob.to_csv(args.csv_dir / "probabilistic_highest_density_selected.csv", index=False)
    selected_components.to_csv(args.csv_dir / "uncertainty_components_highest_density_selected.csv", index=False)
    status.to_csv(args.csv_dir / "concern_resolution_status.csv", index=False)

    if not args.no_png:
        make_plots(trial_summary, footprint, selected_components, args.png_dir)

    print(f"Saved CSV directory: {args.csv_dir}")
    if not args.no_png:
        print(f"Saved PNG directory: {args.png_dir}")


if __name__ == "__main__":
    main()
