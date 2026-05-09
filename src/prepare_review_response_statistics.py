#!/usr/bin/env python3
"""Prepare compact review-response statistics from station-thinning predictions."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_PREDICTIONS = Path(
    "outputs/csv/station_thinning_interpolation_6lower_plus_class/thinning_prediction_errors.csv.gz"
)
DEFAULT_OUT_DIR = Path("data/derived/station_thinning_interpolation_6lower_plus_class")
METHODS = ["idw", "kriging", "gmpe_kriging"]
BINS = [0, 5, 10, 20, 30, 50, 75, np.inf]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260509)
    return parser.parse_args()


def density_bin(series: pd.Series) -> pd.Series:
    return pd.cut(series, bins=BINS, right=False).astype(str)


def rmse(values: pd.Series) -> float:
    arr = values.to_numpy(dtype=float)
    return float(math.sqrt(np.mean(arr**2))) if len(arr) else float("nan")


def metric_row(group: pd.DataFrame) -> dict[str, float | int]:
    return {
        "n": int(len(group)),
        "n_events": int(group["event_id"].nunique()),
        "mae": float(group["abs_error"].mean()),
        "rmse": rmse(group["prediction_error"]),
        "exact_class_accuracy": float(group["class_match"].mean()),
        "underprediction_rate": float((group["class_error"] < 0).mean()),
    }


def event_bootstrap_ci(
    group: pd.DataFrame,
    rng: np.random.Generator,
    n_bootstrap: int,
    metric: str = "exact_class_accuracy",
) -> tuple[float, float]:
    event_ids = group["event_id"].drop_duplicates().to_numpy()
    if len(event_ids) < 2 or n_bootstrap <= 0:
        value = metric_row(group)[metric]
        return float(value), float(value)

    by_event = {event_id: ev for event_id, ev in group.groupby("event_id", sort=False)}
    values = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        sample_ids = rng.choice(event_ids, size=len(event_ids), replace=True)
        sampled = pd.concat([by_event[event_id] for event_id in sample_ids], ignore_index=False)
        values[i] = metric_row(sampled)[metric]
    return (float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)))


def high_density_subset_table(df: pd.DataFrame, rng: np.random.Generator, n_bootstrap: int) -> pd.DataFrame:
    high = df[
        df["method"].isin(METHODS)
        & (df["effective_train_density_per_10000km2"] >= 50)
        & (df["effective_train_density_per_10000km2"] < 75)
    ].copy()
    subsets = [
        ("I < 5-", high["observed_intensity_class_value"] < 5.0),
        ("I >= 5-", high["observed_intensity_class_value"] >= 5.0),
        ("I >= 6-", high["observed_intensity_class_value"] >= 6.0),
    ]

    rows: list[dict[str, object]] = []
    for subset_name, mask in subsets:
        subset = high[mask]
        for method in METHODS:
            group = subset[subset["method"] == method]
            if group.empty:
                continue
            row = {"validation_subset": subset_name, "method": method}
            row.update(metric_row(group))
            lo, hi = event_bootstrap_ci(group, rng, n_bootstrap)
            row["exact_class_accuracy_ci95_low"] = lo
            row["exact_class_accuracy_ci95_high"] = hi
            rows.append(row)
    return pd.DataFrame(rows)


def event_equal_density_table(df: pd.DataFrame, rng: np.random.Generator, n_bootstrap: int) -> pd.DataFrame:
    use = df[df["method"].isin(METHODS)].copy()
    use["density_bin"] = density_bin(use["effective_train_density_per_10000km2"])

    event_rows: list[dict[str, object]] = []
    for (method, bin_label, event_id), group in use.groupby(["method", "density_bin", "event_id"], observed=True):
        if bin_label == "nan" or group.empty:
            continue
        row = {"method": method, "density_bin": bin_label, "event_id": event_id}
        row.update(metric_row(group))
        row["median_density_per_10000km2"] = float(group["effective_train_density_per_10000km2"].median())
        row["median_train_nn_km"] = float(group["effective_train_median_nn_distance_km"].median())
        event_rows.append(row)

    event_df = pd.DataFrame(event_rows)
    rows: list[dict[str, object]] = []
    for (method, bin_label), group in event_df.groupby(["method", "density_bin"], observed=True):
        values = group["exact_class_accuracy"].to_numpy(dtype=float)
        if len(values) >= 2 and n_bootstrap > 0:
            boot = np.empty(n_bootstrap, dtype=float)
            for i in range(n_bootstrap):
                boot[i] = rng.choice(values, size=len(values), replace=True).mean()
            ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
        elif len(values) == 1:
            ci_low = ci_high = values[0]
        else:
            ci_low = ci_high = np.nan

        rows.append(
            {
                "method": method,
                "density_bin": bin_label,
                "n_events": int(group["event_id"].nunique()),
                "n_predictions": int(group["n"].sum()),
                "median_density_per_10000km2": float(group["median_density_per_10000km2"].median()),
                "median_train_nn_km": float(group["median_train_nn_km"].median()),
                "event_equal_mae": float(group["mae"].mean()),
                "event_equal_rmse": float(group["rmse"].mean()),
                "event_equal_exact_class_accuracy": float(values.mean()),
                "event_equal_exact_class_accuracy_ci95_low": float(ci_low),
                "event_equal_exact_class_accuracy_ci95_high": float(ci_high),
                "event_equal_underprediction_rate": float(group["underprediction_rate"].mean()),
            }
        )
    out = pd.DataFrame(rows)
    out["_density_left"] = out["density_bin"].str.extract(r"\[([^,]+),").astype(float)
    out = out.sort_values(["method", "_density_left"]).drop(columns="_density_left")
    return out


def sufficient_density_event_equal(event_equal: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method in METHODS:
        group = event_equal[event_equal["method"] == method].copy()
        group["_density_left"] = group["density_bin"].str.extract(r"\[([^,]+),").astype(float)
        group = group.sort_values("_density_left")
        for target in [0.70, 0.80, 0.90]:
            hit = group[group["event_equal_exact_class_accuracy"] >= target].head(1)
            if hit.empty:
                rows.append(
                    {
                        "method": method,
                        "criterion": f"event_equal_exact>={target:.2f}",
                        "density_per_10000km2": np.nan,
                        "median_nn_km": np.nan,
                        "exact_class_accuracy": np.nan,
                        "status": "not_reached",
                    }
                )
            else:
                r = hit.iloc[0]
                rows.append(
                    {
                        "method": method,
                        "criterion": f"event_equal_exact>={target:.2f}",
                        "density_per_10000km2": r["median_density_per_10000km2"],
                        "median_nn_km": r["median_train_nn_km"],
                        "exact_class_accuracy": r["event_equal_exact_class_accuracy"],
                        "status": "reached",
                    }
                )
    return pd.DataFrame(rows)


def station_region(station_code: object) -> str:
    code = str(station_code).zfill(7)
    try:
        prefix = int(code[:2])
    except ValueError:
        return "Other"
    if 10 <= prefix <= 16:
        return "Hokkaido"
    if 20 <= prefix <= 25:
        return "Tohoku"
    if 30 <= prefix <= 36:
        return "Kanto"
    if 37 <= prefix <= 49:
        return "Chubu-Hokuriku"
    if 50 <= prefix <= 63 or prefix == 70:
        return "Kinki-Chugoku-Shikoku"
    if 71 <= prefix <= 77:
        return "Kyushu"
    return "Other"


def regional_reference_table(df: pd.DataFrame) -> pd.DataFrame:
    region_order = [
        "Hokkaido",
        "Tohoku",
        "Kanto",
        "Chubu-Hokuriku",
        "Kinki-Chugoku-Shikoku",
        "Kyushu",
    ]
    work = df.copy()
    work["station_region"] = work["station_code"].map(station_region)
    raw = work[work["method"] == "gmpe_raw"]
    high = work[
        (work["method"] == "kriging")
        & (work["effective_train_density_per_10000km2"] >= 50)
        & (work["effective_train_density_per_10000km2"] < 75)
    ]

    rows: list[dict[str, object]] = []
    raw_stats: dict[str, dict[str, float | int]] = {}
    for region, group in raw.groupby("station_region", sort=False):
        raw_stats[region] = metric_row(group)

    for region in region_order:
        raw_row = raw_stats.get(region)
        if raw_row is None:
            continue
        rows.append(
            {
                "station_region": region,
                "method": "gmpe_raw",
                "n_predictions": raw_row["n"],
                "exact_class_accuracy": raw_row["exact_class_accuracy"],
                "rmse": raw_row["rmse"],
                "accuracy_gain_points_vs_raw": 0.0,
                "rmse_reduction_fraction_vs_raw": 0.0,
                "comparison_note": "raw_all_validation_points",
            }
        )
        group = high[high["station_region"] == region]
        if group.empty:
            continue
        high_row = metric_row(group)
        rows.append(
            {
                "station_region": region,
                "method": "kriging",
                "n_predictions": high_row["n"],
                "exact_class_accuracy": high_row["exact_class_accuracy"],
                "rmse": high_row["rmse"],
                "accuracy_gain_points_vs_raw": high_row["exact_class_accuracy"] - raw_row["exact_class_accuracy"],
                "rmse_reduction_fraction_vs_raw": (raw_row["rmse"] - high_row["rmse"]) / raw_row["rmse"],
                "comparison_note": "kriging_high_density_vs_raw_region_reference",
            }
        )
    return pd.DataFrame(rows)


def station_region_density_table(df: pd.DataFrame) -> pd.DataFrame:
    work = df[df["method"] == "kriging"].copy()
    work["station_region"] = work["station_code"].map(station_region)
    work["density_bin"] = density_bin(work["effective_train_density_per_10000km2"])

    rows: list[dict[str, object]] = []
    for (region, bin_label), group in work.groupby(["station_region", "density_bin"], observed=True):
        if region == "Other" or bin_label == "nan" or group.empty:
            continue
        rows.append(
            {
                "station_region": region,
                "density_bin": bin_label,
                "n_predictions": int(len(group)),
                "median_density_per_10000km2": float(group["effective_train_density_per_10000km2"].median()),
                "exact_class_accuracy": float(group["class_match"].mean()),
                "rmse": rmse(group["prediction_error"]),
            }
        )
    out = pd.DataFrame(rows)
    out["_density_left"] = out["density_bin"].str.extract(r"\[([^,]+),").astype(float)
    return out.sort_values(["station_region", "_density_left"]).drop(columns="_density_left")


def main() -> None:
    args = parse_args()
    cols = [
        "station_code",
        "event_id",
        "method",
        "effective_train_density_per_10000km2",
        "effective_train_median_nn_distance_km",
        "observed_intensity_class_value",
        "prediction_error",
        "abs_error",
        "class_error",
        "class_match",
    ]
    df = pd.read_csv(args.predictions, usecols=cols)
    rng = np.random.default_rng(args.seed)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    high_density = high_density_subset_table(df, rng, args.bootstrap)
    event_equal = event_equal_density_table(df, rng, args.bootstrap)
    thresholds = sufficient_density_event_equal(event_equal)
    regional = regional_reference_table(df)
    regional_density = station_region_density_table(df)

    high_density.to_csv(args.out_dir / "high_density_intensity_subset_accuracy_ci.csv", index=False)
    event_equal.to_csv(args.out_dir / "event_equal_density_bin_summary.csv", index=False)
    thresholds.to_csv(args.out_dir / "event_equal_sufficient_density_thresholds.csv", index=False)
    regional.to_csv(args.out_dir / "station_region_high_density_reference_comparison.csv", index=False)
    regional_density.to_csv(args.out_dir / "station_region_density_accuracy.csv", index=False)

    print(f"Wrote {args.out_dir / 'high_density_intensity_subset_accuracy_ci.csv'}")
    print(f"Wrote {args.out_dir / 'event_equal_density_bin_summary.csv'}")
    print(f"Wrote {args.out_dir / 'event_equal_sufficient_density_thresholds.csv'}")
    print(f"Wrote {args.out_dir / 'station_region_high_density_reference_comparison.csv'}")
    print(f"Wrote {args.out_dir / 'station_region_density_accuracy.csv'}")


if __name__ == "__main__":
    main()
