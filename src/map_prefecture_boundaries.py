#!/usr/bin/env python3
"""Utilities for drawing Japanese prefectural boundaries on map figures."""

from __future__ import annotations

from pathlib import Path

import numpy as np


DEFAULT_PREFECTURE_BOUNDARY = Path(
    "/Users/mitsuki/02_作業/001_地震/01_地震速報/2025年7月30日_カムチャッカ半島M88/data/gmt_data/Kenkyo.txt"
)


def resolve_prefecture_boundary(path: Path | str | None = None) -> Path | None:
    boundary = Path(path) if path else DEFAULT_PREFECTURE_BOUNDARY
    return boundary if boundary.exists() else None


def read_gmt_multisegment(path: Path | str | None = None) -> list[np.ndarray]:
    """Read a GMT multi-segment lon/lat text file."""

    boundary = resolve_prefecture_boundary(path)
    if boundary is None:
        return []
    segments: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []
    with boundary.open(encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    segments.append(current)
                    current = []
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue
            try:
                lon = float(parts[0])
                lat = float(parts[1])
            except ValueError:
                continue
            if np.isfinite(lon) and np.isfinite(lat):
                current.append((lon, lat))
    if current:
        segments.append(current)
    return [np.asarray(segment, dtype=float) for segment in segments if len(segment) >= 2]


def plot_prefecture_boundaries_pygmt(
    fig,
    path: Path | str | None = None,
    pen: str = "0.18p,#666666",
) -> None:
    boundary = resolve_prefecture_boundary(path)
    if boundary is None:
        return
    fig.plot(data=str(boundary), pen=pen)


def plot_prefecture_boundaries_cartopy(
    ax,
    path: Path | str | None = None,
    color: str = "#666666",
    linewidth: float = 0.32,
    alpha: float = 0.90,
    zorder: int = 2,
) -> None:
    import cartopy.crs as ccrs

    for segment in read_gmt_multisegment(path):
        ax.plot(
            segment[:, 0],
            segment[:, 1],
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            transform=ccrs.PlateCarree(),
            zorder=zorder,
        )
