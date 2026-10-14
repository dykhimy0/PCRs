# -*- coding: utf-8 -*-
"""
_probe_metric_support.py -- does the paper's real-data result survive a
reprojection from the geographic grid to a metric support?

The reviewer's hardest objection is that the eleven windows are read on the
EPSG:4326 grid, where a fixed cell count spans a constant north-south distance
but a shrinking east-west one, so at high latitude a "3 x 3 block" is not 30 m
of ground and the scale geometry is mixed with the spatial dependence.

This probe answers the question that decides the paper's fate:

    recompute the headline statistic (the factor by which independence
    understates the block variance) on a metric support, and on the geographic
    grid, for three windows spanning the latitude range.

Nothing here is written into the manuscript; it is a measurement.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import Resampling, calculate_default_transform, reproject

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
sys.path.insert(0, str(ROOT / "code"))

from windows_v2 import CENTRE, WINDOWS  # noqa: E402

TREE = 10
CELLS_PER_DEG = 12000.0
WINDOW = 3072
BLOCK_M = (4, 25)          # m = 4 (2x2 of 30 m) and m = 25 (5x5 of 30 m)
TARGET_M = 10.0            # metric cell size, matches the product's nominal 10 m

# three windows spanning the latitude range 22-67 N
PROBE = [
    ("N48E006", 26000, 26000, "temperate, 48.7 N"),
    ("N21E099", 24000, 12000, "tropical, 22.5 N"),
    ("N66E024", 26000, 20000, "boreal, 66.8 N"),
]


def load(tile: str, row: int, col: int) -> np.ndarray:
    """Raw WorldCover class codes, recoded to a tree / non-tree indicator."""
    p = DATA / f"worldcover_{tile}_{row}_{col}_{WINDOW}.npy"
    arr = np.load(p)
    return (arr == TREE).astype(np.uint8)


def blocks_var(fine: np.ndarray, m: int) -> tuple[float, float, int]:
    """Block variance for an m-cell block, on a square cell grid.

    `fine` is the 30 m fraction field.  m = 4 -> 2x2, m = 25 -> 5x5.
    Returns (true block variance, independence variance, number of blocks).
    """
    side = int(round(m ** 0.5))
    assert side * side == m, m
    ny, nx = fine.shape
    ny_b, nx_b = ny // side, nx // side
    trimmed = fine[:ny_b * side, :nx_b * side]
    block = trimmed.reshape(ny_b, side, nx_b, side).mean(axis=(1, 3))
    v_true = float(block.var())
    v_ind = float(fine.var()) / m
    return v_true, v_ind, ny_b * nx_b


def analyse(fine: np.ndarray, tag: str) -> dict:
    out = {"tag": tag, "pixel_var": float(fine.var()), "mean": float(fine.mean())}
    for m in BLOCK_M:
        v_true, v_ind, nb = blocks_var(fine, m)
        out[f"m{m}_true"] = v_true
        out[f"m{m}_ind"] = v_ind
        out[f"m{m}_under"] = v_true / v_ind if v_ind > 0 else float("nan")
        out[f"m{m}_deff"] = m * v_true / out["pixel_var"]
        out[f"m{m}_nblocks"] = nb
    return out


def geographic(arr: np.ndarray) -> np.ndarray:
    """The paper's current pipeline: 3x3 cells of the EPSG:4326 grid."""
    ny, nx = arr.shape
    ny3, nx3 = ny // 3, nx // 3
    a = arr[:ny3 * 3, :nx3 * 3].astype(np.float64)
    return a.reshape(ny3, 3, nx3, 3).mean(axis=(1, 3))


def metric(arr: np.ndarray, tile: str, row: int, col: int) -> tuple[np.ndarray, dict]:
    """Reproject the window to a local azimuthal-equidistant metric grid."""
    lat_c, lon_c = CENTRE[f"{tile} ({row // 1000}k,{col // 1000}k)"]
    # source grid: window top-left corner in degrees
    lat0 = float(tile[1:3]) + 3.0 - row / CELLS_PER_DEG
    lon0 = float(tile[4:]) + col / CELLS_PER_DEG
    src_tf = from_origin(lon0, lat0, 1.0 / CELLS_PER_DEG, 1.0 / CELLS_PER_DEG)
    src_crs = "EPSG:4326"
    dst_crs = (f"+proj=aeqd +lat_0={lat_c:.6f} +lon_0={lon_c:.6f} "
               "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

    left, top = lon0, lat0
    right = lon0 + WINDOW / CELLS_PER_DEG
    bottom = lat0 - WINDOW / CELLS_PER_DEG

    dst_tf, dst_w, dst_h = calculate_default_transform(
        src_crs, dst_crs, WINDOW, WINDOW, left, bottom, right, top,
        resolution=TARGET_M)

    dst = np.zeros((dst_h, dst_w), dtype=arr.dtype)
    reproject(
        source=arr, destination=dst,
        src_transform=src_tf, src_crs=src_crs,
        dst_transform=dst_tf, dst_crs=dst_crs,
        resampling=Resampling.nearest)

    # 3x3 metric cells -> 30 m fraction
    ny3, nx3 = dst_h // 3, dst_w // 3
    a = dst[:ny3 * 3, :nx3 * 3].astype(np.float64)
    fine = a.reshape(ny3, 3, nx3, 3).mean(axis=(1, 3))

    # ground geometry, for the record
    a_px = abs(dst_tf.a)
    geo_ew = WINDOW / CELLS_PER_DEG * 111320.0 * np.cos(np.deg2rad(lat_c))
    geo_ns = WINDOW / CELLS_PER_DEG * 110574.0
    return fine, {
        "dst_shape": (dst_h, dst_w),
        "metric_cell_m": (a_px, abs(dst_tf.e)),
        "geographic_window_m": (geo_ew, geo_ns),
        "geographic_block30_m": (30 * geo_ew / (WINDOW * a_px / TARGET_M) if False else
                                 3 * geo_ew / WINDOW * 1.0, 3 * geo_ns / WINDOW * 1.0),
        "geographic_aspect": geo_ns / geo_ew,
    }


def main() -> None:
    print("=" * 78)
    print("PROBE: geographic grid vs metric support, three windows")
    print("=" * 78)

    for tile, row, col, note in PROBE:
        arr = load(tile, row, col)
        fine_g = geographic(arr)
        fine_m, geo = metric(arr, tile, row, col)

        # physical size of one 30 m "pixel" under the geographic grid
        lat_c, lon_c = CENTRE[f"{tile} ({row // 1000}k,{col // 1000}k)"]
        ew_m = 3 * (111320.0 * np.cos(np.deg2rad(lat_c)) / CELLS_PER_DEG)
        ns_m = 3 * (110574.0 / CELLS_PER_DEG)

        print(f"\n--- {tile} ({row//1000}k,{col//1000}k)  {note} ---")
        print(f"    centre lat/lon      : {lat_c:.3f} N, {lon_c:.3f} E")
        print(f"    geographic 30 m px  : {ns_m:6.2f} m N-S x {ew_m:6.2f} m E-W"
              f"   (aspect {ns_m/ew_m:.3f})")
        print(f"    metric  30 m px     : 30.00 m N-S x 30.00 m E-W   (aspect 1.000)")
        print(f"    metric grid shape   : {geo['dst_shape']}, cell "
              f"{geo['metric_cell_m'][0]:.2f} m")

        g = analyse(fine_g, "geographic")
        m = analyse(fine_m, "metric")
        print(f"    tree fraction       : geographic {g['mean']:.4f} | "
              f"metric {m['mean']:.4f}")
        print(f"    pixel-scale var     : geographic {g['pixel_var']:.4f} | "
              f"metric {m['pixel_var']:.4f}")
        print(f"    {'m':>3} {'V_true geo':>11} {'V_true mtr':>11} "
              f"{'indep geo':>10} {'under geo':>10} {'under mtr':>10} {'delta %':>8}")
        for k in BLOCK_M:
            d = 100.0 * (m[f"m{k}_under"] - g[f"m{k}_under"]) / g[f"m{k}_under"]
            print(f"    {k:>3} {g[f'm{k}_true']:>11.4f} {m[f'm{k}_true']:>11.4f} "
                  f"{g[f'm{k}_ind']:>10.4f} {g[f'm{k}_under']:>10.2f}x "
                  f"{m[f'm{k}_under']:>9.2f}x {d:>+7.1f}%")
        print(f"    DEFF geo/mtr m=25   : {g['m25_deff']:.3f} / {m['m25_deff']:.3f}")

    print("\n" + "=" * 78)
    print("Read: 'under' is the factor by which independence understates the block")
    print("variance.  If the two columns agree to within a few per cent, the")
    print("reprojection is a geometry correction and the headline survives; if they")
    print("diverge systematically with latitude, the current numbers are an artefact")
    print("of the geographic grid and must be recomputed.")
    print("=" * 78)


if __name__ == "__main__":
    main()
