# -*- coding: utf-8 -*-
"""
_probe_metric_all.py -- the full 11-window geographic-vs-metric comparison.

Reads the cached WorldCover windows, builds the tree-cover fraction field on
both supports, and reports every quantity the manuscript's real-data section
depends on, so that the revision plan can state measured impacts rather than
guesses.

Writes a CSV next to the script and a human-readable report on stdout.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from rasterio.transform import from_origin
from rasterio.warp import Resampling, calculate_default_transform, reproject

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
sys.path.insert(0, str(ROOT / "code"))

from windows_v2 import (BIOME, CENTRE, LABELS, PLACE, WINDOWS,  # noqa: E402
                        _tile_origin, label_of)

TREE = 10
CELLS_PER_DEG = 12000.0
WINDOW = 3072
BLOCK_M = (4, 9, 16, 25)
TARGET_M = 10.0
OUT = ROOT / "_metric_vs_geographic.csv"

# A reprojection must preserve the class fraction and the footprint area; if
# either check fails the sign or the transform is wrong and every downstream
# number is garbage.  The first version of this script computed the window
# origin by hand and dropped the W sign, so the two western-hemisphere windows
# were projected 15 degrees away from their own projection centre and their
# tree fraction fell by 11 percentage points -- a silent, plausible-looking
# failure that only the full 11-window run exposed.
MAX_MEAN_SHIFT = 0.02        # 2 percentage points
MAX_AREA_ERROR = 0.10        # 10 %


def load(tile: str, row: int, col: int) -> np.ndarray:
    arr = np.load(DATA / f"worldcover_{tile}_{row}_{col}_{WINDOW}.npy")
    return (arr == TREE).astype(np.uint8)


def group3(binary: np.ndarray) -> np.ndarray:
    ny, nx = binary.shape
    a = binary[: ny // 3 * 3, : nx // 3 * 3].astype(np.float64)
    return a.reshape(ny // 3, 3, nx // 3, 3).mean(axis=(1, 3))


def block_var(fine: np.ndarray, m: int) -> tuple[float, int]:
    side = int(round(m ** 0.5))
    ny, nx = fine.shape
    nyb, nxb = ny // side, nx // side
    if nyb < 4 or nxb < 4:
        return float("nan"), 0
    b = fine[: nyb * side, : nxb * side].reshape(nyb, side, nxb, side).mean(axis=(1, 3))
    return float(b.var()), nyb * nxb


def geographic(arr: np.ndarray) -> np.ndarray:
    return group3(arr)


def metric(arr: np.ndarray, tile: str, row: int, col: int) -> tuple[np.ndarray, dict]:
    lab = label_of(tile, row, col)
    lat_c, lon_c = CENTRE[lab]
    # signed tile origin: N48E006 -> (48, 6), N39W009 -> (39, -9)
    lat_tile, lon_tile = _tile_origin(tile)
    lat0 = lat_tile + 3.0 - row / CELLS_PER_DEG
    lon0 = lon_tile + col / CELLS_PER_DEG
    left, top = lon0, lat0
    right, bottom = lon0 + WINDOW / CELLS_PER_DEG, lat0 - WINDOW / CELLS_PER_DEG
    src_tf = from_origin(lon0, lat0, 1.0 / CELLS_PER_DEG, 1.0 / CELLS_PER_DEG)
    dst_crs = (f"+proj=aeqd +lat_0={lat_c:.6f} +lon_0={lon_c:.6f} "
               "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    dst_tf, dw, dh = calculate_default_transform(
        "EPSG:4326", dst_crs, WINDOW, WINDOW, left, bottom, right, top,
        resolution=TARGET_M)
    dst = np.zeros((dh, dw), dtype=arr.dtype)
    reproject(source=arr, destination=dst, src_transform=src_tf, src_crs="EPSG:4326",
              dst_transform=dst_tf, dst_crs=dst_crs, resampling=Resampling.nearest)

    ns_cell = 110574.0 / CELLS_PER_DEG
    ew_cell = 111320.0 * np.cos(np.deg2rad(lat_c)) / CELLS_PER_DEG
    exp_area_m2 = (WINDOW * ns_cell) * (WINDOW * ew_cell)
    got_area_m2 = (dh * TARGET_M) * (dw * TARGET_M)
    area_err = got_area_m2 / exp_area_m2 - 1.0
    assert abs(area_err) < MAX_AREA_ERROR, (
        f"{lab}: destination footprint {got_area_m2/1e6:.1f} Mm2 vs expected "
        f"{exp_area_m2/1e6:.1f} Mm2 ({area_err:+.1%}) -- transform is wrong")

    return group3(dst), {"dst_shape": (dh, dw), "area_err": area_err}


def main() -> None:
    rows = []
    print(f"{'window':<28}{'lat':>7}{'aspect':>8}"
          f"{'pxvar g':>9}{'pxvar m':>9}"
          f"{'m4 g':>7}{'m4 m':>7}{'d%':>7}"
          f"{'m25 g':>8}{'m25 m':>8}{'d%':>7}")
    print("-" * 108)

    for tile, row, col, biome, group, place in WINDOWS:
        lab = label_of(tile, row, col)
        lat_c, lon_c = CENTRE[lab]
        ns = 3 * (110574.0 / CELLS_PER_DEG)
        ew = 3 * (111320.0 * np.cos(np.deg2rad(lat_c)) / CELLS_PER_DEG)

        arr = load(tile, row, col)
        fg, geo = geographic(arr), None
        fm, geo = metric(arr, tile, row, col)

        # the reprojection must preserve the class fraction
        shift = abs(float(fm.mean()) - float(fg.mean()))
        assert shift < MAX_MEAN_SHIFT, (
            f"{lab}: tree fraction moved {shift:.4f} under reprojection "
            f"({fg.mean():.4f} -> {fm.mean():.4f}) -- check the CRS and the sign "
            f"of the longitude")

        rec = {
            "window": lab, "tile": tile, "biome": biome, "place": place,
            "lat": round(lat_c, 4), "lon": round(lon_c, 4),
            "ns_m": round(ns, 3), "ew_m": round(ew, 3), "aspect": round(ns / ew, 4),
            "mean_geo": round(float(fg.mean()), 5), "mean_met": round(float(fm.mean()), 5),
            "pxvar_geo": round(float(fg.var()), 6), "pxvar_met": round(float(fm.var()), 6),
            "met_shape": f"{fm.shape[0]}x{fm.shape[1]}",
        }
        for m in BLOCK_M:
            vt_g, nb_g = block_var(fg, m)
            vt_m, nb_m = block_var(fm, m)
            vi_g = float(fg.var()) / m
            vi_m = float(fm.var()) / m
            rec[f"m{m}_vg"] = round(vt_g, 6)
            rec[f"m{m}_vm"] = round(vt_m, 6)
            rec[f"m{m}_ug"] = round(vt_g / vi_g, 3) if vi_g else None
            rec[f"m{m}_um"] = round(vt_m / vi_m, 3) if vi_m else None
            rec[f"m{m}_d_pct"] = round(100 * (vt_m - vt_g) / vt_g, 2) if vt_g else None
            rec[f"m{m}_nb_geo"] = nb_g
            rec[f"m{m}_nb_met"] = nb_m
        rows.append(rec)

        print(f"{lab:<28}{lat_c:>7.2f}{ns/ew:>8.3f}"
              f"{rec['pxvar_geo']:>9.4f}{rec['pxvar_met']:>9.4f}"
              f"{rec['m4_ug']:>7.2f}{rec['m4_um']:>7.2f}{rec['m4_d_pct']:>7.1f}"
              f"{rec['m25_ug']:>8.2f}{rec['m25_um']:>8.2f}{rec['m25_d_pct']:>7.1f}")

    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("-" * 108)
    d4 = [r["m4_d_pct"] for r in rows]
    d25 = [r["m25_d_pct"] for r in rows]
    print(f"m=4  V_true change : {min(d4):+.1f}% to {max(d4):+.1f}%   "
          f"mean {np.mean(d4):+.1f}%")
    print(f"m=25 V_true change : {min(d25):+.1f}% to {max(d25):+.1f}%   "
          f"mean {np.mean(d25):+.1f}%")
    print(f"pixel-scale var    : "
          f"{min(100*(r['pxvar_met']-r['pxvar_geo'])/r['pxvar_geo'] for r in rows):+.1f}% to "
          f"{max(100*(r['pxvar_met']-r['pxvar_geo'])/r['pxvar_geo'] for r in rows):+.1f}%")
    print(f"tree fraction      : "
          f"{min(100*(r['mean_met']-r['mean_geo']) for r in rows):+.3f} to "
          f"{max(100*(r['mean_met']-r['mean_geo']) for r in rows):+.3f} percentage points")
    ug = [r["m4_ug"] for r in rows]
    um = [r["m4_um"] for r in rows]
    print(f"m=4 understatement spread (max/min): geography {max(ug)/min(ug):.3f} | "
          f"metric {max(um)/min(um):.3f}")
    ug25 = [r["m25_ug"] for r in rows]
    um25 = [r["m25_um"] for r in rows]
    print(f"m=25 understatement spread (max/min): geography {max(ug25)/min(ug25):.3f} | "
          f"metric {max(um25)/min(um25):.3f}")
    print(f"\nCSV: {OUT}")

    # biome ordering check: does the biome ranking by mean bias change?
    print("\nbiome mean of the m=25 understatement factor:")
    for b in sorted({r["biome"] for r in rows}):
        g = [r["m25_ug"] for r in rows if r["biome"] == b]
        m = [r["m25_um"] for r in rows if r["biome"] == b]
        print(f"  {b:<16} geography {np.mean(g):6.2f}   metric {np.mean(m):6.2f}")


if __name__ == "__main__":
    main()
