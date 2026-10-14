# -*- coding: utf-8 -*-
"""
metric_fields.py -- build the metric-support version of every analysis window.

WHY THIS EXISTS
---------------
Every number in the empirical half of the manuscript is currently computed on
the geographic grid of ESA WorldCover: a cell is 1/12000 of a degree, and the
window is read as a 3072 x 3072 block of such cells.  On that grid a "30 m"
pixel (3 x 3 cells) is a rectangle whose ground size depends on latitude:

    north-south   3 * 110574 / 12000            = 27.64 m   (constant)
    east-west     3 * 111320 * cos(lat) / 12000 = 25.83 m at 21.9 N
                                                = 18.37 m at 48.7 N
                                                = 11.01 m at 66.7 N

so the aspect ratio runs from 1.07 to 2.51.  Two consequences follow, and both
touch quantities the manuscript reports:

  (i)  a block size m is not a physical area, so "m = 25" is not the same
       experiment in Lapland and in Yunnan;
  (ii) an empirical variogram binned by INTEGER lag mixes separation vectors
       whose ground lengths differ by up to a factor of 2.5, so the fitted
       range -- which Section 5.6 uses as an explanatory variable -- is not a
       distance at all.

The fix is not to reinterpret the old numbers but to rebuild the fields on a
support where one cell is one cell in both directions.  Each window is
reprojected into its own local azimuthal-equidistant frame at a 10 m grid, so
that

    1 cell  = 10 m x 10 m                 (the product's nominal resolution)
    1 pixel = 3 x 3 cells = 30 m x 30 m   (a square, the same everywhere)
    1 block = b x b pixels = 30b m square

and every lag, range and block size can be quoted in metres.

WHY AZIMUTHAL EQUIDISTANT PER WINDOW
------------------------------------
A window is at most 28 km across.  At that scale the only thing a projection
has to get right is distance, and aeqd centred on the window gives exact
radial distance from the centre and scale error below 1e-5 within 30 km.  A
single continental projection (UTM, Albers) would introduce a scale factor that
varies with the window's position inside its zone, which is precisely the kind
of latitude-dependent artefact this rebuild is meant to remove.

RESAMPLING CHOICE
-----------------
The class array is recoded to the binary tree indicator BEFORE resampling, and
nearest-neighbour is used.  Nearest is the right choice for a categorical
field: averaging class codes is meaningless, and averaging the indicator would
produce a 10 m field of fractions whose variance is not the variance of the
10 m map, which would silently change sigma^2 and therefore every bound.

VALIDATION
----------
A reprojection that silently distorts the data is worse than no reprojection,
because the numbers it produces look plausible.  Three checks are applied per
window and the run aborts if any fails: the tree fraction must not move by more
than 2 percentage points; the destination footprint must match the geodesic
extent to within 10 %; and the destination grid must be large enough for the
block sizes used downstream.  The first draft of this recipe dropped the sign
of western longitudes, which moved the two Iberian windows 15 degrees away
from their own projection centre and cost them 11 percentage points of tree
cover -- a failure that looked entirely plausible in the summary statistics.

Outputs
    data/metric/<tile>_<row>_<col>.npy     uint8 binary tree field, 10 m grid
    data/metric/_metric_fields.csv         one validated record per window
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rasterio.transform import from_origin
from rasterio.warp import Resampling, calculate_default_transform, reproject

ROOT = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
sys.stdout.reconfigure(encoding="utf-8")

from windows_v2 import (BIOME, CENTRE, PLACE, WINDOWS, _tile_origin,  # noqa: E402
                        label_of)

DATA = ROOT / "data"
OUTDIR = DATA / "metric"
META = OUTDIR / "_metric_fields.csv"

TREE = 10
CELLS_PER_DEG = 12000.0
WINDOW = 3072          # window side in geographic cells
TARGET_M = 10.0        # metric grid spacing, metres
PIXEL = 3              # cells per pixel side
MAX_BLOCK_PIXELS = 5   # largest block used downstream

# validation tolerances -- see the module docstring
MAX_MEAN_SHIFT = 0.02      # 2 percentage points of tree cover
MAX_AREA_ERROR = 0.10      # 10 % of the geodesic footprint
MIN_PIXELS_PER_SIDE = 300  # >= 9 km on the short side, so that the per-window
                           # lag window (half the short side) resolves ranges
                           # of a few km; see exp_metric.rmax_for_window


def geographic_binary(tile: str, row: int, col: int) -> np.ndarray:
    """The cached geographic window, recoded to the tree indicator."""
    arr = np.load(DATA / f"worldcover_{tile}_{row}_{col}_{WINDOW}.npy")
    return (arr == TREE).astype(np.uint8)


def window_extent(tile: str, row: int, col: int):
    """(west, south, east, north) of the geographic window, in degrees."""
    lat_tile, lon_tile = _tile_origin(tile)          # signed SW corner
    lat0 = lat_tile + 3.0 - row / CELLS_PER_DEG      # north edge
    lon0 = lon_tile + col / CELLS_PER_DEG            # west edge
    south = lat0 - WINDOW / CELLS_PER_DEG
    east = lon0 + WINDOW / CELLS_PER_DEG
    return lon0, south, east, lat0


def geodesic_extent_m(lat_c: float) -> tuple[float, float]:
    """(north-south, east-west) extent of the window in metres."""
    ns = WINDOW * 110574.0 / CELLS_PER_DEG
    ew = WINDOW * 111320.0 * np.cos(np.deg2rad(lat_c)) / CELLS_PER_DEG
    return ns, ew


def reproject_binary(arr: np.ndarray, tile: str, row: int, col: int):
    """Reproject one binary window onto its local 10 m metric grid."""
    lab = label_of(tile, row, col)
    lat_c, lon_c = CENTRE[lab]
    west, south, east, north = window_extent(tile, row, col)

    src_tf = from_origin(west, north, 1.0 / CELLS_PER_DEG, 1.0 / CELLS_PER_DEG)
    dst_crs = (f"+proj=aeqd +lat_0={lat_c:.6f} +lon_0={lon_c:.6f} "
               "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    dst_tf, dw, dh = calculate_default_transform(
        "EPSG:4326", dst_crs, WINDOW, WINDOW, west, south, east, north,
        resolution=TARGET_M)

    dst = np.zeros((dh, dw), dtype=np.uint8)
    reproject(source=arr, destination=dst,
              src_transform=src_tf, src_crs="EPSG:4326",
              dst_transform=dst_tf, dst_crs=dst_crs,
              resampling=Resampling.nearest)
    return dst, (dst_tf, dst_crs), (dh, dw)


def validate(lab: str, tile: str, row: int, col: int,
             geo: np.ndarray, met: np.ndarray, lat_c: float) -> dict:
    """Abort on any silent distortion; return the measured diagnostics."""
    ns_m, ew_m = geodesic_extent_m(lat_c)
    dh, dw = met.shape
    got = (dh * TARGET_M) * (dw * TARGET_M)
    exp = ns_m * ew_m
    area_err = got / exp - 1.0
    assert abs(area_err) < MAX_AREA_ERROR, (
        f"{lab}: destination footprint {got/1e6:.1f} Mm2 vs geodesic "
        f"{exp/1e6:.1f} Mm2 ({area_err:+.1%}) -- the transform is wrong")

    shift = float(met.mean()) - float(geo.mean())
    assert abs(shift) < MAX_MEAN_SHIFT, (
        f"{lab}: tree fraction moved {shift:+.4f} under reprojection "
        f"({geo.mean():.4f} -> {met.mean():.4f}) -- check the CRS, the sign of "
        f"the longitude, and that the binary recode precedes the resample")

    nx_pix, ny_pix = dw // PIXEL, dh // PIXEL
    assert min(nx_pix, ny_pix) >= MIN_PIXELS_PER_SIDE, (
        f"{lab}: only {ny_pix} x {nx_pix} pixels on the metric grid, below the "
        f"{MIN_PIXELS_PER_SIDE} needed for a resolvable variogram on the short "
        f"side -- widen the geographic window for this site")

    # geographic pixel geometry, for the record: this is the defect being fixed
    ns_pix = PIXEL * 110574.0 / CELLS_PER_DEG
    ew_pix = PIXEL * 111320.0 * np.cos(np.deg2rad(lat_c)) / CELLS_PER_DEG
    return dict(
        window=lab, tile=tile, row=row, col=col,
        biome=BIOME[lab], place=PLACE[lab],
        lat=round(lat_c, 4), lon=round(CENTRE[lab][1], 4),
        geo_shape=f"{WINDOW}x{WINDOW}",
        met_shape=f"{dh}x{dw}",
        met_pixels=f"{ny_pix}x{nx_pix}",
        ns_extent_m=round(ns_m, 1), ew_extent_m=round(ew_m, 1),
        area_err=round(area_err, 4),
        geo_pixel_ns_m=round(ns_pix, 3), geo_pixel_ew_m=round(ew_pix, 3),
        geo_pixel_aspect=round(ns_pix / ew_pix, 4),
        met_pixel_m=PIXEL * TARGET_M,
        tree_geo=round(float(geo.mean()), 5),
        tree_met=round(float(met.mean()), 5),
        tree_shift_pp=round(100.0 * shift, 4),
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    print(f"{'window':18s}{'lat':>7}{'aspect':>8}{'metric grid':>14}"
          f"{'pixels':>12}{'tree geo':>10}{'tree met':>10}{'shift pp':>10}")
    print("-" * 99)

    for tile, row, col, biome, group, place in WINDOWS:
        lab = label_of(tile, row, col)
        t0 = time.time()
        out = OUTDIR / f"{tile}_{row}_{col}.npy"

        geo = geographic_binary(tile, row, col)
        if out.exists():
            met = np.load(out)
        else:
            met, _tf, _shape = reproject_binary(geo, tile, row, col)
            np.save(out, met)

        rec = validate(lab, tile, row, col, geo, met, CENTRE[lab][0])
        rec["group"] = group
        rec["cache"] = out.name
        rows.append(rec)

        print(f"{lab:18s}{rec['lat']:>7.2f}{rec['geo_pixel_aspect']:>8.3f}"
              f"{rec['met_shape']:>14}{rec['met_pixels']:>12}"
              f"{rec['tree_geo']:>10.5f}{rec['tree_met']:>10.5f}"
              f"{rec['tree_shift_pp']:>+10.4f}"
              f"   ({time.time() - t0:4.1f} s)")

    df = pd.DataFrame(rows)
    df.to_csv(META, index=False)
    print("-" * 99)
    print(f"{len(df)} windows reprojected to {TARGET_M:.0f} m "
          f"azimuthal-equidistant grids")
    print(f"tree fraction : max |shift| "
          f"{df['tree_shift_pp'].abs().max():.4f} percentage points")
    print(f"footprint     : max |area error| {df['area_err'].abs().max():.4%}")
    print(f"geographic pixel aspect (N-S / E-W): "
          f"{df['geo_pixel_aspect'].min():.3f} to "
          f"{df['geo_pixel_aspect'].max():.3f}  ->  metric 1.000")
    print(f"metric pixels per window: "
          f"{min(int(r.split('x')[0]) * int(r.split('x')[1]) for r in df['met_pixels']):,}"
          f" to "
          f"{max(int(r.split('x')[0]) * int(r.split('x')[1]) for r in df['met_pixels']):,}")
    print(f"\nCSV: {META}")


if __name__ == "__main__":
    main()
