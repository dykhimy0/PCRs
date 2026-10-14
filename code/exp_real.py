"""
exp_real.py -- the same change-of-support analysis on real land cover.

Reference data: ESA WorldCover v200 (2021), 10 m, EPSG:4326, read directly from
the public COG on S3 (no credentials required).  The 10 m map is treated as the
fine reference.  The raster under study is stored at the pixel level
(3 x 3 cells = 30 m) and the task is change of support to a coarser block.

The stored object is the law of the within-pixel class composition.  Here that
law is estimated *from the reference map itself*, i.e. the marginals are known
exactly -- the most favourable case for the raster.  Any residual uncertainty
in the block-level composition is therefore attributable purely to the missing
cross-pixel dependence.

Outputs
-------
data/worldcover_<tile>_<window>.npy
results/real_blocks.csv
results/real_exceedance.csv
results/real_marginal.csv
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pcp_core import (  # noqa: E402
    block_composition_map,
    comonotonic_variance,
    empirical_law,
    envelope,
    identified_set_bounds,
    independent_aggregate_law,
    independent_variance,
    marginal_moments,
    psd_variance_bounds,
)
from windows import LAT_CENTRE, PLACE, WINDOWS  # noqa: E402

DATA = ROOT / "data"
RES = ROOT / "results"
DATA.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

URL = ("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
       "ESA_WorldCover_10m_2021_v200_{tile}_Map.tif")

PIXEL = 3                                   # 3 x 3 cells = 30 m pixel
WIN = 3072                                  # window side in fine cells
TREE, CROP, GRASS, BUILT, BARE, WATER = 10, 40, 30, 50, 60, 80

BLOCK_PIXELS = (2, 3, 4, 5)


def fetch(tile, row, col, side=WIN):
    """Read a window of the reference map, cached under data/."""
    cache = DATA / f"worldcover_{tile}_{row}_{col}_{side}.npy"
    if cache.exists():
        return np.load(cache)
    import rasterio
    from rasterio.windows import Window
    with rasterio.open(URL.format(tile=tile)) as src:
        arr = src.read(1, window=Window(col, row, side, side))
    np.save(cache, arr)
    return arr


def binary(arr, cls=TREE):
    return (arr == cls).astype(np.int8)


def analyse(arr, label, cls=TREE):
    """Change-of-support analysis at the pixel level for one window."""
    z = binary(arr, cls)
    pix = block_composition_map(z, PIXEL)
    fine_values = np.arange(0, PIXEL ** 2 + 1) / PIXEL ** 2
    vals, probs = empirical_law(pix, values=fine_values)
    _, var_pix = marginal_moments(vals, probs)

    rows, trows = [], []
    for b in BLOCK_PIXELS:
        m = b * b
        blk = block_composition_map(pix, b)
        truth_var = float(blk.var())
        v_ind = independent_variance(vals, probs, m)
        v_plus = comonotonic_variance(vals, probs, m)
        lo_psd, _, _, rho_psd, _ = psd_variance_bounds(vals, probs, m)
        if m <= 4:
            lo_exact = identified_set_bounds(vals, probs, m, kind="variance")["lower"]
            method = "exact LP"
        else:
            lo_exact, method = np.nan, "PSD relaxation"
        v_minus = lo_exact if m <= 4 else lo_psd
        rows.append(dict(
            window=label, block_pixels=b, m=m, method=method,
            support_size=int(len(vals)),
            pixel_marginal_mean=float(probs @ vals),
            pixel_marginal_var=var_pix,
            truth_block_mean=float(blk.mean()),
            truth_block_var=truth_var,
            v_independent=v_ind,
            v_lower=v_minus, v_upper=v_plus,
            implied_reduction_truth=var_pix / truth_var if truth_var > 0 else np.inf,
            independent_understates_by=v_ind / truth_var if truth_var > 0 else np.nan,
            aci=(v_plus - v_minus) / v_plus,
            truth_inside=bool(v_minus - 1e-12 <= truth_var <= v_plus + 1e-12),
        ))
        # lag-1 correlation of the pixel composition (a single scalar summary)
        a = pix - pix.mean()
        rows[-1]["lag1_autocorr_pixel"] = float((a[:, :-1] * a[:, 1:]).sum() / (a * a).sum())

        for tau in (0.25, 0.5, 0.75):
            trows.append(dict(
                window=label, block_pixels=b, m=m, tau=tau,
                truth_tail=float((blk >= tau - 1e-12).mean()),
                independent_tail=float(
                    independent_aggregate_law(vals, probs, m)[1][
                        independent_aggregate_law(vals, probs, m)[0] >= tau - 1e-12].sum()),
                **({"lower_bound": np.nan, "upper_bound": np.nan}
                   if m > 4 else {
                       "lower_bound": envelope(vals, probs, m, [tau])[0]["lower_exact"],
                       "upper_bound": envelope(vals, probs, m, [tau])[0]["upper_exact"]}),
            ))
    return pd.DataFrame(rows), pd.DataFrame(trows), (vals, probs, pix, z)


def main():
    blocks, tails, marg = [], [], []
    for tile, row, col, label, biome, group in WINDOWS:
        arr = fetch(tile, row, col)
        df, tdf, (vals, probs, pix, z) = analyse(arr, label)
        for frame in (df, tdf):
            frame.insert(1, "biome", biome)
            frame.insert(2, "group", group)
        blocks.append(df)
        tails.append(tdf)
        classes = {int(k): int(v) for k, v in zip(*np.unique(arr, return_counts=True))}
        marg.append(dict(window=label, tile=tile, row=row, col=col,
                         biome=biome, group=group,
                         place=PLACE.get(label, ""),
                         lat_centre=LAT_CENTRE.get(label, float("nan")),
                         classes=str(classes),
                         tree_share=float((arr == TREE).mean()),
                         crop_share=float((arr == CROP).mean()),
                         built_share=float((arr == BUILT).mean()),
                         **{f"marg_p_k{int(v * PIXEL ** 2)}": float(p)
                            for v, p in zip(vals, probs)}))
    B = pd.concat(blocks, ignore_index=True)
    T = pd.concat(tails, ignore_index=True)
    M = pd.DataFrame(marg)
    B.to_csv(RES / "real_blocks.csv", index=False)
    T.to_csv(RES / "real_exceedance.csv", index=False)
    M.to_csv(RES / "real_marginal.csv", index=False)

    print(M[["window", "biome", "group", "place", "lat_centre",
             "tree_share", "crop_share", "built_share"]].to_string(index=False))
    print()
    print(B[["window", "biome", "block_pixels", "m", "method", "pixel_marginal_var",
             "truth_block_var", "v_independent", "v_lower", "v_upper",
             "implied_reduction_truth", "independent_understates_by",
             "aci", "truth_inside"]].to_string(
        index=False, float_format=lambda x: f"{x:.5g}"))
    print()
    print(T.to_string(index=False, float_format=lambda x: f"{x:.5g}"))


if __name__ == "__main__":
    main()
