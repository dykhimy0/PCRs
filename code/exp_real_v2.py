"""
exp_real_v2.py -- change-of-support analysis on 11 real ESA WorldCover windows.

Stage 1 of the revision-2 experiment chain.  Reuses the analysis code of
exp_real.py unchanged (so the numbers are computed by exactly the same
functions as in revision 1) and only widens the window table and redirects the
output to results_v2/.

Outputs
    results_v2/real_blocks.csv
    results_v2/real_exceedance.csv
    results_v2/real_marginal.csv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.stdout.reconfigure(encoding="utf-8")

from exp_real import TREE, analyse, fetch          # noqa: E402  (analysis reused)
from windows_v2 import CENTRE, PLACE, WINDOWS, label_of  # noqa: E402

RES = ROOT / "results_v2"
RES.mkdir(exist_ok=True)


def main() -> None:
    blocks, tails, marg = [], [], []
    for tile, row, col, biome, group, _place in WINDOWS:
        label = label_of(tile, row, col)
        t0 = time.time()
        arr = fetch(tile, row, col)
        df, tdf, _extra = analyse(arr, label)
        for frame in (df, tdf):
            frame.insert(1, "biome", biome)
            frame.insert(2, "group", group)
        blocks.append(df)
        tails.append(tdf)
        classes = {int(k): int(v) for k, v in zip(*np.unique(arr, return_counts=True))}
        lat, lon = CENTRE[label]
        marg.append(dict(
            window=label, tile=tile, row=row, col=col, biome=biome, group=group,
            place=PLACE[label], lat_centre=lat, lon_centre=lon,
            classes=str(classes),
            tree_share=float((arr == TREE).mean()),
            crop_share=float((arr == 40).mean()),
            grass_share=float((arr == 30).mean()),
            water_share=float((arr == 80).mean()),
            built_share=float((arr == 50).mean()),
            bare_share=float((arr == 60).mean()),
            wetland_share=float((arr == 90).mean()),
            moss_share=float((arr == 100).mean()),
        ))
        print(f"  {label:18s} {biome:15s} tree {marg[-1]['tree_share']:.3f} "
              f"({time.time() - t0:.1f} s)")

    B = pd.concat(blocks, ignore_index=True)
    T = pd.concat(tails, ignore_index=True)
    M = pd.DataFrame(marg)
    B.to_csv(RES / "real_blocks.csv", index=False)
    T.to_csv(RES / "real_exceedance.csv", index=False)
    M.to_csv(RES / "real_marginal.csv", index=False)

    print("\n" + M[["window", "biome", "place", "lat_centre", "lon_centre",
                     "tree_share"]].to_string(index=False))
    print("\n" + B[["window", "biome", "m", "method", "pixel_marginal_var",
                    "truth_block_var", "v_independent", "v_lower", "v_upper",
                    "independent_understates_by", "aci", "truth_inside"]].to_string(
        index=False, float_format=lambda x: f"{x:.5g}"))
    print(f"\nwindows {B['window'].nunique()}  biomes {B['biome'].nunique()}  "
          f"cases {len(B)}")
    print(f"independent understates by {B['independent_understates_by'].min():.2f}"
          f"-{B['independent_understates_by'].max():.2f}")
    print(f"truth inside the identified band in all cases: {bool(B['truth_inside'].all())}")


if __name__ == "__main__":
    main()
