"""
window_metadata.py -- the descriptive table that every window must carry.

For each window, from the cached 10 m array, record the quantities needed to
argue that the six windows are genuinely different landscapes rather than six
draws from one landscape:

  composition   class shares at 10 m, and the two diversity indices
  contrast      tree-fraction mean and variance at 30 m
  structure     lag-1 autocorrelation of the 30 m tree fraction
  geometry      the ground size of a 3 x 3 cell block, which depends on
                latitude because the map is in EPSG:4326 and is therefore not
                conformal: a degree of longitude shrinks with the cosine of the
                latitude, a degree of latitude does not.

Writes results/window_metadata.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results"

PIXEL = 3
LAT0, LON0 = 111_320.0, 111_320.0     # metres per degree, on the sphere
CELL_DEG = 3.0 / 36000.0              # 10 m cell in degrees of the map

WC = {10: "tree", 20: "shrub", 30: "grass", 40: "crop", 50: "built",
      60: "bare", 70: "snow", 80: "water", 90: "wetland", 95: "mangrove",
      100: "moss"}


def main() -> None:
    marg = pd.read_csv(RES / "real_marginal.csv")
    rows = []
    for r in marg.itertuples():
        f = DATA / f"worldcover_{r.tile}_{int(r.row)}_{int(r.col)}_3072.npy"
        arr = np.load(f)
        vals, cnts = np.unique(arr, return_counts=True)
        shares = {int(v): float(c) / arr.size for v, c in zip(vals, cnts)}

        k = arr.shape[0] // PIXEL
        blk = arr[:k * PIXEL, :k * PIXEL].reshape(k, PIXEL, k, PIXEL)
        tree = (blk == 10).mean(axis=(1, 3))
        major = blk[:, 0, :, 0].reshape(-1)
        _, mc = np.unique(major, return_counts=True)
        p = mc / mc.sum()

        a = tree - tree.mean()
        lag1 = float((a[:, :-1] * a[:, 1:]).sum() / (a * a).sum())

        lat = float(r.lat_centre)
        ew = PIXEL * CELL_DEG * LON0 * np.cos(np.deg2rad(lat))
        ns = PIXEL * CELL_DEG * LAT0
        cell_m = np.sqrt(ew * ns)

        shares_named = {WC.get(c, str(c)): round(v, 4) for c, v in shares.items()}
        rows.append(dict(
            window=r.window, biome=r.biome, group=r.group,
            place=r.place, lat_centre=lat,
            n_classes_ge_1pct=sum(v >= 0.01 for v in shares.values()),
            dominant=WC.get(max(shares, key=shares.get), "?"),
            dominant_share=round(max(shares.values()), 4),
            tree_share_10m=shares.get(10, 0.0),
            tree_mean_30m=float(tree.mean()), sigma2_30m=float(tree.var()),
            shannon_30m=float(-(p * np.log(p)).sum()),
            simpson_30m=float(1.0 - (p ** 2).sum()),
            lag1_tree_30m=lag1,
            block_ew_m=ew, block_ns_m=ns, block_diag_m=float(np.hypot(ew, ns)),
            anisotropy_ns_over_ew=ns / ew,
            class_shares=str(shares_named),
        ))
    M = pd.DataFrame(rows)
    M.to_csv(RES / "window_metadata.csv", index=False)

    show = ["window", "biome", "place", "lat_centre", "n_classes_ge_1pct",
            "dominant", "dominant_share", "tree_mean_30m", "sigma2_30m",
            "shannon_30m", "simpson_30m", "lag1_tree_30m",
            "block_ew_m", "block_ns_m", "anisotropy_ns_over_ew"]
    print(M[show].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("\nclass shares")
    for r in M.itertuples():
        print(f"  {r.window:20s} {r.class_shares}")
    print("\nwrote window_metadata.csv")


if __name__ == "__main__":
    main()
