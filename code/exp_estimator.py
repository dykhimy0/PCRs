"""
exp_estimator.py -- proof of concept: shrinking the identified interval from a
handful of training blocks.

The marginals (level A) are stored and known exactly.  Given K training blocks
with observed composition, the natural plug-in estimate of the block variance
V(m) is the sample variance of the K block means,

    V_hat = S2 = (1/(K-1)) * sum_k ( mean_k - grand_mean )^2 ,

which is unbiased for the true block variance (within-block covariance plus
any large-scale trend).  The implied exchangeable dependence parameter is
reported for reference as  rho_hat = (m*S2/sigma2 - 1)/(m-1), but it is NOT
used to re-enter the variance: on a trended field rho_hat conflates within-block
correlation with large-scale structure, and can exceed 1, so it is reported raw.

We report, over repeated random samples of K blocks, the median bias and the
90th-percentile relative error of V_hat against the enumerated block variance,
and the shrinkage factor: the level-A sharp-band width divided by the
interquartile range of V_hat.  This is a concept check only; no new algorithm.

Output
------
results/estimator_shrinkage.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pcp_core import block_composition_map  # noqa: E402
from windows import WINDOWS  # noqa: E402

DATA = ROOT / "data"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

TREE = 10
PIXEL = 3
BLOCK_PIXELS = (2, 3, 4, 5)          # m = 4, 9, 16, 25
K_GRID = (4, 8, 16, 32, 64, 128, 256)
R_REPS = 400
SEED = 20261013


def load_pixel_field(tile, row, col):
    arr = np.load(DATA / f"worldcover_{tile}_{row}_{col}_3072.npy")
    z = (arr == TREE).astype(np.int8)
    return block_composition_map(z, PIXEL).astype(np.float64)


def block_means(pix, b):
    nb = pix.shape[0] // b
    blk = pix[: nb * b, : nb * b].reshape(nb, b, nb, b).mean(axis=(1, 3))
    return blk.ravel()


def main():
    rng = np.random.default_rng(SEED)
    band = pd.read_csv(RES / "real_blocks.csv")
    band_w = {(r.window, r.m): float(r.v_upper - r.v_lower) for r in band.itertuples()}

    rows = []
    for tile, row, col, label, biome, group in WINDOWS:
        pix = load_pixel_field(tile, row, col)
        sigma2 = float(pix.var())
        for b in BLOCK_PIXELS:
            m = b * b
            means = block_means(pix, b)
            N = means.size
            truth_var = float(means.var())
            rho_oracle = (m * truth_var / sigma2 - 1.0) / (m - 1.0)
            for K in K_GRID:
                if K >= N:
                    continue
                vhats, rhos = [], []
                for _ in range(R_REPS):
                    s = means[rng.choice(N, size=K, replace=False)]
                    S2 = float(((s - s.mean()) ** 2).sum() / (K - 1.0))
                    vhats.append(S2)
                    rhos.append((m * S2 / sigma2 - 1.0) / (m - 1.0))
                vhats = np.asarray(vhats)
                rhos = np.asarray(rhos)
                rel = np.abs(vhats - truth_var) / truth_var * 100.0
                rows.append(dict(
                    window=label, biome=biome, group=group,
                    block_pixels=b, m=m, K=K, N_blocks=N,
                    sigma2=sigma2, truth_block_var=truth_var,
                    rho_oracle=rho_oracle,
                    v_hat_median=float(np.median(vhats)),
                    v_hat_iqr=float(np.percentile(vhats, 75) - np.percentile(vhats, 25)),
                    bias_median_pct=float(np.median(vhats - truth_var) / truth_var * 100.0),
                    err_median_pct=float(np.median(rel)),
                    err_p90_pct=float(np.percentile(rel, 90)),
                    rho_hat_median=float(np.median(rhos)),
                    rho_hat_frac_gt1=float((rhos > 1.0).mean()),
                    levelA_width=band_w.get((label, m), float("nan")),
                ))
    df = pd.DataFrame(rows)
    df["shrinkage_iqr"] = df["levelA_width"] / df["v_hat_iqr"]
    df.to_csv(RES / "estimator_shrinkage.csv", index=False)

    cols = ["window", "biome", "m", "K", "truth_block_var", "rho_oracle",
            "v_hat_median", "v_hat_iqr", "bias_median_pct", "err_p90_pct",
            "rho_hat_median", "levelA_width", "shrinkage_iqr"]
    print(df[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))


if __name__ == "__main__":
    main()
