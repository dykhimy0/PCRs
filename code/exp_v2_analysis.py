"""
exp_v2_analysis.py -- stage 2 of the revision-2 experiment chain.

Three analyses over the 11 windows:

1. geostatistical comparator and decision  (results_v2/baseline_*.csv)
   Reuses the variogram machinery of exp_baseline.py unchanged, but loads each
   window by (tile, row, col) instead of by the fragile index of a sorted glob.

2. training-block estimator                 (results_v2/estimator_shrinkage.csv)
   Reuses the plug-in estimator of exp_estimator.py.

3. NEW: what one stored dependence scalar can and cannot do
   (results_v2/rho_scalar.csv)

   Section 6.7 claims that storing a single exchangeable correlation rho makes
   the block variance point-identified.  That claim is only as good as the
   scalar, so this stage measures the proposition rather than asserting it.

   Under exchangeable coupling with equal weights the block variance is affine
   in rho:

       V(rho) = sigma^2 * [1 + (m - 1) rho] / m,
       DEFF   = m V / sigma^2 = 1 + (m - 1) rho.

   For every window and every m we therefore compute

     rho_exact(m)   the scalar that reproduces the enumerated truth exactly,
                    rho_exact = (m * V_true / sigma^2 - 1) / (m - 1);
     rho_lag1       the lag-1 autocorrelation of the stored pixel field, which
                    is the scalar a producer would most naturally write down;
     V_from_lag1    what the exchangeable formula returns for that scalar;
     spread_exact   the spread of rho_exact across m within a window.

   Two diagnostics follow directly.
     * If exchangeability held, rho_exact would not depend on m.  Its spread is
       therefore a direct measure of how far the fields are from exchangeable,
       and it bounds the accuracy of ANY single stored scalar.
     * The error of the natural scalar is turned into a sampling-cost error
       through DEFF, which is what a user actually pays.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
sys.stdout.reconfigure(encoding="utf-8")

from exp_baseline import (block_distances, fit_exp, fit_nested,  # noqa: E402
                          radial_variogram, var_block_empirical,
                          var_block_model, var_block_model_nested,
                          nested_model, exp_model,
                          RMAX_FIT, RMAX_TRUNC, TARGET_PRECISION)
from pcp_core import block_composition_map  # noqa: E402
from windows_v2 import LABELS, WINDOWS, label_of  # noqa: E402

DATA = ROOT / "data"
RES = ROOT / "results_v2"
RES.mkdir(exist_ok=True)
TREE = 10
PIXEL = 3


def pixel_field(tile: str, row: int, col: int) -> tuple[np.ndarray, np.ndarray]:
    """(30 m tree-fraction field, 10 m binary field) for one window."""
    arr = np.load(DATA / f"worldcover_{tile}_{row}_{col}_3072.npy")
    z = (arr == TREE).astype(np.int8)
    pix = block_composition_map(z, PIXEL).astype(np.float64)
    return pix, z


def block_means(pix: np.ndarray, b: int) -> np.ndarray:
    nb = pix.shape[0] // b
    blk = pix[: nb * b, : nb * b].reshape(nb, b, nb, b).mean(axis=(1, 3))
    return blk.ravel()


def fit_nested_bounded(rs, gs, sig2, rmax):
    """Nested two-structure fit with the long range held inside the data.

    The unbounded fit of exp_baseline can place the second structure far beyond
    the lag window (one of the new windows drives it to a2 ~ 3000 px while the
    data stop at 300 px).  That parameter is then unidentifiable: R^2 is 0.998
    because the extra term is flat over the observed lags, yet the block
    variance it implies is out by a factor of four.  Constraining the long
    range to the fitted window is the standard identifiability fix, and the two
    fits are reported side by side so that the failure is on the record.
    """
    from scipy.optimize import curve_fit

    p0 = [0.0, sig2 * 0.6, 3.0, sig2 * 0.4, min(100.0, rmax)]
    bounds = ([0.0, 0.0, 1e-3, 0.0, 1e-3], [sig2, sig2 * 3.0, 1e3, sig2 * 3.0, rmax])
    popt, _ = curve_fit(nested_model, rs, gs, p0=p0, bounds=bounds, maxfev=200000)
    pred = nested_model(rs, *popt)
    ss_res = float(np.sum((gs - pred) ** 2))
    ss_tot = float(np.sum((gs - gs.mean()) ** 2))
    return tuple(float(v) for v in popt), 1.0 - ss_res / ss_tot


# ==========================================================================
# 1. geostatistical comparator
# ==========================================================================
def run_baseline(blocks: pd.DataFrame) -> None:
    rows, curves, fits, dec = [], [], [], []
    for tile, row, col, biome, group, _p in WINDOWS:
        label = label_of(tile, row, col)
        t0 = time.time()
        pix, _z = pixel_field(tile, row, col)
        sig2 = float(pix.var())

        rs_t, gs_t, _ = radial_variogram(pix, rmax=RMAX_TRUNC)
        (c0t, c1t, at), r2t = fit_exp(rs_t, gs_t, sig2)
        rs, gs, _ = radial_variogram(pix, rmax=RMAX_FIT)
        (c0, c1, a), r2 = fit_exp(rs, gs, sig2)
        (d0, d1, da1, d2, da2), r2n = fit_nested(rs, gs, sig2)
        (b0, b1, ba1, b2, ba2), r2b = fit_nested_bounded(rs, gs, sig2, RMAX_FIT)

        for wl, cc, rl in (("truncated", (c0t, c1t, at), r2t),
                           ("wide", (c0, c1, a), r2),
                           ("nested", (b0, b1, ba1, b2, ba2), r2b),
                           ("nested_unbounded", (d0, d1, da1, d2, da2), r2n)):
            if wl == "nested_unbounded":
                continue
            fits.append(dict(
                window=label, biome=biome, group=group, sigma2=sig2, model=wl,
                rmax=RMAX_TRUNC if wl == "truncated" else RMAX_FIT,
                a1=cc[2], a2=cc[4] if wl == "nested" else np.nan,
                sill=(cc[0] + cc[1]) if wl != "nested" else (cc[0] + cc[1] + cc[3]),
                r2=rl,
                gamma_at_trunc_over_sill=float(np.interp(RMAX_TRUNC, rs, gs) / sig2)))
        fits.append(dict(
            window=label, biome=biome, group=group, sigma2=sig2,
            model="nested_unbounded", rmax=RMAX_FIT, a1=da1, a2=da2,
            sill=d0 + d1 + d2, r2=r2n,
            gamma_at_trunc_over_sill=float(np.interp(RMAX_TRUNC, rs, gs) / sig2)))
        print(f"  baseline {label:18s} nested a2 bounded {ba2:7.1f} / unbounded {da2:8.1f}")

        for lag, gam in zip(rs, gs):
            curves.append(dict(window=label, biome=biome, lag=float(lag),
                               gamma=float(gam), sigma2=sig2,
                               within_trunc=bool(lag <= RMAX_TRUNC)))

        sub = blocks[blocks["window"] == label].sort_values("m")
        for _, r in sub.iterrows():
            m = int(r["m"])
            ms = int(round(np.sqrt(m)))
            if ms ** 2 != m:
                continue
            v_tr = var_block_model(ms, c0t, c1t, at)
            v_wi = var_block_model(ms, c0, c1, a)
            v_ne = var_block_model_nested(ms, b0, b1, ba1, b2, ba2)
            v_em = var_block_empirical(ms, rs, gs, sig2)
            truth = float(r["truth_block_var"])
            lo = max(0.0, float(r["v_lower"]))
            hi = float(r["v_upper"])
            rows.append(dict(
                window=label, biome=biome, group=group, m=m, m_side=ms, sigma2=sig2,
                truth=truth, independence=float(r["v_independent"]),
                identified_lower=lo, identified_upper=hi,
                model_truncated=v_tr, model_wide=v_wi, model_nested=v_ne,
                empirical_variogram=v_em,
                truth_inside=bool(lo - 1e-9 <= truth <= hi + 1e-9),
                err_truncated_pct=100.0 * (v_tr - truth) / truth,
                err_wide_pct=100.0 * (v_wi - truth) / truth,
                err_nested_pct=100.0 * (v_ne - truth) / truth,
                err_empirical_pct=100.0 * (v_em - truth) / truth,
                err_independence_pct=100.0 * (float(r["v_independent"]) - truth) / truth,
            ))
            deff = {"truth": m * truth / sig2, "model_wide": m * v_wi / sig2,
                    "model_truncated": m * v_tr / sig2, "empirical": m * v_em / sig2,
                    "independence": 1.0, "upper": m * hi / sig2, "lower": m * lo / sig2}
            rr = dict(window=label, biome=biome, group=group, m=m)
            for k, v in deff.items():
                rr[f"deff_{k}"] = v
                rr[f"n_{k}"] = sig2 * v / TARGET_PRECISION ** 2
            dec.append(rr)
        print(f"  baseline {label:18s} range_trunc {at:7.1f}  range_wide {a:7.1f}  "
              f"({time.time() - t0:.1f} s)")

    pd.DataFrame(rows).to_csv(RES / "baseline_variogram.csv", index=False)
    pd.DataFrame(curves).to_csv(RES / "baseline_variogram_curves.csv", index=False)
    pd.DataFrame(fits).to_csv(RES / "baseline_fits.csv", index=False)
    pd.DataFrame(dec).to_csv(RES / "baseline_decision.csv", index=False)


# ==========================================================================
# 2. training-block estimator
# ==========================================================================
def run_estimator(blocks: pd.DataFrame, reps: int = 400,
                  seed: int = 20261013) -> None:
    rng = np.random.default_rng(seed)
    band_w = {(r.window, int(r.m)): float(r.v_upper - r.v_lower)
              for r in blocks.itertuples()}
    rows = []
    for tile, row, col, biome, group, _p in WINDOWS:
        label = label_of(tile, row, col)
        pix, _z = pixel_field(tile, row, col)
        sig2 = float(pix.var())
        for b in (2, 3, 4, 5):
            m = b * b
            means = block_means(pix, b)
            N = means.size
            truth = float(means.var())
            rho_oracle = (m * truth / sig2 - 1.0) / (m - 1.0)
            for K in (4, 8, 16, 32, 64, 128, 256):
                if K >= N:
                    continue
                vhat = np.empty(reps)
                for i in range(reps):
                    s = means[rng.choice(N, size=K, replace=False)]
                    vhat[i] = ((s - s.mean()) ** 2).sum() / (K - 1.0)
                rhos = (m * vhat / sig2 - 1.0) / (m - 1.0)
                rel = np.abs(vhat - truth) / truth * 100.0
                rows.append(dict(
                    window=label, biome=biome, group=group, block_pixels=b, m=m,
                    K=K, N_blocks=N, sigma2=sig2, truth_block_var=truth,
                    rho_oracle=rho_oracle,
                    v_hat_median=float(np.median(vhat)),
                    v_hat_iqr=float(np.percentile(vhat, 75) - np.percentile(vhat, 25)),
                    bias_median_pct=float(np.median(vhat - truth) / truth * 100.0),
                    err_median_pct=float(np.median(rel)),
                    err_p90_pct=float(np.percentile(rel, 90)),
                    rho_hat_median=float(np.median(rhos)),
                    rho_hat_frac_gt1=float((rhos > 1.0).mean()),
                    levelA_width=band_w.get((label, m), float("nan")),
                ))
        print(f"  estimator {label}")
    df = pd.DataFrame(rows)
    df["shrinkage_iqr"] = df["levelA_width"] / df["v_hat_iqr"]
    df["rho_hat_frac_outside"] = (
        (df["rho_hat_median"] > 1.0)
        | (df["rho_hat_median"] < -1.0 / (df["m"] - 1.0))).astype(float)
    df.to_csv(RES / "estimator_shrinkage.csv", index=False)


# ==========================================================================
# 3. what one stored dependence scalar can and cannot do
# ==========================================================================
def run_rho_scalar(blocks: pd.DataFrame) -> None:
    rows = []
    for tile, row, col, biome, group, _p in WINDOWS:
        label = label_of(tile, row, col)
        pix, _z = pixel_field(tile, row, col)
        a = pix - pix.mean()
        rho_lag1 = float((a[:, :-1] * a[:, 1:]).sum() / (a * a).sum())
        sig2 = float(pix.var())
        for b in (2, 3, 4, 5):
            m = b * b
            means = block_means(pix, b)
            truth = float(means.var())
            rho_exact = (m * truth / sig2 - 1.0) / (m - 1.0)
            v_lag1 = sig2 * (1.0 + (m - 1.0) * rho_lag1) / m
            deff_true = m * truth / sig2
            deff_lag1 = 1.0 + (m - 1.0) * rho_lag1
            rho_lo = -1.0 / (m - 1.0)
            rows.append(dict(
                window=label, biome=biome, group=group, block_pixels=b, m=m,
                sigma2=sig2, truth_block_var=truth,
                rho_exact=rho_exact, rho_lag1=rho_lag1,
                rho_exact_in_psd=bool(rho_lo - 1e-12 <= rho_exact <= 1.0 + 1e-12),
                rho_lag1_in_psd=bool(rho_lo - 1e-12 <= rho_lag1 <= 1.0 + 1e-12),
                v_from_lag1=v_lag1,
                err_lag1_pct=100.0 * (v_lag1 - truth) / truth,
                deff_true=deff_true, deff_from_lag1=deff_lag1,
                cost_ratio=deff_lag1 / deff_true,
                # how far the natural scalar falls below the exact one, relative
                # to the level-A band width: 1.0 means "no better than no scalar"
                band_frac_explained=float(
                    (v_lag1 - 0.0) / max(truth, 1e-12) if truth > 0 else np.nan),
            ))
    df = pd.DataFrame(rows)
    df["rho_exact_spread"] = df.groupby("window")["rho_exact"].transform(
        lambda s: s.max() - s.min())
    df["rho_exact_ratio"] = df.groupby("window")["rho_exact"].transform(
        lambda s: s.iloc[0] / s.iloc[-1] if s.iloc[-1] != 0 else np.nan)
    df.to_csv(RES / "rho_scalar.csv", index=False)


def main() -> None:
    blocks = pd.read_csv(RES / "real_blocks.csv")
    print("=== 1. geostatistical comparator ===")
    run_baseline(blocks)
    print("=== 2. training-block estimator ===")
    run_estimator(blocks)
    print("=== 3. one stored dependence scalar ===")
    run_rho_scalar(blocks)
    print("\ndone")


if __name__ == "__main__":
    main()
