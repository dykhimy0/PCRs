# -*- coding: utf-8 -*-
"""
exp_metric.py -- recompute the whole empirical half on metric support.

This is the recomputation that metric_fields.py exists to enable.  It re-runs,
on the 10 m azimuthal-equidistant grids, the four analyses whose numbers the
manuscript quotes:

  A. change of support          V_true, the identified band, ACI, exceedance
  B. the geostatistical route   variogram in METRES, fitted range in METRES,
                                model block variance, sampling cost
  C. one stored scalar          rho_exact(m) and its support dependence
  D. the training-block estimator

and it adds one analysis that is only meaningful once the grid is metric:

  E. anisotropy                 gamma along rows (E-W) against gamma along
                                columns (N-S), each in metres

E is not decoration.  The whole reason the geographic numbers cannot be quoted
as distances is that an integer lag bin on that grid mixes separation vectors
of different ground length.  Whether that mixture matters is an empirical
question about the fields, not an argument about the grid: if the two
directional variograms coincide, an omnidirectional read-off is safe and the
old numbers survive reinterpretation; if they differ, the mixture was doing
real damage and the recomputation is mandatory.  Measuring it converts the
reproducibility finding from an assertion into a result.

LAG WINDOW
----------
On the geographic grid the two lag windows were fixed at 60 and 300 pixels,
which is a fixed fraction of the 1024-pixel side.  On the metric grid the same
fraction of the side is no longer the same number of metres, because the
windows are genuinely different shapes once latitude is respected -- Finnish
Lapland is 28 km north-south and 11 km east-west.  A fixed 9 km window would
therefore run past the short side of the high-latitude windows, where the
omnidirectional read-off degenerates into a purely north-south one.

The lag window is consequently defined per window as a fraction of that
window's own short side, capped at 9 km:

    rmax_fit   = min(9000 m, 0.50 * short side in metres)
    rmax_trunc = 0.20 * rmax_fit          (preserves the 60 : 300 = 1 : 5 ratio)

Half the short side is the standard upper limit for an omnidirectional
variogram; one fifth of it is the naive default the manuscript shows to be
biased.  Both are recorded per window so the reader can see the window.

UNITS
-----
Every distance in the outputs is in metres.  Ranges are reported both in
metres and in pixel units (range_m / 30) so that the old numbers can be lined
up against the new ones.

Outputs (results_metric/)
    metric_blocks.csv             A: change of support, one row per window x m
    metric_exceedance.csv         A: P(block >= tau) and its band
    metric_fits.csv               B: fitted parameters per window and model
    metric_variogram.csv          B: block variance vs truth, per window x m
    metric_variogram_curves.csv   B: the binned empirical variograms
    metric_decision.csv           B: sample sizes for a +/- 1 pp target
    metric_rho_scalar.csv         C: rho_exact(m), rho_lag1, the cost ratio
    metric_estimator.csv          D: training-block estimator
    metric_anisotropy.csv         E: directional variograms and the range ratio
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
sys.stdout.reconfigure(encoding="utf-8")

from pcp_core import (block_composition_map, comonotonic_variance,  # noqa: E402
                      empirical_law, envelope, identified_set_bounds,
                      independent_aggregate_law, independent_variance,
                      marginal_moments, psd_variance_bounds)
from windows_v2 import BIOME, CENTRE, PLACE, WINDOWS, label_of  # noqa: E402

METRIC = ROOT / "data" / "metric"
RES = ROOT / "results_metric"
RES.mkdir(exist_ok=True)

TREE = 10
PIXEL = 3              # 10 m cells per pixel side
PIXEL_M = 30.0         # pixel side, metres
BLOCK_PIXELS = (2, 3, 4, 5)
TAUS = (0.25, 0.5, 0.75)
TARGET_PRECISION = 0.01

RMAX_FIT_CAP_M = 9000.0
RMAX_SHORT_FRAC = 0.50
RMAX_TRUNC_FRAC = 0.20
MIN_PAIRS = 8

EXACT_M = 4            # exact LP is feasible up to m = 4


# --------------------------------------------------------------------------
# fields
# --------------------------------------------------------------------------
def metric_binary(tile: str, row: int, col: int) -> np.ndarray:
    return np.load(METRIC / f"{tile}_{row}_{col}.npy")


def pixel_field(binary: np.ndarray) -> np.ndarray:
    """30 m tree-fraction field on the metric grid."""
    return block_composition_map(binary, PIXEL).astype(np.float64)


def block_means(pix: np.ndarray, b: int) -> np.ndarray:
    nb_y, nb_x = pix.shape[0] // b, pix.shape[1] // b
    blk = pix[: nb_y * b, : nb_x * b].reshape(nb_y, b, nb_x, b).mean(axis=(1, 3))
    return blk.ravel()


# --------------------------------------------------------------------------
# variograms in metres
# --------------------------------------------------------------------------
def radial_variogram(z: np.ndarray, rmax_m: float, cell_m: float = PIXEL_M,
                     min_pairs: int = MIN_PAIRS):
    """Omnidirectional gamma(h) from the 2-D FFT autocovariance, h in metres.

    Rectangular input is handled explicitly: the geographic version of this
    function indexed a square lag grid against a possibly non-square
    autocovariance array, which fails outright once a window is not square --
    and on the metric grid no window is square.
    """
    ny, nx = z.shape
    zc = z - z.mean()
    F = np.fft.fft2(zc)
    ac = np.real(np.fft.ifft2(F * np.conj(F))) / zc.size
    ac = np.fft.fftshift(ac)
    oy = (np.arange(ny) - ny // 2)[:, None] * cell_m
    ox = (np.arange(nx) - nx // 2)[None, :] * cell_m
    r = np.sqrt(oy ** 2 + ox ** 2)
    bins = np.floor(r / cell_m).astype(int)
    sig2 = float(zc.var())
    rs, gs, ns = [], [], []
    for k in range(int(rmax_m // cell_m) + 1):
        sel = bins == k
        if sel.sum() < min_pairs:
            continue
        rs.append(k * cell_m)
        gs.append(sig2 - float(ac[sel].mean()))
        ns.append(int(sel.sum()))
    return np.array(rs), np.array(gs), np.array(ns), sig2


def directional_variogram(z: np.ndarray, axis: int, rmax_m: float,
                          cell_m: float = PIXEL_M, min_pairs: int = MIN_PAIRS):
    """gamma(h) along one axis, averaged over the orthogonal direction.

    axis = 1 -> east-west (along rows); axis = 0 -> north-south (columns).
    Returns (h in metres, gamma, number of lines averaged).
    """
    zc = z - z.mean()
    n = z.shape[axis]
    F = np.fft.fft(zc, axis=axis)
    acf = np.real(np.fft.ifft(F * np.conj(F), axis=axis)) / n
    ac = acf.mean(axis=1 - axis)          # average across the other direction
    ac = np.fft.fftshift(ac)
    sig2 = float(zc.var())
    hs, gs = [], []
    for k in range(int(rmax_m // cell_m) + 1):
        if k >= n // 2:
            break
        hs.append(k * cell_m)
        gs.append(sig2 - float(ac[n // 2 + k]))
    return np.array(hs), np.array(gs), int(z.shape[1 - axis])


def exp_model(h, c0, c1, a):
    return c0 + c1 * (1.0 - np.exp(-h / a))


def nested_model(h, c0, c1, a1, c2, a2):
    return (c0 + c1 * (1.0 - np.exp(-h / a1))
            + c2 * (1.0 - np.exp(-h / a2)))


def _r2(gs, pred):
    ss_res = float(np.sum((gs - pred) ** 2))
    ss_tot = float(np.sum((gs - gs.mean()) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def fit_exp(rs, gs, sig2, rmax_m):
    p0 = [0.0, sig2 * 0.9, 300.0]
    bounds = ([0.0, 0.0, 1e-3], [sig2, sig2 * 3.0, rmax_m * 10.0])
    popt, _ = curve_fit(exp_model, rs, gs, p0=p0, bounds=bounds, maxfev=40000)
    return tuple(float(v) for v in popt), _r2(gs, exp_model(rs, *popt))


def fit_nested_bounded(rs, gs, sig2, rmax_m):
    """Two-structure fit with the long range held inside the data.

    The unbounded version can place the second structure far beyond the lag
    window, where it is flat over the observed lags and therefore unidentified:
    R^2 stays at 0.998 while the block variance it implies is out by a factor
    of four.  Constraining it is the standard identifiability fix.
    """
    p0 = [0.0, sig2 * 0.6, 60.0, sig2 * 0.4, min(2000.0, rmax_m)]
    bounds = ([0.0, 0.0, 1e-3, 0.0, 1e-3],
              [sig2, sig2 * 3.0, rmax_m * 10.0, sig2 * 3.0, rmax_m])
    popt, _ = curve_fit(nested_model, rs, gs, p0=p0, bounds=bounds,
                        maxfev=200000)
    return tuple(float(v) for v in popt), _r2(gs, nested_model(rs, *popt))


def fit_nested_unbounded(rs, gs, sig2):
    """The two-structure fit WITHOUT the identifiability constraint.

    Kept on the record because its failure is the point: it can place the long
    structure far beyond the lag window, where it is flat over the observed
    lags, so R^2 stays near 1 while the block variance it implies is badly
    wrong.  The bounded fit is the one quoted in the text; this one is shown
    beside it so the reader can see what was being avoided.
    """
    p0 = [0.0, sig2 * 0.6, 60.0, sig2 * 0.4, 2000.0]
    bounds = ([0.0, 0.0, 1e-3, 0.0, 1e-3],
              [sig2, sig2 * 3.0, 1e5, sig2 * 3.0, 1e5])
    popt, _ = curve_fit(nested_model, rs, gs, p0=p0, bounds=bounds,
                        maxfev=200000)
    return tuple(float(v) for v in popt), _r2(gs, nested_model(rs, *popt))


def block_distances(m_side: int, cell_m: float = PIXEL_M):
    coords = np.stack(np.meshgrid(np.arange(m_side), np.arange(m_side),
                                  indexing="ij"), axis=-1).reshape(-1, 2)
    coords = coords * cell_m
    d = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=-1))
    return d, coords.shape[0]


def var_block_model(m_side: int, c0: float, c1: float, a: float):
    d, m = block_distances(m_side)
    C = c1 * np.exp(-d / a)
    np.fill_diagonal(C, c0 + c1)
    return float(C.sum() / (m * m))


def var_block_model_nested(m_side, c0, c1, a1, c2, a2):
    d, m = block_distances(m_side)
    C = c1 * np.exp(-d / a1) + c2 * np.exp(-d / a2)
    np.fill_diagonal(C, c0 + c1 + c2)
    return float(C.sum() / (m * m))


def var_block_empirical(m_side, rs, gs, sigma2):
    d, m = block_distances(m_side)
    g = np.interp(d, rs, gs)
    np.fill_diagonal(g, 0.0)
    return float(sigma2 - g.sum() / (m * m))


def rmax_for(pix: np.ndarray) -> tuple[float, float]:
    """(fit window, truncated window) in metres, for this window's shape."""
    short_m = min(pix.shape) * PIXEL_M
    fit = min(RMAX_FIT_CAP_M, RMAX_SHORT_FRAC * short_m)
    return fit, RMAX_TRUNC_FRAC * fit


# --------------------------------------------------------------------------
# A. change of support
# --------------------------------------------------------------------------
def run_blocks() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rows, trows, laws = [], [], {}
    for tile, row, col, biome, group, _p in WINDOWS:
        lab = label_of(tile, row, col)
        t0 = time.time()
        z = metric_binary(tile, row, col)
        pix = pixel_field(z)
        laws[lab] = pix
        fine_values = np.arange(0, PIXEL ** 2 + 1) / PIXEL ** 2
        vals, probs = empirical_law(pix, values=fine_values)
        _, var_pix = marginal_moments(vals, probs)

        for b in BLOCK_PIXELS:
            m = b * b
            blk = block_composition_map(pix, b)
            truth_var = float(blk.var())
            v_ind = independent_variance(vals, probs, m)
            v_plus = comonotonic_variance(vals, probs, m)
            lo_psd, _, _, rho_psd, _ = psd_variance_bounds(vals, probs, m)
            if m <= EXACT_M:
                lo_exact = identified_set_bounds(
                    vals, probs, m, kind="variance")["lower"]
                method = "exact LP"
            else:
                lo_exact, method = np.nan, "PSD relaxation"
            v_minus = lo_exact if m <= EXACT_M else lo_psd

            a = pix - pix.mean()
            rows.append(dict(
                window=lab, biome=biome, group=group,
                block_pixels=b, m=m, method=method,
                block_side_m=b * PIXEL_M,
                support_size=int(len(vals)),
                sigma2=var_pix,
                truth_block_mean=float(blk.mean()),
                truth_block_var=truth_var,
                v_independent=v_ind, v_lower=v_minus, v_upper=v_plus,
                independent_understates_by=v_ind / truth_var if truth_var else np.nan,
                aci=(v_plus - v_minus) / v_plus if v_plus else np.nan,
                truth_inside=bool(v_minus - 1e-12 <= truth_var <= v_plus + 1e-12),
                lag1_autocorr_pixel=float((a[:, :-1] * a[:, 1:]).sum()
                                          / (a * a).sum()),
            ))

            for tau in TAUS:
                ind_v, ind_p = independent_aggregate_law(vals, probs, m)
                rec = dict(window=lab, biome=biome, group=group,
                           block_pixels=b, m=m, tau=tau,
                           truth_tail=float((blk >= tau - 1e-12).mean()),
                           independent_tail=float(ind_p[ind_v >= tau - 1e-12].sum()),
                           lower_bound=np.nan, upper_bound=np.nan)
                if m <= EXACT_M:
                    e = envelope(vals, probs, m, [tau])[0]
                    rec["lower_bound"] = e["lower_exact"]
                    rec["upper_bound"] = e["upper_exact"]
                trows.append(rec)
        print(f"  blocks {lab:18s} sigma2 {var_pix:.5f}  ({time.time() - t0:.1f} s)")

    B = pd.DataFrame(rows)
    T = pd.DataFrame(trows)
    B.to_csv(RES / "metric_blocks.csv", index=False)
    T.to_csv(RES / "metric_exceedance.csv", index=False)
    return B, T, laws


# --------------------------------------------------------------------------
# B. geostatistical comparator
# --------------------------------------------------------------------------
def run_baseline(blocks: pd.DataFrame, laws: dict) -> None:
    rows, curves, fits, dec = [], [], [], []
    for tile, row, col, biome, group, _p in WINDOWS:
        lab = label_of(tile, row, col)
        t0 = time.time()
        pix = laws[lab]
        sig2 = float(pix.var())
        rmax_fit, rmax_trunc = rmax_for(pix)

        rs_t, gs_t, _, _ = radial_variogram(pix, rmax_trunc)
        (c0t, c1t, at), r2t = fit_exp(rs_t, gs_t, sig2, rmax_trunc)
        rs, gs, ns, _ = radial_variogram(pix, rmax_fit)
        (c0, c1, a), r2 = fit_exp(rs, gs, sig2, rmax_fit)
        (b0, b1, ba1, b2, ba2), r2b = fit_nested_bounded(
            rs, gs, sig2, rmax_fit)

        (d0, d1, da1, d2, da2), r2n = fit_nested_unbounded(rs, gs, sig2)

        for wl, cc, rl, rm in (("truncated", (c0t, c1t, at), r2t, rmax_trunc),
                               ("wide", (c0, c1, a), r2, rmax_fit),
                               ("nested", (b0, b1, ba1, b2, ba2), r2b, rmax_fit),
                               ("nested_unbounded", (d0, d1, da1, d2, da2),
                                r2n, rmax_fit)):
            nested = wl.startswith("nested")
            sill = cc[0] + cc[1] if not nested else cc[0] + cc[1] + cc[3]
            fits.append(dict(
                window=lab, biome=biome, group=group, sigma2=sig2, model=wl,
                rmax_m=rm,
                a1_m=cc[2], a1_px=cc[2] / PIXEL_M,
                a2_m=cc[4] if nested else np.nan,
                a2_px=cc[4] / PIXEL_M if nested else np.nan,
                sill=sill, r2=rl,
                gamma_at_trunc_over_sill=float(
                    np.interp(rmax_trunc, rs, gs) / sig2),
                range_resolved=bool(cc[2] < rm),
            ))

        for h, gam, npair in zip(rs, gs, ns):
            curves.append(dict(window=lab, biome=biome, group=group,
                               lag_m=float(h), lag_px=float(h / PIXEL_M),
                               gamma=float(gam), n_pairs=int(npair),
                               sigma2=sig2, within_trunc=bool(h <= rmax_trunc)))

        sub = blocks[blocks["window"] == lab].sort_values("m")
        for _, r in sub.iterrows():
            m = int(r["m"])
            ms = int(round(np.sqrt(m)))
            v_tr = var_block_model(ms, c0t, c1t, at)
            v_wi = var_block_model(ms, c0, c1, a)
            v_ne = var_block_model_nested(ms, b0, b1, ba1, b2, ba2)
            v_em = var_block_empirical(ms, rs, gs, sig2)
            truth = float(r["truth_block_var"])
            lo = max(0.0, float(r["v_lower"]))
            hi = float(r["v_upper"])
            rows.append(dict(
                window=lab, biome=biome, group=group, m=m, m_side=ms,
                block_side_m=ms * PIXEL_M, sigma2=sig2,
                truth=truth, independence=float(r["v_independent"]),
                identified_lower=lo, identified_upper=hi,
                model_truncated=v_tr, model_wide=v_wi, model_nested=v_ne,
                empirical_variogram=v_em,
                truth_inside=bool(lo - 1e-9 <= truth <= hi + 1e-9),
                err_truncated_pct=100.0 * (v_tr - truth) / truth,
                err_wide_pct=100.0 * (v_wi - truth) / truth,
                err_nested_pct=100.0 * (v_ne - truth) / truth,
                err_empirical_pct=100.0 * (v_em - truth) / truth,
                err_independence_pct=100.0
                * (float(r["v_independent"]) - truth) / truth,
            ))
            deff = {"truth": m * truth / sig2, "model_wide": m * v_wi / sig2,
                    "model_truncated": m * v_tr / sig2,
                    "empirical": m * v_em / sig2, "independence": 1.0,
                    "upper": m * hi / sig2, "lower": m * lo / sig2}
            rr = dict(window=lab, biome=biome, group=group, m=m,
                      block_side_m=ms * PIXEL_M,
                      n_independence=sig2 / TARGET_PRECISION ** 2)
            for k, v in deff.items():
                rr[f"deff_{k}"] = v
                rr[f"n_{k}"] = sig2 * v / TARGET_PRECISION ** 2
            dec.append(rr)
        print(f"  baseline {lab:18s} range_trunc {at/PIXEL_M:6.1f} px "
              f"({at:7.0f} m)  wide {a/PIXEL_M:6.1f} px ({a:7.0f} m)  "
              f"rmax {rmax_fit:6.0f} m  ({time.time() - t0:.1f} s)")

    pd.DataFrame(rows).to_csv(RES / "metric_variogram.csv", index=False)
    pd.DataFrame(curves).to_csv(RES / "metric_variogram_curves.csv", index=False)
    pd.DataFrame(fits).to_csv(RES / "metric_fits.csv", index=False)
    pd.DataFrame(dec).to_csv(RES / "metric_decision.csv", index=False)


# --------------------------------------------------------------------------
# C. one stored dependence scalar
# --------------------------------------------------------------------------
def run_rho_scalar(blocks: pd.DataFrame, laws: dict) -> None:
    rows = []
    for tile, row, col, biome, group, _p in WINDOWS:
        lab = label_of(tile, row, col)
        pix = laws[lab]
        a = pix - pix.mean()
        rho_lag1 = float((a[:, :-1] * a[:, 1:]).sum() / (a * a).sum())
        sig2 = float(pix.var())
        for b in BLOCK_PIXELS:
            m = b * b
            means = block_means(pix, b)
            truth = float(means.var())
            rho_exact = (m * truth / sig2 - 1.0) / (m - 1.0)
            v_lag1 = sig2 * (1.0 + (m - 1.0) * rho_lag1) / m
            deff_true = m * truth / sig2
            deff_lag1 = 1.0 + (m - 1.0) * rho_lag1
            rho_lo = -1.0 / (m - 1.0)
            rows.append(dict(
                window=lab, biome=biome, group=group, block_pixels=b, m=m,
                block_side_m=b * PIXEL_M, sigma2=sig2, truth_block_var=truth,
                rho_exact=rho_exact, rho_lag1=rho_lag1,
                rho_exact_in_psd=bool(rho_lo - 1e-12 <= rho_exact <= 1.0 + 1e-12),
                rho_lag1_in_psd=bool(rho_lo - 1e-12 <= rho_lag1 <= 1.0 + 1e-12),
                v_from_lag1=v_lag1,
                err_lag1_pct=100.0 * (v_lag1 - truth) / truth if truth else np.nan,
                deff_true=deff_true, deff_from_lag1=deff_lag1,
                cost_ratio=deff_lag1 / deff_true if deff_true else np.nan,
            ))
    df = pd.DataFrame(rows)
    df["rho_exact_spread"] = df.groupby("window")["rho_exact"].transform(
        lambda s: s.max() - s.min())
    df.to_csv(RES / "metric_rho_scalar.csv", index=False)


# --------------------------------------------------------------------------
# D. training-block estimator
# --------------------------------------------------------------------------
def run_estimator(blocks: pd.DataFrame, laws: dict, reps: int = 400,
                  seed: int = 20261013) -> None:
    rng = np.random.default_rng(seed)
    band_w = {(r.window, int(r.m)): float(r.v_upper - r.v_lower)
              for r in blocks.itertuples()}
    rows = []
    for tile, row, col, biome, group, _p in WINDOWS:
        lab = label_of(tile, row, col)
        pix = laws[lab]
        sig2 = float(pix.var())
        for b in BLOCK_PIXELS:
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
                rel = (np.abs(vhat - truth) / truth * 100.0 if truth
                       else np.full(reps, np.nan))
                rows.append(dict(
                    window=lab, biome=biome, group=group, block_pixels=b, m=m,
                    block_side_m=b * PIXEL_M, K=K, N_blocks=N, sigma2=sig2,
                    truth_block_var=truth, rho_oracle=rho_oracle,
                    v_hat_median=float(np.median(vhat)),
                    v_hat_iqr=float(np.percentile(vhat, 75)
                                    - np.percentile(vhat, 25)),
                    bias_median_pct=float(np.median(vhat - truth) / truth * 100.0)
                    if truth else np.nan,
                    err_median_pct=float(np.median(rel)),
                    err_p90_pct=float(np.percentile(rel, 90)),
                    rho_hat_median=float(np.median(rhos)),
                    rho_hat_frac_gt1=float((rhos > 1.0).mean()),
                    levelA_width=band_w.get((lab, m), float("nan")),
                ))
        print(f"  estimator {lab}")
    df = pd.DataFrame(rows)
    df["shrinkage_iqr"] = df["levelA_width"] / df["v_hat_iqr"]
    df.to_csv(RES / "metric_estimator.csv", index=False)


# --------------------------------------------------------------------------
# E. anisotropy -- only meaningful once the grid is metric
# --------------------------------------------------------------------------
def run_anisotropy(laws: dict) -> None:
    rows, curves = [], []
    for tile, row, col, biome, group, _p in WINDOWS:
        lab = label_of(tile, row, col)
        pix = laws[lab]
        sig2 = float(pix.var())
        rmax_fit, _trunc = rmax_for(pix)

        hs, g_ew, n_ew = directional_variogram(pix, axis=1, rmax_m=rmax_fit)
        _, g_ns, n_ns = directional_variogram(pix, axis=0, rmax_m=rmax_fit)
        n = min(len(hs), len(g_ns))
        hs, g_ew, g_ns = hs[:n], g_ew[:n], g_ns[:n]

        (_, _, a_ew), r2_ew = fit_exp(hs, g_ew, sig2, rmax_fit)
        (_, _, a_ns), r2_ns = fit_exp(hs, g_ns, sig2, rmax_fit)

        for h, ge, gn in zip(hs, g_ew, g_ns):
            curves.append(dict(window=lab, biome=biome, group=group,
                               lag_m=float(h), lag_px=float(h / PIXEL_M),
                               gamma_ew=float(ge), gamma_ns=float(gn),
                               sigma2=sig2))

        rows.append(dict(
            window=lab, biome=biome, group=group, sigma2=sig2,
            lat=round(CENTRE[lab][0], 3),
            rmax_m=rmax_fit,
            range_ew_m=a_ew, range_ns_m=a_ns,
            range_ew_px=a_ew / PIXEL_M, range_ns_px=a_ns / PIXEL_M,
            range_ratio_ew_over_ns=a_ew / a_ns if a_ns else np.nan,
            r2_ew=r2_ew, r2_ns=r2_ns,
            lines_ew=n_ew, lines_ns=n_ns,
            gamma_ew_at_rmax=float(g_ew[-1]), gamma_ns_at_rmax=float(g_ns[-1]),
            gamma_ratio_at_rmax=float(g_ew[-1] / g_ns[-1]) if g_ns[-1] else np.nan,
        ))
        print(f"  anisotropy {lab:18s} range E-W {a_ew:7.0f} m | "
              f"N-S {a_ns:7.0f} m  ratio {a_ew / a_ns:5.2f}")

    pd.DataFrame(rows).to_csv(RES / "metric_anisotropy.csv", index=False)
    pd.DataFrame(curves).to_csv(RES / "metric_anisotropy_curves.csv", index=False)


def main() -> None:
    print("=== A. change of support on metric support ===")
    blocks, _tails, laws = run_blocks()
    print(f"truth inside the band in all cases: "
          f"{bool(blocks['truth_inside'].all())}")
    print(f"independence understates by "
          f"{blocks['independent_understates_by'].min():.2f}x - "
          f"{blocks['independent_understates_by'].max():.2f}x")
    print(f"ACI {blocks['aci'].min():.3f} - {blocks['aci'].max():.3f}")

    print("\n=== B. geostatistical comparator, lags in metres ===")
    run_baseline(blocks, laws)

    print("\n=== C. one stored dependence scalar ===")
    run_rho_scalar(blocks, laws)

    print("\n=== D. training-block estimator ===")
    run_estimator(blocks, laws)

    print("\n=== E. anisotropy ===")
    run_anisotropy(laws)

    print(f"\nwrote 10 CSVs to {RES}")


if __name__ == "__main__":
    main()
