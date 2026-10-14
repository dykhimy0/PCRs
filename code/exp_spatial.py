"""
exp_spatial.py -- scale-consistency of a probabilistic-compositional raster
on simulated landscapes with known ground truth.

A fine binary landscape Z is generated from a latent Gaussian field with a
prescribed correlation length.  The raster is stored at the *pixel* level:
each pixel is `pixel` x `pixel` fine cells, and the stored information is the
marginal law Lambda of the pixel composition.  The task is change of support
to a coarser block of `b` x `b` pixels (m = b^2 pixels).  Because the stored
object contains no cross-pixel dependence, the block composition law is only
partially identified; the true law is known here and is used as validation.

Outputs
-------
results/landscape_sweep.csv
results/metadata_ladder.csv
results/coverage_calibration.csv
results/copula_family.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pcp_core import (  # noqa: E402
    block_composition_map,
    comonotonic_variance,
    empirical_law,
    envelope,
    gaussian_landscape,
    identified_set_bounds,
    independent_aggregate_law,
    independent_variance,
    marginal_moments,
    pixel_frame,
    psd_variance_bounds,
)

RES = Path(__file__).resolve().parent.parent / "results"
RES.mkdir(parents=True, exist_ok=True)

N_SIDE = 1024
PIXEL = 2                       # pixel = 2 x 2 fine cells -> q = 5 composition values
RANGES = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
BLOCK_PIXELS = (2, 3, 4)        # m = 4, 9, 16
SEEDS = tuple(range(6))
PROPORTIONS = (0.3, 0.5)
TAUS = (0.25, 0.5, 0.75)
PIXEL_VALUES = np.arange(0, PIXEL ** 2 + 1) / PIXEL ** 2


def stored_marginal(field, pixel=PIXEL):
    """The law of the pixel composition -- the only thing the raster stores."""
    pix = pixel_frame(field, pixel)
    vals, probs = empirical_law(pix, values=PIXEL_VALUES)
    return pix, vals, probs


def spatial_sweep():
    rows = []
    for p in PROPORTIONS:
        for L in RANGES:
            for seed in SEEDS:
                z, _ = gaussian_landscape(N_SIDE, PIXEL, L, seed=seed, proportion=p)
                pix, vals, probs = stored_marginal(z)
                _, var_pix = marginal_moments(vals, probs)
                # lag-1 autocorrelation of the pixel-level composition
                a = pix - pix.mean()
                rho1 = float((a[:, :-1] * a[:, 1:]).sum() / (a * a).sum())

                for b in BLOCK_PIXELS:
                    m = b * b
                    blk = block_composition_map(pix, b)
                    truth_var = float(blk.var())
                    truth_mean = float(blk.mean())

                    v_ind = independent_variance(vals, probs, m)
                    v_plus = comonotonic_variance(vals, probs, m)
                    lo_psd, hi_psd, rho_fh_min, rho_psd_min, _ = psd_variance_bounds(
                        vals, probs, m)
                    if m <= 4:
                        ex = identified_set_bounds(vals, probs, m, kind="variance")
                        v_minus_exact = ex["lower"]
                        method = "exact LP"
                    else:
                        v_minus_exact = np.nan
                        method = "PSD relaxation"

                    rec = dict(
                        proportion=p, correlation_length=L, seed=seed,
                        block_pixels=b, m=m, method=method,
                        pixel_marginal_var=var_pix,
                        lag1_autocorr_pixel=rho1,
                        truth_block_var=truth_var,
                        truth_block_mean=truth_mean,
                        v_independent=v_ind,
                        v_upper_exact=v_plus,
                        v_lower_exact=v_minus_exact,
                        v_lower_psd=lo_psd,
                        rho_psd_min=rho_psd_min, rho_fh_min=rho_fh_min,
                        implied_reduction_truth=var_pix / truth_var if truth_var > 0 else np.inf,
                        implied_reduction_independent=float(m),
                        aci_marginals_only=(v_plus - (v_minus_exact if m <= 4 else lo_psd)) / v_plus,
                        aci_worstcase_bound=1.0 - 1.0 / m,
                        independent_understates_by=v_ind / truth_var if truth_var > 0 else np.nan,
                        truth_inside_bound=bool(
                            (v_minus_exact if m <= 4 else lo_psd) - 1e-12 <= truth_var <= v_plus + 1e-12),
                    )
                    if m == 4:
                        for tau in TAUS:
                            lev = envelope(vals, probs, m, [tau])[0]
                            truth_tail = float((blk >= tau - 1e-12).mean())
                            rec[f"truth_tail_{tau}"] = truth_tail
                            rec[f"lower_tail_{tau}"] = lev["lower_exact"]
                            rec[f"upper_tail_{tau}"] = lev["upper_exact"]
                            rec[f"independent_tail_{tau}"] = lev["independent"]
                            rec[f"comonotonic_tail_{tau}"] = lev["comonotonic"]
                            rec[f"tail_inside_bound_{tau}"] = bool(
                                lev["lower_exact"] - 1e-12 <= truth_tail <= lev["upper_exact"] + 1e-12)
                    rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "landscape_sweep.csv", index=False)
    return df


# --------------------------------------------------------------------------
# One scalar dependence parameter: the exchangeable Gaussian copula
# --------------------------------------------------------------------------

def gaussian_copula_samples(values, probs, m, rho, n=200_000, seed=0):
    """Samples of nu_C under an exchangeable Gaussian copula with correlation rho."""
    rng = np.random.default_rng(seed)
    rho = float(np.clip(rho, -1.0 / (m - 1) + 1e-9, 0.999999))
    L = np.linalg.cholesky((1 - rho) * np.eye(m) + rho * np.ones((m, m)))
    g = rng.standard_normal((n, m)) @ L.T
    from scipy.stats import norm
    u = norm.cdf(g)
    cdf = np.cumsum(probs)
    idx = np.searchsorted(cdf, u, side="left")
    idx = np.clip(idx, 0, len(values) - 1)
    return values[idx].mean(axis=1)


def copula_family(p=0.5, m=4, taus=TAUS, n=200_000):
    """Exceedance probability as a function of the single dependence parameter."""
    vals, probs = gaussian_landscape_law(p)
    rho_min = -1.0 / (m - 1) + 1e-6
    rows = []
    for rho in np.linspace(rho_min, 0.999, 25):
        s = gaussian_copula_samples(vals, probs, m, rho, n=n, seed=11)
        for tau in taus:
            rows.append(dict(p=p, m=m, rho=rho, tau=tau,
                             exceedance=float((s >= tau - 1e-12).mean()),
                             aggregate_var=float(s.var())))
    df = pd.DataFrame(rows)
    df.to_csv(RES / "copula_family.csv", index=False)
    return df


def gaussian_landscape_law(p, pixel=PIXEL, L=4.0, n=N_SIDE, seed=0):
    """Pixel-composition marginal law produced by one simulated landscape."""
    z, _ = gaussian_landscape(n, pixel, L, seed=seed, proportion=p)
    pix = pixel_frame(z, pixel)
    return empirical_law(pix, values=PIXEL_VALUES)


def metadata_ladder():
    """Identified width for variance and exceedance under three metadata levels."""
    rows = []
    for p in PROPORTIONS:
        vals, probs = gaussian_landscape_law(p, L=4.0)
        for b in BLOCK_PIXELS:
            m = b * b
            v_ind = independent_variance(vals, probs, m)
            v_plus = comonotonic_variance(vals, probs, m)
            lo_psd, _, _, _, _ = psd_variance_bounds(vals, probs, m)
            lo_exact = (identified_set_bounds(vals, probs, m, kind="variance")["lower"]
                        if m <= 4 else np.nan)
            v_minus = lo_exact if m <= 4 else lo_psd
            rec = dict(proportion=p, m=m,
                       var_width_M0=v_plus - v_minus,
                       var_width_M1=0.0,          # pairwise covariances pin the variance
                       var_width_M2=0.0,
                       var_value_at_M1=v_ind)
            if m == 4:
                for tau in TAUS:
                    lev = envelope(vals, probs, m, [tau])[0]
                    w0 = lev["upper_exact"] - lev["lower_exact"]
                    # M1: exchangeable Gaussian copula, one scalar rho
                    s_lo = gaussian_copula_samples(vals, probs, m, -1.0 / (m - 1) + 1e-6,
                                                   n=120_000, seed=3)
                    s_hi = gaussian_copula_samples(vals, probs, m, 0.999, n=120_000, seed=4)
                    plo = float((s_lo >= tau - 1e-12).mean())
                    phi = float((s_hi >= tau - 1e-12).mean())
                    w1 = abs(phi - plo)
                    rec[f"tail_tau_{tau}"] = tau
                    rec[f"tail_width_M0_{tau}"] = w0
                    rec[f"tail_width_M1_{tau}"] = w1
                    rec[f"tail_shrinkage_{tau}"] = 1 - w1 / w0 if w0 > 0 else np.nan
            rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "metadata_ladder.csv", index=False)
    return df


# --------------------------------------------------------------------------
# Calibration of the independence assumption
# --------------------------------------------------------------------------

def coverage_calibration(level=0.90):
    """Coverage of a nominal 90% interval built under cross-pixel independence."""
    rows = []
    for p in PROPORTIONS:
        vals, probs = gaussian_landscape_law(p, L=4.0)
        for L in RANGES:
            for seed in SEEDS:
                z, _ = gaussian_landscape(N_SIDE, PIXEL, L, seed=seed, proportion=p)
                pix, v2, p2 = stored_marginal(z)
                for b in BLOCK_PIXELS:
                    m = b * b
                    blk = block_composition_map(pix, b)
                    grid, law = independent_aggregate_law(v2, p2, m)
                    c = np.cumsum(law)
                    lo = grid[np.searchsorted(c, (1 - level) / 2)]
                    hi = grid[np.searchsorted(c, 1 - (1 - level) / 2)]
                    cover = float(((blk >= lo - 1e-12) & (blk <= hi + 1e-12)).mean())
                    rows.append(dict(
                        proportion=p, correlation_length=L, seed=seed, block_pixels=b, m=m,
                        nominal=level, empirical_coverage=cover,
                        nominal_width=float(hi - lo),
                        truth_width=float(np.quantile(blk, 1 - (1 - level) / 2)
                                          - np.quantile(blk, (1 - level) / 2)),
                        truth_var=float(blk.var()),
                        independent_var=independent_variance(v2, p2, m),
                    ))
    df = pd.DataFrame(rows)
    df.to_csv(RES / "coverage_calibration.csv", index=False)
    return df


def main():
    print("== landscape sweep ==", flush=True)
    df = spatial_sweep()
    key = ["proportion", "correlation_length", "m", "truth_block_var", "v_independent",
           "v_lower_exact", "v_upper_exact", "implied_reduction_truth",
           "independent_understates_by", "truth_inside_bound"]
    agg = df.groupby(["proportion", "correlation_length", "m"]).agg(
        truth_block_var=("truth_block_var", "mean"),
        v_independent=("v_independent", "mean"),
        v_lower=("v_lower_exact", "mean"),
        v_lower_psd=("v_lower_psd", "mean"),
        v_upper=("v_upper_exact", "mean"),
        implied_reduction_truth=("implied_reduction_truth", "mean"),
        understatement=("independent_understates_by", "mean"),
        coverage_of_truth_inside=("truth_inside_bound", "mean"),
    ).reset_index()
    print(agg.to_string(index=False, float_format=lambda x: f"{x:.5g}"))

    print("\n== copula family (one scalar dependence parameter) ==", flush=True)
    cf = copula_family()
    print(cf.groupby(["tau", "rho"]).exceedance.mean().unstack(0).to_string(
        float_format=lambda x: f"{x:.4f}"))

    print("\n== metadata ladder ==", flush=True)
    ml = metadata_ladder()
    print(ml.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\n== coverage calibration ==", flush=True)
    cc = coverage_calibration()
    print(cc.groupby(["proportion", "correlation_length", "m"]).empirical_coverage
          .mean().unstack(["proportion", "m"]).to_string(float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
