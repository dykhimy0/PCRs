"""
exp_baseline.py -- geostatistical baseline and a decision-relevant consequence,
on all three real ESA WorldCover windows.

Part A  For each of the three windows, estimate the variogram of the 30 m
        tree-fraction field and turn it into a model-based block variance for
        m = 4, 9, 16, 25, placed against three references:
            independence        sigma^2 / m
            enumerated truth    (every block in the window)
            sharp identified bounds (level A)
        Two lag windows are used deliberately.  A truncated window
        (rmax = 60 pixels) is the naive default and is shown to bias the block
        variance low, because the range is not resolved there.  A wide window
        (rmax = 300 pixels) resolves the range and removes most of the bias.
        Three estimators are compared: a single exponential fitted on the
        truncated window, a single exponential fitted on the wide window, a
        nested two-structure model, and a non-parametric read-off of the binned
        empirical variogram.

Part B  Propagate to a decision: how many sampling units are needed to estimate
        the regional tree-cover fraction to a target precision, under each of
        independence, the fitted variogram model, and the sharp bounds.  The
        design effect DEFF = m * Var_block / sigma^2 is the quantity that turns
        a variance into a sample size.

Outputs
    results/baseline_variogram.csv         per window x m, all estimators
    results/baseline_variogram_curves.csv  the binned empirical variograms
    results/baseline_fits.csv              fitted parameters per window
    results/baseline_decision.csv          per window x m, sample sizes
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from windows import biome_of

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
DATA = ROOT / "data"
RES.mkdir(exist_ok=True)

PIXEL = 3
RMAX_TRUNC = 60              # naive default lag window
RMAX_FIT = 300               # wide enough for the range to be resolved
TARGET_PRECISION = 0.01      # +/- 1 percentage point on the mean share


# --------------------------------------------------------------------------
# empirical variogram via FFT autocovariance
# --------------------------------------------------------------------------
def pixel_field(window_idx: int, pixel: int = PIXEL):
    f = sorted(DATA.glob("worldcover_*.npy"))[window_idx]
    arr = np.load(f)
    n = arr.shape[0]
    k = n // pixel
    blk = arr[:k * pixel, :k * pixel].reshape(k, pixel, k, pixel)
    tree = (blk == 10).mean(axis=(1, 3))
    return tree.astype(float), f.name


def radial_variogram(z: np.ndarray, rmax: int = RMAX_FIT):
    """sigma^2 - C(r), averaged over integer lag radii, from the 2-D FFT."""
    n = z.shape[0]
    zc = z - z.mean()
    F = np.fft.fft2(zc)
    ac = np.real(np.fft.ifft2(F * np.conj(F))) / zc.size
    ac = np.fft.fftshift(ac)
    oy = (np.arange(n) - n // 2)[:, None]
    ox = (np.arange(n) - n // 2)[None, :]
    r = np.sqrt(oy ** 2 + ox ** 2)
    bins = np.floor(r).astype(int)
    sig2 = float(zc.var())
    rs, gs = [], []
    for k in range(rmax + 1):
        sel = bins == k
        if sel.sum() < 8:
            continue
        rs.append(float(k))
        gs.append(sig2 - float(ac[sel].mean()))
    return np.array(rs), np.array(gs), sig2


def exp_model(h, c0, c1, a):
    return c0 + c1 * (1.0 - np.exp(-h / a))


def nested_model(h, c0, c1, a1, c2, a2):
    return (c0 + c1 * (1.0 - np.exp(-h / a1))
            + c2 * (1.0 - np.exp(-h / a2)))


def fit_exp(rs, gs, sig2):
    p0 = [0.0, sig2 * 0.9, 5.0]
    bounds = ([0.0, 0.0, 1e-3], [sig2, sig2 * 3.0, 1e4])
    popt, _ = curve_fit(exp_model, rs, gs, p0=p0, bounds=bounds, maxfev=40000)
    pred = exp_model(rs, *popt)
    ss_res = float(np.sum((gs - pred) ** 2))
    ss_tot = float(np.sum((gs - gs.mean()) ** 2))
    return tuple(float(v) for v in popt), 1.0 - ss_res / ss_tot


def fit_nested(rs, gs, sig2):
    p0 = [0.0, sig2 * 0.6, 3.0, sig2 * 0.4, 100.0]
    bounds = ([0.0, 0.0, 1e-3, 0.0, 1e-3],
              [sig2, sig2 * 3.0, 1e3, sig2 * 3.0, 1e5])
    popt, _ = curve_fit(nested_model, rs, gs, p0=p0, bounds=bounds,
                        maxfev=200000)
    pred = nested_model(rs, *popt)
    ss_res = float(np.sum((gs - pred) ** 2))
    ss_tot = float(np.sum((gs - gs.mean()) ** 2))
    return tuple(float(v) for v in popt), 1.0 - ss_res / ss_tot


# --------------------------------------------------------------------------
# block variance under different estimators
# --------------------------------------------------------------------------
def block_distances(m_side: int):
    coords = np.stack(np.meshgrid(np.arange(m_side), np.arange(m_side),
                                  indexing="ij"), axis=-1).reshape(-1, 2)
    d = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=-1))
    return d, coords.shape[0]


def var_block_model(m_side: int, c0: float, c1: float, a: float):
    """Model-consistent Var(block mean) = (1/m^2) sum_ij C(h_ij)."""
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
    """Var(block mean) = sigma2 - (1/m^2) sum_ij gamma(h_ij), gamma(0) = 0."""
    d, m = block_distances(m_side)
    g = np.interp(d, rs, gs)
    np.fill_diagonal(g, 0.0)
    return float(sigma2 - g.sum() / (m * m))


def main() -> None:
    blocks = pd.read_csv(RES / "real_blocks.csv")
    marg = pd.read_csv(RES / "real_marginal.csv")
    # filename -> analysis label, from the window table's own record
    label_of = {
        f"worldcover_{r.tile}_{int(r.row)}_{int(r.col)}_3072.npy": r.window
        for r in marg.itertuples()
    }
    n_files = len(sorted(DATA.glob("worldcover_*.npy")))

    rows, curves, fits, dec = [], [], [], []

    for idx in range(n_files):
        tree, fname = pixel_field(idx)
        sig2 = float(tree.var())
        # a window is identified by its pixel-scale variance, which is unique
        sub = blocks[np.isclose(blocks["pixel_marginal_var"], sig2, atol=1e-9)]
        sub = sub.sort_values("m")
        if sub.empty:
            print(f"[skip] no enumerated truth for {fname} (var {sig2:.6f})")
            continue
        label = label_of.get(fname)
        if label is None:
            print(f"[skip] {fname} is not in the window table")
            continue
        biome = biome_of(label)

        # ---- truncated fit (the failure mode) ---------------------------
        rs_t, gs_t, _ = radial_variogram(tree, rmax=RMAX_TRUNC)
        (c0t, c1t, at), r2t = fit_exp(rs_t, gs_t, sig2)

        # ---- wide fit: single and nested -------------------------------
        rs, gs, _ = radial_variogram(tree, rmax=RMAX_FIT)
        (c0, c1, a), r2 = fit_exp(rs, gs, sig2)
        (d0, d1, da1, d2, da2), r2n = fit_nested(rs, gs, sig2)

        print(f"\n===== {label}   [{fname}] =====")
        print(f"  mean {tree.mean():.5f}  sigma^2 {sig2:.5f}")
        print(f"  truncated rmax={RMAX_TRUNC}: c0={c0t:.5f} c1={c1t:.5f} "
              f"a={at:.2f}px sill={c0t + c1t:.5f} "
              f"(sill-sigma^2 {c0t + c1t - sig2:+.5f}) R^2={r2t:.4f}")
        print(f"  wide      rmax={RMAX_FIT}: c0={c0:.5f} c1={c1:.5f} "
              f"a={a:.2f}px sill={c0 + c1:.5f} "
              f"(sill-sigma^2 {c0 + c1 - sig2:+.5f}) R^2={r2:.4f}")
        print(f"  nested                : a1={da1:.2f} a2={da2:.1f} "
              f"sill={d0 + d1 + d2:.5f} R^2={r2n:.4f}")

        for wl, cc, rl in (("truncated", (c0t, c1t, at), r2t),
                           ("wide", (c0, c1, a), r2),
                           ("nested", (d0, d1, da1, d2, da2), r2n)):
            fits.append(dict(
                window=label, biome=biome, window_idx=idx, sigma2=sig2, model=wl,
                rmax=RMAX_TRUNC if wl == "truncated" else RMAX_FIT,
                c0=cc[0], c1=cc[1] if wl != "nested" else cc[1],
                a1=cc[2] if wl != "nested" else cc[2],
                c2=cc[3] if wl == "nested" else np.nan,
                a2=cc[4] if wl == "nested" else np.nan,
                sill=(cc[0] + cc[1]) if wl != "nested"
                else (cc[0] + cc[1] + cc[3]),
                sill_minus_sigma2=((cc[0] + cc[1]) if wl != "nested"
                                   else (cc[0] + cc[1] + cc[3])) - sig2,
                r2=rl,
                gamma_at_trunc=float(np.interp(RMAX_TRUNC, rs, gs)),
                gamma_at_trunc_over_sill=float(
                    np.interp(RMAX_TRUNC, rs, gs) / sig2)))

        # variogram curves for the figure
        for lag, gam in zip(rs, gs):
            curves.append(dict(window=label, biome=biome, window_idx=idx,
                               lag=float(lag),
                               gamma=float(gam), sigma2=sig2,
                               within_trunc=bool(lag <= RMAX_TRUNC)))

        # ---- block variance vs references ------------------------------
        for _, r in sub.iterrows():
            m = int(r["m"])
            m_side = int(round(np.sqrt(m)))
            if m_side ** 2 != m:
                continue
            v_tr = var_block_model(m_side, c0t, c1t, at)
            v_wi = var_block_model(m_side, c0, c1, a)
            v_ne = var_block_model_nested(m_side, d0, d1, da1, d2, da2)
            v_em = var_block_empirical(m_side, rs, gs, sig2)
            truth = float(r["truth_block_var"])
            lo = max(0.0, float(r["v_lower"]))
            hi = float(r["v_upper"])
            rows.append(dict(
                window=label, biome=biome, window_idx=idx, m=m,
                m_side=m_side, sigma2=sig2,
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
            deff = {
                "truth": m * truth / sig2,
                "model_wide": m * v_wi / sig2,
                "model_truncated": m * v_tr / sig2,
                "empirical": m * v_em / sig2,
                "independence": 1.0,
                "upper": m * hi / sig2,
                "lower": m * lo / sig2,
            }
            row = dict(window=label, biome=biome, window_idx=idx, m=m,
                       n_independence=sig2 / TARGET_PRECISION ** 2)
            for k, v in deff.items():
                row[f"deff_{k}"] = v
                row[f"n_{k}"] = sig2 * v / TARGET_PRECISION ** 2
            dec.append(row)

    df = pd.DataFrame(rows)
    cv = pd.DataFrame(curves)
    ft = pd.DataFrame(fits)
    dd = pd.DataFrame(dec)
    df.to_csv(RES / "baseline_variogram.csv", index=False)
    cv.to_csv(RES / "baseline_variogram_curves.csv", index=False)
    ft.to_csv(RES / "baseline_fits.csv", index=False)
    dd.to_csv(RES / "baseline_decision.csv", index=False)

    # ---- summary across all windows -------------------------------------
    print("\n\n################ SUMMARY ACROSS ALL WINDOWS ################")
    n_win = df["window"].nunique()
    print(f"cases: {len(df)} ({n_win} windows x "
          f"{df['m'].nunique()} block sizes, "
          f"{df['biome'].nunique()} biomes)")
    print("\n== error against enumerated truth, by estimator (%) ==")
    errs = ["err_truncated_pct", "err_wide_pct", "err_nested_pct",
            "err_empirical_pct", "err_independence_pct"]
    print(df[errs].describe().loc[["min", "max", "mean"]].to_string())
    print("\n== mean |error| by window and estimator ==")
    absdf = df.assign(**{c: df[c].abs() for c in errs})
    print(absdf.groupby("window")[errs].mean().to_string())
    print("\n== is the truth inside the identified band? ==")
    print(df.groupby("window")["truth_inside"].all().to_string())
    print("\n== sample size for +/-1 pp on the regional mean ==")
    print(dd[["window", "m", "n_independence", "n_truth", "n_model_wide",
              "n_model_truncated", "n_empirical", "n_upper"]].to_string(
                  index=False, float_format=lambda v: f"{v:,.0f}"))
    infl = (dd["n_truth"] / dd["n_independence"])
    print(f"\ninflation of the sample size vs independence: "
          f"min {infl.min():.1f}x  max {infl.max():.1f}x")

    print("\n== cross-biome: is the understatement a function of m alone? ==")
    per = df.assign(under=df["independence"] / df["truth"],
                    norm=df["truth"] / df["sigma2"])
    print(per.groupby(["m"])[["under", "norm"]]
          .agg(["min", "max", "mean"]).to_string(float_format=lambda v: f"{v:.4f}"))
    print("\n== same two quantities, split by biome (m = 25) ==")
    m25 = per[per["m"] == 25]
    print(m25.groupby("biome")[["under", "norm"]]
          .agg(["min", "max"]).to_string(float_format=lambda v: f"{v:.4f}"))
    print("\n== per-window fit summary ==")
    print(ft[ft["model"] != "nested"][
        ["window", "biome", "model", "sigma2", "sill", "sill_minus_sigma2",
         "a1", "r2", "gamma_at_trunc_over_sill"]].to_string(
             index=False, float_format=lambda v: f"{v:.4f}"))
    print("wrote baseline_variogram.csv baseline_variogram_curves.csv "
          "baseline_fits.csv baseline_decision.csv")


if __name__ == "__main__":
    main()
