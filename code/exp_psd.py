"""
exp_psd.py -- closed-form sharp bounds for the exchangeable two-state case,
and the asymptotic sharpness of the positive-semidefinite relaxation.

Theory being tested
-------------------
m exchangeable units, common two-state marginal Bernoulli(p), equal area
weights w_j = 1/m.  Write pi = frac(m p), sigma^2 = p(1-p), v = sigma^2.

  (i)   V+ = comonotonic variance = sigma^2, always attainable.
  (ii)  V- (sharp) = pi (1 - pi) / m^2, attained by the *cyclic-rotation*
        coupling  nu_j = 1{ (U + j p) mod 1 < p },  U ~ Uniform[0,1),
        whose sum S takes only the two values floor(m p), ceil(m p).
        Hence the level-A interval for the variance is closed-form:
            [ pi (1-pi)/m^2 ,  p (1-p) ].
  (iii) For a two-state marginal the sharp minimum pairwise correlation is
        rho_FH = -min(p,1-p)/max(p,1-p).  The PSD relaxation replaces it by
        rho_PSD = max(rho_FH, -1/(m-1)), giving
            V-_PSD = v [ 1/m + ((m-1)/m) rho_PSD ].
        Two regimes:
          * m p <= 1  ->  the pairwise FH bound binds, V-_PSD = V-  (gap 0);
          * m p >  1  ->  the PSD corner rho = -1/(m-1) binds, V-_PSD = 0,
                          and the joint-mixability gap is pi (1-pi)/m^2.
        So gap(m,p) > 0  <=>  m p > 1 and m p not an integer.
  (iv)  For m p > 1 the relative gap obeys
            gap / v  <=  1 / (4 m^2 p (1-p))  =  O(m^-2),
        with the envelope 1/(4 m^2) attained when pi = 1/2.  The PSD
        relaxation is therefore asymptotically exact and needs no exact LP.

Checks: (i)-(iii) against the exact LP of pcp_core.identified_set_bounds and
pcp_core.psd_variance_bounds for m <= 8; (ii) additionally by simulating the
cyclic-rotation coupling directly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pcp_core import (
    identified_set_bounds,
    psd_variance_bounds,
    transport_cov_bounds,
)

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

VALUES = np.array([0.0, 1.0])


def frac(x: float) -> float:
    return x - np.floor(x)


def closed_lower(p: float, m: int) -> float:
    """Sharp minimum of Var(nu_C): pi(1-pi)/m^2, pi = frac(m p)."""
    pi = frac(m * p)
    return pi * (1.0 - pi) / (m * m)


def closed_upper(p: float) -> float:
    return p * (1.0 - p)


def rho_fh_min(p: float) -> float:
    """Sharp minimum pairwise correlation of a two-state marginal."""
    return -min(p, 1.0 - p) / max(p, 1.0 - p)


def psd_lower_closed(p: float, m: int) -> float:
    if m == 1:
        return p * (1.0 - p)
    v = p * (1.0 - p)
    rho_psd = max(rho_fh_min(p), -1.0 / (m - 1))
    return v * (1.0 / m + ((m - 1.0) / m) * rho_psd)


def gap_closed(p: float, m: int) -> float:
    """Joint-mixability gap = exact sharp minimum minus the PSD minimum."""
    return float(max(0.0, closed_lower(p, m) - psd_lower_closed(p, m)))


def gap_envelope(m: int) -> float:
    return 1.0 / (4.0 * m * m)


def rel_gap_envelope(p: float, m: int) -> float:
    return 1.0 / (4.0 * m * m * p * (1.0 - p))


def rotation_sum(u: float, p: float, m: int) -> int:
    """Number of j < m with (u + j p) mod 1 < p -- the cyclic-rotation sum."""
    return int(np.sum(((u + np.arange(m) * p) % 1.0) < p))


# --------------------------------------------------------------------------
# 1. Closed form vs exact LP, and vs the PSD relaxation, for small m
# --------------------------------------------------------------------------
def verify_against_lp():
    probs = [0.2, 0.3, 1.0 / np.pi, 0.5, 0.7, 0.8]
    rows = []
    for p in probs:
        for m in range(2, 9):
            ex = identified_set_bounds(VALUES, np.array([1 - p, p]), m,
                                       kind="variance")
            lo_psd, _, _, _, _ = psd_variance_bounds(
                VALUES, np.array([1 - p, p]), m)
            rows.append(dict(
                p=p, m=m, mp=m * p,
                lp_lower=ex["lower"], closed_lower=closed_lower(p, m),
                lower_error=closed_lower(p, m) - ex["lower"],
                lp_upper=ex["upper"], closed_upper=closed_upper(p),
                lp_psd_lower=max(0.0, lo_psd),
                closed_psd_lower=psd_lower_closed(p, m),
                psd_lower_error=psd_lower_closed(p, m) - max(0.0, lo_psd),
                lp_gap=ex["lower"] - max(0.0, lo_psd),
                closed_gap=gap_closed(p, m),
                gap_error=(ex["lower"] - max(0.0, lo_psd))
                          - gap_closed(p, m),
                regime=("PSD-corner" if m * min(p, 1 - p) > 1
                        else "pairwise-FH"),
                fh_lower=(1.0 / m) * p * (1 - p)
                         + ((m * m - m) / (m * m)) * transport_cov_bounds(
                             VALUES, np.array([1 - p, p]))[0],
                mp_is_integer=bool(frac(m * p) < 1e-12
                                   or 1 - frac(m * p) < 1e-12),
            ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 2. Achievability of the cyclic-rotation coupling
# --------------------------------------------------------------------------
def verify_rotation():
    rows = []
    for p in [0.2, 0.3, 1.0 / np.pi, 0.5, 0.7]:
        for m in range(2, 13):
            us = (np.arange(2000) + 0.5) / 2000.0
            mult = np.array([rotation_sum(float(u), p, m) for u in us])
            lo, hi = int(np.floor(m * p)), int(np.ceil(m * p))
            rows.append(dict(
                p=p, m=m, floor_mp=lo, ceil_mp=hi,
                observed_min=int(mult.min()), observed_max=int(mult.max()),
                support_matches=bool(mult.min() >= lo and mult.max() <= hi),
                mean_sum=float(mult.mean()), target_mean=m * p,
                mean_error=float(abs(mult.mean() - m * p)),
                var_sum=float(mult.var()),
                closed_var_sum=frac(m * p) * (1 - frac(m * p)),
                var_error=float(abs(mult.var()
                                    - frac(m * p) * (1 - frac(m * p)))),
            ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 3. Sweep in m: gap, envelope, relative gap
# --------------------------------------------------------------------------
def sweep_m():
    ps = [0.2, 0.3, 0.5, 0.7, 0.8, 1.0 / np.pi]
    rows = []
    for m in range(2, 65):
        env = gap_envelope(m)
        for p in ps:
            v = p * (1.0 - p)
            g = gap_closed(p, m)
            rel = g / v
            rel_env = rel_gap_envelope(p, m)
            rows.append(dict(
                m=m, p=p, mp=m * p,
                gap=g, gap_envelope=env,
                gap_over_envelope=g / env if env else np.nan,
                rel_gap=rel, rel_gap_envelope=rel_env,
                rel_over_envelope=rel / rel_env if rel_env else np.nan,
                regime=("PSD-corner" if m * min(p, 1 - p) > 1
                        else "pairwise-FH"),
                psd_exact=bool(g <= 1e-15),
            ))
    return pd.DataFrame(rows)


def sweep_p():
    """For fixed m, maximise the *true* gap over p; confirm sup = 1/(4m^2)."""
    rows = []
    for m in [2, 3, 4, 5, 8, 16, 32, 64]:
        ps = np.linspace(0.02, 0.98, 1921)
        gaps = np.array([gap_closed(float(p), m) for p in ps])
        j = int(gaps.argmax())
        rows.append(dict(
            m=m,
            gap_sup_observed=float(gaps[j]),
            gap_envelope=gap_envelope(m),
            ratio=float(gaps[j] / gap_envelope(m)),
            argmax_p=float(ps[j]), argmax_mp=float(ps[j] * m),
        ))
    return pd.DataFrame(rows)


def fh_table():
    """Normalised lower endpoints: exact, PSD, and Frechet-Hoeffding."""
    rows = []
    for p in [0.1, 0.3, 0.5]:
        v = p * (1 - p)
        for m in range(2, 17):
            cmin, cmax = transport_cov_bounds(VALUES, np.array([1 - p, p]))
            quad_cross = (m * m - m) / (m * m)
            fh_lower = (1.0 / m) * v + quad_cross * cmin
            rows.append(dict(
                p=p, m=m, mp=m * p,
                exact_norm=closed_lower(p, m) / v,
                psd_norm=psd_lower_closed(p, m) / v,
                fh_norm=fh_lower / v,
                fh_is_negative=bool(fh_lower < 0),
            ))
    return pd.DataFrame(rows)


def main():
    df_lp = verify_against_lp()
    df_lp.to_csv(RES / "psd_verify_lp.csv", index=False)
    df_rot = verify_rotation()
    df_rot.to_csv(RES / "psd_rotation.csv", index=False)
    df_m = sweep_m()
    df_m.to_csv(RES / "psd_sweep_m.csv", index=False)
    df_p = sweep_p()
    df_p.to_csv(RES / "psd_sweep_p.csv", index=False)
    df_fh = fh_table()
    df_fh.to_csv(RES / "psd_fh_table.csv", index=False)

    print("== closed form vs exact LP ==")
    print(f"  max |closed_lower - LP lower|  : "
          f"{np.abs(df_lp['lower_error']).max():.3e}")
    print(f"  max |closed_upper - LP upper|  : "
          f"{np.abs(df_lp['lp_upper'] - df_lp['closed_upper']).max():.3e}")
    print(f"  max |closed_psd - LP psd lower|: "
          f"{np.abs(df_lp['psd_lower_error']).max():.3e}")
    print(f"  max |closed_gap - LP gap|      : "
          f"{np.abs(df_lp['gap_error']).max():.3e}")
    print(f"  rows in PSD-corner regime       : "
          f"{int((df_lp['regime'] == 'PSD-corner').sum())} / {len(df_lp)}")
    sparse = (df_lp['m'] * np.minimum(df_lp['p'], 1 - df_lp['p']) <= 1)
    print("  gap > 0 exactly when m*min(p,1-p) > 1 and mp not integer:",
          bool(((df_lp['closed_gap'] > 1e-15)
                == (~sparse & ~df_lp['mp_is_integer'])).all()))
    print(f"  sparse regime (FH exact, gap 0): {int(sparse.sum())}/{len(df_lp)}"
          f"  -> all gaps zero: {bool((df_lp.loc[sparse,'closed_gap']<=1e-15).all())}")
    print("== rotation coupling realises {floor(mp), ceil(mp)} ==")
    print(f"  support_matches all : {bool(df_rot['support_matches'].all())}")
    print(f"  max |mean - mp|     : {df_rot['mean_error'].max():.3e}")
    print(f"  max |var(S)-pi(1-pi)|: {df_rot['var_error'].max():.3e}")
    print("== sup over p of the true gap / (1/4m^2) ==")
    print(df_p.to_string(index=False))
    print("== sweep m: gap/envelope and rel/envelope over m = 2..64 ==")
    print(f"  gap/envelope in [{df_m['gap_over_envelope'].min():.4f}, "
          f"{df_m['gap_over_envelope'].max():.4f}]")
    print(f"  rel/envelope in [{df_m['rel_over_envelope'].min():.4f}, "
          f"{df_m['rel_over_envelope'].max():.4f}]")
    print("== worst-case relative gap (envelope) by m, p = 0.5 ==")
    for m in [4, 10, 16, 25, 64]:
        print(f"  m={m:3d}  1/(4 m^2 p(1-p)) = {rel_gap_envelope(0.5, m):.5f}")
    nz = df_m[(df_m.p == 0.5) & (df_m.gap > 1e-15)].sort_values("m")
    print(f"  p=0.5: nonzero gap only at odd m; largest plot gap "
          f"{nz['gap'].max():.5f} at m={int(nz.loc[nz['gap'].idxmax(), 'm'])}")
    print("== pairwise FH lower endpoint goes negative for m >= 3 ==")
    print(df_fh[df_fh.p == 0.5][["m", "exact_norm", "psd_norm", "fh_norm",
                                 "fh_is_negative"]].head(6).to_string(index=False))
    print("wrote psd_verify_lp.csv psd_rotation.csv psd_sweep_m.csv "
          "psd_sweep_p.csv psd_fh_table.csv")


if __name__ == "__main__":
    main()
