"""
exp_theory.py -- exact finite-dimensional checks of the change-of-support
identified set, and the scale-amplification result.

Outputs
-------
results/scale_amplification.csv
results/exceedance_envelope.csv
results/pairwise_relaxation_gap.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pcp_core import (  # noqa: E402
    comonotonic_variance,
    envelope,
    identified_set_bounds,
    independent_variance,
    marginal_moments,
    pair_variance_bounds,
    psd_variance_bounds,
)

RES = Path(__file__).resolve().parent.parent / "results"
RES.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------- experiment 1
def scale_amplification(ms=(2, 3, 4, 6, 8, 10, 12, 16), ps=(0.2, 0.5, 0.8)):
    """V+ , V- , V_ind and the ambiguity index as the aggregation ratio grows."""
    rows = []
    for p in ps:
        values = np.array([0.0, 1.0])
        probs = np.array([1 - p, p])
        mu, var = marginal_moments(values, probs)
        for m in ms:
            exact = identified_set_bounds(values, probs, m, kind="variance")
            v_plus = exact["upper"]
            v_minus = exact["lower"]
            v_ind = independent_variance(values, probs, m)
            v_como = comonotonic_variance(values, probs, m)
            rows.append(
                dict(
                    p=p,
                    m=m,
                    marginal_mean=mu,
                    marginal_var=var,
                    v_independent=v_ind,
                    v_comonotonic_closed=v_como,
                    v_upper_exact=v_plus,
                    v_lower_exact=v_minus,
                    aci_exact=(v_plus - v_minus) / v_plus if v_plus > 0 else np.nan,
                    aci_worstcase_bound=1.0 - 1.0 / m,
                    ratio_upper_over_independent=v_plus / v_ind,
                    ratio_m=m,
                    # validity flags for the analytic statements
                    chk_upper_equals_marginal_var=abs(v_plus - var) < 1e-9,
                    chk_lower_below_independent=v_minus <= v_ind + 1e-9,
                    chk_aci_bound=(v_plus - v_minus) / v_plus >= 1.0 - 1.0 / m - 1e-9,
                )
            )
    df = pd.DataFrame(rows)
    df.to_csv(RES / "scale_amplification.csv", index=False)
    return df


# ---------------------------------------------------------------- experiment 2
def exceedance_envelope(p=0.5, m=4, taus=None):
    """Sharp upper/lower P(nu_C >= tau) against the named couplings.

    Shows that no single coupling attains both envelopes: the ordering of
    comonotonic vs countermonotonic flips with the threshold.
    """
    values = np.array([0.0, 1.0])
    probs = np.array([1 - p, p])
    w = np.full(m, 1.0 / m)
    if taus is None:
        taus = np.linspace(1.0 / (2 * m), 1.0, 40)

    # comonotonic coupling: all components equal
    idx = np.arange(len(values))
    como_grid = np.zeros((len(values),) * m)
    for i in range(len(values)):
        como_grid[(i,) * m] = probs[i]
    # countermonotonic coupling for even m with p=0.5: exactly m/2 occupied
    if abs(p - 0.5) < 1e-12 and m % 2 == 0:
        ctr_values = [0.5]
        ctr_probs = [1.0]
    else:
        ctr_values, ctr_probs = None, None

    rows = []
    for tau in taus:
        ex = identified_set_bounds(values, probs, m, kind="exceedance", tau=tau)
        # independent
        from numpy.polynomial import polynomial as _p
        # distribution of sum of m iid Bernoulli(p) -> binomial
        from math import comb
        s = np.array([k / m for k in range(m + 1)])
        pk = np.array([comb(m, k) * p ** k * (1 - p) ** (m - k) for k in range(m + 1)])
        ind = float(pk[s >= tau - 1e-12].sum())
        como = float(sum(pr for v, pr in zip(values, probs) if v >= tau - 1e-12))
        ctr = np.nan
        if ctr_values is not None:
            ctr = float(sum(pr for v, pr in zip(ctr_values, ctr_probs) if v >= tau - 1e-12))
        rows.append(
            dict(
                p=p, m=m, tau=tau,
                lower_exact=ex["lower"], upper_exact=ex["upper"],
                independent=ind, comonotonic=como, countermonotonic=ctr,
            )
        )
    df = pd.DataFrame(rows)
    df.to_csv(RES / "exceedance_envelope.csv", index=False)
    return df


# ---------------------------------------------------------------- experiment 3
def pairwise_relaxation_gap(ms=(2, 3, 4, 5, 6, 8), ps=(0.3, 0.5, 0.7)):
    """Cost of replacing the exact identified set by two tractable relaxations."""
    rows = []
    for p in ps:
        values = np.array([0.0, 1.0])
        probs = np.array([1 - p, p])
        for m in ms:
            exact = identified_set_bounds(values, probs, m, kind="variance")
            lo_fh, hi_fh, cmin, cmax = pair_variance_bounds(values, probs, m)
            lo_psd, hi_psd, rho_fh, rho_psd, var = psd_variance_bounds(values, probs, m)
            w_exact = exact["upper"] - exact["lower"]
            rows.append(
                dict(
                    p=p, m=m,
                    lower_exact=exact["lower"], upper_exact=exact["upper"],
                    lower_fh=lo_fh, upper_fh=hi_fh,
                    lower_psd=lo_psd, upper_psd=hi_psd,
                    rho_fh_min=rho_fh, rho_psd_min=rho_psd, rho_fh_max=cmax / var,
                    width_exact=w_exact,
                    width_fh=hi_fh - lo_fh,
                    width_psd=hi_psd - lo_psd,
                    fh_is_vacuous=lo_fh < 0,
                    excess_width_fh_pct=100.0 * ((hi_fh - lo_fh) - w_exact) / max(1e-12, w_exact),
                    excess_width_psd_pct=100.0 * ((hi_psd - lo_psd) - w_exact) / max(1e-12, w_exact),
                    joint_mixability_gap=exact["lower"] - lo_psd,
                )
            )
    df = pd.DataFrame(rows)
    df.to_csv(RES / "pairwise_relaxation_gap.csv", index=False)
    return df


# ---------------------------------------------------------------- experiment 4
def compositional_case(m=2, n_side=2, q=(0.3, 0.5, 0.2)):
    """K = 3 pixel compositions: a pixel is an n_side x n_side block of cells.

    The stored object is the law Lambda of the pixel composition vector on the
    simplex.  The target is the class-1 (say forest) share of the coarse
    support formed from m such pixels.  A linear functional of nu_C depends only
    on the projected marginals, so the scalar bounds applied to the projected
    law are valid outer bounds on the sharp set; we also compute the exact LP
    over the full composition support for small m.
    """
    from itertools import product as iproduct
    from math import factorial

    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix

    q = np.asarray(q, float)
    q = q / q.sum()
    n_cells = n_side ** 2
    counts = [c for c in iproduct(range(n_cells + 1), repeat=3) if sum(c) == n_cells]
    support = np.array(counts, float) / n_cells
    probs = np.array([
        factorial(n_cells) / np.prod([factorial(int(x)) for x in c]) * np.prod(q ** np.array(c))
        for c in counts
    ])
    probs = probs / probs.sum()

    S = len(support)
    states = np.array(list(iproduct(range(S), repeat=m)))
    comps = support[states]                                  # (S^m, m, 3)
    w = np.full(m, 1.0 / m)
    nu_c = np.einsum("j,ljc->lc", w, comps)
    target = np.array([1.0, 0.0, 0.0])                       # class-1 share
    x = nu_c @ target
    mu = float((support @ target) @ probs)

    S_tot = len(states)
    rows, cols, data, r = [], [], [], 0
    for j in range(m):
        for i in range(S):
            sel = np.flatnonzero(states[:, j] == i)
            rows.append(np.full(sel.size, r)); cols.append(sel)
            data.append(np.ones(sel.size)); r += 1
    rows.append(np.full(S_tot, r)); cols.append(np.arange(S_tot)); data.append(np.ones(S_tot))
    A = csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(r + 1, S_tot))
    b = np.concatenate([np.tile(probs, m), [1.0]])

    out = {"n_composition_states": S, "n_joint_states": S_tot, "mean_target": mu}
    for name, cvec in (("var", (x - mu) ** 2), ("tail_at_mean", (x >= mu).astype(float))):
        lo = linprog(cvec, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
        hi = linprog(-cvec, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
        out[f"{name}_lower"] = float(lo.fun)
        out[f"{name}_upper"] = float(-hi.fun)

    # marginal law of the target functional (the projection actually stored)
    proj = np.round(support @ target, 12)
    vals, inv = np.unique(proj, return_inverse=True)
    pp = np.bincount(inv, weights=probs)
    pp = pp / pp.sum()
    _, var = marginal_moments(vals, pp)
    out["projection_marginal_var"] = var
    out["independent_var"] = independent_variance(vals, pp, m)
    out["comonotonic_var"] = var
    out["projection_support_size"] = int(len(vals))
    pd.DataFrame([out]).to_csv(RES / "compositional_case.csv", index=False)
    return out


# ---------------------------------------------------------------- experiment 5
def _empirical_block_cov(pix, mu):
    """Pairwise covariance matrix of the four pixels of a 2 x 2 pixel block."""
    a = pix - mu
    ch = float((a[:, :-1] * a[:, 1:]).mean())
    cv = float((a[:-1, :] * a[1:, :]).mean())
    cd = float((a[:-1, :-1] * a[1:, 1:]).mean())
    var = float((a * a).mean())
    C = np.array([[var, ch, cv, cd],
                  [ch, var, cd, cv],
                  [cv, cd, var, ch],
                  [cd, cv, ch, var]])
    return C, dict(var=var, lag1_h=ch, lag1_v=cv, diag=cd)


def identification_ladder(m=4, tau_values=(0.25, 0.5, 0.75), L=4.0):
    """Information levels: marginals / + true pairwise covariance / + full law."""
    from pcp_core import (  # local import to keep the top of the file tidy
        block_composition_map,
        bounds_with_pairwise_dependence,
        exchangeable_cov_matrix,
        gaussian_landscape,
        pixel_frame,
    )

    z, _ = gaussian_landscape(1024, 2, L, seed=0, proportion=0.5)
    pix = pixel_frame(z, 2)
    mu = float(pix.mean())
    vals, counts = np.unique(pix.ravel(), return_counts=True)
    vals = vals.astype(float)
    probs = counts / counts.sum()
    _, var = marginal_moments(vals, probs)
    C_true, cparts = _empirical_block_cov(pix, mu)
    rho_equiv = float(C_true[0, 1] / var)

    blk4 = block_composition_map(pix, 2)
    truth_tail = {t: float((blk4 >= t - 1e-12).mean()) for t in tau_values}
    truth_var = float(blk4.var())

    rows = []
    for tau in tau_values:
        lev = envelope(vals, probs, m, [tau])[0]
        rec = dict(
            tau=tau,
            level_A_lower=lev["lower_exact"], level_A_upper=lev["upper_exact"],
            level_A_width=lev["upper_exact"] - lev["lower_exact"],
            independent=lev["independent"], comonotonic=lev["comonotonic"],
            truth_block_tail=truth_tail[tau],
            truth_block_var=truth_var,
            pair_cov_var=cparts["var"], pair_cov_lag1=cparts["lag1_h"],
            pair_cov_diag=cparts["diag"], rho_equiv=rho_equiv,
        )
        lo_B, hi_B, resid = bounds_with_pairwise_dependence(
            vals, probs, m, C_true, tau)
        rec.update(level_B_true_lower=lo_B, level_B_true_upper=hi_B,
                   level_B_true_width=hi_B - lo_B,
                   level_B_true_contains_truth=bool(lo_B - 1e-12 <= truth_tail[tau]
                                                    <= hi_B + 1e-12),
                   level_B_max_residual=float(resid))
        for rho in (-1.0 / (m - 1) + 1e-6, 0.0, rho_equiv, 0.5, 1.0):
            C = exchangeable_cov_matrix(vals, probs, m, rho)
            try:
                lo, hi, _ = bounds_with_pairwise_dependence(vals, probs, m, C, tau)
            except RuntimeError:
                continue
            rec[f"level_B_rho{rho:+.4f}_width"] = hi - lo
            rec[f"level_B_rho{rho:+.4f}_contains_truth"] = bool(
                lo - 1e-12 <= truth_tail[tau] <= hi + 1e-12)
        vlo, vhi, _ = bounds_with_pairwise_dependence(
            vals, probs, m, C_true, tau, kind="variance")
        rec["level_B_variance_width"] = vhi - vlo
        rec["level_B_variance_value"] = vlo
        rec["level_A_contains_truth"] = bool(
            rec["level_A_lower"] - 1e-12 <= truth_tail[tau] <= rec["level_A_upper"] + 1e-12)
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "identification_ladder.csv", index=False)
    return df


def main():
    print("== scale amplification ==")
    df1 = scale_amplification()
    cols = ["p", "m", "v_independent", "v_upper_exact", "v_lower_exact",
            "aci_exact", "aci_worstcase_bound", "ratio_upper_over_independent"]
    print(df1[cols].to_string(index=False))
    print("validity:", df1[["chk_upper_equals_marginal_var", "chk_lower_below_independent",
                            "chk_aci_bound"]].all().to_dict())

    print("\n== exceedance envelope ==")
    df2 = exceedance_envelope()
    print(df2.to_string(index=False))

    print("\n== pairwise relaxation gap ==")
    df3 = pairwise_relaxation_gap()
    print(df3[["p", "m", "width_exact", "width_fh", "width_psd",
               "excess_width_fh_pct", "excess_width_psd_pct",
               "joint_mixability_gap", "fh_is_vacuous"]].to_string(index=False))

    print("\n== compositional case (K=3) ==")
    print(compositional_case())

    print("\n== identification ladder ==")
    dl = identification_ladder()
    show = ["tau", "level_A_lower", "level_A_upper", "level_A_width",
            "level_B_true_lower", "level_B_true_upper", "level_B_true_width",
            "level_B_variance_width", "level_B_variance_value",
            "truth_block_tail", "truth_block_var", "rho_equiv",
            "level_A_contains_truth", "level_B_true_contains_truth",
            "level_B_max_residual"]
    print(dl[show].to_string(index=False, float_format=lambda x: f"{x:.5g}"))


if __name__ == "__main__":
    main()
