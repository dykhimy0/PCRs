"""
pcp_core.py -- core machinery for the scale-consistency study of
probabilistic-compositional pixels (PCPs).

Objects
-------
A fine categorical field Z on a lattice.  A *pixel* j carries a composition
nu_j (a vector of class area shares).  A *probabilistic-compositional pixel*
(PCP) stores the law Lambda_j = Law(nu_j | O).  A *probabilistic-compositional
raster* (PCR) at one support level stores {Lambda_j}.  Aggregating m pixels of
area weight w_j into a coarser support C gives

    nu_C = sum_j w_j nu_j = g_w(nu_1,...,nu_m),
    Lambda_C = (g_w)_# Gamma,      Gamma = Law((nu_1,...,nu_m) | O).

Gamma is never stored, so Lambda_C is only partially identified by {Lambda_j}.

Everything below works in the *scalar* (binary occupancy) reduction: the class
of interest is a single class, nu_j in [0,1].  For a multi-class composition and
a linear functional a'nu_C the same code applies to the projection a'nu_j; the
resulting bounds are valid outer bounds on the sharp set.

All optimisation is exact linear programming over the finite-dimensional
polytope of couplings with prescribed marginals.
"""

from __future__ import annotations

from itertools import product
from typing import Sequence

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, vstack


# --------------------------------------------------------------------------
# 1. Exact identified set for m exchangeable components
# --------------------------------------------------------------------------

def _state_index(m: int, q: int) -> np.ndarray:
    """Return an (m, q**m) int array: component j, state s -> value index."""
    return np.indices((q,) * m, dtype=np.int32).reshape(m, -1)


def identified_set_bounds(
    values: np.ndarray,
    probs: np.ndarray,
    m: int,
    weights: np.ndarray | None = None,
    kind: str = "variance",
    tau: float | None = None,
) -> dict:
    """Sharp lower/upper bound of a functional of nu_C over all couplings.

    Parameters
    ----------
    values : support of the (common) marginal law, length q.
    probs  : marginal probabilities, length q, summing to 1.
    m      : number of pixels aggregated.
    weights: area weights, length m (default equal weights).
    kind   : 'variance' -> Var(nu_C); 'exceedance' -> P(nu_C >= tau).
    tau    : threshold, required for kind='exceedance'.

    Returns lower, upper and the identity of the optimising couplings is not
    tracked; the LP objective values are exact up to solver tolerance.
    """
    values = np.asarray(values, dtype=float)
    probs = np.asarray(probs, dtype=float)
    probs = probs / probs.sum()
    q = len(values)
    if weights is None:
        weights = np.full(m, 1.0 / m)
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()

    idx = _state_index(m, q)                       # (m, S)
    vals = values[idx]                             # (m, S)
    s = weights @ vals                             # (S,) aggregate composition
    mu = float(probs @ values)                     # weights sum to 1

    if kind == "variance":
        c = (s - mu) ** 2
    elif kind == "exceedance":
        if tau is None:
            raise ValueError("tau required for kind='exceedance'")
        c = (s >= tau - 1e-12).astype(float)
    else:
        raise ValueError(kind)

    S = s.shape[0]
    # Equality constraints: marginal of every component equals `probs`, and sum = 1.
    rows, cols, data = [], [], []
    r = 0
    for j in range(m):
        aj = idx[j]
        for i in range(q):
            sel = np.flatnonzero(aj == i)
            rows.append(np.full(sel.size, r))
            cols.append(sel)
            data.append(np.ones(sel.size))
            r += 1
    rows.append(np.full(S, r))
    cols.append(np.arange(S))
    data.append(np.ones(S))
    A = csr_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(r + 1, S),
    )
    b = np.concatenate([np.tile(probs, m), [1.0]])

    lo = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    hi = linprog(-c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if not (lo.success and hi.success):
        raise RuntimeError(f"LP failed: {lo.message} / {hi.message}")
    return {
        "m": m,
        "kind": kind,
        "tau": tau,
        "lower": float(max(0.0, lo.fun if kind == "variance" else lo.fun)),
        "upper": float(min(1.0, -hi.fun) if kind == "exceedance" else -hi.fun),
        "n_states": int(S),
        "mean": mu,
    }


# --------------------------------------------------------------------------
# 2. Closed forms: independence, comonotonic, countermonotonic
# --------------------------------------------------------------------------

def marginal_moments(values, probs):
    values = np.asarray(values, float)
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    mu = float(probs @ values)
    var = float(probs @ (values - mu) ** 2)
    return mu, var


def independent_variance(values, probs, m, weights=None):
    """Var of the area-weighted mean under cross-pixel independence."""
    if weights is None:
        weights = np.full(m, 1.0 / m)
    weights = np.asarray(weights, float)
    weights = weights / weights.sum()
    _, var = marginal_moments(values, probs)
    return float(np.sum(weights ** 2) * var)


def comonotonic_variance(values, probs, m, weights=None):
    """Var of nu_C when the m components are comonotonic (identical).

    For exchangeable components the comonotonic coupling sets all components
    equal, so nu_C = nu_1 and the variance equals the marginal variance.
    This is the largest attainable value for *either* weighting.
    """
    _, var = marginal_moments(values, probs)
    return float(var)


def transport_cov_bounds(values, probs):
    """Sharp bounds on Cov(X, Y) over couplings of two identical marginals.

    Solves  min / max  sum_ij (v_i - mu)(v_j - mu) pi_ij  s.t. marginals = probs.
    The minimiser is the opposite-monotone rearrangement coupling, the
    maximiser the comonotonic coupling (Hoeffding 1940; Frechet 1951).
    """
    values = np.asarray(values, float)
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    q = len(values)
    mu = float(probs @ values)
    dev = values - mu
    cost = np.outer(dev, dev).ravel()
    A_eq = np.zeros((2 * q, q * q))
    b_eq = np.concatenate([probs, probs])
    for i in range(q):
        A_eq[i, i * q:(i + 1) * q] = 1.0            # row sums
        A_eq[q + i, i::q] = 1.0                     # column sums
    lo = linprog(cost, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    hi = linprog(-cost, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    if not (lo.success and hi.success):
        raise RuntimeError("transport LP failed")
    return float(lo.fun), float(-hi.fun)


def pair_variance_bounds(values, probs, m, weights=None):
    """Variance bounds from the pairwise Hoeffding-Frechet bounds.

    Sharp for m = 2; a valid outer bound (superset of the identified set) for
    m >= 3, because pairwise extreme covariances need not be jointly attainable.
    """
    if weights is None:
        weights = np.full(m, 1.0 / m)
    weights = np.asarray(weights, float)
    weights = weights / weights.sum()
    _, var = marginal_moments(values, probs)
    cmin, cmax = transport_cov_bounds(values, probs)
    quad = (weights[:, None] * weights[None, :])
    np.fill_diagonal(quad, 0.0)
    base = float(np.sum(weights ** 2) * var)
    cross = float(quad.sum())
    return base + cross * cmin, base + cross * cmax, cmin, cmax


def psd_variance_bounds(values, probs, m, weights=None):
    """Pairwise Frechet bounds further restricted by covariance-matrix PSDness.

    For exchangeable marginals the correlation matrix ((1-rho) I + rho J) is
    positive semidefinite iff rho >= -1/(m-1).  Intersecting this with the
    Hoeffding-Frechet interval gives a strictly tighter outer bound than
    pair_variance_bounds whenever the Frechet minimum is below -1/(m-1).
    The relaxation is exact exactly when rho = -1/(m-1) is jointly attainable
    (the joint-mixability condition); otherwise the residual gap is the
    joint-mixability gap of the marginal law at that m.
    """
    if weights is None:
        weights = np.full(m, 1.0 / m)
    weights = np.asarray(weights, float)
    weights = weights / weights.sum()
    _, var = marginal_moments(values, probs)
    cmin, cmax = transport_cov_bounds(values, probs)
    rho_fh_min = cmin / var
    rho_fh_max = cmax / var
    rho_psd_min = max(rho_fh_min, -1.0 / (m - 1)) if m > 1 else 1.0
    quad = (weights[:, None] * weights[None, :])
    np.fill_diagonal(quad, 0.0)
    base = float(np.sum(weights ** 2) * var)
    cross = float(quad.sum())
    return (base + cross * rho_psd_min * var,
            base + cross * rho_fh_max * var,
            rho_fh_min, rho_psd_min, var)


def joint_mixability_gap(values, probs, m):
    """Exact LP minimum variance minus the PSD-relaxation minimum variance."""
    ex = identified_set_bounds(values, probs, m, kind="variance")
    lo_psd, _, _, _, _ = psd_variance_bounds(values, probs, m)
    return float(ex["lower"] - lo_psd), float(ex["lower"]), float(lo_psd)


# --------------------------------------------------------------------------
# 3. Landscape simulation (latent Gaussian field, thresholded)
# --------------------------------------------------------------------------

def gaussian_landscape(n, cell, cell_range, seed=0, proportion=0.5):
    """Binary categorical field on an n x n lattice.

    A stationary latent Gaussian field with an isotropic Gaussian covariance of
    standard deviation `cell_range` (in fine cells) is thresholded at the
    quantile that reproduces the requested class proportion.
    """
    rng = np.random.default_rng(seed)
    white = rng.standard_normal((n, n))
    fy = np.fft.fftfreq(n)[:, None]
    fx = np.fft.fftfreq(n)[None, :]
    radius2 = fx ** 2 + fy ** 2
    if cell_range <= 0:
        field = white
    else:
        # Gaussian spectral filter with correlation length `cell_range` cells.
        spectral = np.exp(-2.0 * (np.pi ** 2) * radius2 * (cell_range ** 2))
        field = np.real(np.fft.ifft2(np.fft.fft2(white) * spectral))
    field = (field - field.mean()) / field.std()
    thr = float(np.quantile(field, 1.0 - proportion))
    z = (field > thr).astype(np.int8)
    return z, field


def block_composition_map(z, block):
    """Non-overlapping block composition (occupied fraction) map.

    Rectangular input is handled: each axis is cropped to its own multiple of
    `block`.  On the metric grids built by metric_fields.py no window is square
    (a 28 km north-south window is 11 km east-west at 66.7 N), and the previous
    version cropped both axes by the same count derived from the first axis
    alone, which fails outright on such input.  For a square array the result
    is unchanged.
    """
    ky, kx = z.shape[0] // block, z.shape[1] // block
    zz = z[:ky * block, :kx * block].reshape(ky, block, kx, block)
    return zz.mean(axis=(1, 3))


def pixel_frame(z, cell):
    """Reduce a fine binary field to the pixel level (`cell` x `cell` cells)."""
    return block_composition_map(z, cell)


# --------------------------------------------------------------------------
# 4. Small utilities
# --------------------------------------------------------------------------

def independent_aggregate_law(values, probs, m, weights=None, decimals=12):
    """Law of nu_C when the m pixels are mutually independent.

    Implemented as a dictionary dynamic program over the (small) lattice of
    attainable aggregates rather than a dense q**m tensor.
    """
    values = np.asarray(values, float)
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    if weights is None:
        weights = np.full(m, 1.0 / m)
    weights = np.asarray(weights, float)
    weights = weights / weights.sum()
    state = {0.0: 1.0}
    for j in range(m):
        nxt: dict[float, float] = {}
        for agg, p in state.items():
            for v, pv in zip(values, probs):
                key = round(agg + weights[j] * v, decimals)
                nxt[key] = nxt.get(key, 0.0) + p * pv
        state = nxt
    grid = np.array(sorted(state))
    law = np.array([state[g] for g in grid])
    law = law / law.sum()
    return grid, law


def comonotonic_aggregate_law(values, probs, m, weights=None):
    """Law of nu_C under the comonotonic coupling (all pixels identical)."""
    values = np.asarray(values, float)
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    if weights is None:
        weights = np.full(m, 1.0 / m)
    weights = np.asarray(weights, float)
    weights = weights / weights.sum()
    return values.copy(), probs.copy()          # nu_C = nu_1 exactly


def countermonotonic_pair_law(values, probs):
    """Opposite-monotone coupling of two identical marginals (m = 2)."""
    values = np.asarray(values, float)
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    q = len(values)
    mu = float(probs @ values)
    dev = values - mu
    cost = np.outer(dev, dev).ravel()
    A_eq = np.zeros((2 * q, q * q))
    b_eq = np.concatenate([probs, probs])
    for i in range(q):
        A_eq[i, i * q:(i + 1) * q] = 1.0
        A_eq[q + i, i::q] = 1.0
    res = linprog(cost, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    pi = res.x.reshape(q, q)
    vals = (values[:, None] + values[None, :]).ravel() / 2.0
    p = pi.ravel()
    keep = p > 1e-14
    order = np.argsort(vals[keep])
    return vals[keep][order], p[keep][order] / p[keep].sum()


def exceedance_from_law(law_grid, law_probs, tau):
    return float(law_probs[law_grid >= tau - 1e-12].sum())


def envelope(values, probs, m, taus, weights=None):
    """Sharp bounds on P(nu_C >= tau) plus the named couplings, over a tau grid."""
    rows = []
    ci_grid, ci_p = comonotonic_aggregate_law(values, probs, m, weights)
    ii_grid, ii_p = independent_aggregate_law(values, probs, m, weights)
    if m == 2:
        cc_grid, cc_p = countermonotonic_pair_law(values, probs)
    else:
        cc_grid, cc_p = None, None
    for tau in taus:
        ex = identified_set_bounds(values, probs, m, weights=weights,
                                   kind="exceedance", tau=tau)
        rows.append(dict(
            m=m, tau=float(tau),
            lower_exact=ex["lower"], upper_exact=ex["upper"],
            comonotonic=exceedance_from_law(ci_grid, ci_p, tau),
            independent=exceedance_from_law(ii_grid, ii_p, tau),
            countermonotonic=(
                exceedance_from_law(cc_grid, cc_p, tau) if cc_grid is not None else np.nan),
        ))
    return rows


def bounds_with_exchangeable_dependence(values, probs, m, rho, tau,
                                        weights=None, kind="exceedance"):
    """Sharp bounds when, in addition to the marginals, the pairwise
    cross-pixel correlation rho of the component compositions is known.

    For a linear functional of nu_C the variance is a function of the pairwise
    covariances alone, so at this information level the variance is point
    identified.  Exceedance probabilities are not: this LP returns their sharp
    range over couplings that also reproduce rho.
    """
    values = np.asarray(values, float)
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    q = len(values)
    if weights is None:
        weights = np.full(m, 1.0 / m)
    weights = np.asarray(weights, float)
    weights = weights / weights.sum()
    mu, var = marginal_moments(values, probs)

    idx = _state_index(m, q)
    vals = values[idx]
    s = weights @ vals
    target_cov = rho * var

    rows, cols, data, r = [], [], [], 0
    for j in range(m):
        for i in range(q):
            sel = np.flatnonzero(idx[j] == i)
            rows.append(np.full(sel.size, r)); cols.append(sel)
            data.append(np.ones(sel.size)); r += 1
    b = list(np.tile(probs, m))
    for j in range(m):
        for k in range(j + 1, m):
            rows.append(np.full(s.size, r)); cols.append(np.arange(s.size))
            data.append((vals[j] - mu) * (vals[k] - mu)); r += 1
            b.append(target_cov)
    rows.append(np.full(s.size, r)); cols.append(np.arange(s.size))
    data.append(np.ones(s.size)); r += 1
    b.append(1.0)

    S = s.size
    A = csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(r, S))
    b = np.asarray(b, float)
    if kind == "exceedance":
        c = (s >= tau - 1e-12).astype(float)
    elif kind == "variance":
        c = (s - mu) ** 2
    else:
        raise ValueError(kind)
    lo = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    hi = linprog(-c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if not (lo.success and hi.success):
        raise RuntimeError("LP infeasible or failed for this rho")
    return float(lo.fun), float(-hi.fun)


def bounds_with_pairwise_dependence(values, probs, m, cov_matrix, tau,
                                    weights=None, kind="exceedance"):
    """Sharp bounds when, in addition to the marginals, the full pairwise
    cross-pixel covariance matrix of the component compositions is known.

    Because any linear functional's variance is a function of the marginal
    variances and the pairwise covariances alone (Proposition 2), `kind =
    'variance'` returns a degenerate interval at this information level.  For
    `kind = 'exceedance'` the LP returns the sharp range over all couplings
    that also reproduce `cov_matrix`; a strictly positive width proves that
    pairwise dependence information does not identify threshold events.
    """
    values = np.asarray(values, float)
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    q = len(values)
    if weights is None:
        weights = np.full(m, 1.0 / m)
    weights = np.asarray(weights, float)
    weights = weights / weights.sum()
    mu, var = marginal_moments(values, probs)

    idx = _state_index(m, q)
    vals = values[idx]
    s = weights @ vals

    rows, cols, data, r = [], [], [], 0
    for j in range(m):
        for i in range(q):
            sel = np.flatnonzero(idx[j] == i)
            rows.append(np.full(sel.size, r)); cols.append(sel)
            data.append(np.ones(sel.size)); r += 1
    b = list(np.tile(probs, m))
    for j in range(m):
        for k in range(j + 1, m):
            rows.append(np.full(s.size, r)); cols.append(np.arange(s.size))
            data.append((vals[j] - mu) * (vals[k] - mu)); r += 1
            b.append(float(cov_matrix[j, k]))
    rows.append(np.full(s.size, r)); cols.append(np.arange(s.size))
    data.append(np.ones(s.size)); r += 1
    b.append(1.0)

    S = s.size
    A = csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(r, S))
    b = np.asarray(b, float)
    if kind == "exceedance":
        c = (s >= tau - 1e-12).astype(float)
    elif kind == "variance":
        c = (s - mu) ** 2
    else:
        raise ValueError(kind)
    lo = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    hi = linprog(-c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if not (lo.success and hi.success):
        raise RuntimeError("LP infeasible or failed for this covariance matrix")
    return float(lo.fun), float(-hi.fun), np.abs(A @ lo.x - b).max()


def exchangeable_cov_matrix(values, probs, m, rho):
    _, var = marginal_moments(values, probs)
    C = np.full((m, m), rho * var)
    np.fill_diagonal(C, var)
    return C


def empirical_law(samples, values=None, tol=1e-9):
    """Discretised empirical law of a sample, snapped to `values` if given."""
    samples = np.asarray(samples, float).ravel()
    if values is None:
        values = np.unique(np.round(samples, 12))
    values = np.asarray(values, float)
    idx = np.argmin(np.abs(samples[:, None] - values[None, :]), axis=1)
    probs = np.bincount(idx, minlength=len(values)).astype(float)
    return values, probs / probs.sum()


def aci(lower, upper, point):
    """Aggregation-ambiguity index: share of the range not pinned down."""
    if upper <= 0:
        return np.nan
    return float((upper - lower) / upper)
