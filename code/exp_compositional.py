"""
exp_compositional.py -- genuine K >= 3 compositional experiments.

Why this file exists
--------------------
The spatial and real-data experiments in exp_spatial.py / exp_real.py reduce the
raster to a single binary indicator.  The paper, however, claims to be about
*compositional* rasters, so that reduction has to be paid back.  Here the stored
object is a full composition nu_j on the simplex Delta^{K-1} with K = 3.

Design
------
For a *linear* functional a'nu_C the aggregate is a weighted sum of scalar
random variables with known marginals, so the scalar machinery applies verbatim
with a discrete (multi-state) marginal in place of the Bernoulli one.  For a
*nonlinear* functional (a diversity index) the multivariate structure is what
matters, and the coupling LP must be solved on the composition lattice directly.

Both are done:

  Part 1  controlled K = 3 lattice (3 levels per class -> 10 composition states),
          m = 2..5.  Functionals: a class share (linear), the Simpson diversity
          1 - sum_k s_k^2 (nonlinear), and the majority share max_k s_k
          (nonlinear).  Sharp LP bounds on E[phi(nu_C)] and on P(phi >= tau),
          against the independent and comonotonic couplings.

  Part 2  real ESA WorldCover, K = 3 (tree / crop / other) at 30 m, on the same
          3-level lattice, m = 2..4.  The truth is enumerated over every block in
          the window, so the comparison is against ground truth, not a model.

Cost note: the coupling polytope has q**m entries.  With q = 10 the exact LP is
used up to m = 5 (10^5); the guard skips anything larger.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results_v2"
DATA = ROOT / "data"
RES.mkdir(exist_ok=True)

MAX_STATES = 500_000          # guard on q**m


# --------------------------------------------------------------------------
# composition lattice
# --------------------------------------------------------------------------
def compose_states(K: int, levels: int) -> np.ndarray:
    """All compositions of `levels` into K parts, scaled to the simplex."""
    out: list[tuple[int, ...]] = []

    def rec(k: int, rem: int, cur: list[int]) -> None:
        if k == K - 1:
            out.append(tuple(cur + [rem]))
            return
        for i in range(rem + 1):
            rec(k + 1, rem - i, cur + [i])

    rec(0, levels, [])
    return np.array(out, float) / levels


def phi_simpson(s: np.ndarray) -> np.ndarray:
    """1 - sum_k s_k^2."""
    return 1.0 - (s ** 2).sum(axis=-1)


def phi_majority(s: np.ndarray) -> np.ndarray:
    """max_k s_k."""
    return s.max(axis=-1)


def phi_share(s: np.ndarray, k: int = 0) -> np.ndarray:
    return s[..., k]


# --------------------------------------------------------------------------
# exact coupling LP for E[phi(nu_C)] and P(phi(nu_C) >= tau)
# --------------------------------------------------------------------------
def coupling_lp(states: np.ndarray, probs: np.ndarray, m: int,
                objective: np.ndarray):
    """Sharp min/max of sum_states obj(state) * P(state) over couplings.

    `objective` is indexed by the flattened m-tuple of unit states.
    """
    q = len(states)
    idx = np.indices((q,) * m, dtype=np.int32).reshape(m, -1)   # (m, q**m)
    S = idx.shape[1]

    rows, cols, data, r = [], [], [], 0
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
    A = csr_matrix((np.concatenate(data),
                    (np.concatenate(rows), np.concatenate(cols))),
                   shape=(r + 1, S))
    b = np.concatenate([np.tile(probs, m), [1.0]])

    lo = linprog(objective, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    hi = linprog(-objective, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if not (lo.success and hi.success):
        raise RuntimeError(f"LP failed: {lo.message} / {hi.message}")
    return float(lo.fun), float(-hi.fun)


def aggregate_compositions(states: np.ndarray, m: int) -> np.ndarray:
    """Per-tuple weighted aggregate for equal weights (m, q**m) -> (q**m, K)."""
    idx = np.indices((len(states),) * m, dtype=np.int32).reshape(m, -1)
    return states[idx].mean(axis=0)


# --------------------------------------------------------------------------
# reference couplings
# --------------------------------------------------------------------------
def independent_expectation(states: np.ndarray, probs: np.ndarray, m: int,
                            phi) -> float:
    """E[phi(nu_C)] under the independent coupling, exactly, by DP on counts.

    Every unit state lies on the grid {0,1/levels,...}; the sum of the m unit
    states is therefore a lattice point n/L, and the DP runs over count vectors
    n with sum(n) = m * L.  Exact and cheap.
    """
    K = states.shape[1]
    L = int(round(1.0 / np.min(np.diff(np.unique(states[:, 0]))))) \
        if states.shape[1] == 1 else None
    # infer the grid denominator from the states themselves
    denom = 0
    for s in states:
        for v in s:
            denom = max(denom, _rational_denominator(v))
    L = denom
    counts = np.round(states * L).astype(int)

    dist: dict[tuple[int, ...], float] = {tuple(np.zeros(K, int)): 1.0}
    for _ in range(m):
        nxt: dict[tuple[int, ...], float] = {}
        for key, p in dist.items():
            for c, pv in zip(counts, probs):
                nk = tuple(key[k] + c[k] for k in range(K))
                nxt[nk] = nxt.get(nk, 0.0) + p * pv
        dist = nxt
    tot = 0.0
    for key, p in dist.items():
        s = np.array(key, float) / (m * L)
        tot += p * float(phi(s[None, :])[0])
    return tot


def _rational_denominator(v: float, maxden: int = 64) -> int:
    for d in range(1, maxden + 1):
        if abs(v * d - round(v * d)) < 1e-9:
            return d
    return maxden


def comonotonic_expectation(states: np.ndarray, probs: np.ndarray, phi) -> float:
    return float((probs * phi(states)) .sum())


# --------------------------------------------------------------------------
# Part 1 -- controlled K = 3 lattice
# --------------------------------------------------------------------------
def controlled(K: int = 3, L: int = 3, m_values=(2, 3, 4, 5)):
    states = compose_states(K, L)
    # a moderately spread marginal on the lattice: Dirichlet(1.4, 1.0, 0.7)
    # evaluated on the lattice points and renormalised (deterministic, no RNG).
    alpha = np.array([1.4, 1.0, 0.7])
    logp = ((alpha - 1.0)[None, :] * np.log(np.clip(states, 1e-12, None))).sum(axis=1)
    probs = np.exp(logp - logp.max())
    probs = probs / probs.sum()

    rows = []
    for m in m_values:
        if len(states) ** m > MAX_STATES:
            print(f"  skip K={K} L={L} m={m}: {len(states)**m:.1e} states")
            continue
        agg = aggregate_compositions(states, m)
        for name, phi in (("simpson", phi_simpson),
                          ("majority", phi_majority)):
            vals = phi(agg)
            lo, hi = coupling_lp(states, probs, m, vals)
            ind = independent_expectation(states, probs, m, phi)
            com = comonotonic_expectation(states, probs, phi)
            rows.append(dict(
                K=K, levels=L, n_states=len(states), m=m,
                functional=name, n_coupling_states=len(states) ** m,
                mean_lower=lo, mean_upper=hi, mean_width=hi - lo,
                independent=ind, comonotonic=com,
                rel_width=(hi - lo) / hi if hi > 0 else np.nan,
                independent_inside=bool(lo - 1e-9 <= ind <= hi + 1e-9),
                comonotonic_inside=bool(lo - 1e-9 <= com <= hi + 1e-9),
            ))
        # exceedance of the diversity index
        agg_s = phi_simpson(agg)
        for tau in (0.4, 0.5, 0.6):
            obj = (agg_s >= tau - 1e-12).astype(float)
            lo, hi = coupling_lp(states, probs, m, obj)
            rows.append(dict(
                K=K, levels=L, n_states=len(states), m=m,
                functional=f"P(simpson>={tau})", n_coupling_states=len(states) ** m,
                mean_lower=lo, mean_upper=hi, mean_width=hi - lo,
                independent=np.nan, comonotonic=np.nan,
                rel_width=(hi - lo) / hi if hi > 0 else np.nan,
                independent_inside=np.nan, comonotonic_inside=np.nan,
            ))

    # class share: linear functional, variance bounds on the projected marginal
    share_states = np.round(states[:, 0] * L).astype(int) / L
    uniq, cnt = np.unique(share_states, return_counts=True)
    p_share = np.array([probs[share_states == u].sum() for u in uniq])
    for m in m_values:
        q = len(uniq)
        if q ** m > MAX_STATES:
            continue
        # Var of the aggregate share: objective (s - E s)^2, E s identified
        Es = float(p_share @ uniq)
        idx = np.indices((q,) * m, dtype=np.int32).reshape(m, -1)
        s = uniq[idx].mean(axis=0)
        obj = (s - Es) ** 2
        lo, hi = coupling_lp(uniq[:, None], p_share, m, obj)
        rows.append(dict(
            K=K, levels=L, n_states=len(states), m=m,
            functional="Var(share k=0)", n_coupling_states=q ** m,
            mean_lower=lo, mean_upper=hi, mean_width=hi - lo,
            independent=float(np.sum(np.full(m, 1.0 / m) ** 2)
                              * float(p_share @ (uniq - Es) ** 2)),
            comonotonic=float(p_share @ (uniq - Es) ** 2),
            rel_width=(hi - lo) / hi if hi > 0 else np.nan,
            independent_inside=bool(lo - 1e-9
                                    <= np.sum(np.full(m, 1.0 / m) ** 2)
                                    * float(p_share @ (uniq - Es) ** 2)
                                    <= hi + 1e-9),
            comonotonic_inside=True,
        ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Part 2 -- real WorldCover, K = 3
# --------------------------------------------------------------------------
def _lattice_counts(n1: int, n2: int, n3: int, L: int):
    """Largest-remainder projection of a cell-count triple onto the 1/L lattice.

    The triple sums to 9 (fine cells); the result sums to L, so the projected
    composition lies exactly on the composition lattice.
    """
    total = n1 + n2 + n3
    x = np.array([n1, n2, n3], float) * (L / total)
    fl = np.floor(x).astype(int)
    rem = L - int(fl.sum())
    if rem > 0:
        order = np.argsort(-(x - fl))
        for i in range(rem):
            fl[order[i % 3]] += 1
    elif rem < 0:
        order = np.argsort(x - fl)
        for i in range(-rem):
            fl[order[i % 3]] -= 1
    return int(fl[0]), int(fl[1]), int(fl[2])


def real_compositional(window_name: str = "worldcover_N48E006_26000_26000_3072.npy",
                       L: int = 3, m_values=(2, 3, 4, 5), pixel: int = 3,
                       max_states: int = MAX_STATES):
    """Real WorldCover, K=3, on one named window.

    The window is named rather than indexed.  v9 selected it as
    ``sorted(DATA.glob("worldcover_*.npy"))[2]``, which was the Black Forest
    while the reference set had six windows; a seventh window added in the
    eleven-window revision would have silently moved index 2 to a different
    biome and changed the figure without a word.
    """
    path = DATA / window_name
    if not path.exists():
        raise FileNotFoundError(path)
    arr = np.load(path)
    n = arr.shape[0]
    k = n // pixel
    blocks = arr[:k * pixel, :k * pixel].reshape(k, pixel, k, pixel)
    n_tree = (blocks == 10).sum(axis=(1, 3)).astype(int).ravel()
    n_crop = (blocks == 40).sum(axis=(1, 3)).astype(int).ravel()
    n_cell = pixel * pixel
    tree = (n_tree.reshape(k, k) / n_cell).astype(float)
    crop = (n_crop.reshape(k, k) / n_cell).astype(float)

    # project every *distinct* cell-count triple once, then weight by its count
    trip = n_tree * 100 + n_crop
    uniq_trip, trip_cnt = np.unique(trip, return_counts=True)
    states = compose_states(3, L)
    key = {tuple(np.round(v, 9)): i for i, v in enumerate(states)}
    counts = np.zeros(len(states))
    for t, c in zip(uniq_trip, trip_cnt):
        a, b = divmod(int(t), 100)
        c1, c2, c3 = _lattice_counts(a, b, n_cell - a - b, L)
        counts[key[tuple(np.round(np.array([c1, c2, c3]) / L, 9))]] += c
    probs = counts / counts.sum()

    def block_truth(m):
        bb = k // m
        if bb < 2:
            return np.nan, np.nan, 0
        t = tree[:bb * m, :bb * m].reshape(bb, m, bb, m).mean(axis=(1, 3))
        c = crop[:bb * m, :bb * m].reshape(bb, m, bb, m).mean(axis=(1, 3))
        cb = np.stack([t, c, 1.0 - t - c], axis=-1).reshape(-1, 3)
        d = phi_simpson(cb)
        return float(d.mean()), float(d.std()), int(d.size)

    rows = []
    for m in m_values:
        if len(states) ** m > max_states:
            print(f"  skip real L={L} m={m}: {len(states)**m:.1e} states")
            continue
        agg = aggregate_compositions(states, m)
        vals = phi_simpson(agg)
        lo, hi = coupling_lp(states, probs, m, vals)
        ind = independent_expectation(states, probs, m, phi_simpson)
        com = comonotonic_expectation(states, probs, phi_simpson)
        tmean, tsd, tn = block_truth(m)
        rows.append(dict(
            window=path.name, K=3, levels=L,
            n_states_used=int((counts > 0).sum()), n_states_lattice=len(states),
            m=m, blocks_compared=tn, functional="simpson",
            truth=tmean, truth_sd=tsd, lower=lo, upper=hi, width=hi - lo,
            independent=ind, comonotonic=com,
            truth_inside=bool(lo - 1e-6 <= tmean <= hi + 1e-6),
            width_over_span=((hi - lo) / hi if hi > 0 else np.nan),
            independent_gap=abs(ind - tmean),
        ))
    return pd.DataFrame(rows), probs, states, tree, crop


def main():
    print("== Part 1: controlled K=3 lattice ==")
    df1 = controlled()
    df1.to_csv(RES / "compositional_controlled.csv", index=False)
    print(df1.to_string(index=False))

    frames = []
    for L, ms in ((9, (2, 3)), (3, (2, 3, 4, 5))):
        print(f"\n== Part 2: real WorldCover, K=3, 30 m, L={L} ==")
        df2, probs, states, tree, crop = real_compositional(L=L, m_values=ms)
        df2.to_csv(RES / f"compositional_real_L{L}.csv", index=False)
        print(df2.to_string(index=False))
        frames.append(df2)
        print(f"  lattice states occupied: {int((probs > 0).sum())} "
              f"of {len(states)}")
    pd.concat(frames).to_csv(RES / "compositional_real.csv", index=False)

    print("\n== real pixel composition moments (K=3, 30 m) ==")
    print(f"  tree  mean {tree.mean():.4f}  var {tree.var():.5f}")
    print(f"  crop  mean {crop.mean():.4f}  var {crop.var():.5f}")
    print(f"  other mean {1 - tree.mean() - crop.mean():.4f}  "
          f"var {(1 - tree - crop).var():.5f}")
    print("wrote compositional_controlled.csv compositional_real_L9.csv "
          "compositional_real_L3.csv compositional_real.csv")


if __name__ == "__main__":
    main()
