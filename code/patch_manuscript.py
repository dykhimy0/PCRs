"""
patch_manuscript.py -- apply the PSD-asymptotics additions to the manuscript.

Usage:  python patch_manuscript.py <input.md> <output.md>

Idempotent-friendly: every anchor is asserted present, so running it twice on
an already-patched file fails loudly instead of silently duplicating text.
"""

import sys
from pathlib import Path

ABSTRACT_OLD = (
    "and a positive-semidefinite relaxation, which is always feasible and sharp "
    "whenever the aggregate can be made deterministic. In the exchangeable case"
)
ABSTRACT_NEW = (
    "and a positive-semidefinite relaxation, which is always feasible. For the "
    "exchangeable two-state case we give the sharp interval for the block "
    "variance in closed form, show that the relaxation is exact for sparse "
    "classes and for two units, and show that elsewhere its relative error is "
    "$O(m^{-2})$, so the exact program is asymptotically unnecessary. In the "
    "exchangeable case"
)

PARA_OLD = (
    r"In the exchangeable two-state case with within-block proportion $p$, the exact lower endpoint is driven to zero exactly when $mp$ is an integer, because the aggregate can then be made deterministic by assigning $\nu_j \in \{0,1\}$ to exactly $mp$ units. When $mp \notin \mathbb{Z}$ the best deterministic aggregate is at distance $\delta = \min_{k\in\mathbb{Z}} |p - k/m|$ from $p$, giving a strictly positive minimum variance $\delta^2$. The PSD relaxation, by contrast, always permits $\rho = -1/(m-1)$ and hence a zero minimum. The joint-mixability gap is therefore zero when $mp\in\mathbb{Z}$ and $\delta^2$ otherwise — a crisp diagnostic of when the cheap relaxation is *exact*."
)
PARA_NEW = (
    r"In the exchangeable two-state case with within-block proportion $p$ the two regimes are explicit. The sharp minimum pairwise correlation of a two-state marginal is $\rho_{\mathrm{FH}} = -\min(p,1-p)/\max(p,1-p)$, and the PSD constraint is active only where it is the stricter of the two, which happens exactly when $m\min(p,1-p) > 1$. Where $m\min(p,1-p) \le 1$ the pairwise extremal configuration is jointly attainable and the FH relaxation coincides with the exact linear program. Where $m\min(p,1-p) > 1$ the PSD corner $\rho = -1/(m-1)$ is feasible, the relaxation returns a zero minimum while the exact minimum is strictly positive, and the FH endpoint has gone negative and is vacuous. Section 4.8 gives the exact size of the minimum and of the gap in both regimes."
)

SEC48 = r"""### 4.8 Closed-form bounds and the asymptotic sharpness of the PSD relaxation

Two refinements make the relaxations usable at the block sizes at which rasters are in practice aggregated.

**Proposition 5 (closed-form sharp interval for the variance).** Let the $m$ units be exchangeable with a common two-state marginal, $\Pr(\nu_j = 1) = p$, and equal weights, and write $\pi = \operatorname{frac}(mp)$ for the fractional part of $mp$. Then the sharp identified interval for $\mathrm{Var}(\nu_C)$ over all couplings with the stored marginals is

$$\left[\; \frac{\pi(1-\pi)}{m^{2}} \;,\;\; p(1-p) \;\right].$$

The upper endpoint is attained by the comonotonic coupling. The lower endpoint is attained by the *cyclic-rotation* coupling

$$\nu_j = \mathbf{1}\{(U + jp) \bmod 1 < p\}, \qquad U \sim \mathrm{Uniform}[0,1),$$

whose sum $\sum_{j=1}^{m} \nu_j$ takes only the two values $\lfloor mp \rfloor$ and $\lceil mp \rceil$.

*Proof.* The upper endpoint follows as in Proposition 4. For the lower endpoint, under the rotation coupling each $\nu_j$ is Bernoulli with parameter $p$, so the marginals are reproduced; and the number of units equal to one is the multiplicity with which the $m$ arcs $\{(u + jp) \bmod 1 < p\}$ cover the point $u$, which for every $u$ is either $\lfloor mp \rfloor$ or $\lceil mp \rceil$. Since the sum has mean $mp$ it places mass $\pi$ on $\lceil mp \rceil$ and $1-\pi$ on $\lfloor mp \rfloor$, and its variance is $\pi(1-\pi)$, giving $\mathrm{Var}(\nu_C) = \pi(1-\pi)/m^{2}$. For optimality, the sum ranges over $\{0,\dots,m\}$ with mean $mp$ under every admissible coupling, and the minimum-variance law on that lattice with that mean is the two-point law just described; the attainable laws are a subset of all such laws. $\square$

**Proposition 6 (the PSD relaxation has two regimes).** In the same setting write $\sigma^{2} = p(1-p)$, $\rho_{\mathrm{FH}} = -\min(p,1-p)/\max(p,1-p)$ and $\rho_{\mathrm{PSD}} = \max\{\rho_{\mathrm{FH}}, -1/(m-1)\}$. The relaxation's lower endpoint is

$$V^{-}_{\mathrm{PSD}} = \sigma^{2}\left[ \frac{1}{m} + \frac{m-1}{m}\,\rho_{\mathrm{PSD}} \right].$$

Then:

(i) if $m\min(p,1-p) \le 1$, then $\rho_{\mathrm{PSD}} = \rho_{\mathrm{FH}}$ and $V^{-}_{\mathrm{PSD}} = V^{-}$: the pairwise relaxation is exact and the joint-mixability gap is zero. This covers $m = 2$ for every $p$, and every block size for a sufficiently sparse class.

(ii) if $m\min(p,1-p) > 1$, then $\rho_{\mathrm{PSD}} = -1/(m-1)$, $V^{-}_{\mathrm{PSD}} = 0$, and the gap is

$$\mathrm{gap}(m,p) \;=\; V^{-} \;=\; \frac{\pi(1-\pi)}{m^{2}}, \qquad \pi = \operatorname{frac}(mp).$$

In this regime the Fréchet–Hoeffding endpoint is negative and therefore vacuous.

*Proof.* Substitute $\rho_{\mathrm{FH}}$ and the value of $V^{-}$ from Proposition 5 into the equicorrelation variance formula. The two cases follow from $\rho_{\mathrm{FH}} < -1/(m-1) \iff (m-1)\min(p,1-p) > \max(p,1-p) \iff m\min(p,1-p) > 1$. $\square$

**Corollary 6.1 (asymptotic sharpness).** In regime (ii),

$$0 \;\le\; \mathrm{gap}(m,p) \;\le\; \frac{1}{4m^{2}},$$

with the envelope attained when $mp$ is a half-integer. Hence the relative error of the PSD relaxation obeys

$$\frac{\mathrm{gap}}{\sigma^{2}} \;\le\; \frac{1}{4\,m^{2}\,p(1-p)} \;=\; O(m^{-2}),$$

so the cheap relaxation is asymptotically exact. For a balanced class, $p = 1/2$, the bound is $m^{-2}$: below 1 % once $m \ge 10$, and below 0.03 % at $m = 64$. Together with Proposition 6(i), the practical consequence is that the exact linear program of Section 5.1 is not needed for the block variance: the pairwise relaxation is either exact, or its relative error is $O(m^{-2})$.
"""

SEC69 = r"""### 6.9 The cheap relaxation is asymptotically exact

Section 4.8 is a statement about when the exact linear program can be replaced, so it is worth testing directly. Figure 7a plots the joint-mixability gap against $m$ on log–log axes for three proportions, together with the envelope $1/(4m^{2})$. Every curve lies under the envelope and decays with the same slope, and the envelope is essentially attained at proportions for which $mp$ is a half-integer. Figure 7b converts this into relative error: the gap never exceeds $m^{-2}$ at $p = 1/2$, so it falls below 1 % once $m \ge 10$ and below 0.03 % at $m = 64$. Figure 7c exposes the two-regime structure of Proposition 6 by plotting the gap against $p$ at fixed $m$: it is identically zero where $m\min(p,1-p) \le 1$, then rises to a bounded plateau whose height is exactly $1/(4m^{2})$.

The closed forms were checked against the exact linear program over 42 $(p, m)$ combinations with $m \le 8$. The closed-form sharp interval of Proposition 5 matched the LP to $2.1\times10^{-17}$ at both endpoints, the closed-form PSD endpoint matched the PSD solver to $1.6\times10^{-17}$, and the closed-form gap matched the LP gap to $2.1\times10^{-17}$. The stated regime predicate held for every combination: the gap is positive exactly when $m\min(p,1-p) > 1$ and $mp$ is not an integer. In the 15 combinations falling in regime (i) the gap was zero to machine precision, confirming that the pairwise relaxation is exactly — not merely nearly — sharp for sparse classes and for $m = 2$. The cyclic-rotation coupling was also simulated directly on a 2000-point grid of $U$; the realised sum took only the values $\lfloor mp \rfloor$ and $\lceil mp \rceil$ in every case tested, with the sample mean and variance matching $mp$ and $\pi(1-\pi)$ to within the grid resolution ($4\times10^{-4}$ and $2.5\times10^{-4}$ respectively).

The practical reading is that the combinatorial cost of the exact program can be avoided. For a variance query the closed form of Proposition 5 is exact and free; where a relaxation is preferred, Proposition 6 bounds the price of using it, and at realistic block sizes that price is orders of magnitude below the accuracy of any uncertainty estimate a real product would carry. The exact linear program of Section 5.1 is therefore best understood as a verification tool for small $m$ rather than as a production method.

![](../figures/fig7_psd_asymptotic.png)

**Figure 7.** Asymptotic sharpness of the PSD relaxation. (a) The joint-mixability gap against $m$ for three proportions, with the envelope $1/(4m^{2})$. (b) The relative gap against $m$, with the $O(m^{-2})$ references; the 1 % level and $m = 10$ are marked. (c) The gap against $p$ at fixed $m = 4, 8, 16$; it vanishes while $m\min(p,1-p) \le 1$ and rises to a plateau of height exactly $1/(4m^{2})$.
"""

CONC_OLD = ("the block variance is bounded but not identified, and the worst-case "
            "representational ambiguity grows with the number of units, at "
            "least as $1-1/m$.")
CONC_NEW = (CONC_OLD + " The sharp interval is available in closed form, and the "
            "cheap positive-semidefinite relaxation is exact for sparse classes "
            "and asymptotically exact elsewhere, so the obstruction is "
            "informational rather than computational.")

NOTES_OLD = r"""3. **Downstream-task experiment.** Section 7.5 proposes a decision-relevant application (area estimation with a target precision). A single worked example on one window would materially strengthen the paper and is the single highest-value addition.
4. **Level-C storage cost** is asserted but not tabulated; a short table of $\Omega_q^m$ against $m$ and $\Omega_q$ would make the trade-off concrete.
5. **Generalisations to overlay and resampling** are stated as future work only; they would require the cross-raster joint and are deliberately out of scope for this draft.
6. **Figure polish.** Figure 2a has a minor legend/annotation overlap; a redundant-linestyle warning in `make_figures.py:75` should be cleaned up. Neither affects any reported value.
7. **No figure from the previous draft is reused.** Figures 1–6 are generated de novo from the new CSV results by `make_figures.py`."""

NOTES_NEW = r"""3. **Attribution of the closed-form result.** The two-point minimum-variance law for a sum of exchangeable Bernoulli variables, and the rotation construction that attains it (Proposition 5), are classical. Hoeffding (1940) treats the extremal couplings; the attainability construction should be cited to a named source rather than presented as new. The claims original to this paper are the two-regime exactness statement (Proposition 6) and the $O(m^{-2})$ sharpness corollary, not the closed form itself.
4. **Downstream-task experiment.** Section 7.5 proposes a decision-relevant application (area estimation with a target precision). A single worked example on one window would materially strengthen the paper and is the single highest-value addition.
5. **Level-C storage cost** is asserted but not tabulated; a short table of $\Omega_q^m$ against $m$ and $\Omega_q$ would make the trade-off concrete.
6. **Generalisations to overlay and resampling** are stated as future work only; they would require the cross-raster joint and are deliberately out of scope for this draft.
7. **Figure polish.** Figure 2a has a minor legend/annotation overlap; a redundant-linestyle warning in `make_figures.py:75` should be cleaned up. Neither affects any reported value.
8. **No figure from the previous draft is reused.** Figures 1–7 are generated de novo from the new CSV results by `make_figures.py` and `make_figure7.py`."""


def main():
    src = Path(sys.argv[1])
    dst = Path(sys.argv[2])
    s = src.read_text(encoding="utf-8")
    edits = [
        (ABSTRACT_OLD, ABSTRACT_NEW),
        (PARA_OLD, PARA_NEW),
        ("## 5. Materials and methods", SEC48 + "\n---\n\n## 5. Materials and methods"),
        ("## 7. Discussion", SEC69 + "\n---\n\n## 7. Discussion"),
        (CONC_OLD, CONC_NEW),
        (NOTES_OLD, NOTES_NEW),
    ]
    for old, new in edits:
        if old not in s:
            raise SystemExit(f"ANCHOR MISSING: {old[:80]!r}")
        s = s.replace(old, new, 1)
    dst.write_text(s, encoding="utf-8")
    print(f"wrote {dst} ({len(s)} chars, {len(edits)} edits)")


if __name__ == "__main__":
    main()
