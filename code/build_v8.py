"""
build_v8.py -- build manuscript_draft_v8.md from v7.

Applies the reviewer checklist in one deterministic pass:
  * IJGIS conventions (5 keywords, Data and codes availability statement,
    Disclosure / Funding / Notes on contributors, explicit contribution list)
  * research-gap paragraph, related-work comparison table, plain-language
    sentence after every key proposition, discussion split, conclusion rewrite
  * PCR abbreviation, glosses for identified set / joint-mixability,
    ACI recap at the head of the results
  * new Section 6.8 (estimator proof of concept) and its figure
  * captions rewritten so the first sentence states the result
  * figure renumbering for the two new front-matter figures
  * DOI verification fixes and DOI completion for every verified entry
  * Appendix A (notation) and Appendix B (proofs); draft notes removed
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "manuscript" / "manuscript_draft_v7.md"
DST = ROOT / "manuscript" / "manuscript_draft_v8.md"

text = SRC.read_text(encoding="utf-8")
_log: list[str] = []


def swap(old: str, new: str, count: int = 1) -> None:
    global text
    n = text.count(old)
    if n < count:
        raise SystemExit(f"ANCHOR FAIL ({n} found, need {count}):\n---\n{old[:300]}\n---")
    text = text.replace(old, new, count)
    _log.append(f"ok  swap  {old[:64]!r}")


def splice_after(anchor: str, block: str) -> None:
    global text
    if text.count(anchor) != 1:
        raise SystemExit(f"ANCHOR FAIL (splice_after, {text.count(anchor)} found):\n{anchor[:200]}")
    text = text.replace(anchor, anchor + block, 1)
    _log.append(f"ok  after {anchor[:64]!r}")


# ---------------------------------------------------------------- 0. draft notes
draft_start = text.index("## Draft notes for the authors")
text = text[:draft_start].rstrip() + "\n"
_log.append("ok  removed draft-notes section")

# ---------------------------------------------------------------- 1. preamble
swap("**Draft v0.4 — working manuscript, not for circulation**",
     "**Draft v0.5 — working manuscript, not for circulation**")

ABS_START = "Geographic data models increasingly attach"
ABS_END = "shrinks but does not close the gap."
i0 = text.index(ABS_START)
i1 = text.index(ABS_END) + len(ABS_END)

NEW_ABSTRACT = r"""Geographic data models increasingly attach a distribution, rather than a single value, to each pixel: a pixel carries an uncertain *composition*, not merely an uncertain class. Yet the operations that make a raster useful — aggregation to coarser support, overlay, resampling, and area estimation — are defined on compositions, not on the distributions that encode them. We ask a question that current probabilistic sub-pixel models do not answer: **given the information a probabilistic-compositional raster (PCR) actually stores, what can and cannot be identified about a coarser unit derived from it?** We answer it as a problem of partial identification, not of estimation. We formalise the stored object as the marginal law of each fine unit's composition, the derived object as the law of an area-weighted aggregate, and three stored levels: (A) marginals only, (B) marginals plus the pairwise covariance matrix, (C) the full joint law. The aggregate mean is always identified. The aggregate variance is identified exactly at level B and not at level A, and exact identification needs only pair information when two units are combined but is genuinely partial for three or more. Exceedance queries are sharper still in their demands: their extremal dependence structure changes with the threshold, so no single coupling stands in for the unknown joint. The sharp interval is available in closed form for the two-state case, and a cheap positive-semidefinite relaxation is exact for sparse classes and asymptotically exact elsewhere, with relative error $O(m^{-2})$. We verify every closed form against an exact linear program, and we test the theory on six windows of the ESA WorldCover 10 m product spanning temperate, Mediterranean and boreal biomes. A per-pixel-only representation understates block variance by a factor of 3.5 to 19.8 and degrades a nominal 90 % interval to 3–7 % empirical coverage. A variogram resolved over an adequate lag window recovers the block variance to within 3.6 % in all 24 cases, and the residual bias of a truncated lag window is set by the fitted range rather than by the biome, so a $\pm1$ percentage-point target on a regional mean costs 3.5 to 19.8 times more samples than independence implies. A few dozen training blocks contract the identified interval to a usable estimate. The remedy is therefore representational rather than algorithmic: a PCR is only as useful as the cross-unit dependence information it stores, so store that dependence, or store the sharp bounds and say that they are bounds."""

text = text[:i0] + NEW_ABSTRACT + text[i1:]
_log.append("ok  replaced abstract (tightened, PCR abbreviation introduced)")

swap("**Keywords:** probabilistic composition; change of support; partial identification; compositional data; uncertainty representation; land cover; Fréchet–Hoeffding bounds; copula",
     "**Keywords:** change of support; partial identification; uncertainty representation; compositional data; land cover")
_log.append("ok  keywords trimmed to 5 (IJGIS convention)")

# ---------------------------------------------------------------- 2. introduction
swap("The paper is organised as follows.",
     r"""The contributions are threefold. First, we pose change of support as an identifiability problem for a PCR and characterise exactly what each stored level determines (Propositions 1–3 and Corollary 2.1). Second, we give sharp and tractable answers — a closed form in the two-state case, an asymptotically exact relaxation in general — and show that the cheap relaxation fails in a specific, diagnosable way (Propositions 4–6 and Corollary 6.1). Third, we test the theory against exact linear programs, synthetic landscapes with known truth, and six real reference windows across three biomes, and we translate the result into the quantity a practitioner actually pays: samples.

Figure 1 summarises the route the paper takes from that question to those findings.

![](../figures_v2/fig0_roadmap.png)

**Figure 1.** The paper reduces scale-consistency to an identifiability question and then answers it in four steps: formalisation, theory, experiments, and a decision-level consequence.

The paper is organised as follows.""")

swap("We are careful about novelty. Bayesian spectral unmixing (Eches et al., 2010), spatial compositional models (Pirzamanbein et al., 2018), and change-of-support frameworks in geostatistics all address *estimation* of sub-pixel composition or its uncertainty.",
     r"""**Research gap.** Existing work supplies either a better estimate or a better propagation machinery. Bayesian spectral unmixing and spatial compositional models estimate the per-pixel composition and its uncertainty (Eches et al., 2010; Pirzamanbein et al., 2018), and geostatistical change of support propagates that uncertainty through a fitted variogram (Heuvelink & Pebesma, 1999; Gotway & Young, 2002; Wadoux & Heuvelink, 2023). Both lines assume that the dependence between units is available, and both therefore answer a question about a *model*. Neither asks the prior question that a data product must answer: given only what the raster actually stores, what can be said at all about a coarser unit? That is a question of identifiability, and it is the question this paper answers.

We are careful about novelty. Bayesian spectral unmixing (Eches et al., 2010), spatial compositional models (Pirzamanbein et al., 2018), and change-of-support frameworks in geostatistics all address *estimation* of sub-pixel composition or its uncertainty.""")

# ---------------------------------------------------------------- 3. related work
swap("We import this stance into change of support and treat the aggregate law as partially identified.",
     r"""We import this stance into change of support and treat the aggregate law as partially identified. The identified set of a functional is the range of values it can take under every coupling consistent with the stored information; it is the honest object when the data do not pin the functional down.

### 2.5 Why the geographic literature has not asked this question

The four literatures we build on ask different questions, and the identification question falls between them. Sub-pixel estimation asks what the composition *is* and delivers a per-pixel posterior. Geostatistical change of support asks what the block mean *is* and delivers a model-based variance. Area-to-point kriging (Kyriakakis, 2004) asks how to move a coarse map to a fine support and delivers a fine map. None of these asks whether the stored object determines the answer, because none of them treats the stored object as fixed: in each, the dependence between units is supplied by the modeller, so the question of what happens without it never arises. The question becomes visible only when the raster is read as a data product whose metadata are already written, which is precisely the partial-identification stance imported from econometrics (Molinari, 2020; Kline & Tamer, 2023). Table 2.1 makes the contrast explicit.

| Literature | What it estimates | What it assumes | What the raster stores | What it can say at a coarser support |
|---|---|---|---|---|
| Sub-pixel composition (Eches et al., 2010; Pirzamanbein et al., 2018) | per-pixel composition and its posterior | a spectral or spatial model | per-pixel marginal law | the aggregate mean, under an assumed dependence |
| Geostatistical change of support (Matheron, 1963; Wadoux & Heuvelink, 2023) | block mean and its variance | a variogram | a fitted model, not marginals | a point estimate with a model-based variance |
| Area-to-point kriging (Kyriakidis, 2004) | a fine map from coarse data | a variogram and an area-to-point operator | a coarse map plus a variogram | a fine map, model-based |
| **This paper** | what the stored information *determines* | none beyond the stored marginals | the per-pixel marginal law | a sharp identified set, with no model |

*Table 2.1. The four literatures side by side. The first three supply a dependence structure and then estimate; this paper takes the dependence as absent and characterises what remains.*""")

swap("The essential fact we exploit — that the Fréchet–Hoeffding pairwise endpoints cannot all hold at once in general, so the pairwise relaxation overshoots — is the \"joint mixability\" phenomenon.",
     r"""The essential fact we exploit — that the Fréchet–Hoeffding pairwise endpoints cannot all hold at once in general, so the pairwise relaxation overshoots — is the *joint-mixability* phenomenon: the pairwise extremal dependence structures are not simultaneously realisable by any single coupling.""")

# ---------------------------------------------------------------- 4. problem formulation
swap("Level C is the ideal that a probabilistic-compositional raster should, in principle, approximate. The identifiability question is:",
     r"""Level C is the ideal that a PCR should, in principle, approximate. The identifiability question is:""")

splice_after("**which functionals of $\\Lambda_C$ are the same for every $\\Gamma$ consistent with the stored information at each level?**",
             r"""

The set of values a functional can take under all such couplings is its *identified set*: the range the stored data leave possible, as opposed to the single point an estimator returns. Figure 2 draws the gap and the ladder.

![](../figures_v2/fig0_schematic.png)

**Figure 2.** Storing marginals is not the same as storing a joint. (a) Fine units carry marginal laws; the aggregate law needs the joint, which the raster does not store. (b) The metadata ladder: what each stored level lets one answer for the mean, the variance and an exceedance probability.""")

# ---------------------------------------------------------------- 5. theory glosses
swap("Everything above specialises to a two-class indicator only for exposition.",
     r"""In words: a per-pixel distribution, however rich, is still a per-pixel statement, and no amount of per-pixel richness substitutes for the cross-unit dependence.

Everything above specialises to a two-class indicator only for exposition.""")

swap("**Proposition 5 (closed-form sharp interval for the variance).** Let the $m$ units be exchangeable",
     r"""**Proposition 5 (closed-form sharp interval for the variance; Hoeffding, 1940; Rüschendorf, 1984).** Let the $m$ units be exchangeable""")

splice_after("whose sum $\\sum_{j=1}^{m} \\nu_j$ takes only the two values $\\lfloor mp \\rfloor$ and $\\lceil mp \\rceil$.",
             r"""

In words: the sharp interval is two explicit numbers, and the lower one is attained by spreading the ones as evenly as the marginals allow. The proof is deferred to Appendix B.""")

swap("*Proof.* The upper endpoint follows as in Proposition 4. For the lower endpoint, under the rotation coupling each $\\nu_j$ is Bernoulli with parameter $p$, so the marginals are reproduced; and the number of units equal to one is the multiplicity with which the $m$ arcs $\\{(u + jp) \\bmod 1 < p\\}$ cover the point $u$, which for every $u$ is either $\\lfloor mp \\rfloor$ or $\\lceil mp \\rceil$. Since the sum has mean $mp$ it places mass $\\pi$ on $\\lceil mp \\rceil$ and $1-\\pi$ on $\\lfloor mp \\rfloor$, and its variance is $\\pi(1-\\pi)$, giving $\\mathrm{Var}(\\nu_C) = \\pi(1-\\pi)/m^{2}$. For optimality, the sum ranges over $\\{0,\\dots,m\\}$ with mean $mp$ under every admissible coupling, and the minimum-variance law on that lattice with that mean is the two-point law just described; the laws the coupling can actually produce form a subset of all such laws. $\\square$",
     "*Proof.* See Appendix B. $\\square$")

splice_after("In this regime the Fréchet–Hoeffding endpoint is negative and therefore unusable.",
             r"""

In words: either the pairwise relaxation is exactly right, or it reports zero where the truth is a small positive number — and the crossover between the two regimes is a single inequality, $m\min(p,1-p) > 1$, that a practitioner can evaluate before using it.""")

swap("*Proof.* Substitute $\\rho_{\\mathrm{FH}}$ and the value of $V^{-}$ from Proposition 5 into the equicorrelation variance formula. The two cases follow from $\\rho_{\\mathrm{FH}} < -1/(m-1) \\iff (m-1)\\min(p,1-p) > \\max(p,1-p) \\iff m\\min(p,1-p) > 1$. $\\square$",
     "*Proof.* See Appendix B. $\\square$")

splice_after("so the cheap relaxation is asymptotically exact.",
             r"""

In words: the price of avoiding the exact program is a relative error that decays like the square of the block size, far below the accuracy of any real product's per-pixel uncertainty.""")

# ---------------------------------------------------------------- 6. results
swap("## 6. Results\n",
     r"""## 6. Results

Throughout this section, **ACI(m)** denotes the *aggregation-ambiguity index* defined in Section 4.5: the width of the identified interval for the block variance, normalised by its sharp upper endpoint. ACI = 0 means the variance is identified; ACI = 1 means the stored marginals leave it completely undetermined. Widths of exceedance intervals are reported in absolute probability units.
""")

# --- 6.7 split: retitle and promote the exactness result
swap("### 6.7 One stored dependence parameter is worth a lot, and still not enough",
     r"""### 6.7 One stored dependence parameter identifies the block variance exactly""")

swap("A natural mitigation is to store a single dependence parameter alongside the marginals: an exchangeable copula with correlation $\\rho$, a $\\mathcal{O}(1)$-size addition to the metadata.",
     r"""The most important positive result of the paper is easy to state. At level A the identified interval for the block variance has width 0.184–0.226 on the synthetic landscapes, that is, essentially the whole range the variance could take. Adding a *single* dependence parameter — the exchangeable correlation $\rho$ — reduces that width to exactly zero: the block variance becomes point-identified, as Corollary 2.1 requires, because level B is exactly one stored scalar in the exchangeable case. This is the concrete metadata contract a data producer can offer.

A natural mitigation is therefore to store a single dependence parameter alongside the marginals: an exchangeable copula with correlation $\rho$, a $\mathcal{O}(1)$-size addition to the metadata.""")

splice_after("A single parameter is therefore useful only if it is estimated honestly; a wrong one produces false confidence.",
             r"""

What the scalar does **not** buy is any threshold query. Exceedance intervals narrow by 54–74 % but never close, and the residual width is not a modelling residual: it is the higher-order structure that one pair parameter cannot see. Section 6.8 asks whether that parameter can at least be *estimated* from data, which is the missing practical half of this result.

![](../figures_v2/fig11_estimator.png)

**Figure 9.** A few training blocks turn the identified interval into a usable estimate. (a) The 90th-percentile relative error of the plug-in block variance against the number of training blocks K, at m = 16, for the six reference windows; the dotted line is the level-A band width and the dashed line the mean independence error. (b) The contraction factor, that is, the level-A band width divided by the interquartile range of the estimate.

### 6.8 Shrinking the interval from a few training blocks

The bounds are honest but wide, and a reader may ask whether the missing dependence can be recovered from data rather than declared. This subsection is a proof of concept, not a new method: it asks how many observed blocks are needed to replace the identified interval with a usable point estimate.

The estimator is the obvious one. Given K training blocks with observed composition, the plug-in estimate of the block variance is the sample variance of the K block means,

$$\widehat{V}(m) \;=\; \frac{1}{K-1}\sum_{k=1}^{K}\big(\bar{\nu}_k-\bar{\nu}\big)^{2},$$

which is unbiased for the true block variance, large-scale trend included. The implied exchangeable dependence parameter is $\hat\rho = (m\widehat V/\sigma^{2}-1)/(m-1)$; it is reported for reference but deliberately not re-entered into the variance, because on a trended field it conflates within-block correlation with regional structure and can exceed one. It reads as a total dependence inflation factor rather than as a pure correlation.

Figure 9 reports the result on all six windows at $m=16$, over 400 random samples at each K. Three facts matter. First, the estimator is essentially unbiased from about thirty-two blocks: the median bias stays within $\pm3$ % for every window and every block size from $K=32$ onward, and within $\pm1.5$ % from $K=128$. Second, the residual uncertainty contracts steadily: the interquartile range of $\widehat V$ is 2.9 to 14 times narrower than the level-A band at $K=64$, 4.1 to 22 times at $K=128$, and 6.7 to 29 times at $K=256$. Third, the estimator is worse than the band when K is very small: at $K=4$ its 90th-percentile relative error is 83–184 %, at or above the width of the level-A band itself, so a handful of blocks should be reported as a bound rather than as an estimate.

The comparison with independence is the practical one. Independence is off by 72–95 % at every K; the plug-in estimator reaches a 90th-percentile error below 32 % by $K=128$ and below 23 % by $K=256$ on every window. The level-A band therefore sits between two useful things: wider than a modest training effort can deliver, and the only honest statement when no training data exist. The caveat is the one Section 6.7 already makes for a stored parameter — the exchangeable model is a summary, and on a landscape with a strong regional trend the single number it produces summarises everything the field does rather than a separable within-block correlation.""")

# ---------------------------------------------------------------- 7. discussion
swap("### 7.3 Practical implications",
     r"""### 7.3 Practical recommendations

The results condense into four rules, two for whoever produces a PCR and two for whoever uses one.""")

swap("Three design rules follow for anyone building or consuming a probabilistic-compositional raster.",
     "")

swap("First, **do not report a variance for any aggregate finer than the whole map unless you also report the dependence structure**. The mean is safe; the variance is not. If dependence is unknown, report the identified interval, not a point.",
     r"""*For the data producer.* First, **do not report a variance for any aggregate finer than the whole map unless you also report the dependence structure**. The mean is safe; the variance is not. If dependence is unknown, report the identified interval, not a point. Second, **store pair structure, and store it honestly**: pair covariance is exactly what the variance needs and captures most of the exceedance gap, and a single exchangeable parameter is enough to identify the variance outright.""")

swap("Second, **store pair structure, and store it honestly**. Pair covariance is exactly what is needed for the variance and captures most of the gap for exceedance. If a single exchangeable parameter is used, its estimation error propagates directly into confidence, and a wrong parameter can make the interval spuriously narrow.",
     r"""*For the data user.* Third, **expect queries, not scenes**. Because the worst case for an exceedance depends on the threshold, a data product cannot pre-compute one adversarial realisation; either the joint is stored, or the query is solved as a program over the coupling polytope. Fourth, **treat the short-range variogram as the part that matters most, and the identified set as the fallback**. If a variogram is available and its range is resolved over a lag window several times longer than the practical range, it will estimate a block variance to a few per cent and should be used; if it is not, the sharp bounds require no model, are never wrong, and can be reported alongside the model estimate, since Section 6.13 shows they bracket the decision-relevant quantity from opposite sides.""")

swap("Third, **expect queries, not scenes**. Because the worst case for an exceedance depends on the threshold, a data product cannot pre-compute one adversarial realisation. Either the joint is stored, or the query is solved as a program over the coupling polytope.",
     "")

swap("Fourth, **treat the short-range variogram as the part that matters most, and the identified set as the fallback**. If a variogram is available and its range is resolved over a lag window several times longer than the practical range, it will estimate a block variance to a few per cent and should be used. If it is not available, the sharp bounds require no model, are never wrong, and are the honest alternative to a wrong number that looks precise; the two can be reported together, since Section 6.12 shows they bracket the decision-relevant quantity from opposite sides.",
     "")

# --- PSD trust guidance + area-to-point, into limitations & comparator
splice_after("The reported interval widths were, however, insensitive to the grid, changing by less than $0.01$ when $L$ was raised from 3 to 9.",
             r""" A practical consequence for the reader is *when to trust the PSD relaxation*. Proposition 6 gives the exact rule: the relaxation is exact whenever $m\min(p,1-p) \le 1$ — always true for $m=2$, and for any block size when the class of interest is sparse enough — and otherwise overstates the lower endpoint by at most $1/(4m^{2})$, which is below 1 % of the marginal variance once $m \ge 10$ for a balanced class. The relaxation is therefore safe for the variance at every realistic block size, and the place to be wary is not the variance but the *exceedance* width, where no relaxation is sharp and the bound should be reported as a bound.""")

# ---------------------------------------------------------------- 8. conclusion
i0 = text.index("## 8. Conclusion")
i1 = text.index("## Data and code availability")
NEW_CONCL = r"""## 8. Conclusion

A probabilistic-compositional raster claims to know more than a value: it claims to know a distribution. We have shown that the claim is only partially transferable across scales, and we have located exactly where it fails.

For the **producer** of such a raster the implication is a metadata contract. The aggregate mean is identified from what is stored, so it can be advertised. The aggregate variance is not: it is identified exactly by one additional scalar in the exchangeable case, and by the pairwise covariance matrix in general, so storing that scalar or that matrix is the difference between a variance that means something and one that does not. If neither is stored, the honest product is the sharp bound, which requires no model and is never wrong. A wrong dependence parameter is worse than none, because it closes the interval while excluding the truth.

For the **user** the implication is a sampling cost. On real ESA WorldCover land cover a per-pixel-only representation understates block variance by a factor of 3.5 to 19.8, and a nominal 90 % interval covers 3–7 % of the time at realistic autocorrelation. A $\pm1$ percentage-point target on a regional mean therefore needs 3.5 to 19.8 times more pixels than independence implies. Where a variogram is available and its range is resolved, it recovers the block variance to within 3.6 %; where it is not, the identified set is the fallback, and the two bracket the answer from opposite sides. And where a few dozen training blocks can be observed, they contract the interval to a usable estimate — the bridge between a bound and a number.

The obstruction is representational rather than computational. Scale-consistency is not something a probabilistic raster gets for free by attaching distributions to pixels; it is a property of what those distributions are allowed to carry.

"""

text = text[:i0] + NEW_CONCL + text[i1:]
_log.append("ok  conclusion rewritten (producer / user split, symbols compressed)")

# ---------------------------------------------------------------- 9. back matter
swap("## Data and code availability",
     r"""## Data and codes availability statement

Every number in this paper is produced by a script and read from a CSV; nothing is transcribed by hand. The closed forms and the relaxations are verified against the exact linear program in `exp_psd.py` and `pcp_core.py`; the spatial and real-data results come from `exp_spatial.py` and `exp_real.py`; the multi-class compositional experiments from `exp_compositional.py`; the geostatistical comparator and the decision analysis from `exp_baseline.py`, which reads the window definitions from a single table in `windows.py`; and the dependence-parameter proof of concept from `exp_estimator.py`. Figures are produced by `make_figures.py`, `make_figure7.py`–`make_figure10.py`, `make_figure_roadmap.py`, `make_figure_schematic.py` and `make_figure_estimator.py`, all of which share one style module (`figstyle.py`, Times New Roman) and read only their own CSV output.

**Reproducibility conditions.** The synthetic experiments (`exp_theory.py`, `exp_spatial.py`, `exp_compositional.py`, `exp_psd.py`) are fully offline: seeds are fixed and the variogram and block-variance computations are deterministic, so every reported number is reproducible bit-for-bit with only NumPy, SciPy and pandas. The real-data experiments additionally require network access to the ESA WorldCover 2021 v200 public cloud-optimised GeoTIFF on the AWS Open Data endpoint (`esa-worldcover.s3.eu-central-1.amazonaws.com`); no credentials are needed, and each window is read once and cached under `data/`, after which the analysis is offline. The six windows are named in `windows.py`, which is the single source of truth for the study design. Scripts, configuration, the cached window definitions and the CSV outputs underlying every figure are provided in the supplementary repository accompanying this manuscript.

## Acknowledgements

We thank the editors and the anonymous reviewers for comments that improved the manuscript.

## Disclosure statement

No potential conflict of interest was reported by the author(s).

## Funding

To be completed before submission.

## Notes on contributors

To be completed before submission: author names, affiliations, ORCID and short biographies in IJGIS house style.""")

# ---------------------------------------------------------------- 10. appendices
appendix = r"""

## Appendix A. Notation

| Symbol | Meaning |
|---|---|
| $m$ | number of fine units aggregated into one block |
| $w_j$ | area weight of fine unit $j$, $\sum_j w_j = 1$ |
| $K$ | number of composition classes |
| $\nu_j$ | composition of fine unit $j$, a point of the simplex $\Delta^{K-1}$ |
| $\nu_C$ | composition of the block, $\nu_C = \sum_j w_j \nu_j$ |
| $\Lambda_j$ | stored marginal law of $\nu_j$ (level A) |
| $\Gamma$ | joint law of $(\nu_1,\dots,\nu_m)$; never stored |
| $\Lambda_C$ | law of $\nu_C$, the push-forward of $\Gamma$ |
| $\sigma^{2}$ | common fine-unit variance |
| $V^{-}, V^{+}$ | sharp lower and upper endpoint of the identified set for $\mathrm{Var}(\nu_C)$ |
| $\mathrm{ACI}(m)$ | aggregation-ambiguity index, $(V^{+}-V^{-})/V^{+}$ |
| $\rho$ | exchangeable pairwise correlation |
| $\rho_{\mathrm{FH}}$ | sharp minimum pairwise correlation of a two-state marginal |
| $\rho_{\mathrm{PSD}}$ | PSD-feasible pairwise correlation, $\max\{\rho_{\mathrm{FH}},-1/(m-1)\}$ |
| $\Omega_q$ | number of lattice states of one fine composition |
| $L$ | number of grid levels per class in the compositional experiments |
| $\gamma(h)$ | variogram; $c_0, c_1, a$ its nugget, partial sill and range |
| $\mathrm{DEFF}$ | design effect, $m\,\mathrm{Var}(\bar{\nu}_B)/\sigma^{2}$ |
| $R$ | achieved variance-reduction factor, $\sigma^{2}/\mathrm{Var}(\nu_C)$ |

## Appendix B. Proofs deferred from Section 4

**Proof of Proposition 5.** The upper endpoint is $p(1-p)$, attained by the comonotonic coupling exactly as in Proposition 4. For the lower endpoint, under the rotation coupling each $\nu_j = \mathbf{1}\{(U+jp)\bmod 1 < p\}$ is Bernoulli with parameter $p$, so the marginals are reproduced. The number of units equal to one is the multiplicity with which the $m$ arcs $\{(u+jp)\bmod 1 < p\}$ cover the point $u$, and because the arc endpoints are spaced exactly $p$ apart on the circle, that multiplicity is either $\lfloor mp\rfloor$ or $\lceil mp\rceil$ for every $u$. The sum therefore takes only those two values; since its mean is $mp$ it places mass $\pi=\operatorname{frac}(mp)$ on $\lceil mp\rceil$ and $1-\pi$ on $\lfloor mp\rfloor$, and its variance is $\pi(1-\pi)$, giving $\mathrm{Var}(\nu_C)=\pi(1-\pi)/m^{2}$.

For optimality, under any admissible coupling the sum $S=\sum_j \nu_j$ takes values in $\{0,\dots,m\}$ with mean $mp$. Among all laws on that lattice with that mean, the minimum variance is attained by the two-point law concentrated on $\lfloor mp\rfloor$ and $\lceil mp\rceil$, whose variance is $\pi(1-\pi)$; every other admissible law has variance at least as large. The couplings the marginals admit produce a subset of these laws, so no coupling does better, and the rotation coupling attains the bound. $\square$

**Proof of Proposition 6.** Substituting $\rho_{\mathrm{FH}}=-\min(p,1-p)/\max(p,1-p)$ and the value of $V^{-}$ from Proposition 5 into the equicorrelation formula $\mathrm{Var}(\nu_C)=\sigma^{2}\big[\frac1m+\frac{m-1}{m}\rho\big]$ gives the two regimes directly, because
$\rho_{\mathrm{FH}} < -\frac{1}{m-1} \iff (m-1)\min(p,1-p) > \max(p,1-p) \iff m\min(p,1-p) > 1.$
In regime (i) the PSD constraint is slack, so the relaxation returns the exact pairwise endpoint $V^{-}$. In regime (ii) the constraint binds at $\rho=-\frac1{m-1}$, where the equicorrelation formula gives exactly zero; the exact minimum is $\pi(1-\pi)/m^{2}>0$, and the Fréchet–Hoeffding endpoint is negative and therefore infeasible. $\square$

**Proof of Corollary 6.1.** In regime (ii) the gap is $\mathrm{gap}(m,p)=\pi(1-\pi)/m^{2}$ with $\pi=\operatorname{frac}(mp)\in[0,1)$, so $\mathrm{gap}\le\frac1{4m^{2}}$, with equality when $mp$ is a half-integer. Dividing by $\sigma^{2}=p(1-p)$ gives the relative error $\frac{\mathrm{gap}}{\sigma^{2}}\le\frac{1}{4m^{2}p(1-p)}$, which is $O(m^{-2})$ and equals $m^{-2}$ at $p=\frac12$. $\square$
"""

splice_after("To be completed before submission: author names, affiliations, ORCID and short biographies in IJGIS house style.",
             appendix)

# ---------------------------------------------------------------- 11. figure paths + renumber
text = text.replace("../figures/", "../figures_v2/")

FIG_MAP = {1: 3, 2: 4, 3: 5, 4: 6, 5: 7, 6: 8, 7: 10, 8: 11, 9: 12, 10: 13}
for n in sorted(FIG_MAP, reverse=True):
    text = re.sub(rf"\bFigure {n}\b", f"@@F{n}@@", text)
    text = re.sub(rf"\bPanel {n}\b", f"@@P{n}@@", text)
    text = re.sub(rf"\bFigures {n}\b", f"@@F{n}@@", text)
for n, new in FIG_MAP.items():
    text = text.replace(f"@@F{n}@@", f"Figure {new}")
    text = text.replace(f"@@P{n}@@", f"Panel {new}")
_log.append("ok  figure numbers remapped for the two new front-matter figures")

# ---------------------------------------------------------------- 12. captions
CAPTIONS = [
 (r"""**Figure 1.** Aggregation widens the identified set. (a) Sharp bounds on the block variance $\mathrm{Var}(\nu_C)$ at level A as a function of the number of aggregated pixels $m$, for $p=0.5$ and $p=0.3$. The upper endpoint is the single-pixel variance $V^{+}=V_1$; the independent coupling sits at $\sigma^2/m$ at the bottom of the interval. (b) The aggregation-ambiguity index $\mathrm{ACI}(m)$ computed by exact linear programming, against the worst-case bound $1-1/m$.""",
  r"""**Figure 1.** Aggregating more pixels widens, not narrows, the identified set. (a) Sharp level-A bounds on $\mathrm{Var}(\nu_C)$ against $m$; the upper endpoint is the single-pixel variance and the independent coupling sits at $\sigma^2/m$. (b) The exact $\mathrm{ACI}(m)$ against the worst-case bound $1-1/m$."""),
 (r"""**Figure 2.** (a) Sharp bounds on the exceedance $\Pr(\nu_C\ge\tau)$ for $m=4$, $p=0.5$, together with the comonotonic, countermonotonic and independent couplings. Four regimes are visible; no coupling is extremal at all thresholds. (b) The identification ladder at $m=4$: level A (marginals only), level B (marginals plus the true pairwise covariance matrix) and level C (full joint law), with the true value marked.""",
  r"""**Figure 2.** No coupling is extremal at every threshold. (a) Sharp exceedance bounds for $m=4$, $p=0.5$, with the comonotonic, countermonotonic and independent couplings; four regimes are visible. (b) The identification ladder at $m=4$: levels A, B and C, with the true value marked."""),
 (r"""**Figure 3.** Simulated landscapes with known ground truth, $p=0.5$, mean of six seeds. (a) True block variance and the level-A identified set against correlation length. (b) The variance-reduction factor actually achieved against the $R=m$ promised by independence, for $m=4,9,16$. (c) Exceedance at $\tau=0.5$, $m=4$, against the sharp envelope.""",
  r"""**Figure 3.** The truth is always inside the bounds; independence is far outside. (a) True block variance and the level-A identified set against correlation length. (b) The variance-reduction factor achieved against the $R=m$ promised by independence. (c) Exceedance at $\tau=0.5$ against the sharp envelope."""),
 (r"""**Figure 4.** Calibration of the independence assumption. (a) Empirical coverage of a nominal 90 % interval built from the stored marginals, as a function of correlation length, for $m=4,9,16$ and $p=0.3,0.5$. (b) The ratio of the nominal interval width to the true interval width; the interval is the wrong shape as well as the wrong size.""",
  r"""**Figure 4.** A nominal 90 % interval built under independence covers as little as 3 %. (a) Empirical coverage against correlation length. (b) The interval is the wrong shape as well as the wrong size."""),
 (r"""**Figure 5.** Real land cover, ESA WorldCover 10 m. (a) The 10 m reference with the 30 m pixel grid. (b) The stored marginal law of the within-pixel tree fraction. (c) True block variance against the identified bounds and the independent estimate for $m=4,9,16,25$. (d) The variance-reduction factor actually achieved against the factor implied by independence, with the understatement factor labelled.""",
  r"""**Figure 5.** The gap is not an artefact of the simulation. (a) ESA WorldCover 10 m reference with the 30 m pixel grid. (b) The stored marginal law of the within-pixel tree fraction. (c) True block variance against the identified bounds. (d) The variance-reduction factor achieved against the one independence implies."""),
 (r"""**Figure 6.** The value of one stored dependence parameter. (a) Exceedance as a function of the exchangeable correlation $\rho$ for $\tau=0.25,0.5,0.75$, $m=4$, $p=0.5$; the PSD bound is marked. (b) Interval width at level A and with the true $\rho$ added, showing shrinkage without closure.""",
  r"""**Figure 6.** One stored scalar identifies the variance but not the exceedance. (a) Exceedance against the exchangeable correlation $\rho$; the PSD bound is marked. (b) Interval width at level A and with the true $\rho$ added."""),
 (r"""**Figure 7.** Asymptotic sharpness of the PSD relaxation. (a) The joint-mixability gap against $m$ for three proportions, with the envelope $1/(4m^{2})$. (b) The relative gap against $m$, with the $O(m^{-2})$ references; the 1 % level and $m = 10$ are marked. (c) The gap against $p$ at fixed $m = 4, 8, 16$; it vanishes while $m\min(p,1-p) \le 1$ and rises to a plateau of height exactly $1/(4m^{2})$.""",
  r"""**Figure 7.** The cheap relaxation is asymptotically exact. (a) The joint-mixability gap against $m$ with the envelope $1/(4m^{2})$. (b) The relative gap, with the $O(m^{-2})$ references. (c) The gap against $p$: zero while $m\min(p,1-p)\le1$, then a plateau of height exactly $1/(4m^{2})$."""),
 (r"""**Figure 8.** Ambiguity of a genuine simplex composition ($K=3$). (a) The sharp
interval for $\Pr(\text{diversity} \ge 0.4)$ on the controlled $L=3$ lattice; the
labelled width grows with $m$. (b) The sharp interval for mean Simpson diversity
itself, with the independent and comonotonic couplings. (c) Real ESA
WorldCover, $K=3$: enumerated truth, the sharp band, and the independent and
comonotonic couplings; the labels give the factor by which independence
overstates the truth. The band width is 89–92 % of the range it could take and the
truth is inside it at every $m$.""",
  r"""**Figure 8.** The ambiguity survives a genuine simplex and a nonlinear query. (a) The sharp interval for $\Pr(\text{diversity} \ge 0.4)$ on the $L=3$ lattice; the labelled width grows with $m$. (b) The sharp interval for mean Simpson diversity, with the independent and comonotonic couplings. (c) Real ESA WorldCover, $K=3$: enumerated truth inside a band that covers 89–92 % of its range; the labels give the factor by which independence overstates the truth."""),
 (r"""**Figure 9.** The identified set against the geostatistical route, and the
decision consequence, on all six reference windows, grouped by biome. (a) The
normalised empirical variogram of each window, with the 60-pixel fit window
shaded; triangles mark the first lag at which $\gamma$ reaches
$0.95\sigma^{2}$. The boreal curve is already at the sill inside the shaded
window and the Mediterranean curves are the furthest from it, which is why the
truncated fit is least damaging in the boreal window. (b) Error against the
enumerated block variance for five estimators, over all 24 cases; horizontal
bars are medians. (c) Pixel samples needed for $\pm1$ percentage point on the
regional mean, as a multiple of the independence requirement; solid lines use
the variogram estimate and dotted lines the sharp upper bound. The six curves
are nearly superimposed, which is the visual form of the biome-invariance of the
design effect. (d) Mean absolute error by biome for the four variogram-based
estimators. The ordering boreal $<$ temperate $<$ Mediterranean holds for every
estimator and mirrors the ordering of the fitted range, not of the sample size.""",
  r"""**Figure 9.** The identified set against the geostatistical route, on all six windows. (a) Normalised empirical variograms with the 60-pixel fit window shaded; the boreal curve is already at the sill inside it and the Mediterranean curves are furthest from it. (b) Error against the enumerated block variance for five estimators over all 24 cases. (c) Pixel samples needed for a $\pm1$ percentage-point target, as a multiple of the independence requirement. (d) Mean absolute error by biome, which orders by fitted range rather than by sample size."""),
 (r"""**Figure 10.** The mechanism behind the biome grading of the truncated-window bias. (a) The mean absolute bias of the 60-pixel fit within each window against the fraction of the sill the empirical variogram has reached at 60 pixels. The line is a least-squares fit on that fraction alone, with no biome term, and accounts for 92 % of the variance in the bias. (b) The same resolved fraction against the fitted range, with an empirical fit in $1/a$. A short range saturates inside the 60-pixel window and a long one does not, which is the whole of the biome effect: the boreal window reaches 0.95 of its sill before the window closes, and the Mediterranean windows reach only about 0.78.""",
  r"""**Figure 10.** The truncated-window bias is set by the fitted range, not by the biome. (a) Mean absolute bias of the 60-pixel fit against the resolved fraction of the sill, a least-squares fit with no biome term ($R^{2}=0.92$). (b) The resolved fraction against the fitted range, with an empirical fit in $1/a$ ($R^{2}=0.87$)."""),
]
for old, new in CAPTIONS:
    swap(old, new)
_log.append(f"ok  rewrote {len(CAPTIONS)} captions (first sentence states the result)")

# ---------------------------------------------------------------- 13. PCR pass
for probe in ("a probabilistic-compositional raster (PCR) actually stores",
              "A *probabilistic-compositional raster* (PCR) over the block stores"):
    if probe not in text:
        raise SystemExit(f"PCR definition anchor missing: {probe[:60]}")
parts = text.split("(PCR)")
head = parts[0] + "(PCR)"
tail = "(PCR)".join(parts[1:])
tail = tail.replace("probabilistic-compositional raster", "PCR")
tail = tail.replace("probabilistic-compositional pixel", "probabilistic-compositional pixel")
text = head + tail
_log.append("ok  PCR abbreviation applied after first definition")

# ---------------------------------------------------------------- 14. references
REF_FIX = [
 ("Embrechts, P., & Puccetti, G. (2009). Bounds for the sum of dependent risks having overlapping marginals.",
  "Embrechts, P., & Puccetti, G. (2010). Bounds for the sum of dependent risks having overlapping marginals."),
 ("Puccetti, G., & Rüschendorf, L. (2011). Sharp bounds on the distribution of a function of dependent risks. *Journal of Computational and Applied Mathematics*, 236(7), 1833–1837.",
  "Puccetti, G., & Rüschendorf, L. (2012). Computation of sharp bounds on the distribution of a function of dependent risks. *Journal of Computational and Applied Mathematics*, 236(7), 1833–1840. doi:10.1016/j.cam.2011.10.015"),
]
for old, new in REF_FIX:
    swap(old, new)

DOI_ADD = {
 "Cova, T. J., & Goodchild, M. F. (2002). Extending geographical representation to include fields of spatial objects. *International Journal of Geographical Information Science*, 16(6), 509–532.":
 " doi:10.1080/13658810210137040",
 "Eches, O., Dobigeon, N., Mailhes, C., & Tourneret, J.-Y. (2010). Bayesian estimation of linear mixtures using the normal compositional model: application to hyperspectral imagery. *IEEE Transactions on Image Processing*, 19(6), 1403–1413.":
 " doi:10.1109/TIP.2010.2042993",
 "Goodchild, M. F., Guo, H., Annoni, A., Bian, L., de Bie, K., Campbell, F., … & Zhang, C. (2012). Next-generation digital earth. *Proceedings of the National Academy of Sciences*, 109(28), 11088–11094.":
 " doi:10.1073/pnas.1202383109",
 "Gotway, C. A., & Young, L. J. (2002). Combining incompatible spatial data. *Journal of the American Statistical Association*, 97(458), 632–648.":
 " doi:10.1198/016214502760047140",
 "Heuvelink, G. B. M., Brown, J. D., & van Loon, E. E. (2007). A probabilistic framework for representing and simulating uncertain environmental variables. *International Journal of Geographical Information Science*, 21(5), 497–513.":
 " doi:10.1080/13658810601063951",
 "Heuvelink, G. B. M., & Pebesma, E. J. (1999). Spatial aggregation and soil process modelling. *Geoderma*, 89(1–2), 47–65.":
 " doi:10.1016/S0016-7061(98)00077-9",
 "Kaas, R., Dhaene, J., & Vyncke, D. (2002). A simple geometric proof that comonotonic risks have the convex-largest sum. *ASTIN Bulletin*, 32(1), 71–80.":
 " doi:10.2143/AST.32.1.1015",
 "Kline, B., & Tamer, E. (2023). Recent developments in partial identification. *Annual Review of Economics*, 15, 125–150.":
 " doi:10.1146/annurev-economics-051520-021124",
 "Molinari, F. (2020). Microeconometrics with partial identification. In *Handbook of Econometrics*, Vol. 7A (pp. 355–486). Elsevier.":
 " doi:10.1016/bs.hoe.2020.05.002",
 "Nelsen, R. B. (2006). *An Introduction to Copulas* (2nd ed.). Springer.":
 " doi:10.1007/0-387-28678-0",
 "Olofsson, P., Foody, G. M., Herold, M., Stehman, S. V., Woodcock, C. E., & Wulder, M. A. (2014). Good practices for estimating area and assessing accuracy of land change. *Remote Sensing of Environment*, 148, 42–57.":
 " doi:10.1016/j.rse.2014.02.015",
 "Prékopa, A. (1990). Sharp bounds on probabilities using linear programming. *Operations Research*, 38(2), 227–239.":
 " doi:10.1287/opre.38.2.227",
 "Chiles, J.-P., & Delfiner, P. (2012). *Geostatistics: Modeling Spatial Uncertainty* (2nd ed.). Wiley.":
 " doi:10.1002/9781118136188",
 "Cressie, N. (1993). *Statistics for Spatial Data* (rev. ed.). Wiley.":
 " doi:10.1002/9781119115151",
 "Goovaerts, P. (1997). *Geostatistics for Natural Resources Evaluation.* Oxford University Press.":
 " doi:10.1093/oso/9780195115383.001.0001",
 "Webster, R., & Oliver, M. A. (2007). *Geostatistics for Environmental Scientists* (2nd ed.). Wiley.":
 " doi:10.1002/9780470517277",
}
for ref, doi in DOI_ADD.items():
    swap(ref, ref + doi)
_log.append(f"ok  added {len(DOI_ADD)} Crossref-verified DOIs")

for old, new in (("Puccetti & Rüschendorf, 2011", "Puccetti & Rüschendorf, 2012"),
                 ("Puccetti and Rüschendorf (2011)", "Puccetti and Rüschendorf (2012)"),
                 ("Embrechts & Puccetti, 2009", "Embrechts & Puccetti, 2010"),
                 ("Embrechts and Puccetti (2009)", "Embrechts and Puccetti (2010)")):
    text = text.replace(old, new)
_log.append("ok  in-text citation years aligned with Crossref")

# ---------------------------------------------------------------- write
DST.write_text(text, encoding="utf-8")
print("\n".join(_log))
print(f"\nwrote {DST}  ({len(text.splitlines())} lines, {len(text.split())} words)")
