# Scale-consistency of probabilistic-compositional rasters: what a per-pixel distribution does and does not determine

*Submitted to the International Journal of Geographical Information Science (IJGIS) as a research paper. Author names, affiliations and contact details appear on the separate title page and are withheld here to preserve double-blind review.*

---

## Abstract

Geographic data models increasingly attach a distribution, rather than a
single value, to each pixel. Yet the operations that make a raster useful
(aggregation to a coarser support, overlay, resampling, area estimation) are
defined on the compositions, not on the distributions that encode them. We ask
what current probabilistic sub-pixel models do not answer:
**given what a probabilistic-compositional raster (PCR) actually stores, what can and cannot be identified about a coarser unit derived from it?**
We answer it as a problem of partial identification, not of estimation. The
stored object is the marginal law of each fine unit's composition; the derived
object is the law of an area-weighted aggregate; and the three levels a
product might store are (A) marginals only, (B) marginals plus the pairwise
covariance matrix, (C) the full joint law. The aggregate mean is always
identified. The aggregate variance is identified exactly at level B and not at
level A, and exact identification needs only pair information when two units
are combined but is genuinely partial for three or more. Threshold queries are
more demanding still: the extremal dependence structure changes with the
threshold, so no single worst-case scene stands in for the unknown joint. For
the two-state case the interval is available in closed form, and a cheap
positive-semidefinite relaxation is exact for sparse classes and
asymptotically exact elsewhere, its error measured against the single-pixel
variance decaying as $O(m^{-2})$. We verify every closed form against an exact
linear program and test the theory on eleven windows of the ESA WorldCover 10
m product spanning seven biomes, four of them in mainland China. A
per-pixel-only representation understates block variance by a factor of 3.4 to
21.2 on this reference field and degrades a nominal 90 % interval to 3–7 %
empirical coverage. A variogram read non-parametrically off an adequate lag
window recovers the block variance to within 4.9 % in all 44 cases, while a
single exponential fitted over the common 60-pixel default is biased low by
10–43 %; a $\pm1$ percentage-point target on a regional mean therefore costs
3.4 to 21.2 times more samples than independence implies. One stored
dependence scalar identifies the variance at one block size but does not
transfer to another, overstating it by 8.6–34.1 % when exported from four
pixels to twenty-five. We do not argue against Bayesian sub-pixel models or
against geostatistics: a model that stores or outputs the full joint is
unaffected. A variogram is an effective remedy, but it is an external
modelling assumption, not a consequence of the stored marginals, and the
metadata must say so.

The remedy is therefore representational rather than algorithmic: a PCR is
only as useful as the cross-unit dependence information it stores, so store
that dependence, and state the block size at which it holds, or store the
sharp bounds and say that they are bounds.

**Keywords:** change of support; partial identification; raster metadata; uncertainty representation; land cover

---

## 1. Introduction

Geographical information science has long treated the pixel as an elementary
value: a cell holds one label or one measurement, and the map is the
collection of those values. Two decades of work on uncertainty and
representation have loosened that assumption without abandoning it. Extended
object-field models let a cell carry field-like and object-like content at
once (Cova & Goodchild, 2002). Probabilistic frameworks attach a distribution
to an uncertain environmental variable and propagate it through GIS operations
(Heuvelink et al., 2007). Continuous sub-pixel models estimate the *fraction*
of each class inside a cell rather than a single winner (Eches et al., 2010).
Each move enriches what a cell *is*. None settles the question that becomes
unavoidable as soon as a cell holds a composition rather than a value:
**how does the information stored at fine support determine, or fail to determine, the information required at coarse support?**

The question is not academic. Almost every quantity a land-cover or
land-surface map is asked to support is defined on a support coarser than the
one the map is stored at: area estimates, biome fractions, deforestation rates
and national greenhouse-gas inventories are all aggregates. Change of support
is therefore not an optional post-processing step but the point at which a
raster becomes evidence about the world. If the stored cells are uncertain,
the honest question is whether their uncertainty can be combined correctly,
and if not, how badly the combination can go wrong.

We approach this through a specific, deliberately narrow object. A
*probabilistic-compositional pixel* stores, for each fine unit $j$, the law of
its composition over the simplex of $K$ class fractions,
$$\Lambda_j = \mathrm{Law}(\nu_j \mid \mathcal{O}) \in \mathcal{P}(\Delta^{K-1}), \qquad \nu_j \in \Delta^{K-1},$$
where $\mathcal{O}$ is the observation. This sits in a nested hierarchy of
increasing expressive power: an ordinary *categorical pixel* records one
class; a *deterministic compositional pixel* records one vector of fractions
(as in spectral unmixing without uncertainty); a *probabilistic-compositional
pixel* records a distribution over vectors. Moving down the hierarchy discards
information, and each step is a strict generalisation of the one below it.

What is not natural, and what this paper is about, is that the hierarchy does
**not** close under aggregation. The physically meaningful object for a coarse
unit $C$ formed from fine units $1,\dots,m$ with area weights $w_j$ is the
joint law of the fine compositions,
$$\Gamma = \mathrm{Law}\big((\nu_1,\dots,\nu_m)\mid\mathcal{O}\big),$$ because
the aggregate composition is the deterministic function
$$\nu_C = g_w(\nu_1,\dots,\nu_m) = \sum_{j=1}^{m} w_j\,\nu_j,$$ and the law of
$\nu_C$ is the push-forward $\Lambda_C = (g_w)_\#\,\Gamma$. A PCR that stores
only the marginals stores everything needed to form $g_w$ but **not**
everything needed to form $\Gamma$, so the aggregate law is not determined by
the stored raster. That is the gap we formalise, quantify and test.

**Research gap.** Existing work supplies either a better estimate or a better
propagation machinery. Bayesian spectral unmixing and spatial compositional
models estimate the per-pixel composition and its uncertainty (Eches et al.,
2010; Pirzamanbein et al., 2018), and geostatistical change of support
propagates that uncertainty through a fitted variogram (Heuvelink & Pebesma,
1999; Gotway & Young, 2002; Wadoux & Heuvelink, 2023). Both assume the
dependence between units is available and both therefore answer a question
about a *model*. Neither asks the prior question a data product must answer:
given only what the raster stores, what can be said at all about a coarser
unit? That is a question of identifiability, and it is the question this paper
answers.

The machinery is borrowed; the question is new. We pose change of support as
an identifiability problem for a PCR, characterise exactly what each stored
level determines (Propositions 1–3 and Corollary 2.1), give sharp and
tractable answers (Propositions 4–6 and Corollary 6.1), and test them against
exact linear programs, synthetic landscapes with known truth and eleven real
reference windows across seven biomes, four of them in mainland China. Figure
1 summarises the route from that question to those findings.

![](../figures_v10/fig0_roadmap.png)

**Figure 1.** Scale-consistency is an identifiability problem, not a
computational one. *1 Question*: a stored marginal law is mapped by the
aggregation operator $g_w$ to an aggregate $\nu_C$ whose law is not thereby
determined. *2 Formalisation*: the three rungs of the metadata ladder (A
marginals, B marginals plus pair covariance, C full joint). *3 Theory*: the
identified interval $[V^-,V^+]$, solid where the stored information pins the
variance and open where it leaves it free. *4 Experiments*: eleven reference
windows over seven biomes. *5 Findings*: the variance independence implies
against the enumerated truth.

**What is new here, for GIS and remote sensing specifically.** The
mathematical tools are borrowed, so we say plainly which parts of the paper
are specific to geospatial data products. First, the three-level ladder (A, B,
C) is a statement about what a raster *product* can ship. The
dependence-modelling literature takes the joint law as given and has no
counterpart for it. Second, the aggregation-ambiguity index turns that ladder
into a single number a producer can compute from the stored marginals alone
and print in the metadata. Third, the finding that the reproducing dependence
parameter falls with block size is a change-of-support result rather than a
probabilistic one: it says a dependence parameter is a right answer at one
support and a wrong answer at another, so a metadata contract has to carry a
support label. Fourth, both the two-regime rule $m\min(p,1-p)>1$ and the
variogram-lag-window diagnosis of Section 5.6 can be checked before any model
is fitted. None of the four is a new theorem; all four are new obligations on
a raster that advertises per-pixel uncertainty.

The paper is organised as follows. Section 2 positions the work against
sub-pixel estimation and change of support, Section 3 formalises the raster
and states what each level determines, Section 4 describes the data and the
design, Section 5 reports the experiments with their reading, implications and
limits, and Section 6 concludes.

---

## 2. Background and related work

### 2.1 Sub-pixel composition and probabilistic sub-pixel models

The idea that a cell contains a mixture rather than a member is as old as
spectral unmixing. Linear mixture models express a pixel spectrum as a convex
combination of endmembers, and the normal compositional model recasts these
fractions as random variables with a tractable posterior (Eches et al., 2010).
Compositional-data theory supplies the geometry: fractions live on a simplex,
not in a Euclidean box, and the log-ratio transform provides an inner product
in which the usual multivariate methods are admissible (Aitchison, 1986).
Pirzamanbein et al. (2018) bring this geometry into a spatial setting,
modelling the composition field jointly so that neighbouring cells share
information.

This literature is about *recovering* composition and its uncertainty. It is
largely silent about what happens once the recovered object is aggregated,
apart from noting that the aggregate is a weighted mean. The distinction
matters because the uncertainty attached to a composition is a full
distribution on a simplex, and the operation that produces an aggregate is
nonlinear in that distribution. The mean of a sum is the sum of means; almost
nothing else about a sum is the sum of the corresponding summaries.

### 2.2 Change of support and aggregation of uncertainty

Geostatistics has a mature theory of change of support. Block kriging gives
the variance of a block average given a variogram, and the regularisation of a
variogram to a block support is a standard operation. Heuvelink and Pebesma
(1999) and Gotway and Young (2002) analyse the more general problem of
combining data on incompatible supports. Heuvelink et al. (2007) provide a
probabilistic framework for uncertain environmental variables and their
propagation. Most directly, Wadoux and Heuvelink (2023) quantify the
uncertainty of spatial averages and totals of natural resource maps, the exact
operation we consider, and show that ignoring spatial correlation corrupts the
interval for a spatial mean. The same dependence corrupts validation as well
as aggregation, since spatially autocorrelated samples inflate apparent
accuracy (Segurado et al., 2006; Kattenborn et al., 2022).

This body of work assumes that the dependence between units is known, or
estimable from a model fitted to observations, and therefore answers a
question about a *model* rather than about the metadata a raster stores. The
two questions are not in competition; one comes before the other. Under a
given covariance or variogram model, block kriging and area-to-point kriging
return a block value and its variance. This paper asks what is determined
*before* any model is chosen: the identified set is the range every such
estimate must fall in, so it is the answer available when no model is
justified, and it is the check that any model's number has to survive.

### 2.3 Uncertainty representation, partial identification, and sharp bounds

Three strands inform the framing. First, the GIS literature on uncertainty
representation has repeatedly argued that a raster should expose the semantics
of its uncertainty, not just a numeric error field (Goodchild et al., 2012;
Cova & Goodchild, 2002). A per-pixel variance without a cross-pixel dependence
structure is exactly the kind of under-specified metadata that this literature
warns about. Second, the econometrics of partial identification supplies the
right vocabulary: when the data do not pin down a parameter, the honest object
is the identified set, not a point estimate, and the useful question is how
the set shrinks as assumptions are added (Molinari, 2020; Kline & Tamer,
2023). We import this stance into change of support and treat the aggregate
law as partially identified. The identified set of a functional is the range
of values it can take under every coupling consistent with the stored
information; it is the honest object when the data do not pin the functional
down.

Third, the technical core is borrowed from actuarial and operations-research
work on dependent risks. Hoeffding (1940) and Fréchet (1951) characterise the
extremal couplings with fixed marginals (Nelsen, 2006, gives the modern
account); Prékopa (1990) shows that sharp bounds on probabilities of events
defined by linear constraints are obtained by linear programming over the
coupling polytope; Puccetti and Rüschendorf (2012) give sharp bounds on the
distribution of a function of dependent risks; Embrechts and Puccetti (2010)
handle overlapping marginals; Wang and Wang (2015) develop extreme negative
dependence. Kaas et al. (2002) give the geometric proof that comonotonic risks
have the convex-largest sum. The essential fact we exploit (that the
Fréchet–Hoeffding pairwise endpoints cannot all hold at once in general, so
the pairwise relaxation overshoots) is the *joint-mixability* phenomenon: the
pairwise extremal dependence structures are not simultaneously realisable by
any single coupling. We are not aware of it being used to critique the
metadata of a geospatial raster, which is the contribution here.

---

## 3. Methodology

This section formalises the raster, the three levels of stored information,
and the aggregation problem, and then states what each level determines.
Nothing here is fitted: every result is a statement about the information a
product stores, not about a model of the landscape.

### 3.1 Support, composition, and the information a raster stores

Let a study region be partitioned into fine cells and grouped into coarse
units, or *blocks*, $C$. Block $C$ contains $m$ fine units $1,\dots,m$ with
areas $a_j > 0$; write $w_j = a_j / \sum_k a_k$ so that $\sum_j w_j = 1$. Each
fine unit has a composition
$$\nu_j = (\nu_{j1},\dots,\nu_{jK}) \in \Delta^{K-1} = \Big\{ v \in \mathbb{R}^{K}_{\ge 0} : \textstyle\sum_{k=1}^{K} v_k = 1 \Big\}.$$
The composition of $C$ is the area-weighted average
$$\nu_C = \sum_{j=1}^{m} w_j\,\nu_j.$$ A *probabilistic-compositional raster*
(PCR) over the block stores a distribution for each fine composition. The
strongest object it could store is the joint law
$\Gamma = \mathrm{Law}((\nu_1,\dots,\nu_m)\mid\mathcal{O})$; in practice it
stores the marginals $\Lambda_j = \mathrm{Law}(\nu_j\mid\mathcal{O})$.
Aggregating is a fixed operation on the joint law, so the law of the block
composition follows from it:
$$\Lambda_C = (g_w)_\#\,\Gamma, \qquad g_w(\nu_1,\dots,\nu_m) = \sum_{j=1}^{m} w_j \nu_j.$$
We distinguish three levels by how much of $\Gamma$ is retained. -
**Level A (marginals).** The stored object is $\{\Lambda_j\}_{j\le m}$. The
joint $\Gamma$ can be any *coupling*: any joint distribution whose per-unit
marginals are the $\Lambda_j$. -
**Level B (marginals plus pairwise covariance).** In addition to
$\{\Lambda_j\}$, the raster stores $\mathrm{Cov}(\nu_j,\nu_k)$ for every pair
$j\ne k$. The pair marginals of $\Gamma$ are fixed; higher-order structure is
free. - **Level C (full joint).** The raster stores $\Gamma$. The aggregate
law is then available by push-forward. These are the natural rungs of a
metadata ladder. Level A is the default produced by any method that estimates
uncertainty pixel by pixel; level B is what a raster achieves if it also ships
a covariance summary, a realistic and cheap extension; level C is the ideal.
The identifiability question is:

The set of values a functional can take under all such couplings is its
*identified set*: the range the stored data leave possible, as opposed to the
single point an estimator returns. "Sharp" is used below in its plain sense: a
sharp bound is one that some admissible coupling actually attains, so it
cannot be improved without more information. Figure 2 draws the gap and the
ladder.

![](../figures_v10/fig0_schematic.png)

**Figure 2.** Storing marginals is not the same as storing a joint. (a) Fine units carry marginal laws; the aggregate law needs the joint, which the raster does not store. (b) The metadata ladder: what each stored level lets one answer for the mean, the variance and an exceedance probability.

### 3.2 The generating process used in the experiments

For most of the paper $K$ is restricted to two classes and $\nu_j$ reduces to
an indicator, so that $\nu_j$ is Bernoulli with within-block marginal $p$ and
$\Lambda_j$ is a two-point law. This is the simplest setting in which the gap
appears and the only one in which the exact coupling polytope is small enough
to enumerate. Section 5.5 confirms that the effect survives genuine
three-class simplex compositions.

**Scope of this section.** Sections 3.1 and 3.2 set up the object and the
experiments. Sections 3.3 to 3.7 state what the stored information determines.
The assumptions below govern 3.3 to 3.7 only.

**Standing assumptions (SA).**

**Standing assumptions (SA).** Throughout Section 3 the $m$ fine units are
exchangeable: they share a common marginal mean $\mu$ and a common variance
$\sigma^{2}$, the weights are equal, $w_j = 1/m$, and the underlying field is
second-order stationary, so that a covariance depends only on the separation
vector. Every proposition below states the premise under which it holds. These
assumptions make the mechanism visible; the mechanism itself is not special to
them.

One distinction inside SA matters more than it looks, and it is the one
Section 5.4 turns on. Second-order stationarity only says that the covariance
between two units is a function of their separation. It does **not** say that
every pair shares the same covariance. Equal correlation for every pair
(equicorrelation) is strictly stronger, and it is what the closed forms of
Sections 3.4 to 3.7 need, because it is what reduces a covariance matrix to
the single number $\rho$. A stationary field with a distance-decaying
covariance is not equicorrelated at any block size, so a single $\rho$ cannot
reproduce it exactly. The simulated and reference fields of Section 4 are not
assumed to be equicorrelated: their covariance varies with lag, so the closed
forms are used as analytical benchmarks and are not applied to heterogeneous
windows as if they were exact.

**What is new here, and what is not.** The mathematics is classical.
Propositions 1 and 2 are elementary identities; Proposition 5 is the
closed-form sharp interval of Hoeffding (1940) and Rüschendorf (1984); the
Fréchet–Hoeffding relaxation it sharpens is due to Prékopa (1990) and Puccetti
and Rüschendorf (2012). Only Propositions 3, 4 and 6 and Corollary 6.1 are
new, and they are new as consequences of a question rather than as machinery.
The question is the one a probabilistic raster forces: not how to propagate a
dependence structure one has *assumed*, but what the dependence structure a
raster *stores* determines about a coarser unit. Asked that way,
joint-mixability stops being a curiosity about artificial marginals and
becomes an engineering defect, and it forces three checkable diagnostics that
no earlier statement of the phenomenon supplies: the aggregation-ambiguity
index (Section 3.5), the two-regime rule $m\min(p,1-p) > 1$ (Section 3.7), and
the support-specific scalar $\rho_{\mathrm{exact}}(m)$ that must be stored
with the block size at which it was measured (Section 5.4). The contribution
is representational and identification-theoretic rather than mathematical: the
theorems are used, not developed. Proofs of Propositions 1–6 and Corollary 6.1
are in Appendix B.

### 3.3 What the mean and the variance require

**Proposition 1 (mean identification).** For any coupling $\Gamma$ with the given marginals,

$$\mathbb{E}[\nu_C] = \sum_{j=1}^{m} w_j\, \mathbb{E}[\nu_j].$$

Proposition 1 says the aggregate mean is exactly recoverable from level A.
This is why the representational gap is easy to miss: the most commonly
reported aggregate quantity, the mean, is unaffected. The gap is entirely in
the second and higher moments.

**Proposition 2 (variance decomposition).** For any coupling $\Gamma$,

$$\mathrm{Var}(\nu_C) = \sum_{j=1}^{m} w_j^2 \,\mathrm{Var}(\nu_j) \;+\; \sum_{j\ne k} w_j w_k\, \mathrm{Cov}(\nu_j,\nu_k).$$

**Corollary 2.1.** At level B the aggregate variance is exactly identified. At level A it is not, and the identified set at level A is the interval swept by the right-hand side as the pairwise covariances range over those compatible with the marginals.

Corollary 2.1 is the pivot of the paper. It localises the missing information
precisely: a PCR needs pair covariance, nothing more, to make the aggregate
variance well defined. It also explains why a per-pixel variance field is not
enough: the diagonal terms $\mathrm{Var}(\nu_j)$ are exactly the part that
*is* stored and exactly the part that does *not* suffice.

### 3.4 Sharp bounds and two cheaper relaxations

The set of all couplings with the prescribed marginals is a convex polytope: a
flat-sided, many-cornered set of joint distributions. Any functional that is a
linear combination of the coupling's entries has sharp bounds given by a
linear program, and because $\mathrm{Var}(\nu_C)$ and $\Pr(\nu_C \ge \tau)$
are such functionals when the states are finite, their sharp bounds are
computable exactly (Prékopa, 1990; Puccetti & Rüschendorf, 2012).

*Premise (SA): exchangeable units with equal weights; the pairwise structure
is otherwise unrestricted.*

**Proposition 3 (exactness for two units, genuine partiality for three or more).**

(i) For $m=2$ the sharp bounds on $\mathrm{Var}(\nu_C)$ over couplings with
fixed marginals are attained, and they coincide with the classical bounds

$$(w_1\sigma_1 - w_2\sigma_2)^2 \;\le\; \mathrm{Var}(\nu_C) \;\le\; (w_1\sigma_1 + w_2\sigma_2)^2.$$

The upper endpoint is attained by the comonotonic coupling (all units moving
together) and the lower by the countermonotonic coupling (units moving in
exact opposition) (Hoeffding, 1940; Kaas et al., 2002). This
standard-deviation form is a statement about two units only. It is not a
general sharp bound: once $m\ge3$, or once the marginals are general discrete
laws rather than two-point laws, the sharp bounds are those of the coupling
program of Section 4.1, and the expression above is neither attainable nor
sharp.

(ii) For $m\ge 3$ the pairwise Fréchet–Hoeffding relaxation, which sets every
pair to its own extremal value, cannot in general hold all at once: the
pairwise endpoint configurations are not jointly mixable. The relaxation
therefore yields a strictly wider interval than the exact linear program, and
its lower endpoint can be negative, which is not a feasible variance.

Proposition 3 explains a numerical symptom that is easy to misread as a bug:
the naïve Fréchet–Hoeffding lower bound on a block variance is frequently
negative. A negative bound is not a tight statement about uncertainty; it is
evidence that the relaxation has left the set of distributions that exist.

Because the exact linear program scales with the number of states, we also
study two tractable outer sets. The **Fréchet–Hoeffding (FH) relaxation**
plugs the pairwise extremal covariances into Proposition 2; it is valid but,
by Proposition 3, loose and often unusable. The
**positive-semidefinite (PSD) relaxation** requires the matrix of pairwise
correlations to be a valid correlation matrix, which in the exchangeable case
with correlation $\rho$ for every pair is the single constraint
$$\rho \;\ge\; -\frac{1}{m-1},$$ because
$(1-\rho)I + \rho \mathbf{1}\mathbf{1}^\top \succeq 0$ for the equicorrelation
matrix. The PSD set is strictly tighter than FH and is never unusable. The
*joint-mixability gap* is the amount by which the PSD relaxation is wider than
the exact linear program at the lower endpoint.

In the exchangeable two-state case with within-block proportion $p$ the two
regimes are explicit. The sharp minimum pairwise correlation of a two-state
marginal is $\rho_{\mathrm{FH}} = -\min(p,1-p)/\max(p,1-p)$, and the PSD
constraint is active only where it is the stricter of the two, which happens
exactly when $m\min(p,1-p) > 1$. Below that threshold the pairwise extremal
configuration can hold all at once and the FH relaxation coincides with the
exact linear program; above it the PSD corner $\rho = -1/(m-1)$ is feasible,
the relaxation returns a zero minimum while the exact minimum is strictly
positive, and the FH endpoint has gone negative and is unusable. Section 3.7
gives the exact size of both.

Figure 3 shows the two relaxations against the exact program, and the single
inequality that decides when the cheap one may be used.

![](../figures_v10/fig12_identified_set_geometry.png)

**Figure 3.** The cheap relaxation is safe, but it is exact only in one of the two regimes, and which one is decided by a single inequality. (a) The three lower endpoints for $m=4$ and $p=0.30$. The Fréchet–Hoeffding relaxation returns $-0.015$, below zero and therefore not a value a variance can take; the positive-semidefinite relaxation repairs the sign and returns $0$, but the sharp minimum is $0.0100$. In this regime the cheap relaxation is therefore usable but *not* exact, and it is inexact by exactly the joint-mixability gap, $0.0100 \le 1/(4m^{2}) = 0.0156$. (b) The boundary between the regimes is the single inequality $m\min(p,1-p)=1$: below it the pairwise relaxation is exactly sharp (for every $m=2$ block and for any block size when the class is sparse enough) and above it the relaxation returns zero where the truth is a small positive number. Only $p\le 1/2$ is drawn, because the two-state problem is symmetric under $p\mapsto 1-p$.

### 3.5 The aggregation-ambiguity index

To report the size of the gap on a scale that has meaning for applications,
define the *aggregation-ambiguity index*

$$\mathrm{ACI}(m) = \frac{V^{+}(m) - V^{-}(m)}{V^{+}(m)},$$

where $V^{+}$ and $V^{-}$ are the sharp maximum and minimum of
$\mathrm{Var}(\nu_C)$ over couplings with the stored marginals.
$\mathrm{ACI}=0$ means the aggregate variance is identified; $\mathrm{ACI}=1$
means the stored information leaves it completely undetermined in relative
terms.

*Premise (SA): exchangeable two-state units with equal weights.*

**Proposition 4 (worst-case aggregation ambiguity).** In the exchangeable two-state case, for every $m\ge 2$,

$$\mathrm{ACI}(m) \;\ge\; 1 - \frac{1}{m}.$$

The content of Proposition 4 is a factor-of- $m$ statement. Independence
asserts that aggregating $m$ pixels divides the variance by $m$; the sharp
upper bound asserts that the aggregate can carry the full single-pixel
variance. Both are consistent with the same stored marginals. The stored
information cannot choose between them.

### 3.6 Threshold queries and compositional states

Variance is not the only functional that matters; threshold and exceedance
queries (the fraction of a block above a tree-cover threshold, say) are
equally common. For those, no closed form replaces the linear program, and the
structure of the sharp bounds is qualitatively different.

**Observation 5.** As the threshold $\tau$ varies, the identity of the extremal coupling changes. Comonotonic dependence is extremal at small thresholds and countermonotonic dependence is extremal at large thresholds, but there is a band of intermediate thresholds at which neither is extremal, and the sharp bound is attained by an interior coupling.

Observation 5 is established numerically in Section 5.1. Its practical message
is simple: one cannot ship a single "worst-case" dependence scenario and use
it at every threshold. The worst case is query-dependent, so the aggregate
must either be stored jointly or the query must be solved as a program.

A per-pixel distribution, however rich, remains a per-pixel statement:
per-pixel richness does not substitute for cross-unit dependence.

Everything above specialises to a two-class indicator only for exposition. Let
$K\ge2$ and let each $\nu_j$ range over a finite grid on $\Delta^{K-1}$; let
$\Omega_q$ be the number of grid states. The coupling polytope now has
$\Omega_q^m$ entries, but it is still a transportation polytope and
Propositions 1–3 hold verbatim, with $\Lambda_j$ an arbitrary law on the
simplex and the covariance in Proposition 2 read componentwise. The exceedance
$\Pr(\nu_C \ge \tau)$ for a scalar functional of the composition is again a
linear functional of the coupling, so its sharp bounds are again a linear
program. The compository structure changes the geometry, not the
identification status.

**Observation 6.** The ambiguity is a property of the *marginal* information, not of the binary special case. Storing a full composition distribution per pixel does not identify the aggregate law any more than storing a Bernoulli per pixel does; only the joint does.

### 3.7 Closed-form bounds and how good the cheap relaxation is

The closed forms below are exact under SA and only under SA. A real land-cover
field meets none of the three assumptions, because class proportions drift
across the scene, blocks straddle boundaries, and equal weights hold only on a
regular grid. For such a field there is no closed form, and two routes replace
it: the exact linear program of Section 4.1, which assumes nothing beyond the
marginals but grows as $\Omega_q^{m}$ and is practical only for small $m$; and
the positive-semidefinite relaxation, which is cheap at every $m$ and, by
Proposition 6, either exactly right or wrong by a bounded and checkable
amount. What the relaxation buys is bounded for the variance and for nothing
else: it constrains second moments, it is
**not available for threshold probabilities**, which depend on the whole joint
and have no single extremal coupling (Section 3.6). Supplementary S3 gives the
construction behind the closed forms and the reason the relaxation is a
statement about second moments only.

*Premise (SA): exchangeable two-state units with equal weights and a common
marginal $p$.*

**Proposition 5 (closed-form sharp interval for the variance; Hoeffding, 1940; Rüschendorf, 1984).** Let the $m$ units be exchangeable with a common two-state marginal, $\Pr(\nu_j = 1) = p$, and equal weights, and write $\pi = \operatorname{frac}(mp)$ for the fractional part of $mp$. Then the sharp identified interval for $\mathrm{Var}(\nu_C)$ over all couplings with the stored marginals is

$$\left[\; \frac{\pi(1-\pi)}{m^{2}} \;,\;\; p(1-p) \;\right].$$

The sharp interval is therefore two explicit numbers, and the lower one is
attained by spreading the ones as evenly as the marginals allow (the
construction is in Supplementary S3).

*Premise (SA): as in Proposition 5; the relaxation constrains only the correlation matrix.*

**Proposition 6 (the PSD relaxation has two regimes).** In the same setting write $\sigma^{2} = p(1-p)$, $\rho_{\mathrm{FH}} = -\min(p,1-p)/\max(p,1-p)$ and $\rho_{\mathrm{PSD}} = \max\{\rho_{\mathrm{FH}}, -1/(m-1)\}$. The relaxation's lower endpoint is

$$V^{-}_{\mathrm{PSD}} = \sigma^{2}\left[ \frac{1}{m} + \frac{m-1}{m}\,\rho_{\mathrm{PSD}} \right].$$

Then:

(i) if $m\min(p,1-p) \le 1$, then $\rho_{\mathrm{PSD}} = \rho_{\mathrm{FH}}$
and $V^{-}_{\mathrm{PSD}} = V^{-}$: the pairwise relaxation is exact and the
joint-mixability gap is zero. This covers $m = 2$ for every $p$, and every
block size for a sufficiently sparse class.

(ii) if $m\min(p,1-p) > 1$, then $\rho_{\mathrm{PSD}} = -1/(m-1)$,
$V^{-}_{\mathrm{PSD}} = 0$, and the gap is

$$\mathrm{gap}(m,p) \;=\; V^{-} \;=\; \frac{\pi(1-\pi)}{m^{2}}, \qquad \pi = \operatorname{frac}(mp).$$

In this regime the Fréchet–Hoeffding endpoint is negative and therefore
unusable.

Equivalently, either the pairwise relaxation is exactly right, or it reports
zero where the truth is a small positive number; and the crossover between the
two regimes is a single inequality, $m\min(p,1-p) > 1$, that a practitioner
can evaluate before using it.

**Corollary 6.1 (asymptotic sharpness).** In regime (ii),

$$0 \;\le\; \mathrm{gap}(m,p) \;\le\; \frac{1}{4m^{2}},$$

with the envelope attained when $mp$ is a half-integer. Hence the gap,
measured against the single-unit marginal variance, obeys

$$\frac{\mathrm{gap}}{\sigma^{2}} \;\le\; \frac{1}{4\,m^{2}\,p(1-p)} \;=\; O(m^{-2}),$$

so the cheap relaxation is asymptotically exact.

This quantity is **not** a relative error. In regime (ii) the relaxation
reports a lower endpoint of zero while the exact lower endpoint is strictly
positive, so no ratio of the two is defined. What is bounded is an absolute
difference, expressed in units of the single-unit marginal variance
$\sigma^{2}$.

The price of avoiding the exact program is therefore an absolute difference,
measured against the single-pixel variance, that decays like the square of the
block size. For a balanced class, $p = 1/2$, the bound is $m^{-2}$: below 1 %
once $m \ge 10$, and below 0.03 % at $m = 64$. Together with Proposition 6(i),
the practical consequence is that the exact linear program of Section 4.1 is
not needed for the block variance: the pairwise relaxation is either exact, or
its error against the single-pixel variance is $O(m^{-2})$.

---

## 4. Data and experimental design

This section describes the data, the exact solver used to check the closed
forms, the simulated fields with known truth, and the measures reported in
Section 5.

### 4.1 Exact coupling bounds and synthetic landscapes

For a fine marginal with $q$ states and values $x_1,\dots,x_q$ and
probabilities $p_1,\dots,p_q$, the coupling of $m$ units is a non-negative
tensor $P \in \mathbb{R}_{\ge0}^{q\times\dots\times q}$ whose one-dimensional
marginals equal $p$. Its entries are vectorised, and the marginal constraints
form a sparse linear system. The aggregate $\nu_C = \sum_j w_j X_j$ takes
values in a finite set, and $\Pr(\nu_C \ge \tau)$ and $\mathrm{Var}(\nu_C)$
are linear in $P$, so their extremes are linear programs solved with the HiGHS
simplex/interior-point solver via SciPy's `linprog`. The implementation uses a
sparse constraint matrix and solves both the minimisation and maximisation.
Comonotonic, countermonotonic, and independent couplings are available in
closed form and serve as sanity checks: the exact upper endpoint of variance
equals the marginal variance, and the exact lower endpoint never exceeds the
independent variance.

For the simulated and real data, $q=\Omega_q$ grows and the exact program
becomes infeasible beyond $m\approx8$ for the spatial maps; there we report
the PSD relaxation and, where a lower endpoint exists, the exact value for
small $m$ alongside it.

To obtain ground truth for both the marginals and the joint, we generated
latent Gaussian random fields on a $1024\times1024$ grid of fine cells with
isotropic correlation length $L$ controlled by spectral filtering of white
noise, thresholded at the quantile corresponding to a target proportion $p$ to
produce a binary field. Fine cells are grouped into $b\times b$ blocks ($b=2$
for $m=4$, $b=3$ for $m=9$, $b=4$ for $m=16$), and all block-level quantities
(marginals, pairwise covariances, true block variance, true exceedance) are
computed by exhaustive enumeration over blocks. Correlation lengths ranged
over $L\in\{0,0.5,1,2,4,8,16,32\}$ fine cells, and six independent seeds were
averaged. Because the field is generated, the true block variance and the true
exceedance at any threshold are known exactly, and the claim "the truth lies
inside the identified bounds" is a falsifiable check, not an assertion.

### 4.2 Real reference data

We used the ESA WorldCover 2021 v200 product at 10 m resolution as the fine
support. The producer's own global validation reports an overall accuracy of
76.7 % for this version; that figure is a global aggregate over all classes
and all validation sites, not a per-class accuracy and not an accuracy for any
window used here, and we do not use it as an independent regional validation
of anything this paper measures. Eleven $3072\times3072$ -cell windows were
read directly from the public cloud-optimised GeoTIFF on Amazon S3 with
windowed reads, so that no full-tile download was needed. The windows span
seven biomes over eight tiles (Table 1).

Throughout, *enumerated truth* is the reference field's own block value,
computed by exhaustive enumeration rather than estimated. It is a property of
the reference field and not a ground measurement, because the field is a
classified product with its own error; it is used only as the value against
which the identified bounds and the geostatistical estimators are compared,
and Section 5.9 states what that restricts.

The design was widened twice, and the reason matters. Three windows are the
original temperate-forest sites of southern Germany (N48E006, N48E009). A
first revision added two Mediterranean windows (a Portuguese *montado* on
N39W009, a Spanish *dehesa* on N36W006) and one boreal window in Finnish
Lapland on N66E024, a biome in which land-cover mapping is reported to be
hardest (Bartsch et al., 2016). That left an objection open: three biomes
cannot support a claim about biomes in general. A second boreal window in
Kainuu (N63E027) gives boreal an estimable within-biome spread, and four
windows in mainland China (tropical Xishuangbanna N21E099, subtropical Wuyi
Shan N24E117, alpine Bomi N27E093, cool-temperate Xiaoxing'anling N45E126)
supply three belts the earlier design lacked and a fourth it sampled only
indirectly. The set now spans $22^\circ$ – $67^\circ$ N, so the geographic
distortion of a fixed cell count is exercised over a wider span (Section 4.3),
and its pixel-scale variance range extends an order of magnitude further down,
because Xiaoxing'anling is a near-monospecific stand.

The windows were selected to be different landscapes rather than different
addresses on a map, and the selection was made before any result was
inspected. Three criteria were applied. *Compositional contrast*: the
Mediterranean windows are a grass–tree mosaic dominated by neither class; the
boreal windows carry open water and wetland the temperate windows lack; and
the Chinese windows run from a closed subtropical canopy at 83 % tree cover
through a tropical mosaic at 71 % to an almost pure conifer stand at 98 %, the
opposite end of the range from the Swabian Jura window at 25 % tree cover.
*Structural contrast*: the fitted variogram range over a common lag window
must differ by at least a factor of three across the set, which it does by a
wide margin (10 to 88 pixels, Section 5.6). *An external accuracy benchmark*:
the Portuguese window lies in the country for which WorldCover has been
benchmarked against national and European products by Duarte et al. (2023) and
Fonte et al. (2024), using the standard machinery of Olofsson et al. (2014),
Boschetti et al. (2019) and Tsendbazar et al. (2021), so our statements can be
read against an independent accuracy estimate.

Where a biome contributes more than one window the replicate is deliberate, so
that a biome effect can be separated from a site effect inside one biome:
temperate contributes three windows, Mediterranean two and boreal two. The
four Chinese windows contribute one each, so those belts are represented as
sites and not as biome means; Section 5.6 states what the design can and
cannot attribute to a biome, and Section 5.10 explains why the design was
widened in this direction rather than by adding a second product.

For the main experiments the class of interest is tree cover (class 10) mapped
to a binary indicator; all other classes are pooled. Table 1 records the class
mix and the fine-scale variance. The 10 m cells are grouped into $30$ m pixels
(nine cells, matching the scale of a Sentinel-2 band) and those pixels into
square blocks of $m \in \{4,9,16,25\}$, corresponding to $60$ to $150$ m.
Ground truth is computed by exhaustive enumeration over all blocks in each
window. WorldCover's own classification error is not modelled: the point of
the experiment is the geometry of change of support given a reference field,
not the absolute accuracy of the product.

| Window | Biome | Site | Dominant class (share) | Tree cover | Shannon entropy | Pixel-scale variance |
|---|---|---|---|---|---|---|
| N48E006 (26 000, 26 000) | temperate | Black Forest, Germany | tree (0.861) | 0.860 | 0.542 | 0.1066 |
| N48E006 (28 000, 20 000) | temperate | Kraichgau, Germany | tree (0.315) | 0.315 | 1.488 | 0.1887 |
| N48E009 (26 000, 26 000) | temperate | Swabian Jura, Germany | crop (0.491) | 0.246 | 1.273 | 0.1694 |
| N39W009 (18 000, 18 000) | Mediterranean | Alentejo, Portugal | grass (0.477) | 0.379 | 1.143 | 0.2051 |
| N36W006 (15 000, 15 000) | Mediterranean | Andalusia, Spain | tree (0.545) | 0.541 | 1.158 | 0.2240 |
| N66E024 (26 000, 20 000) | boreal | Finnish Lapland | tree (0.760) | 0.758 | 0.764 | 0.1684 |
| N63E027 (14 000, 6 000) | boreal | Kainuu, Finland | tree (0.799) | 0.797 | 0.661 | 0.1454 |
| N21E099 (24 000, 12 000) | tropical | Xishuangbanna, China | tree (0.708) | 0.705 | 0.853 | 0.1848 |
| N24E117 (14 000, 6 000) | subtropical | Wuyi Shan, China | tree (0.838) | 0.834 | 0.599 | 0.1164 |
| N27E093 (10 000, 16 000) | alpine | Bomi, Tibet, China | tree (0.609) | 0.607 | 1.085 | 0.2193 |
| N45E126 (2 000, 25 000) | cool-temperate | Xiaoxing'anling, China | tree (0.983) | 0.983 | 0.102 | 0.0157 |

*Table 1. The eleven reference windows, grouped by biome. "Tree cover" is the
mean 30 m tree-cover fraction; "Shannon entropy" is computed over the 30 m
majority class; "Pixel-scale variance" is the variance of the 30 m tree-cover
fraction, the quantity that sets the scale of every variance below. The
cool-temperate window is a near-monospecific conifer stand whose pixel-scale
variance is an order of magnitude below the rest of the set; that is a property
of the landscape rather than of the product, and it is the one window in which
the training-block estimator of Section 5.4 fails (Section 5.9).*

Figure 4 shows the eleven windows that make up the study design, drawn in one
common style so that the panels are directly comparable.

![](../figures_v10/fig13_reference_windows.png)

**Figure 4.** The eleven reference windows are eleven different landscapes,
not eleven addresses in one biome. Within-pixel tree-cover fraction at 30 m,
one panel per window, with a shared colour scale and a biome key inset. *Top
row, temperate*: Black Forest is the closed-canopy extreme at 86 % tree cover,
Kraichgau the agricultural low at 31 %, Swabian Jura the crop–forest mixture
dominated by neither class. *Middle row, dry and high-latitude*: the two
Mediterranean windows exercise a grass–tree mosaic, and the two boreal windows
a short growing season with open water; having two of each gives those biomes
an estimable within-biome spread. *Bottom row, China*: the four Chinese
windows supply the belts the earlier design lacked and the two ends of the
compositional range, from the 83 % subtropical canopy to the 98 %
cool-temperate stand that sets the variance floor. The set spans $22^\circ$ –
$67^\circ$ N, 0.016 to 0.224 in pixel-scale variance and 0.10 to 1.49 in
Shannon entropy.

### 4.3 Metrics, the geostatistical comparator, and the design effect

We report four quantities. The **identified width** is the length of the sharp
interval, at level A or B, for the functional of interest. The
**aggregation-ambiguity index** normalises that width by the sharp upper
endpoint (Section 3.5). The **achieved variance-reduction factor** $R$ is the
ratio of the fine-unit variance to the true block variance, to be compared
against the $R=m$ promised by independence. **Empirical coverage** is the
fraction of blocks whose true 90 % interval contains the truth, where the
interval is built under the independence assumption from the stored marginals;
a nominal 90 % interval with 5 % empirical coverage is a quantitative
statement that the metadata are misleading.

All experiments are scripted and every figure is generated from the resulting
CSV files, so no number in this paper is transcribed by hand.

The bounds of this paper need only the marginals. The estimator a practitioner
would reach for instead is a geostatistical model (Matheron, 1963; Cressie,
1993; Goovaerts, 1997; Webster & Oliver, 2007), which needs the marginals
*and* a variogram, so the two routes must be compared on the same field. We
estimate the empirical variogram of the $30$ m tree-fraction field from its
two- dimensional FFT autocovariance and average it into integer lag bins. An
isotropic exponential model $\gamma(h) = c_0 + c_1(1 - e^{-h/a})$ is then
fitted by non-linear least squares, and a nested two-structure model
$\gamma(h) = c_0 + c_1(1 - e^{-h/a_1}) + c_2(1 - e^{-h/a_2})$ is fitted
alongside it for comparison. The model-based variance of a block mean is

$$\mathrm{Var}(\bar{\nu}_B) \;=\; \frac{1}{m^{2}}\sum_{i,j} C(h_{ij}), \qquad C(h) = c_0 + c_1 - \gamma(h)\ \ (h>0),\quad C(0) = c_0 + c_1.$$

The area-to-point framework of Kyriakidis (2004) makes the same
change-of-support step at the level of the variogram, and the two routes are
complementary rather than competing: area-to-point kriging imposes a model and
returns a point value, while the identified set imposes nothing beyond the
stored marginals and returns the range the stored information leaves possible.
A model converts the set into a number and buys precision at the cost of an
assumption; the set says what the model's number would have to beat. The
comparison is run on the same fields, with the same block sizes and the same
target, in Section 5.6.

The variance is then carried into a decision. For a regional mean estimated
from $n$ pixels of block size $m$, the design effect is

$$\mathrm{DEFF} \;=\; \frac{m\,\mathrm{Var}(\bar{\nu}_B)}{\sigma^{2}}, \qquad n \;=\; \frac{\sigma^{2}\,\mathrm{DEFF}}{\varepsilon^{2}}$$

for a target precision $\pm\varepsilon$ on the regional share. DEFF equals $1$
under independence and $m$ under perfect within-block correlation, so it is
exactly the quantity that turns a variance statement into a sampling cost.

Because the eleven windows differ in biome as well as in place, the pipeline
above is run unchanged on all of them: the same window size, the same
$3\times3$ aggregation, the same block sizes, the same four estimators and the
same precision target, which turns the biome into a controlled factor rather
than a narrative device. One geometric consequence of using a geographic
reference grid is reported rather than corrected: in EPSG:4326 a fixed number
of cells spans a constant north–south distance but a shrinking east–west one,
so a $3\times3$ block is 27.8 m north–south everywhere but 25.8 m east–west at
21.9°N and only 11.0 m at 66.7°N, and the block aspect ratio runs from 1.08 in
the tropical window to 2.53 in the Lapland window. Reporting the ground
dimensions lets the reader judge whether any conclusion is an artefact of
block geometry.

---

## 5. Results and discussion

Results and their reading are reported together, because each subsection
answers one question and the interpretation is what the numbers are for.
Subsections 5.1 to 5.6 report the experiments; 5.7 to 5.10 relate them to
prior work, to practice, and to the limits of the design.

Two terms recur below and are defined once. An *identified interval* is the
range a functional can take across every coupling consistent with the stored
marginals (Section 2.3); the *joint-mixability gap* is the amount by which the
PSD relaxation is wider than the exact program (Section 3.4). Throughout this
section, **ACI(m)** denotes the *aggregation-ambiguity index* defined in
Section 3.5: the width of the identified interval for the block variance,
normalised by its sharp upper endpoint. ACI = 0 means the variance is
identified; ACI = 1 means the stored marginals leave it completely
undetermined. Widths of exceedance intervals are reported in absolute
probability units.

### 5.1 What the stored marginals identify, and what they leave open

The first result is the size of the level-A identified set for the aggregate
variance. Figure 5a shows, for two-state marginals with $p=0.5$, the sharp
interval $[V^{-},V^{+}]$ of $\mathrm{Var}(\nu_C)$ as the number of aggregated
pixels grows. The upper endpoint is constant at $\sigma^2$, the single-pixel
variance; the lower endpoint falls to zero; and the independent coupling sits
at the bottom of the interval, at $\sigma^2/m$, and is therefore the most
optimistic admissible scenario. The ratio of the upper endpoint to the
independent variance is exactly $m$ for every $m$ and every $p$ tested:
Proposition 4 restated numerically.

Figure 5b makes the growth of ambiguity explicit: the exact ACI reaches 1.0
for even $m$ at $p=0.5$, because $mp\in\mathbb{Z}$ allows a deterministic
aggregate, and tracks the worst-case bound $1-1/m$ from below for unequal
proportions (at $p=0.2$ it is 0.625, 0.938 and 0.996 at $m=2,4,16$, against
0.5, 0.75 and 0.9375). The bound is tight in the sense that it is attained in
the limit, but it is a lower bound and not an estimate. The practical reading
is that a raster shipping only per-pixel variances, with no cross-pixel
dependence, cannot bound the variance of a $4\times4$ block more tightly than
the interval from zero to the full single-pixel variance: the interval does
not shrink as the raster stores more pixels, it grows in relative terms.

The practical reading is that a raster which ships only per-pixel variances,
with no cross-pixel dependence, cannot bound the variance of a $4\times4$
block more tightly than the interval from zero to the full single-pixel
variance. The interval does not shrink as the raster stores more pixels; it
grows in relative terms.

![](../figures_v10/fig1_scale_amplification.png)

**Figure 5.** Aggregating more pixels widens, not narrows, the identified set. (a) Sharp level-A bounds on $\mathrm{Var}(\nu_C)$ against $m$; the upper endpoint is the single-pixel variance and the independent coupling sits at $\sigma^2/m$. (b) The exact $\mathrm{ACI}(m)$ against the worst-case bound $1-1/m$.

Exceedance probabilities behave differently from variance. Figure 6a plots the
sharp bounds on $\Pr(\nu_C \ge \tau)$ for $m=4$, $p=0.5$ as $\tau$ sweeps the
unit interval, together with the comonotonic, countermonotonic and independent
laws.

Four regimes are visible in Figure 6a, and they are confirmed by the
underlying enumeration. For $\tau \le 0.24$ the sharp lower bound equals the
comonotonic law (0.5) and the upper bound the countermonotonic law (1.0). For
$0.26 \le \tau \le 0.48$ the lower bound is an interior value (1/3) that
neither coupling produces, while the countermonotonic law remains the upper
extremal. For $0.51 \le \tau \le 0.73$ the roles swap: the countermonotonic
law is the lower extremal (0) and the upper bound is interior (2/3). For
$\tau \ge 0.75$ the comonotonic law is the upper extremal (0.5) and the
countermonotonic law the lower extremal (0).

Two consequences follow. The independent coupling is never extremal: at
$\tau=0.5$ it predicts 0.3125 while the truth can lie anywhere in $[0, 0.667]$.
And the comonotonic and countermonotonic couplings, the extreme dependence
structures in the usual sense, are extremal only at the ends of the threshold
range and strictly interior in the middle, so a practitioner cannot
pre-compute one adversarial scene and reuse it: the adversarial scene depends
on the query (Observation 5).

Figure 6b quantifies the value of moving up the ladder for the exceedance
functional at $m=4$. At level A the identified widths are 0.556, 0.742 and
0.667 for $\tau=0.25,0.5,0.75$. Supplying the true pairwise covariance matrix
(level B) reduces these to 0.208, 0.211 and 0.217, reductions of 63 %, 72 %
and 67 %. The true value always lies inside the level-A and level-B intervals.
At level B the *variance* is not merely narrowed but identified exactly: the
width of the variance interval is zero to machine precision.

The ladder therefore has a sharp shape. Pairwise covariance is worth roughly a
two-thirds reduction in interval width for exceedance, and is sufficient in
principle for exact variance; the residual level-B width for exceedance
reflects the higher-order structure that pair information cannot see. Full
joint information (level C) closes the interval, but at the cost of storing a
$\Omega_q^m$ -dimensional object, which is precisely what a per-pixel raster
does not do.

![](../figures_v10/fig2_envelope_and_ladder.png)

**Figure 6.** No coupling is extremal at every threshold. (a) Sharp exceedance bounds for $m=4$, $p=0.5$, with the comonotonic, countermonotonic and independent couplings; four regimes are visible. (b) The identification ladder at $m=4$: levels A, B and C, with the true value marked.

### 5.2 Simulated landscapes: the truth stays inside, independence does not

These landscapes are a controlled surrogate, not a model of real land cover;
Section 5.9 states where that idealisation stops applying. Two things are
checked against them, because here the truth is known by construction rather
than measured: whether the truth stays inside the identified bounds, and how
badly an interval built under independence fails.

The first check is the one that must hold and could have failed: the truth
lies inside the identified bounds in every configuration, and at $m=4$, $L=16$
it is 0.236 against an upper endpoint of 0.244, so the sharp bound becomes
tight under strong autocorrelation (Supplementary Figure 7a). The second is
the uncomfortable one. Independence promises a variance-reduction factor $R=m$,
and the landscape delivers it only when the field is white: at $L=0$ the
achieved $R$ is 4.01, 8.95 and 16.06 for $m=4,9,16$, and at $L=32$ it is 1.02,
1.03 and 1.05, so aggregating sixteen pixels removes almost no variance when
the field is smooth. Supplementary Figure 7b and 7c show the same failure for
a threshold query, and the sweep of correlation lengths and seeds is in
Supplementary S4.

If a practitioner ignores dependence and builds an interval from the stored
per-pixel marginals, how often does that interval contain the truth? Figure 8a
answers by computing the empirical coverage of a nominal 90 % interval as the
correlation length grows. At $L=0$ the coverage is at its nominal level
(0.900, 0.940 for $m=9,16$ at $p=0.5$). It degrades monotonically: at $L=1$ it
falls to 0.71–0.74, at $L=4$ to 0.24–0.33, and at $L=32$ to 0.032 and 0.036
for $m=16$ and $m=9$. A nominal 90 % interval that covers 3 % of the time is
not merely conservative or anti-conservative in a mild sense; it is the wrong
object.

Figure 8b shows that the interval is also the wrong *shape*, not just the
wrong width. The ratio of the nominal interval width to the true interval
width falls to 0.31–0.44 in the strongly correlated regime. Independence
simultaneously centres the interval in the wrong place and gives it the wrong
length, so no simple rescaling repairs it. This is the calibration analogue of
the covariance-decomposition result: the failure is structural, not a matter
of choosing a better constant.

![](../figures_v10/fig4_calibration.png)

**Figure 8.** A nominal 90 % interval built under independence covers as little as 3 %. (a) Empirical coverage against correlation length. (b) The interval is the wrong shape as well as the wrong size.

### 5.3 A classified reference field: the gap is not an artefact of the simulation

Every number in this section is a property of one reference field and eleven
windows, not a universal constant. The identifiability results of Section 3
are theorems and hold for any product and any landscape; the empirical
quantities that follow (the $3.4$ to $21.2\times$ understatement of block
variance, the $3$ to $7$ \% coverage of a nominal $90$ \% interval of Section
5.2) are measured on ESA WorldCover 2021 v200 alone and move with the product,
the class definition and the landscape. One further limit on reading them:
WorldCover pixels are class labels, not samples from a per-pixel posterior
distribution of within-pixel composition, so this section illustrates the
spatial-dependence mechanism on a classified reference field and does **not**
validate the posterior composition distributions a probabilistic-compositional
raster would carry.

Figure 9 repeats the analysis on ESA WorldCover: panel (a) the 10 m reference
window with the 30 m pixel grid, panel (b) the stored marginal law of the
within-pixel tree fraction, panel (c) the true block variance against the
identified bounds and the independent estimate. For the first window (86.0 %
tree cover, pixel variance 0.1066) the true block variance falls from 0.0957
at $m=4$ to 0.0794 at $m=25$, while independence predicts 0.0266 and 0.0043,
so independence understates it by a factor of 3.59 at four pixels and 18.63 at
twenty-five, and the level-A ACI is 0.889 at $m=4$ and 1.000 from $m=9$
onward. Over all eleven windows the same quantity runs from $\times3.38$ to
$\times3.76$ at four pixels and from $\times14.95$ to $\times21.20$ at
twenty-five, growing steadily with the block as Proposition 4 requires. The
spread across seven biomes is narrow (the ratio of largest to smallest is only
1.11 at $m=4$, widening to 1.42 at $m=25$), so the *relative* size of the
error is a property of the aggregation geometry rather than of the landscape,
even though absolute block variance varies by a factor of 14. Here and
throughout, *enumerated truth* is the block variance of the aggregated
reference field, not ground-measured truth.[^enum] It lies inside the
identified bounds in all 44 cases, and the level-A ACI is 0.764–0.999 at $m=4$
and never below 0.906 from $m=9$ onward.

The exceedance query behaves the same way on real data, and here the direction
of the error is not even consistent. At $\tau=0.5$ and $m=4$, the Black Forest
window's true exceedance is 0.863 against an independent estimate of 0.985,
with sharp bounds $[0.761,1.0]$; Swabian Jura's is 0.241 against 0.196, with
bounds $[0,0.492]$; and Xishuangbanna's is 0.711 against 0.886, with bounds
$[0.520,1.0]$. Across all eleven windows the independent estimate runs high by
up to 0.336 at $\tau=0.25$ and low by up to 0.216 at $\tau=0.75$, on the wrong
side of the truth in one regime of thresholds and the right side in another,
the signature of using one coupling to stand in for an unknown joint. The
sharp band is wide at every threshold (mean width 0.410, 0.444 and 0.478 at
$\tau=0.25,0.5,0.75$) and the independent estimate lies inside it in all 33
window–threshold combinations: the problem is not that independence is
impossible but that it is confidently wrong by an amount the metadata cannot
bound.

Panel 9d states the result in the units a practitioner cares about: the
variance-reduction factor actually achieved on real land cover is 1.109,
1.205, 1.290 and 1.368 at $m=4,9,16,25$ averaged over the eleven windows,
against the 4, 9, 16 and 25 promised by independence. The reference field has
a pixel-scale lag-1 autocorrelation between 0.824 and 0.956, so this is not an
extreme scenario.

![](../figures_v10/fig5_real_landcover.png)

**Figure 9.** The gap is not an artefact of the simulation. (a) ESA WorldCover
10 m reference with the 30 m pixel grid, for the Black Forest window. (b) The
stored marginal law of the within-pixel tree fraction, the object a per-pixel probabilistic raster would ship, with its mean and variance annotated. (c) True
block variance against the identified bounds and the independent estimate, for
that window. (d) The variance-reduction factor achieved against the one
independence implies; the labels give the factor by which independence
understates the truth.

[^enum]: It is computed by exhaustive enumeration over ESA WorldCover 2021 v200 (Section 4.2), so it inherits the product's classification error and is a property of the reference field rather than an estimate of on-the-ground composition (Section 5.9).

### 5.4 One stored dependence scalar: identified at one block size, wrong at another

At level A the identified interval for the block variance has width
0.015–0.224 on the real windows, the whole range the variance could take.
Adding a *single* dependence parameter, the exchangeable correlation $\rho$,
reduces that width to exactly zero: the block variance becomes
point-identified, as Corollary 2.1 requires, because level B is exactly one
stored scalar in the exchangeable case. Supplementary Figure 10 shows both
sides of the resulting contract. The parameter matters — across the
PSD-feasible range it moves an exceedance prediction across a band of width up
to 0.46 — so one scalar is an attractive proposal (panel a). But it does not
close the interval: adding the *true* $\rho$ narrows the level-A widths by 54
%, 74 % and 66 % at $\tau=0.25,0.5,0.75$ without closing them, and the
PSD-feasible value $\rho=-1/3$ drives the width to approximately zero while
*excluding* the true value (panel b). A single parameter is therefore useful
only if it is estimated honestly; a wrong one produces false confidence. What
it does **not** buy is any threshold query at all, because the residual width
is not a modelling residual but the higher-order structure one pair parameter
cannot see.

A second half is missing as well, and it is the one that mattered most when we
tested the contract rather than asserting it. Corollary 2.1 is stated for a
fixed block size: pair covariance identifies the block variance *at that $m$
*. The exchangeable model additionally implies that a single $\rho$ reproduces
the variance at every $m$, because in that model
$\mathrm{Var}(\nu_C)=\sigma^{2}[1+(m-1)\rho]/m$ and hence
$\mathrm{DEFF}=m\,\mathrm{Var}(\nu_C)/\sigma^{2}=1+(m-1)\rho$. Under
exchangeability the parameter that reproduces the enumerated truth,
$$\rho_{\mathrm{exact}}(m) \;=\; \frac{m\,V_{\mathrm{true}}(m)/\sigma^{2}-1}{m-1},$$
must therefore be **the same at every block size**. It is not: Figure 11a
plots $\rho_{\mathrm{exact}}(m)$ for the eleven windows, and in 11 of 11 it
falls monotonically, from 0.794–0.918 at $m=4$ to 0.581–0.842 at $m=25$, a
spread of 0.076 to 0.212 inside a single window.

Under exchangeability every line in Figure 11a would be flat; every line
slopes down, and that is the point. Under a stationary equicorrelated model
the parameter that reproduces the saturated block variance is fixed by the
covariance function and cannot move with $m$, so a monotone decline is direct
evidence that the exchangeable model is misspecified for these fields, and
what it measures is the non-Gaussian, non-stationary character of real land
cover: a $150$ m block straddles field boundaries, forest edges and riparian
strips that a $60$ m block does not, so the effective dependence among its
pixels is weaker. That reading is an interpretation rather than a measurement
— non-stationary class proportions, a non-Gaussian patch geometry and a
covariance that is not equicorrelated at any scale cannot be separated with
eleven windows, and do not need to be, because each produces the same
practical consequence, which Figure 11b measures.

The consequence is quantitative, and it is a bias rather than an
infeasibility. Figure 11b shows the error at $m=25$ from two plausible
scalars. A scalar stored at $m=4$ and applied at $m=25$, the reading the
exchangeable contract invites since it says nothing about support, overstates
the block variance by **8.6 % to 34.1 %** (median 19.0 %), because
$\rho_{\mathrm{exact}}(4)$ exceeds $\rho_{\mathrm{exact}}(25)$ by a factor of
1.09 to 1.37. The fine-scale lag-1 autocorrelation, the scalar a producer
would reach for when no exchangeable fit is available, is worse, overstating
by **10.2 % to 38.9 %** (median 23.2 %). Both lie inside the PSD-feasible
range in 11 of 11 windows, so the defect is not that the stored number is
impossible but that it is *biased*: no feasibility check fires, and the
resulting interval is not wide but displaced. Figure 11c translates that into
the currency of Section 5.6: the lag-1 scalar inflates the sample size by 1.10
to 1.39 times, so a survey sized on it is under-sized by up to 39 %.

The exchangeable scalar is therefore best read as a cheap upgrade with a known
failure mode: the cheapest sufficient metadata extension we know of, but a
contract with two clauses. The scalar must be *estimated, never declared*; it
must be *stored with the block size at which it was measured*, because it does
not transfer; and it should be checked against a resolved variogram before it
is exported. Storing it without a stated support converts a wide honest
interval into a narrow dishonest one, the same failure mode the paper
diagnoses for independence but easier to miss, because the interval looks
correctly narrow.

![](../figures_v10/fig_rho_transfer.png)

**Figure 11.** A stored dependence scalar identifies the variance at one block
size but does not transfer to another. (a) The exchangeable correlation that
reproduces the enumerated block variance, $\rho_{\mathrm{exact}}(m)$, for each
of the eleven windows, coloured by biome; under exchangeability every line
would be flat, and all eleven fall monotonically, with a within-window spread
of 0.076 to 0.212 between $m=4$ and $m=25$. (b) Error in block variance at
$m=25$ from reusing the $m=4$ scalar (filled bars, $+8.6$ % to $+34.1$ %) and
from the fine-scale lag-1 autocorrelation (hatched bars, $+10.2$ % to $+38.9$
%). (c) The sample size implied by the lag-1 scalar as a multiple of the true
requirement, $1.10\times$ to $1.39\times$. Both candidate scalars are
PSD-feasible in every window, so the defect is bias and not infeasibility.

A stored scalar would still be useful if it could be estimated rather than
declared, and the obvious estimator is the sample variance of the block means
of a small number of observed coarse units. Supplementary Figure 12 asks how
many training blocks that takes: a few dozen contract the identified band to a
usable estimate on all but one window (median 90th-percentile relative error
falls from 100 % at $K=4$ to 13 % at $K=256$, against the 70.4–95.3 % that
independence carries at every $K$). The exception is the one window whose band
is already narrow, the near-monospecific cool-temperate stand at
Xiaoxing'anling, where no affordable number of blocks beats the bound. The
boundary of the estimator is therefore the width of the band it must beat
rather than the sample size.

### 5.5 Compositions on a simplex, and the cost of the cheap relaxation

Section 3.7 states when the exact linear program can be replaced, and
Supplementary Figure 13 tests that statement: the joint-mixability gap lies
under the envelope $1/(4m^{2})$ and decays with the same slope (panel a);
expressed against the single-pixel variance it never exceeds $m^{-2}$ at
$p=1/2$, so it falls below 1 % once $m\ge10$ and below 0.03 % at $m=64$ (panel
b); and it is identically zero where $m\min(p,1-p)\le1$, then rises to a
bounded plateau of height $1/(4m^{2})$ (panel c).

The closed forms were checked against the exact linear program over 42 $(p,m)$
combinations with $m\le8$; every endpoint matched to machine precision
($2.1\times10^{-17}$ or better), and the stated regime predicate held in every
case (Supplementary S5).

In practice, the combinatorial cost of the exact program can be avoided: for a
variance query the closed form of Proposition 5 is exact and free, and where a
relaxation is preferred Proposition 6 bounds its price at orders of magnitude
below the accuracy of any uncertainty estimate a real product would carry. The
exact linear program of Section 4.1 is therefore best understood as a
verification tool for small $m$ rather than as a production method.

A natural objection is that the ambiguity comes from collapsing a composition
to one number. Section 3.6 says otherwise, and the simplex experiments of
Supplementary S5 test it on three genuine classes: the sharp interval for a
threshold on Simpson diversity widens with $m$ rather than narrowing, reaching
81 % of the whole probability range at $m=5$, and on real WorldCover
compositions the truth lies inside a band covering 83–89 % of its range while
independence overstates block diversity by a factor of 2.3 to 2.6. Observation
6 is therefore not a consequence of the binary case: a raster that stores
per-pixel composition posteriors and assumes independence gives a block
diversity that is systematically and substantially too high.

### 5.6 The geostatistical route, and what it costs a decision

The bounds are assumption-light: they need only the marginals. The estimator a
practitioner would use instead is a geostatistical model (Cressie, 1993;
Goovaerts, 1997; Chiles & Delfiner, 2012), which needs the marginals and a
variogram. The two routes must therefore be compared on the same fields,
because the answer decides whether the identification result is a practical
improvement or merely a restatement. We therefore repeat the comparison on all
eleven reference windows, giving 44 cases (11 windows $\times$ 4 block sizes)
drawn from seven biomes.

Figure 15a isolates the mechanism. The normalised empirical variograms neither
coincide nor order themselves by biome; what governs the damage a truncated
fit does is the fitted range. The two boreal windows decorrelate fastest (10
and 16 pixels), so at 60 pixels their variograms have reached 0.99 and 0.95 of
the sill and a truncated fit has almost nothing left to miss; the alpine
window decorrelates slowest (88 pixels), so at 60 pixels only 0.60 of the sill
is visible, with the Mediterranean windows near it at 68 and 43 pixels. In all
eleven windows a single exponential fitted over the truncated window lands its
sill below the field variance, by 3.4 % to 26.3 %, and because the variance of
a block mean is built almost entirely from separations below
$\sqrt{2(m-1)} \approx 5$ pixels, that error falls exactly where it hurts.

Figure 15b collects all 44 cases. The truncated fit is biased low by 9.9 % to
42.8 % (mean $-26.0$ %). Widening the lag window to 300 pixels brings the
fitted sill to within $-6.7$ % to $+2.6$ % of the field variance and reduces
the bias to 7.4–39.7 % (mean $-20.1$ %), but a single exponential still
mis-fits the origin. A nested two-structure model (short ranges 4–39 pixels,
long ranges 18–300 pixels, $R^{2}=0.993$ – $0.999$) straddles zero, with
errors from $-7.8$ % to $+18.4$ %. Reading the variance directly off the
binned empirical variogram, with no functional form imposed, is the most
accurate and by far the most stable: all 44 cases land between $-4.9$ % and
$-1.3$ %, with a mean of $-3.0$ %. For comparison, independence is off by
70.4–95.3 % (mean $-86.3$ %) and its error grows with block size. In all 44
cases the enumerated truth lies inside the sharp identified band.

One numerical hazard is invisible to the fit statistic. Fitting the nested
model without constraining the long range lets it escape the data: in the
alpine window the unconstrained fit places $a_2$ at 2979 pixels, an order of
magnitude beyond the 300-pixel window it was fitted to, while still returning
$R^{2}=0.998$, and its block variance at $m=16$ is out by 303 %. The repair is
to bound $a_2$ by the maximum fitted lag, which is what the reported nested
fit does ($a_2=300$ pixels, error $+16.6$ % at $m=16$): a range larger than
the window it was estimated over is not a range estimate, and an $R^{2}$
computed over that window cannot detect the problem. The comparison is
therefore not "geostatistics versus bounds" but a statement about what each
route needs. The identified set requires nothing beyond the stored marginals
and is never wrong, but it is wide; the geostatistical route is
point-identified and, once the range is resolved over an adequate lag window,
accurate to about 3 %.

A variance is an intermediate quantity; the decision it feeds is usually a
sample size. With a target of $\pm1$ percentage point on the regional mean,
the independence-based requirement is 157 to 2240 pixel observations depending
on the window, the same for every block size because independence makes the
design effect a constant equal to 1. Using the enumerated truth instead
multiplies these by 3.4 to 3.8 at four pixels and by 14.9 to 21.2 at
twenty-five (Figure 15c; estimated from a variogram fitted over the wide
window the multiples are 2.9 to 15.1, and over the truncated window 2.7 to
13.7, so the truncated fit understates the requirement by roughly a quarter).
A survey designed on the independence assumption is therefore under-sized by a
factor of up to 21, and the bound and the variogram bracket the answer from
opposite sides.

![](../figures_v10/fig9_baseline_decision.png)

**Figure 15.** A variogram resolved over an adequate lag window recovers the
block variance; one truncated to 60 pixels does not. (a) Normalised empirical
variograms for the eleven windows, with the 60-pixel fit window shaded; the
two boreal curves are already at the sill inside it and the alpine curve is
furthest from it. (b) Error against the enumerated block variance for five
estimators over all 44 cases, with the bar marking the median. (c) Pixel
samples needed for a $\pm1$ percentage-point target, as a multiple of the
independence requirement. (d) Mean absolute error by biome, which orders by
fitted range rather than by sample size; the number in brackets is the count
of windows contributing to each biome. Since four of the seven biomes
contribute one window, the ordering is a site ordering that covaries with
biome.

Do the failure modes belong to the product, the place, or the mechanism? With
the caveat that four of the seven biomes contribute a single window each, so
nothing here separates a biome effect from a site effect: the bound carries no
biome information at all, the enumerated truth lying inside the sharp band in
all 44 cases across a fourteenfold change in field variance (0.0157 to 0.2240)
and a change in class mix from 25 % to 98 % tree cover; and the relative
independence error is likewise flat, 84.8 % to 87.3 %, because it is set by
the ratio of block variance to field variance. The truncated-window error is
the one quantity that is graded, and it is graded by the fitted range: the
fastest-decorrelating sites carry a mean absolute bias of 13.3 % and the two
Mediterranean sites 38.2 %, while the nested model (2.9 % to 15.6 %) and the
non-parametric estimator (1.5 % to 4.5 %) nearly remove that spread. Because
fitted range and biome label covary across these eleven windows, that ordering
is correlated with biome in this sample without being caused by it; what the
data do show is that the grading variable is continuous rather than
categorical (Supplementary Figure 16: the bias of the 60-pixel fit is monotone
in the fraction of the sill resolved at 60 pixels, $R^{2}=0.72$), so the
plainest reading is that neither fit uses the biome label as such. The full
argument, including the sample-size penalty and the Mediterranean exception,
is in Supplementary S6.

Replication is nonetheless uneven, and Section 5.9 sets out the consequences:
tropical, subtropical, alpine and cool-temperate are represented by one site
each, so for them nothing here separates a biome effect from the particular
site, and no statement of the form "biome X has property Y" is supported by
this design, only statements about the sites measured.

### 5.7 What the results mean, and how they relate to prior work

The experiments give one diagnosis from several angles. The mean of a coarse
unit is safe: it is a linear functional of what is stored (Proposition 1).
Everything else is not. The variance is identified exactly by pair covariances
(Proposition 2, Corollary 2.1) but not by the per-pixel variances a
probabilistic raster typically stores, and the shortfall is quantified by the
aggregation-ambiguity index, bounded below by $1-1/m$ (Proposition 4).
Threshold queries are worse than variance, because their sharp bounds have no
single extremal coupling and the worst case depends on the threshold
(Observation 5, Section 5.1). All of this survives the move from an indicator
to a simplex composition (Observation 6, Section 5.5). In one sentence,
**scale-consistency is a representational requirement, not a computational one**:
the mathematics of combining uncertain units under change of support is
classical and settled, and what is unsettled is whether a raster that
advertises per-pixel uncertainty carries the information that makes the
combination well posed. On real land cover it does not, and the discrepancy is
large — a factor of 3.4 to 21.2 in block variance over eleven windows in seven
biomes, and a nominal 90 % interval that covers 3–7 % of the time at realistic
correlation lengths. Section 5.6 sharpens the diagnosis: read
non-parametrically off a resolved lag window, the standard route recovers the
block variance to within 4.9 %. What decides the outcome is therefore not the
mathematics of change of support but the assumption that comes with it, an
unresolved variogram or an independence that was never stated; and Section 5.6
adds that the one quantity which is graded, the bias of a truncated lag
window, is graded by the fitted range and is therefore correlated with, but
not caused by, biome in this sample. Section 5.4 adds a second, less obvious
version of the same failure: the one-parameter fix is real, but the parameter
is support-specific, so a single stored scalar is a correct answer about one
block size and a wrong answer about another, and the cost is paid in samples —
a $\pm1$ percentage-point target needs 3.4 to 21.2 times more pixels than
independence implies.

In one sentence,
**scale-consistency is a representational requirement, not a computational one**.
The mathematics of combining uncertain units under change of support is
classical and settled; what is unsettled is whether a raster that advertises
per-pixel uncertainty actually carries the information that makes the
combination well posed. On real land cover it does not, and the discrepancy is
large, a factor of 3.4 to 21.2 in block variance over eleven windows in seven
biomes, and a nominal 90 % interval that covers 3–7 % of the time at realistic
correlation lengths. The geostatistical comparator of Section 5.6 sharpens the
diagnosis: read non-parametrically off a resolved lag window, the standard
route recovers the block variance to within 4.9 %. What decides the outcome is
therefore not the mathematics of change of support but the assumption that
comes with it: an unresolved variogram, or an independence that was never
stated.

Section 5.6 adds that this diagnosis is stable across the sites sampled here
except for one quantity: the bias of a truncated lag window is graded by the
fitted range, and is therefore correlated with, but not caused by, biome in
this sample, so it is largest where the landscape is coarsest and where the
sill shortfall is largest in absolute terms.

Section 5.4 adds a second, less obvious version of the same failure: the
one-parameter fix is real, but the parameter is support-specific, so a single
stored scalar is a correct answer to a question about one block size and a
wrong answer to the same question at another. And the cost is paid in samples
rather than in principle: a $\pm1$ percentage-point target needs 3.4 to 21.2
times more pixels than independence implies.

The result that ignoring spatial correlation corrupts the uncertainty of a
spatial mean is well documented, and Wadoux and Heuvelink (2023) make the case
forcefully. Our contribution is to locate the corruption in the
*representation* rather than in the *model*: we do not assume a variogram is
unavailable and then propagate a wrong one; we ask what follows from the
metadata a probabilistic raster actually stores, and the answer is a wide
identified set, not a wrong point. The partial-identification stance, borrowed
from econometrics (Molinari, 2020; Kline & Tamer, 2023), is what turns "your
interval is wrong" into "your interval is one point in a set we can compute
exactly". This is not an argument against either of the two literatures it
sits between: a model that stores or outputs the full joint is not affected by
the gap identified here, and a variogram is an effective remedy, but an
external modelling assumption rather than a consequence of the stored
marginals, and the metadata must say so.

The joint-mixability phenomenon is known in the dependence-modelling
literature (Embrechts & Puccetti, 2010; Wang & Wang, 2015) but is usually a
curiosity about artificial marginals. Here it has an everyday engineering
consequence: the cheap relaxation of a per-pixel uncertainty model is not just
loose but infeasible, producing negative variances, and the PSD relaxation
repairs feasibility at negligible cost and is exact precisely when
$mp\in\mathbb{Z}$, giving a checkable diagnostic. A second comparison isolates
the same point: what decides whether a block variance is estimated well is the
behaviour of the variogram *near the origin*, because only separations up to
$\sqrt{2(m-1)}$ pixels enter a block. A lag window too short to show the rise
to the sill mis-fits exactly that region while still returning an $R^{2}$ of
0.95–0.99, and widening it to 300 pixels removes most of the bias, so a fit
diagnostic computed over the whole curve is not diagnostic of the short-range
behaviour that matters and reporting $R^{2}$ without the lag window is
uninformative.

Finally, this work supports the long-standing argument that uncertainty
metadata must expose semantics, not just magnitudes (Cova & Goodchild, 2002;
Goodchild et al., 2012). A per-pixel variance field is a semantics-poor
object; the minimal upgrade our results support is a pairwise covariance or
one exchangeable dependence parameter — cheap, but, as Section 5.4 shows, not
sufficient alone, since the parameter must be estimated rather than declared
and stored with the block size at which it was measured, because it does not
transfer.

### 5.8 Practical recommendations

The results condense into four rules, two for whoever produces a PCR and two
for whoever uses one.

*For the data producer.* First, **do not report a variance for any aggregate finer than the whole map unless you also report the dependence structure**. The mean is safe; the variance is not. If dependence is unknown, report the identified interval, not a point. Second, **store pair structure, and store it with the support it was measured at**: pair covariance is exactly what the variance needs and captures most of the exceedance gap, and a single exchangeable parameter is enough to identify the variance outright, but only at one block size, because the reproducing parameter falls monotonically with $m$ and a scalar exported from four pixels to twenty-five overstates the block variance by 8.6 % to 34.1 % (Section 5.4). A dependence scalar without an advertised support is not metadata; it is a number carrying a hidden assumption.

*For the data user.* Third, **expect queries, not scenes**. Because the worst case for an exceedance depends on the threshold, a data product cannot pre-compute one adversarial realisation; either the joint is stored, or the query is solved as a program over the coupling polytope. Fourth, **treat the short-range variogram as the part that matters most, and the identified set as the fallback**. If a variogram is available and its range is resolved over a lag window several times longer than the practical range, it will estimate a block variance to a few per cent and should be used. If it is not, the sharp bounds require no model and are never wrong. Report them alongside the model estimate, since Section 5.6 shows the two bracket the decision-relevant quantity from opposite sides.

### 5.9 Limitations

The exact linear programs become infeasible beyond a handful of units because
the coupling polytope grows as $\Omega_q^m$, so the spatial and real-data
results report the PSD relaxation alongside exact values where $m$ is small,
and the PSD relaxation carries its own (small) joint-mixability gap when
$mp\notin\mathbb{Z}$. The spatial experiments use a latent-Gaussian-threshold
landscape, which reproduces tuneable autocorrelation but is not a calibrated
model of any particular land surface: real land cover does not strictly obey
it, its patches are not Gaussian and its class proportions are not stationary.
The real experiments treat tree cover as a binary class of interest and pool
the rest, and the multi-class compositional experiments of Section 5.5 place
each class on an $L$ -level lattice, so the aggregate law is computed on a
grid rather than over all compositions; with $K=3$ and $L=3$ the exact program
is feasible to $m=5$, and larger $K$, finer $L$ and larger $m$ are out of
reach of the exact method. The reported interval widths were, however,
insensitive to the grid, changing by at most $0.0034$ when $L$ was raised from
3 to 9. For the reader the practical consequence is *when to trust the PSD
relaxation*, and Proposition 6 gives the exact rule: it is exact whenever
$m\min(p,1-p) \le 1$ (always for $m=2$, and for any block size when the class
of interest is sparse enough) and otherwise it overstates the lower endpoint
by at most $1/(4m^{2})$, below 1 % of the marginal variance once $m \ge 10$
for a balanced class. The place to be wary is not the variance but the
*exceedance* width, where no relaxation is sharp: the relaxation is a
statement about second moments and nothing more, and it is
**not applicable to exceedance probabilities**, because $\Pr(\nu_C \ge \tau)$
depends on the whole joint and has no single extremal coupling (Section 3.6),
so an exceedance query still requires the exact linear program or a simulation
over the coupling polytope, and its answer should be reported as a bound.

The reference data carry their own limits, and they are the binding ones.
Coverage is not global: tundra, arid and wetland biomes are unsampled, and
Venter et al. (2022) report the lowest per-class accuracy of the 10 m products
precisely in tundra, so the biome most likely to strain the near-invariance of
Section 5.6 is one this design does not reach. Four of the seven biomes
contribute a single window, so for those the biome factor is varied but not
replicated and no within-biome variance component can be estimated. Every
window reads one product, so the product factor is held fixed throughout;
Section 5.10 weighs the alternative and states the audit that would bound it.
And because the reference field is the product under study, the comparison
speaks to change of support and not to map accuracy: a disagreement between
our block variance and a variogram estimate is a statement about the
estimator, not about the map. The producer's global overall accuracy of 76.7 %
(Section 4.2) is not a substitute for regional validation of any of the eleven
windows.

The residual bias depends on the fitted functional form: the non-parametric
estimator was stable to within 4.9 % across all 44 cases, but the single
exponential's undershoot and the nested model's overshoot are properties of
these fields and of those two forms rather than general constants, and the
unbounded nested fit of Section 5.6 shows that a two-structure model can
report an excellent $R^{2}$ while placing its long range outside the window it
was fitted to. Two further limits attach to the variogram route as we ran it.
The fits are isotropic, since the empirical variogram is averaged over
directions; if the range differs by direction, that average truncates the
longest direction first, and surface anisotropy would amplify the bias
precisely where the 60-pixel default cannot see it. The training-block
estimator of Section 5.4 is bounded on both sides: its boundary is the width
of the band it must beat rather than the sample size, so it is useful on a
wide-banded landscape and useless on a near-homogeneous one.

Finally, two modelling assumptions are inherited, not tested. The reference
data carry their own classification error, which is not modelled, so the
reported variances are properties of the reference field and not estimates of
true on-the-ground composition. And the analyses assume a fixed, agreed set of
composition classes; ambiguity in the class taxonomy itself is a separate and
equally real source of change-of-support error that this framework does not
address.

### 5.10 A second product, a targeted audit, and open questions

The eleven windows span seven biomes but read one product, ESA WorldCover 2021
v200, so every empirical claim is conditional on that reference field, whereas
the identifiability results (Propositions 1–6) are theorems about the stored
marginals and hold for any product. Because the driver of the empirical error
is the fitted range rather than the labelling, adding windows in unsampled
belts buys more independent variation in that quantity than adding a second
product over the same windows, which would hold the landscape fixed and vary
only the labelling. A second product remains the right next step, but as a
targeted audit of two or three windows rather than a redesign; Table 2 sets
the two options side by side, and the full accounting, with the open questions
it leaves, is in Supplementary S2.

| | Adding a second 10 m product (over the same windows) | Adding windows in unsampled belts (adopted here) |
|---|---|---|
| *What varies* | labelling, at fixed landscape and range | landscape, and with it the range and the class mix |
| Biomes covered | unchanged at seven | seven, up from three |
| Windows | unchanged at eleven | eleven, up from six |
| *Range span of the fitted variogram* | unchanged, 10–88 px | 10–88 px, up from 16–68 px |
| Effect on the theorems | none | none |
| Effect on the empirical claims | tests whether the mechanism is product-stable | widens the range and class-mix span the mechanism is tested over |
| Effect on the reference field | neither product is the reference, so the ground-truth check becomes conditional | unchanged: the truth check of Section 5.2 stays falsifiable |
| Work | re-define classes, re-project, re-read, re-enumerate for every window | the same read–cache–enumerate pipeline, already scripted and offline after caching |
| Verdict | right as a targeted audit of two or three windows | adopted in this revision |

*Table 2. The two ways to strengthen the empirical half of the paper, and why the second was chosen. The rows in italics are the decisive ones: the theorems are unaffected by either choice, so the question is only which option buys more independent variation in the quantity that drives the mechanism.*

---

## 6. Conclusion

A PCR claims to know more than a value: it claims to know a distribution. We
have shown that the claim is only partially transferable across scales, and
have located exactly where it fails. For the **producer** the implication is a
metadata contract with two clauses. The aggregate mean is identified from what
is stored, so it can be advertised; the aggregate variance is not. It is
identified exactly by one additional scalar in the exchangeable case and by
the pairwise covariance matrix in general, so storing that scalar or that
matrix is the difference between a variance that means something and one that
does not. But the scalar is support-specific, overstating the block variance
by 8.6 % to 34.1 % when exported from four pixels to twenty-five, so it must
be stored with the block size at which it holds and must be measured rather
than declared. If neither a covariance nor an honest scalar is stored, the
honest product is the sharp bound, which requires no model and is never wrong:
a wrong dependence parameter is worse than none, because it closes the
interval while excluding the truth.

For the **user** the implication is a sampling cost. On real ESA WorldCover
land cover, across eleven windows in seven biomes, a per-pixel-only
representation understates block variance by a factor of 3.4 to 21.2, and a
nominal 90 % interval covers 3–7 % of the time at realistic autocorrelation,
so a $\pm1$ percentage-point target on a regional mean needs 3.4 to 21.2 times
more pixels than independence implies. Where a variogram is available and its
range is resolved, reading the variance off the binned curve recovers the
block variance to within 4.9 %; where the lag window is the common 60-pixel
default, the same fit is biased low by 10–43 % while still reporting an
$R^{2}$ of 0.95–0.99. Where neither is available, the identified set is the
fallback, and the two bracket the answer from opposite sides; where a few
dozen training blocks can be observed, they contract the interval to a usable
estimate, except on a landscape whose block variance is already so small that
the bound is the best available statement.

The obstruction is representational rather than computational.
Scale-consistency is not something a probabilistic raster gets for free by
attaching distributions to pixels; it is a property of what those
distributions are allowed to carry, and of whether the dependence they carry
is labelled with the support at which it was measured.

## Data and codes availability statement

Every number in this paper is produced by a script and read from a CSV;
nothing is transcribed by hand. The closed forms and the relaxations are
verified against the exact linear program in `exp_psd.py` and `pcp_core.py`.
The spatial and real-data results come from `exp_spatial.py` and `exp_real.py`,
the multi-class compositional experiments from `exp_compositional.py`, the
geostatistical comparator and the decision analysis from `exp_baseline.py`,
which reads the window definitions from a single table in `windows_v2.py`, and
the training-block proof of concept from `exp_estimator.py`. Figures come from
the scripts in `figscripts_v10/`, which share one style module (`figstyle.py`)
and read only their own CSV output.

**Reproducibility conditions.** The synthetic experiments (`exp_theory.py`, `exp_spatial.py`, `exp_compositional.py`, `exp_psd.py`) are fully offline: seeds are fixed and the variogram and block-variance computations are deterministic, so every reported number is reproducible bit-for-bit with only NumPy, SciPy and pandas. The real-data experiments additionally require network access to the ESA WorldCover 2021 v200 public cloud-optimised GeoTIFF on the AWS Open Data endpoint (`esa-worldcover.s3.eu-central-1.amazonaws.com`); no credentials are needed, and each of the eleven windows is read once and cached under `data/` (159 MB in total), after which the analysis is offline. The eleven windows are named in `windows_v2.py`, which is the single source of truth for the study design, including the tile transform used to place each window's centre. All results underlying the figures are in `results_v2/`. Scripts, configuration, the cached window definitions and the CSV outputs underlying every figure are provided in the supplementary repository accompanying this manuscript.

## Acknowledgements

We thank the editors and the anonymous reviewers for comments that improved
the manuscript.

## Disclosure statement

No potential conflict of interest was reported by the author(s).

## Funding

This research received no specific grant from any funding agency in the
public, commercial, or not-for-profit sectors.

## Notes on contributors

Author names, affiliations, ORCID identifiers and short biographies are
supplied on the separate title page and are withheld from this file to
preserve double-blind review.

## Appendix A. Notation

| Symbol | Meaning |
|---|---|
| $m$ | number of fine units aggregated into one block |
| $w_j$ | area weight of fine unit $j$; the weights sum to 1 |
| $K$ | number of composition classes; in Section 5.4, number of training blocks |
| $\nu_j$ | composition of fine unit $j$, a point of the simplex $\Delta^{K-1}$ |
| $\nu_C$ | composition of the block, the area-weighted average of the fine compositions |
| $\Lambda_j$ | stored marginal law of $\nu_j$ (level A) |
| $\Gamma$ | joint law of $(\nu_1,\dots,\nu_m)$; never stored |
| $\Lambda_C$ | law of $\nu_C$, the push-forward of $\Gamma$ |
| $\sigma^{2}$ | common fine-unit variance |
| $V^{-}, V^{+}$ | sharp lower and upper endpoint of the identified set for $\mathrm{Var}(\nu_C)$ |
| $V^{-}_{\mathrm{PSD}}$ | lower endpoint of the PSD relaxation |
| $\mathrm{ACI}(m)$ | aggregation-ambiguity index, $(V^{+}-V^{-})/V^{+}$ |
| $\rho$ | exchangeable pairwise correlation |
| $\rho_{\mathrm{FH}}$ | sharp minimum pairwise correlation of a two-state marginal |
| $\rho_{\mathrm{PSD}}$ | PSD-feasible pairwise correlation, $\max\{\rho_{\mathrm{FH}},-1/(m-1)\}$ |
| $\rho_{\mathrm{exact}}(m)$ | exchangeable correlation that reproduces the block variance at block size $m$ |
| $\rho_{\mathrm{lag1}}$ | fine-scale lag-1 autocorrelation, used as a stand-in scalar |
| $\Omega_q$ | number of lattice states of one fine composition |
| $L$ | number of grid levels per class in the compositional experiments |
| $\gamma(h)$ | variogram; $c_0, c_1, a$ its nugget, partial sill and range |
| $\mathrm{DEFF}$ | design effect, $m\,\mathrm{Var}(\bar{\nu}_B)/\sigma^{2}$ |
| $R$ | achieved variance-reduction factor, $\sigma^{2}/\mathrm{Var}(\nu_C)$ |
| $\pi$ | $\operatorname{frac}(mp)$, the fractional part of the expected block sum |
| $q$ | number of states of a single fine marginal in the linear program |

## Appendix B. Proofs deferred from Section 3

**Proof of Proposition 1.** Linearity of expectation applied to $g_w$. $\square$

**Proof of Proposition 2.** Bilinearity of the covariance under a linear combination. $\square$

**Proof of Proposition 3.** For $m=2$ the transport polytope of a two-state marginal is a line segment, so its endpoints are the two monotone rearrangements. For $m\ge3$, attaining the pairwise minimum covariance for every pair simultaneously requires a correlation matrix that is not positive semidefinite; the exact optimum solves the same objective subject to the coupling polytope, whose extreme points need not realise all pairwise extrema at once. $\square$

**Proof of Proposition 4.** The comonotonic coupling sets $\nu_1=\dots=\nu_m$, so $\nu_C=\nu_1$ and $V^{+}=\sigma^2$. The independent coupling is admissible, so $V^{-} \le \mathrm{Var}_{\text{ind}}(\nu_C) = \sigma^2/m$. Hence $\mathrm{ACI} = 1 - V^{-}/V^{+} \ge 1 - (\sigma^2/m)/\sigma^2 = 1 - 1/m$. $\square$

**Proof of Proposition 5.** The upper endpoint is $p(1-p)$, attained by the comonotonic coupling exactly as in Proposition 4. For the lower endpoint, under the rotation coupling each $\nu_j = \mathbf{1}\{(U+jp)\bmod 1 < p\}$ is Bernoulli with parameter $p$, so the marginals are reproduced. The number of units equal to one is the multiplicity with which the $m$ arcs $\{(u+jp)\bmod 1 < p\}$ cover the point $u$, and because the arc endpoints are spaced exactly $p$ apart on the circle, that multiplicity is either $\lfloor mp\rfloor$ or $\lceil mp\rceil$ for every $u$. The sum therefore takes only those two values; since its mean is $mp$ it places mass $\pi=\operatorname{frac}(mp)$ on $\lceil mp\rceil$ and $1-\pi$ on $\lfloor mp\rfloor$, and its variance is $\pi(1-\pi)$, giving $\mathrm{Var}(\nu_C)=\pi(1-\pi)/m^{2}$.

For optimality, under any admissible coupling the sum $S=\sum_j \nu_j$ takes
values in $\{0,\dots,m\}$ with mean $mp$. Among all laws on that lattice with
that mean, the minimum variance is attained by the two-point law concentrated
on $\lfloor mp\rfloor$ and $\lceil mp\rceil$, whose variance is $\pi(1-\pi)$;
every other admissible law has variance at least as large. The couplings the
marginals admit produce a subset of these laws, so no coupling does better,
and the rotation coupling attains the bound. $\square$

**Proof of Proposition 6.** Substituting $\rho_{\mathrm{FH}}=-\min(p,1-p)/\max(p,1-p)$ and the value of $V^{-}$ from Proposition 5 into the equicorrelation formula $\mathrm{Var}(\nu_C)=\sigma^{2}\big[\frac1m+\frac{m-1}{m}\rho\big]$ gives the two regimes directly, because
$\rho_{\mathrm{FH}} < -\frac{1}{m-1} \iff (m-1)\min(p,1-p) > \max(p,1-p) \iff m\min(p,1-p) > 1.$
In regime (i) the PSD constraint is slack, so the relaxation returns the exact pairwise endpoint $V^{-}$. In regime (ii) the constraint binds at $\rho=-\frac1{m-1}$, where the equicorrelation formula gives exactly zero; the exact minimum is $\pi(1-\pi)/m^{2}>0$, and the Fréchet–Hoeffding endpoint is negative and therefore infeasible. $\square$

**Proof of Corollary 6.1.** In regime (ii) the gap is $\mathrm{gap}(m,p)=\pi(1-\pi)/m^{2}$ with $\pi=\operatorname{frac}(mp)\in[0,1)$, so $\mathrm{gap}\le\frac1{4m^{2}}$, with equality when $mp$ is a half-integer. Dividing by $\sigma^{2}=p(1-p)$ expresses the gap in units of the single-unit marginal variance, $\frac{\mathrm{gap}}{\sigma^{2}}\le\frac{1}{4m^{2}p(1-p)}$, which is $O(m^{-2})$ and equals $m^{-2}$ at $p=\frac12$. This is an absolute difference scaled by $\sigma^{2}$, not a relative error: in regime (ii) the relaxation reports zero while the exact endpoint is positive, so the ratio of the two is undefined. $\square$

---

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data.*
Chapman & Hall, London.

Bartsch, A., Höfler, A., Kroisleitner, C., & Trofaier, A. M. (2016). Land
cover mapping in northern high latitude permafrost regions with satellite
data: achievements and remaining challenges. *Remote Sensing*, 8(12), 979.
doi:10.3390/rs8120979

Boschetti, L., Roy, D. P., Giglio, L., Huang, H., & Zubkova, M. (2019). Global
validation of the collection 6 MODIS burned area product. *Remote Sensing of
Environment*, 235, 111490. doi:10.1016/j.rse.2019.111490

Chiles, J.-P., & Delfiner, P. (2012). *Geostatistics: Modeling Spatial
Uncertainty* (2nd ed.). Wiley. doi:10.1002/9781118136188

Cova, T. J., & Goodchild, M. F. (2002). Extending geographical representation
to include fields of spatial objects. *International Journal of Geographical
Information Science*, 16(6), 509–532. doi:10.1080/13658810210137040

Cressie, N. (1993). *Statistics for Spatial Data* (rev. ed.). Wiley.
doi:10.1002/9781119115151

Duarte, D., Fonte, C. C., Costa, H., & Caetano, M. R. (2023). Thematic
comparison between ESA WorldCover 2020 land cover product and a national land
use land cover map. *Land*, 12(2), 490. doi:10.3390/land12020490

Eches, O., Dobigeon, N., Mailhes, C., & Tourneret, J.-Y. (2010). Bayesian
estimation of linear mixtures using the normal compositional model:
application to hyperspectral imagery. *IEEE Transactions on Image Processing*,
19(6), 1403–1413. doi:10.1109/TIP.2010.2042993

Embrechts, P., & Puccetti, G. (2010). Bounds for the sum of dependent risks
having overlapping marginals. *Journal of Multivariate Analysis*, 101(1),
120–134. doi:10.1016/j.jmva.2009.07.004

Embrechts, P., Puccetti, G., & Rüschendorf, L. (2013). Model uncertainty and
VaR aggregation. *Journal of Banking & Finance*, 37(8), 2750–2764.
doi:10.1016/j.jbankfin.2013.03.014

Fonte, C. C., Duarte, D., Jesus, I., Costa, H., Benevides, P. J., et al.
(2024). Accuracy assessment and comparison of national, European and global
land use land cover maps at the national scale — case study: Portugal. *Remote
Sensing*, 16(9), 1504. doi:10.3390/rs16091504

Fréchet, M. (1951). Sur les tableaux de corrélation dont les marges sont
données. *Annales de l'Université de Lyon, Section A*, 9, 53–77.

Ge, Y., Jin, Y., Stein, A., Chen, Y., Wang, J., et al. (2019). Principles and
methods of scaling geospatial Earth science data. *Earth-Science Reviews*,
196, 102897. doi:10.1016/j.earscirev.2019.102897

Goodchild, M. F., Guo, H., Annoni, A., Bian, L., de Bie, K., Campbell, F., … &
Zhang, C. (2012). Next-generation digital earth. *Proceedings of the National
Academy of Sciences*, 109(28), 11088–11094. doi:10.1073/pnas.1202383109

Goovaerts, P. (1997). *Geostatistics for Natural Resources Evaluation.* Oxford
University Press. doi:10.1093/oso/9780195115383.001.0001

Gotway, C. A., & Young, L. J. (2002). Combining incompatible spatial data.
*Journal of the American Statistical Association*, 97(458), 632–648.
doi:10.1198/016214502760047140

Heuvelink, G. B. M., & Pebesma, E. J. (1999). Spatial aggregation and soil
process modelling. *Geoderma*, 89(1–2), 47–65.
doi:10.1016/S0016-7061(98)00077-9

Heuvelink, G. B. M., Brown, J. D., & van Loon, E. E. (2007). A probabilistic
framework for representing and simulating uncertain environmental variables.
*International Journal of Geographical Information Science*, 21(5), 497–513.
doi:10.1080/13658810601063951

Hoeffding, W. (1940). Masstabinvariante Korrelationstheorie. *Schriften des
Mathematischen Instituts und des Instituts für Angewandte Mathematik der
Universität Berlin*, 5, 181–233.

Kaas, R., Dhaene, J., & Vyncke, D. (2002). A simple geometric proof that
comonotonic risks have the convex-largest sum. *ASTIN Bulletin*, 32(1), 71–80.
doi:10.2143/AST.32.1.1015

Kattenborn, T., Schiefer, F., Frey, J., Feilhauer, H., & Mahecha, M. D.
(2022). Spatially autocorrelated training and validation samples inflate
performance assessment of convolutional neural networks. *ISPRS Open Journal
of Photogrammetry and Remote Sensing*, 5, 100018.
doi:10.1016/j.ophoto.2022.100018

Kline, B., & Tamer, E. (2023). Recent developments in partial identification.
*Annual Review of Economics*, 15, 125–150.
doi:10.1146/annurev-economics-051520-021124

Kyriakidis, P. C. (2004). A geostatistical framework for area-to-point spatial
interpolation. *Geographical Analysis*, 36(3), 259–289.
doi:10.1111/j.1538-4632.2004.tb01135.x

Matheron, G. (1963). Principles of geostatistics. *Economic Geology*, 58(8),
1246–1266. doi:10.2113/gsecongeo.58.8.1246

Molinari, F. (2020). Microeconometrics with partial identification. In
*Handbook of Econometrics*, Vol. 7A (pp. 355–486). Elsevier.
doi:10.1016/bs.hoe.2020.05.002

Nelsen, R. B. (2006). *An Introduction to Copulas* (2nd ed.). Springer.
doi:10.1007/0-387-28678-0

Olofsson, P., Foody, G. M., Herold, M., Stehman, S. V., Woodcock, C. E., &
Wulder, M. A. (2014). Good practices for estimating area and assessing
accuracy of land change. *Remote Sensing of Environment*, 148, 42–57.
doi:10.1016/j.rse.2014.02.015

Pirzamanbein, B., Lindström, J., Poska, A., & Gaillard, M.-J. (2018).
Modelling spatial compositional data: reconstructions of past land cover and
uncertainties. *Spatial Statistics*, 24, 14–31.
doi:10.1016/j.spasta.2018.03.005

Prékopa, A. (1990). Sharp bounds on probabilities using linear programming.
*Operations Research*, 38(2), 227–239. doi:10.1287/opre.38.2.227

Puccetti, G., & Rüschendorf, L. (2012). Computation of sharp bounds on the
distribution of a function of dependent risks. *Journal of Computational and
Applied Mathematics*, 236(7), 1833–1840. doi:10.1016/j.cam.2011.10.015
doi:10.1016/j.cam.2011.10.015

Rüschendorf, L. (1984). On the minimum discrimination information theorem.
*Statistics & Decisions*, Supplement 1, 263–283.

Segurado, P., Araújo, M. B., & Kunin, W. E. (2006). Consequences of spatial
autocorrelation for niche-based models. *Journal of Applied Ecology*, 43(3),
433–444. doi:10.1111/j.1365-2664.2006.01162.x

Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163(4148), 688.
doi:10.1038/163688a0

Tsendbazar, N.-E., Herold, M., Li, L., Tarko, A., de Bruin, S., et al. (2021).
Towards operational validation of annual global land cover maps. *Remote
Sensing of Environment*, 266, 112686. doi:10.1016/j.rse.2021.112686

Venter, Z. S., Barton, D. N., Chakraborty, T., Simensen, T., & Singh, G.
(2022). Global 10 m land use land cover datasets: a comparison of Dynamic
World, World Cover and Esri Land Cover. *Remote Sensing*, 14(16), 4101.
doi:10.3390/rs14164101

Wadoux, A. M. J.-C., & Heuvelink, G. B. M. (2023). Uncertainty of spatial
averages and totals of natural resource maps. *Methods in Ecology and
Evolution*, 14(5), 1320–1332. doi:10.1111/2041-210X.14106

Wang, R., & Wang, L. (2015). Extreme negative dependence and risk aggregation.
*Journal of Multivariate Analysis*, 136, 12–25. doi:10.1016/j.jmva.2015.01.006

Webster, R., & Oliver, M. A. (2007). *Geostatistics for Environmental
Scientists* (2nd ed.). Wiley. doi:10.1002/9780470517277