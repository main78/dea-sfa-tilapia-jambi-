# Breaking the Frontier: A Learning Guide to Statistical Model Reliability

> **Note:** This document is a pedagogical companion to the paper
> *"When and Why Do DEA and SFA Diverge? Evidence from Spatially
> Heterogeneous Tilapia Aquaculture"* (submitted to the *Journal of
> Productivity Analysis*). It explains the core methodological findings
> in accessible language, with reference to the exact tables and figures
> in the manuscript. All statistics cited here are reproducible from the
> scripts in this repository.

---

In the advanced study of productivity, we employ mathematical "frontiers"
to benchmark the performance of decision-making units (DMUs) against an
idealized production limit. However, the reliability of these
estimators — specifically their ability to distinguish between managerial
failure and mere statistical noise — is heavily contingent upon sample
characteristics. This guide explores the mechanical divergence between
the two primary traditions of efficiency analysis under conditions of
**spatial-ecological heterogeneity**.

---

## 1. The Two Philosophies of Efficiency: DEA vs. SFA

Efficiency measurement is bifurcated into two methodological traditions:
**Data Envelopment Analysis (DEA)**, a non-parametric approach, and
**Stochastic Frontier Analysis (SFA)**, a parametric approach. While
both seek to identify the frontier of maximum possible output, their
internal "engines" process data through fundamentally different logics.

| Dimension | Data Envelopment Analysis (DEA) | Stochastic Frontier Analysis (SFA) |
|---|---|---|
| **Mathematical Logic** | **Deterministic Envelopment:** Employs linear programming to "wrap" a geometric boundary around the best performers. | **Maximum Likelihood:** Employs statistical estimation to find a functional form that best represents the production process. |
| **Handling of Noise** | **Deterministic:** Attributes all deviations from the frontier to inefficiency. | **Stochastic:** Separates random noise (*v_i*) from actual technical inefficiency (*u_i*). |
| **Structural Insulation from Sample Size** | **High:** The boundary remains geometrically well-defined even with small group-specific sample sizes. | **Low:** Highly sensitive to sample size; requires sufficient *n* to identify the variance parameters. |
| **Efficiency Scores** | Often yields higher, more optimistic scores unless **bootstrap bias correction** is applied (Simar & Wilson, 2007). | Generally reports lower mean efficiency by filtering out statistical luck/noise. |

> **The "So What?" for the Learner**
>
> DEA is a deterministic envelopment — it creates a local reference
> point regardless of sample density. Conversely, SFA relies on
> numerical identification on a likelihood surface. If the sample size
> is insufficient, the SFA engine cannot reliably separate luck from
> skill, leading to **boundary-degenerate** results. This becomes
> critically problematic in the presence of **spatial-ecological
> heterogeneity**, where different environments (such as varying
> altitudes) necessitate cluster-specific analysis.

---

## 2. The K2 Cluster: A Case Study in Estimation Fragility

The structural limits of SFA are vividly illustrated by a study of 380
Nile tilapia (*Oreochromis niloticus*) farms in Jambi Province,
Indonesia. Farms were stratified into three elevation clusters:
**K1** (0–100 m a.s.l., *n* = 193), **K2** (>100–500 m a.s.l.,
*n* = 51), and **K3** (>500 m a.s.l., *n* = 136).

Although the total sample was robust, the **K2 cluster (*n* = 51)**
served as the epicenter of model failure.

### The Boundary Problem and the Variance-Ratio Parameter

In SFA, the variance-ratio parameter γ is defined as:

```
γ = σ²_u / (σ²_u + σ²_v)
```

This parameter determines how much of the composite error is attributed
to inefficiency (σ²_u) versus noise (σ²_v). A **Boundary Problem**
occurs when γ hits the limit of its [0, 1] parameter space.

**Evidence of Fragility** (*Manuscript Table 2 — Bootstrap Stability,
B = 500 replications*):

| Cluster | *n* | % Bootstrap Replications at Boundary | Assessment |
|---|---|---|---|
| **K2** | **51** | **80.8%** | **Fragile** |
| K1 | 193 | 26.8% | More stable |
| K3 | 136 | 15.6% | Most stable |

> **Caution on K1:** While K1 and K3 are clearly more stable than K2,
> K1's 26.8% boundary rate is not negligible in absolute terms — it
> simply reflects the same sample-size mechanism operating at a less
> severe level.

### The Flat Likelihood Mechanism

Small sample sizes result in a **flat likelihood surface**. In this
state, the optimizer cannot precisely identify the noise-inefficiency
split. Because the Battese & Coelli (1988) estimator is a non-linear
function of γ, arbitrary choices on a flat surface lead to large,
unreliable swings in farm-level efficiency scores. This failure is not
a localized fluke — it is a predictable statistical breakdown caused by
applying a parametric tool on a sample size insufficient to support its
complexity.

---

## 3. The Monte Carlo Simulation: Proof of a General Rule

*(Manuscript Sections 3.10 and 4.6; script: `J2_08_MonteCarlo_SmallSample.R`)*

To confirm that the K2 failure reflects a structural property of the
**estimator** rather than a Jambi-specific artifact, the authors
performed a Monte Carlo simulation using a **known, fixed
data-generating process** (true γ = 0.80; Cobb-Douglas, 4 inputs).
Sample size *n* was the only variable.

**Results across 1,000 replications per sample size** (*Manuscript Table 14*):

| *n* | Non-Convergence Rate | % Boundary-Degenerate γ | ρ (Estimated vs. True TE) |
|---|---|---|---|
| 50 (~K2) | **25.6%** | **5.2%** | 0.688 |
| 100 | 7.5% | 1.1% | 0.715 |
| 200 (~K1/K3) | 4.8% | 0.0% | 0.721 |

All three diagnostics move **monotonically** in the direction predicted
by the theory — confirming that sample size alone is sufficient to
degrade SFA reliability.

> **Nuance:** The gap in rank-correlation accuracy across sample sizes
> in the simulation (0.688 → 0.721) is considerably narrower than the
> gap observed in the real data between K2 and K1/K3. This suggests
> that sample size is a genuine but **partial** explanation of K2's
> empirical fragility — ecological or measurement features specific to
> that cluster compound the effect.

---

## 4. The Anatomy of Disagreement: Predicting Model Divergence

*(Manuscript Section 4.5; script: `J2_04_Divergence_Predictors_RQ4.R`)*

"Absolute Divergence" |TE_DEA − TE_SFA| is **not random**. An OLS
regression on the 380 farm-level pairs explains 45.5% of the variance
(cross-validated R² = 0.43, confirming genuine out-of-sample predictive
structure; script: `J2_07_CrossValidation_Divergence_Model.R`):

| Predictor | Effect on Divergence | Hypothesized Mechanism |
|---|---|---|
| **Cluster Sample Size (*n*)** | ↓ Decreases | H1: Larger samples provide the numerical identification SFA needs to converge with the stable DEA envelope. |
| **Water Quality (AWQI)** | ↑ Increases | H2: High water quality raises the ecological ceiling. Flexible DEA envelops this ceiling locally; restrictive SFA struggles to approximate it via a global functional form — so the two estimators disagree most where realized output diverges most from ecological potential. |
| **Floating Net Cage System (WAD = 1)** | ↓ Decreases | H3: Standardized technology is easier for a parametric model to approximate than heterogeneous earthen ponds. |
| **Feed-Conversion Ratio (FCR)** | ↑ Increases | High feed intensity introduces technical complexity that the two models interpret differently. |

### The FCR Paradox: Agreement in Disagreement

Despite their structural disagreements on rankings, both paradigms
**converge robustly** on the feed-conversion-ratio paradox: high FCR
is identified as the primary driver of technical inefficiency
(*p* < 0.001 in both DEA two-stage and SFA-BC95 bootstrap models;
3 of 4 determinants show directional concordance). This demonstrates
that even when tools disagree on efficiency *rankings*, they can agree
on major management *failures*.

---

## 5. Synthesis: Choosing the Right Tool for the Data

*(Manuscript Section 6; script: `J2_06_Paired_Bootstrap_Compare_Rho.R`)*

A common mistake is assuming that **architectural symmetry** —
making SFA look like DEA by adding cluster-specific slope interactions
— guarantees reliability. In reality, this "starves" the SFA model of
the data required to identify its parameters.

### The Learner's Decision Matrix

| Condition | Recommendation | Why |
|---|---|---|
| **Small or unbalanced samples** (any cluster approaching *n* ≈ 50) | Prefer DEA with bootstrap bias correction as the primary estimator | DEA's envelopment logic is structurally insulated from group size; SFA's MLE is not |
| **Must use SFA with small groups** | Prefer **Model B** (intercept shifters only, not full slope interactions) | Pools all 380 farms to estimate variance parameters, providing stability that complex interaction models sacrifice |
| **Checking for parameter failure** | Always inspect γ after fitting | If γ = 1.0 or 0.0, the model has failed to identify the noise-inefficiency split — resulting scores are mathematical artifacts, not measures of managerial skill |
| **Comparing DEA and SFA results** | Use a **paired bootstrap** to test whether rank-correlation differences are statistically meaningful — NOT a standard Fisher r-to-z test | Fisher r-to-z assumes independent samples; DEA and SFA scores computed on the same *n* farms are paired by construction |

> **Note on the "*n* < 100" threshold:** Some practitioners use
> *n* < 100 as a rule of thumb for "small sample" in SFA. While
> consistent with the simulation results here (which show clear
> degradation at *n* = 50 and meaningful improvement by *n* = 200),
> this specific threshold was not formally tested in this study.
> Treat it as a heuristic, not a precise cutoff.

---

## Final Summary

**Architectural symmetry between SFA and DEA does not guarantee
reliability.** SFA is a parametric tool constrained by the structural
limits of sample size and the necessity of numerical identification on
a likelihood surface. To produce valid research, the choice of tool
must respect the data's structure:

> *Use complexity only where the sample size provides the "fuel"
> to sustain it.*

---

## References

- Aigner, D., Lovell, C. A. K., & Schmidt, P. (1977). Formulation and
  estimation of stochastic frontier production function models.
  *Journal of Econometrics*, 6(1), 21–37.
- Battese, G. E., & Coelli, T. J. (1988). Prediction of firm-level
  technical efficiencies with a generalized frontier production
  function and panel data. *Journal of Econometrics*, 38(3), 387–399.
- Simar, L., & Wilson, P. W. (2007). Estimation and inference in
  two-stage, semi-parametric models of production processes.
  *Journal of Econometrics*, 136(1), 31–64.
- Wang, H.-J., & Schmidt, P. (2002). One-step and two-step estimation
  of the effects of exogenous variables on technical efficiency levels.
  *Journal of Productivity Analysis*, 18(2), 129–144.

---

*For questions about the methodology, see the manuscript or contact
the corresponding author via the details on the manuscript title page.*
