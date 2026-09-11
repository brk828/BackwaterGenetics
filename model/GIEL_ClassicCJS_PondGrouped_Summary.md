# Bonytail (GIEL) Classic Annualized CJS Model — Pond-Grouped Extension

**Purpose:** extends `model/GIEL_ClassicCJS_Summary.md` (the pooled, no-pond-effect classic CJS
model) by adding IPCA pond identity (IP2/IP5/IP6) as a grouping factor, per user request
2026-09-11: independent Phi and p by pond, built from the previous top-2 pooled models, AIC
comparison across a pond-grouped model set, and a parametric bootstrap goodness-of-fit (GOF)
on the top model itself (not just the null/pooled model). Code: `model/GIEL_ClassicCJS_PondGrouped.R`;
output: `data/GIEL_ClassicCJS_PondGrouped.RData`, `model/GIEL_ClassicCJS_PondGrouped_bootGOF.png`.

## 1. Pond-exclusion check

The 300 GIEL fish dropped for missing length-at-tagging (see `GIEL_ClassicCJS_Summary.md` §6)
are **all stocked fish at IP2**. Checked whether this makes IP2's pond-specific parameters
inestimable: after dropping them, IP2 still has 1,156 tagged fish with a fully populated
location × origin × size cross-tab (137 wild-adult, 720 wild-young/netting-tagged, 291
stocked-adult, 8 stocked-young — no empty cells at any of the three ponds). **No pond was
excluded** — IP2, IP5, and IP6 are all retained with pond-specific Phi and p.

## 2. Model set and AIC comparison

Twelve models were fit: the previous top-2 pooled models as a baseline, plus ten pond-grouped
variants ranging from fully independent Phi/p by pond down to an additive pond effect only.
A first pass (8 models) showed 4 candidates within 20 AIC of the leader; a second pass (4 more
models varying which `location:*` interaction terms were dropped from Phi) was run to resolve
the ambiguity, per the standard practice of exploring further when AIC results aren't clear.

| Model | Phi | p | k | AIC | ΔAIC | weight |
|---|---|---|---|---|---|---|
| **M11 (top)** | `~location*time + location*sizeclass + origin` | `~time` | 40 | 4958.8 | 0.00 | 0.604 |
| M9 | `~location*time + sizeclass + origin` | `~time` | 38 | 4961.8 | 3.00 | 0.135 |
| M8 (fully independent Phi) | `~location*time + location*sizeclass + location*origin` | `~time` | 42 | 4961.8 | 3.04 | 0.132 |
| M10 (fully independent Phi & p, additive) | `~location*time + location*sizeclass + location*origin` | `~location + time` | 44 | 4963.1 | 4.25 | 0.072 |
| M12 | `~location*time + sizeclass + location*origin` | `~time` | 40 | 4963.8 | 5.03 | 0.049 |
| M4 (fully independent Phi, pond-constant p) | `~location*time + location*sizeclass + location*origin` | `~location` | 36 | 4967.4 | 8.62 | 0.008 |
| M3 (fully independent Phi & p) | same as M8 | `~location*time` | 60 | 4980.1 | 21.25 | <0.001 |
| M5 | `~location*time + sizeclass + origin` | `~location*time` | 56 | 4981.8 | 22.94 | <0.001 |
| M7 (pooled Phi, pond-specific p) | `~time + sizeclass + origin` | `~location*time` | 38 | 5351.3 | 392.7 | ~0 |
| M6 (additive pond effect only) | `~location + time + sizeclass + origin` | `~location` | 16 | 5461.0 | 502.2 | ~0 |
| M1 (previous pooled top model) | `~time + sizeclass + origin` | `~time` | 20 | 5876.8 | 918.2 | ~0 |
| M2 (previous pooled 2nd model) | `~time + sizeclass + origin` | `~1` | 12 | 5880.8 | 922.2 | ~0 |

**Adding pond identity to Phi is an enormous improvement** (ΔAIC ≈ 900+ vs. the pooled models) —
by far the largest single structural change tested in either the pooled or pond-grouped model
sets. **Detection (p) does not need a pond effect**: every model with `p(~location*time)` or
`p(~location)` fits worse (by AIC) than the matched model with pooled `p(~time)`, so a single,
pond-shared detection-by-year curve is adequate. **Within Phi, the full three-way pond
interaction (pond × time, pond × sizeclass, pond × origin) is not needed either** — the top
model (M11) drops the `location:origin` interaction (origin's stocked-vs-wild survival gap is
shared across ponds) while keeping `location:time` and `location:sizeclass`; M9 (drops
`location:sizeclass` too) and M8 (keeps all three) are statistically indistinguishable from
M11 (ΔAIC ≈ 3, i.e., within the range where model averaging would ordinarily be considered).
**The requested "independent calculations for all phi and p" per pond (M8) is essentially tied
with the leaner M11** — the AIC evidence says pond-specific *time trajectories and size effects*
in survival are clearly needed, but the origin effect need not be estimated separately per
pond. Given the fine margin between M8/M9/M10/M11/M12, M11 (strict AIC winner) is used as *the*
top model below, but this is not a decisively unique choice among these five.

## 3. A caution on the top model's realized estimates: die-off separation

Fitting `location*time` interactions lets the model directly "discover" both documented GIEL
die-off events (see `GIEL_ModelResults_Summary.md` §3) without any explicit event term — but at
the cost of **quasi-complete separation** in exactly those cells, exactly as would be expected
from an unregularized MLE:

- `Phi[IP5, 2020, wild]` (the IP5 Jan–Sep 2020 crash) is estimated at ~3×10⁻⁸–3×10⁻⁹ with SEs
  in the hundreds of logit units.
- `Phi[IP2, 2022, wild]` (the IP2 Jan–Jul 2022 crash) is similarly estimated near 0 (~10⁻⁶–10⁻⁷)
  with comparably enormous SEs.
- `p[2022]` (pooled across ponds) is pinned at exactly 1.0 with a degenerate CI — the same
  "boundary estimate" phenomenon already documented for the pooled model
  (`GIEL_ClassicCJS_Summary.md` §4) and for the Chapman estimator elsewhere in this project.

**These are not convergence failures** (the optimizer reports `convergence = 0` and the
log-likelihood surface is smooth near the estimate) — they are genuine boundary/separation
estimates driven by near-total mortality in a single pond-year, the same phenomenon that
motivated adding an explicit, prior-regularized `dE` die-off offset in the Bayesian
robust-design model (`GIEL_ModelMethodology.md` §3) rather than letting a free time×pond
survival parameter absorb it. The point estimates and any derived annual-survival figure for
IP5-2020 and IP2-2022 specifically should be treated as "collapsed to near zero, magnitude not
meaningfully estimable from this design" rather than as precise numbers. All other pond-years
have plausible, reasonably-precise estimates (illustrative IP2/wild examples: 2018 adult 0.24,
2019 adult 0.27, 2021 adult 0.14, 2023 adult 0.26, 2024 adult 0.08, 2025 adult 0.07 — note IP2's
generally low, and declining, non-crash-year survival, consistent with the "chronically low
IP2 survival" flag already raised in `GIEL_ModelResults_Summary.md` §4).

## 4. Parametric bootstrap GOF on the top model (M11)

Rather than relying only on the pooled/null-model R2ucare test (`GIEL_ClassicCJS_Summary.md`
§5, which tests the fully time-dependent CJS model with no covariates), a proper ĉ was computed
**for M11 itself**:

1. **Simulate.** For every one of the 2,637 actual tagged fish (its real location, origin,
   entry occasion, and initial size all held fixed), draw a stochastic encounter history
   forward from its own entry occasion using M11's fitted real Phi/p estimates as the true
   generating probabilities (survive → detect, Markov chain, one draw per remaining interval).
2. **Refit.** Re-fit the identical M11 model structure to each simulated dataset (warm-started
   from the observed fit's coefficients for speed; no Hessian needed since only the
   log-likelihood is used).
3. **Deviance.** For both the observed data and each simulated dataset, compute deviance
   relative to a *cohort-saturated* reference model — cohort = (pond, origin, initial size
   class, entry occasion), i.e., every piece of covariate information available to
   differentiate individuals at release, with cohort-specific cell probabilities set equal to
   observed pattern frequencies (the standard saturated-model definition of deviance,
   analogous to Program MARK's bootstrap GOF procedure).
4. **c-hat** = observed deviance / mean(bootstrap deviance).

**Results (300 replicates, all converged):**

| Quantity | Value |
|---|---|
| Observed deviance | 189.0 |
| Bootstrap deviance: mean / median / SD | 138.2 / 138.6 / 17.8 |
| Bootstrap deviance range | 92.3 – 193.4 |
| **ĉ (mean-based)** | **1.37** |
| ĉ (median-based) | 1.36 |
| Bootstrap p-value, P(deviance_boot ≥ deviance_obs) | 0.007 |

See `model/GIEL_ClassicCJS_PondGrouped_bootGOF.png` for the bootstrap distribution with the
observed deviance marked.

**Interpretation.** ĉ ≈ 1.37 indicates **mild overdispersion** even in this fully pond-,
time-, size-, and origin-structured model — the observed deviance sits at roughly the 99th
percentile of what the fitted model itself would produce by chance (p = 0.007), so some
extra-binomial variation remains unexplained. This is a modest amount by CJS standards (ĉ < 2
is generally considered a minor correction, typically handled via QAICc rather than indicating
a fundamentally misspecified model) and is much smaller than the crude null-model ĉ ≈ 8.4
reported for the pooled model's overall RELEASE-style GOF test (`GIEL_ClassicCJS_Summary.md`
§5) — consistent with pond, size, and origin structure absorbing most, but not quite all, of
the heterogeneity in this population. Likely residual sources (not further decomposed here):
the boundary/separation cells noted in §3 contribute some of this excess deviance since a
saturated-cohort reference does not "know" that those cells are structurally near-zero: and/or
mild individual heterogeneity in survival/detection beyond what pond×time×size×origin capture.
**Recommendation: use ĉ ≈ 1.37 (not 1.0) to inflate standard errors / compute QAICc if this
pond-grouped model's parameter estimates are used for inference**, and treat the IP5-2020 and
IP2-2022 point estimates per the §3 caution regardless of ĉ adjustment.

## 5. Caveats (carried over / new)

- Everything in `GIEL_ClassicCJS_Summary.md` §6 still applies (no untagged/recruit pool, no
  netting join, no explicit die-off term, size class purely time-since-tagging based, 300
  excluded fish).
- **The near-tied AIC cluster (M8/M9/M10/M11/M12, §2) means the choice of "which pond
  interactions Phi needs" is somewhat data-thin** — a modeler preferring the literal
  "independent phi and p by pond" specification requested (M8) rather than the strict AIC
  winner (M11) would lose essentially nothing (ΔAIC = 3.04, weight 13% vs. 60%). The two
  practical conclusions (pond effects needed; origin effect can be pooled across ponds) are
  robust to this choice.
- **The bootstrap GOF assumes the top model's fitted real parameters are the true data-generating
  process** — this is a standard limitation of parametric bootstrap GOF (it tests "if this model
  is exactly right, is the observed level of extra variability plausible," not "is this model
  the right structure").

## 6. Key files

| File | Contents |
|---|---|
| `model/GIEL_ClassicCJS_PondGrouped.R` | Build + fit script (data prep with pond as a `marked` group, 12-model AIC comparison, parametric bootstrap GOF) |
| `data/GIEL_ClassicCJS_PondGrouped.RData` | `aic_final` (12-model AIC table), `fits` (all 12 fitted `crm` objects), `top_name`/`top_ponded` (AIC-best model, M11), `preds` (real-scale Phi/p estimates), `deviance_obs`/`dev_boot`/`c_hat_mean`/`c_hat_median`/`p_val` (bootstrap GOF results), `giel_data`/`dp2`/`ddl2`/`sim_base` (processed data, design data, and the individual-level template used for simulation) |
| `model/GIEL_ClassicCJS_PondGrouped_bootGOF.png` | Bootstrap deviance histogram with the observed deviance marked |
