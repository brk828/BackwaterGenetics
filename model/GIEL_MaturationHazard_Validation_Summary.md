# GIEL Maturation-Hazard Sensitivity — Validation Checklist (F1, F2, A1)

**Date:** 2026-09-15
**Script:** `model/GIEL_MaturationHazard_Validation.R` (output: `data/GIEL_MaturationHazard_Validation.RData`)

Validates the three size-tied maturation-hazard sensitivity fits of the joint
2-stage robust-design model, GIEL-only (`RUN_PONDS = c(3, 6, 7)`; pond_new
1 = IP2, 2 = IP5, 3 = IP6), against the on-completion checklist in
`AGENTS.md`'s "GIEL Size-Tied Maturation Hazard" section:

| Track | Growth-model source | File | Run time |
|---|---|---|---|
| F1 | Fabens increment, IPCA-only | `data/BWRobustDesign_MCMC_giel_F1.RData` | 69.1 min (completed 2026-09-14 16:26) |
| F2 | Fabens increment, IPCA + Cibola | `data/BWRobustDesign_MCMC_giel_F2.RData` | 100.0 min (completed 2026-09-14 18:17) |
| A1 | Age-based VBGF, IPCA-only | `data/BWRobustDesign_MCMC_giel_A1.RData` | 100.3 min (completed 2026-09-14 18:19) |
| baseline | scalar `lpsi[GIEL]`, no hazard table | `data/BWRobustDesign_MCMC_giel.RData` | (2026-09-10 fit, pre-dates this sensitivity work) |

All three sensitivity fits use identical settings (3 chains, 30k iter / 10k
burn-in / thin 10, seeds 4127/8351/2960) and differ only in which
`data/GIEL_MaturationHazard_<track>.RData` lookup table was loaded during the
data build feeding each fit.

**Correction to prior belief:** F1 had already been run and was mistakenly
believed not to have been (see conversation of 2026-09-15). All three tracks
are validated together below rather than launching a new F1 run.

## 1. Known-alive constraint (`N_tag >= known alive`)

For every pond x month with at least one fish known alive (`entry_t <= t <=
end_t`, `last_alive_t >= t`, `count_from <= t`), the posterior 2.5% quantile
of `sum(N_tag[pond, t, 1:2])` (tagged juvenile + adult) must be at or above
the directly-counted known-alive total.

| Track | Pond-months checked | Violations |
|---|---|---|
| F1 | 183 | **0** |
| F2 | 183 | **0** |
| A1 | 183 | **0** |

All three tracks pass with no violations.

## 2. Convergence (Rhat)

Excluding the known chronically-poor-mixing nuisance families (`N_tag`,
`Tal`, `U`, `U0`, `RecJ`, `r`, `Atot` — recruitment/pool quantities that mix
poorly in every prior GIEL/XYTE fit in this project), the 58 core
hyperparameters (survival, detection, maturation, offsets, hyper-priors)
converge acceptably in all three tracks:

| Track | Core: max Rhat | Core: n > 1.1 | Nuisance: max Rhat | Nuisance: n > 1.1 (of 616) |
|---|---|---|---|---|
| F1 | 1.046 | 0 | 1.294 | 64 |
| F2 | 1.027 | 0 | 1.289 | 41 |
| A1 | 1.053 | 0 | 1.291 | 32 |

Worst core-parameter Rhats are `lphiA[·,·]` (season-adult-survival) and
`dE[·]` cells in the low single digits above 1.0 (max 1.053, A1) — consistent
with, and no worse than, prior GIEL fits in this project. The elevated
nuisance-family Rhats are the same floating-point/mixing artifact documented
throughout `AGENTS.md` for `N_tag`/`Tal` (pinned by the known-alive
constraint) and for recruitment nuisance parameters, not new problems
introduced by the hazard mechanism.

## 3. `lpsi_adj[GIEL]` plausibility

`lpsi_adj[1]` is the logit-scale correction the model applies on top of each
track's raw growth-model-derived hazard lookup table. A value near 0 means
the robust-design contact/stage data agree with that track's schedule as-is.

| Track | mean | 95% CI |
|---|---|---|
| F1 | **+0.227** | (0.044, 0.402) |
| F2 | **−0.645** | (−0.852, −0.430) |
| A1 | **−0.643** | (−1.081, −0.118) |

All three are credibly nonzero (none of the CIs include 0), but the
*direction* differs by track and is informative: F1's positive offset means
the contact data want maturation **faster** than the raw Fabens-IPCA-only
schedule predicts (consistent with the 2026-09-14 finding in
`AGENTS.md`/`GIEL_VBGF_AgeLength_Summary.md` that Track F1 gives the
*lowest* P(TL >= 250mm) curve at low TL1 of the three tracks — the model
compensates by speeding up maturation). F2 and A1's negative offsets, of
similar magnitude to each other, mean their schedules predict maturation
**too fast** and the data want a slower correction. This is a genuine,
non-trivial difference between tracks, not just noise — see the plan file's
"Open risks" bullet about the smooth-curve divergence between the age-based
and increment growth models, which shows up here as opposite-signed
corrections.

## 4. Core hyperparameter comparison (F1 / F2 / A1 / baseline)

| Parameter | F1 | F2 | A1 | baseline |
|---|---|---|---|---|
| `mu_phi[1]` | 2.71 (2.40, 3.04) | 2.72 (2.40, 3.05) | 2.73 (2.42, 3.05) | 2.71 (2.38, 3.03) |
| `sigma_phi[1]` | 0.848 | 0.848 | 0.847 | 0.849 |
| `dJ[1]` | 0.413 (0.233, 0.592) | 0.361 (0.196, 0.524) | 0.313 (0.151, 0.472) | 0.453 (0.284, 0.634) |
| `dE[1]` (IP2) | −1.73 | −1.78 | −1.74 | −1.74 |
| `dE[2]` (IP5) | −1.42 | −1.46 | −1.54 | −1.45 |
| `dE[3]` (IP6, prior-only) | −2.37 | −2.36 | −2.39 | −2.42 |
| `lpsi[1]` (untagged-pool annual maturation) | −2.23 | −2.16 | −2.03 | **−1.70** |
| `lp[1,1]` (IP2 juvenile detection) | −2.03 | −1.75 | −1.46 | −2.20 |
| `lp[2,1]` (IP5 juvenile detection) | 0.158 | −0.038 | −0.369 | 0.009 |
| `g_eff[1]` | 0.916 | 0.909 | 0.902 | 0.912 |

Survival hyperparameters (`mu_phi`, `sigma_phi`, `dE[1:3]`) are essentially
identical across all three hazard tracks and the baseline (differences <
0.1 logit) — the size-tied maturation submodel remains largely separable
from survival/die-off estimation, matching the pattern already documented
for the XYTE `a1` size-split submodel and for the original F1-only
validation. The real movement is concentrated in `lpsi` (untagged-pool
maturation, which loses information as the tagged-fish pathway becomes
size-tied) and juvenile detection (`lp[·,1]`), both of which shift modestly
and monotonically from F1 -> F2 -> A1 as the underlying schedule's implied
maturation speed changes. `dJ[1]` also drifts down (0.413 -> 0.361 -> 0.313)
as more juvenile time gets reallocated by the faster A1/F2 schedules — still
positive in all three (juvenile survival modeled above adult on the logit
scale), consistent with the pre-existing, separately-investigated IP2
detection-asymmetry finding, not a new artifact of this sensitivity sweep.

## 5. FY2022 (IP2) / FY2020 (IP5) event-adjusted realized annual survival

Using each track's own `lphiA[pond,season]` and `dE[pond]` posteriors, with
the confirmed die-off windows (IP2: 7 of 12 FY2022 months, `nEvt[IP2,FY2022]
= 7`; IP5: 9 of 12 FY2020 months, `nEvt[IP5,FY2020] = 9`, both from
`data/BWRobustDesign_data.RData`):

| Track | Pond | Season | Event-adjusted realized survival (mean, 95% CI) | No-event-adjustment comparison |
|---|---|---|---|---|
| F1 | IP2 | FY2022 | 0.090 (0.045, 0.147) | 0.474 |
| F2 | IP2 | FY2022 | 0.094 (0.045, 0.156) | 0.489 |
| A1 | IP2 | FY2022 | 0.097 (0.048, 0.157) | 0.488 |
| F1 | IP5 | FY2020 | 0.046 (0.020, 0.085) | 0.353 |
| F2 | IP5 | FY2020 | 0.044 (0.019, 0.080) | 0.363 |
| F2 | IP5 | FY2020 | 0.043 (0.019, 0.081) | 0.383 |

All three tracks give essentially the same event-adjusted realized survival
for both die-off seasons (differences <= 0.01), confirming the maturation-
hazard track choice does not materially affect the die-off (`dE`) estimates
or their downstream realized-survival implications — consistent with §4's
finding that the tracks differ mainly in maturation-timing quantities, not
survival.

## Overall conclusion

All three tracks (F1, F2, A1) pass the known-alive constraint with zero
violations and converge acceptably (max core Rhat 1.05). Survival and
die-off estimates are robust to the choice of maturation-hazard track. The
substantive difference between tracks is the sign and magnitude of
`lpsi_adj[GIEL]` (F1 wants faster maturation than its schedule predicts; F2
and A1 want slower), reflecting the same underlying divergence between the
Fabens-increment and age-based VBGF growth-model curves already flagged as
an open reconciliation item in the plan file and in `AGENTS.md`. No track is
disqualified by this validation pass; choosing among them (or reporting all
three as a sensitivity range) remains a modeling decision for BK.
