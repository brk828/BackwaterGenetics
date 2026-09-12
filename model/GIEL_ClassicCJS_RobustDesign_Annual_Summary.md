# Bonytail (GIEL) Classic Annualized CJS Model — IP5-Only, Robust-Design-Motivated

**Purpose:** an IP5-only classic (frequentist, MLE) CJS cross-check on the joint hierarchical
Bayesian robust-design model (`data/BWRobustDesign_MCMC_giel.RData`), built specifically to
target the FY2020 IP5 die-off and structured, where possible, like the Bayesian model's own
survival/gap logic. Built and fit 2026-09-11 with the `marked` R package. Code:
`model/GIEL_ClassicCJS_RobustDesign_Annual.R`; output:
`data/GIEL_ClassicCJS_RobustDesign_Annual.RData`,
`model/GIEL_ClassicCJS_RobustDesign_Annual_bootGOF.png`.

## 0. Why this script exists (collapsed from a failed monthly attempt)

An earlier script, `model/GIEL_ClassicCJS_RobustDesign.R` (not separately documented), tried to
reproduce the Bayesian joint model's *exact* monthly robust-design time structure for IP5 stocked
fish: 70 monthly occasions, 60 of 69 within-season transitions fixed to Phi = 1 via `marked`'s
`fix = 1` mechanism, with all real survival information concentrated into 9 annual Apr→Oct gaps.
`crm()` reproducibly hard-crashed (non-catchable, even under `callr` process isolation) on a large
fraction of bootstrap-simulated capture histories under that design — traced to the within-season
closure assumption producing degenerate/boundary outcomes in the small (~18–20 fish) "young"
stratum. **Fix:** collapse to 10 ANNUAL occasions (FY2017–FY2026, the 9 between-year intervals
only) with no within-season transitions and no `fix` mechanism at all — every one of the 9
transitions is a real, freely-estimated annual survival parameter. Each bootstrap replicate's
refit is additionally isolated in its own `callr` subprocess (protects the job even if `crm()`
crashes again).

## 1. Design

- **Population:** IP5 **stocked-only** GIEL, 600 fish, two cohorts — FY2017 (280 tagged as
  "adult"-size, 20 as "young") and FY2022 (300, all adult-size; `pond_check` in the saved `.RData`
  shows all 600 have a recorded tagging length, unlike IP2's 300 length-less stocked fish). No
  netting-tagged ("wild") fish and no other pond are included — this is a narrower population than
  `GIEL_ClassicCJS.R`/`GIEL_ClassicCJS_PondGrouped.R`, chosen to focus specifically on IP5's FY2020
  die-off.
- **Occasions:** annual, FY2017–FY2026 (10 occasions, 9 survival intervals). Detection at an
  occasion = any PIT contact anywhere in that FY's Oct–Apr window (May–Sep contacts ignored, same
  convention as the pooled/pond-grouped classic CJS models). 0 contact-occasions were suppressed
  as pre-entry (`n_suppressed = 0`).
- **Size class ("young"/"adult"):** identical logic to the other classic CJS scripts — a fish
  tagged <250 mm TL is "young" for exactly the one annual interval immediately following its own
  tagging occasion, "adult" for every interval after that (and immediately, if tagged ≥250 mm).
  Only 3 of 126 Phi design rows are "young" (the FY2017 cohort's 20 young-tagged fish collapse to 3
  unique covariate combinations after `marked` deduplicates identical histories).
- **Phi covariates tested (in place of the pooled models' `time`/`origin`, since origin is
  constant here — all fish are stocked):**
  - `gap` — a **9-level year factor** (FY2017…FY2025, one level per interval-start year); despite
    the name, this is full annual time-variation, not a two-level dummy.
  - `dieoff` — a 2-level factor (`event` = the FY2020 interval only, `normal` = every other
    interval), a more parsimonious representation of "one die-off year, otherwise pooled."
  - `const` — sizeclass only, no time variation.
- **p covariates tested:** `time` (one estimate per year) or constant (`.`).
- **Model set:** 3 Phi structures × 2 p structures = 6 models (a much smaller set than the
  pooled/pond-grouped scripts since origin is dropped and pond identity is fixed at IP5).

## 2. Model selection

| Model | k | −2lnL | AIC | ΔAIC | weight |
|---|---|---|---|---|---|
| **gap_pdot (top)** | 11 | 1304.63 | 1326.63 | 0.00 | 0.997 |
| gap_ptime | 19 | 1300.39 | 1338.39 | 11.76 | 0.003 |
| const_ptime | 11 | 1592.01 | 1614.01 | 287.4 | ~0 |
| dieoff_pdot | 4 | 1890.48 | 1898.48 | 571.9 | ~0 |
| dieoff_ptime | 12 | 1886.03 | 1910.03 | 583.4 | ~0 |
| const_pdot | 3 | 1916.42 | 1922.42 | 595.8 | ~0 |

(full table in `data/GIEL_ClassicCJS_RobustDesign_Annual.RData$aic_ip5`). **`gap` (full annual
time variation) crushes `dieoff` (single event dummy) and `const` by 500+ AIC units** — a single
"crash year" indicator is far too coarse for IP5; survival genuinely varies year to year beyond
just the FY2020 collapse. **Detection does not need to vary by year**: `gap_pdot` beats
`gap_ptime` by 11.8 AIC (weight 0.997 vs 0.003) even though `gap_ptime` has 8 more free detection
parameters — a single shared detection probability across all years fits this smaller, IP5-only
dataset better once the extra parameters are penalized. Results below are from `gap_pdot`.

## 3. Survival (Phi) estimates, top model (`gap_pdot`)

| Interval (FY → FY+1) | Sizeclass | Survival | 95% CI |
|---|---|---|---|
| 2017 → 2018 | adult | 0.897 | 0.852–0.929 |
| 2017 → 2018 | young | 0.958 | 0.672–0.996 |
| 2018 → 2019 | adult | 0.056 | 0.034–0.090 |
| 2019 → 2020 | adult | 0.941 | 0.610–0.994 |
| **2020 → 2021** | adult | **8.3×10⁻⁶** | ~0–1 (degenerate) |
| 2021 → 2022 | adult | 0.888 | fixed, SE = 0 |
| 2022 → 2023 | adult | 0.846 | 0.800–0.883 |
| 2023 → 2024 | adult | 0.657 | 0.596–0.713 |
| 2024 → 2025 | adult | 0.579 | 0.502–0.652 |
| 2025 → 2026 | adult | 0.460 | 0.363–0.561 |

**The FY2020 interval (the documented Jan–Sep 2020 die-off window) collapses to an essentially
zero, boundary/degenerate estimate** (8.3×10⁻⁶, SE 7.7×10⁻⁴, CI spanning the entire [0,1] range) —
the same quasi-complete-separation phenomenon already documented for the die-off cells in
`GIEL_ClassicCJS_PondGrouped_Summary.md` §3, arising naturally here because `gap` gives every year
its own free parameter. This matches the independent empirical floor exactly: **0 of the 300 IP5
FY2017-cohort fish known alive at the FY2020 occasion were ever detected again at FY2021 or later**
(`floor_2020`, §7 below).

**The following FY2021 → FY2022 estimate (0.888) has SE = 0 and should not be trusted as a real,
data-informed number.** Because essentially none of the original FY2017 cohort survived past the
FY2020 crash, there is almost no fish actually at risk during the FY2021 interval to inform this
parameter — the optimizer reports a value with zero estimated uncertainty, a classic sign of a
parameter pinned at (or never moved from) its starting value rather than a genuinely
likelihood-informed estimate. Treat FY2021 → FY2022 survival as **not identifiable from this
IP5-only, 2-cohort design**, not as "recovery to 89% survival."

All other intervals (FY2017–FY2019 and FY2022–FY2025) have plausible point estimates with
finite, reasonable SEs.

## 4. Detection (p) estimate, top model

| p (constant across years) | 95% CI |
|---|---|
| 0.991 | 0.973–0.997 |

Detection is very high and estimated with a single parameter (`p ~ 1` beat `p ~ time` on AIC,
§2) — consistent with the high detection rates found in every other GIEL analysis in this
project (pooled and pond-grouped classic CJS, Bayesian joint model).

## 5. Parametric bootstrap GOF on the top model

Same recipe as `GIEL_ClassicCJS_PondGrouped_Summary.md` §4 (simulate each of the 600 real fish's
history forward under the fitted model's real Phi/p, refit the identical model, compare deviance
to a cohort-saturated [entry-FY × initial sizeclass] reference), but with each bootstrap
replicate's refit isolated in its own `callr` subprocess (see §0).

| Quantity | Value |
|---|---|
| Observed deviance | 6.98 |
| Bootstrap deviance: mean / SD | 7.65 / 3.37 |
| **ĉ (mean-based)** | **0.912** |
| ĉ (median-based) | 0.971 |
| Bootstrap p-value, P(deviance_boot ≥ deviance_obs) | 0.53 |
| Replicates crashed/failed | **0 of 300** |

**Clean fit, no overdispersion** (ĉ < 1, bootstrap p = 0.53 — the observed deviance sits right in
the middle of what the fitted model itself produces by chance). This is a much better-behaved GOF
result than either the pooled all-3-pond model (crude null-model ĉ ≈ 8.4) or the pond-grouped
model (ĉ ≈ 1.37) — expected, since restricting to one pond and letting survival vary completely
freely by year leaves very little unexplained heterogeneity for a saturated-cohort reference to
detect. See `model/GIEL_ClassicCJS_RobustDesign_Annual_bootGOF.png` for the bootstrap
distribution with the observed deviance marked. The **0-of-300 crash/failure rate** confirms the
annual-occasion redesign (§0) resolved the earlier monthly-design crashes.

## 6. Empirical floor check (FY2020 die-off)

Because occasions are annual, the "known alive going into the crash, ever seen again after" floor
is direct: of the 300 IP5 stocked fish known alive at the FY2020 occasion (entered on or before
FY2020), **0 were detected at the FY2021 occasion or later** (`floor_2020`:
`alive = 300, seen_after = 0`). This is the same total-loss signature documented from the raw
known-alive series in `AGENTS.md` (IP5 known-alive fell from a post-tagging peak of 145 in Dec
2019 to 0 by Sep 2020) and independently confirms the near-zero Phi estimate in §3.

## 7. Comparison against the Bayesian GIEL-only robust-design fit (IP5 = pond 2)

The script computes a quick, **simplified** annualized-survival companion from
`data/BWRobustDesign_MCMC_giel.RData` for reference: `plogis(lphiA[2, y] + 1{y=2020}·dE[2])^12`
(posterior mean), i.e. it applies the pond's `dE` die-off offset across **all 12 months** of the
FY2020 season as a simplification, rather than only the calendar months actually flagged by the
event window (`nEvt[IP5, FY2020] = 9` of 12 — see the more precise decomposition used in
`AGENTS.md`'s 2026-09-11 cross-check entry, which instead computes
`plogis(lphiA)^(12-9) · plogis(lphiA+dE)^9 ≈ 0.041` for FY2020, vs. this script's cruder
`plogis(lphiA+dE)^12 ≈ 0.021`; both are near-zero and consistent with the 0/300 floor in §6).

| FY (season) | Bayesian, naive 12-month dE (this script) |
|---|---|
| 2017 | 0.814 |
| 2018 | 0.154 |
| 2019 | 0.745 |
| 2020 | 0.021 (or 0.041, precise nEvt-weighted version — see above) |
| 2021 | 0.136 |
| 2022 | 0.783 |
| 2023 | 0.627 |
| 2024 | 0.542 |
| 2025 | 0.419 |
| 2026 | 0.463 |

**Important indexing caveat — this is not an apples-to-apples comparison with §3.** The Bayesian
column above is a **pond-level, all-fish (stocked + netting-tagged) seasonal survival rate**
derived from `lphiA[IP5, y]`, whereas §3's classic CJS Phi is fit to **stocked-only** fish across
just 2 cohorts and lets every year float completely freely.

### 7.1 Stocked-only Bayesian re-run (fairer comparison, 2026-09-11)

To test whether restricting the Bayesian comparison to stocked-only fish narrows the gap with the
classic CJS, the same posterior (`data/BWRobustDesign_MCMC_giel.RData`) was re-summarized using
only IP5's 600 stocked GIEL (`rd_Fish |> filter(pond_new==2, origin==1)`) — exactly the classic
CJS's 2 cohorts (300 tagged 2017-03-21, 300 tagged 2022-03-09), confirming the two models' fish
sets line up as expected. The **whole-season `b_post` approximation** is used (apply the
stocked-origin post-release offset `b_post[1,1]` to the entire stocking-FY season, mirroring the
existing whole-season `dE` convention already used for the die-off year); a calendar-precise
version that tracks each cohort's actual 12 months post-stocking was also tried but rejected as
unstable (it mixes two different seasons' `lphiA` depending on the exact tagging date, and gives
very different first-year survival for the two cohorts — 0.41 vs 0.71 — driven by an arbitrary
choice of calendar reference point that annual CJS occasions don't have).

| FY (season) | Bayesian, all IP5 GIEL (§7) | **Bayesian, stocked-only (n=600)** | Classic CJS Phi (§3, adult) |
|---|---|---|---|
| 2017 | 0.814 | **0.809 (0.748–0.863)** | 0.897 |
| 2018 | 0.154 | **0.154 (0.120–0.191)** | 0.056 |
| 2019 | 0.745 | **0.745 (0.617–0.856)** | 0.941 |
| 2020 | 0.021 (naive) / 0.041 (precise) | **0.041 (0.017–0.077)** | 8.3×10⁻⁶ (degenerate) |
| 2021 | 0.136 | **0.136 (0.005–0.528)** | 0.888 (SE=0, unidentifiable) |
| 2022 | 0.783 | **0.777 (0.713–0.832)** | 0.846 |
| 2023 | 0.627 | **0.627 (0.567–0.685)** | 0.657 |
| 2024 | 0.542 | **0.542 (0.467–0.616)** | 0.579 |
| 2025 | 0.419 | **0.419 (0.327–0.517)** | 0.460 |
| 2026 | 0.463 | **0.463 (0.255–0.691)** | — (no 10th interval in classic CJS) |

**Conclusion: restricting to stocked-only fish moves the Bayesian trajectory by ≤ 0.006 in every
season** (compare the "all IP5 GIEL" and "stocked-only" columns) — because `b_post[1,1]`
(stocked-origin post-release survival offset) has posterior mean only **−0.031 logit units**
(95% CI −0.110 to −0.001), a small and only marginally-distinguishable-from-zero effect. **This is
a structural, not a data-subsetting, limitation:** the joint Bayesian model has no mechanism for
stocked- vs. netting-tagged/established-origin fish to have different *long-run* annual survival
at IP5 — origin only enters as a small, temporary, imprecisely-estimated 6-month post-release
offset, so pooling in netting-tagged fish barely perturbs the trajectory in the first place.

**The year-by-year disagreement with the classic CJS therefore persists after this fix** — e.g.
FY2018 (classic 0.056 vs Bayesian stocked-only 0.154), FY2019 (0.941 vs 0.745), FY2021 (0.888,
itself unreliable per §3, vs 0.136) — and traces to **model structure** (the Bayesian model's
single shared `lphiA[IP5, y]` season effect vs. the classic CJS's fully free, unconstrained annual
Phi for a 2-cohort population), not to which fish are included in the Bayesian posterior summary.
FY2022–FY2025, where both cohorts are simultaneously at risk in the classic CJS and the Bayesian
season effect is presumably better informed by more contemporaneous data, show the closest
agreement (differences of 0.02–0.30, versus 0.09–0.75 in FY2018/FY2019/FY2021).

**Where the two models remain most directly comparable — the FY2020 crash itself and the 0/300
empirical floor — they still agree closely**: classic CJS ≈ 8×10⁻⁶ (degenerate low), Bayesian
stocked-only ≈ 0.041, both effectively "collapsed to near zero." This remains the extent of the
genuine cross-validation between the two models; outside the crash year, and now independent of
which fish subset is used, the two are not expected to match season-by-season.

### 7.2 Is the disagreement just the crash years, and is the classic CJS more accurate? (2026-09-11)

Two follow-up questions worth recording. **First: does real year-to-year variability persist even
excluding the two uninformative rows (FY2020's true near-zero die-off and FY2021's SE = 0 boundary
artifact)?** Yes. The remaining classic CJS estimates (§3) are all well-identified — no boundary
values, all with finite, reasonably tight SEs — and still swing substantially: FY2017 0.897
(SE 0.019), **FY2018 0.056 (SE 0.014, CI 0.034–0.090)**, FY2019 0.941 (SE 0.066), FY2022 0.846
(SE 0.021), FY2023 0.657 (SE 0.030), FY2024 0.579 (SE 0.039), FY2025 0.460 (SE 0.051). The FY2018
value is the most notable: a **precisely-estimated near-total loss in a year with no documented
die-off event** in `data/BWDieOffEvents.csv` (only the FY2020 IP5 and FY2022 IP2 crashes are
recorded) — either a real undocumented mortality event at IP5, or an artifact of the classic CJS's
own assumptions (e.g., size-class/detection structure, or thin cohort information once the FY2017
cohort has mostly aged into "adult"). This has **not yet been cross-checked against the raw
known-alive/contact series** for IP5 around FY2018 — a natural next step. Given this, the
Bayesian model's `lphiA[IP5, y] ~ N(mu_phi, sigma_phi)` hierarchical structure — a single, fairly
modest `sigma_phi` shared across seasons — will structurally shrink genuine sharp,
regime-like drops (like FY2018) toward the pond-level mean; it has no mechanism (beyond the
event-specific `dE` term for documented die-offs) to represent an occasional, precise, large
deviation.

**Second: does this mean the classic CJS's (MLE) annual estimates are simply more accurate?** Not
unconditionally — the two estimators fail in different ways, and this data set shows both failure
modes:

- The classic CJS has no shrinkage, so it can track a real regime shift precisely *when there is
  enough data in that interval* — but that is also exactly where it broke down for FY2020/FY2021
  (boundary/degenerate estimates, §3), and the pond-grouped model (`GIEL_ClassicCJS_PondGrouped_Summary.md`
  §3) showed the same quasi-complete-separation-sized SEs are a general risk in crash cells.
- The Bayesian model's partial pooling stabilizes estimates and avoids boundary degeneracy, but at
  the cost of underrepresenting real episodic shocks — it treats survival as varying smoothly
  around a pond-level mean rather than allowing occasional sharp drops.
- The parametric bootstrap GOF for this IP5-only classic CJS was clean (ĉ ≈ 0.91 mean-based,
  bootstrap p = 0.53, §5) — some evidence the freely-varying annual Phi is tracking real signal
  rather than overfitting noise, at least away from the crash-year boundary cases.

**Working conclusion:** the classic CJS's freely-varying annual Phi is probably a better
description of IP5's actual year-to-year survival dynamics, because those dynamics apparently
include real episodic shocks (documented FY2020, and possibly an undocumented FY2018 event) that
the Bayesian model's hierarchical structure is not designed to capture. But this is not a general
claim that MLE point estimates are "more accurate" — the classic CJS's own crash-year estimates
(FY2020, FY2021) are themselves unreliable at the boundary, and its narrower stocked-only,
single-pond design trades the Bayesian model's richer structure (untagged pool, recruitment,
netting joins) for that additional survival flexibility. Any future reconciliation should probably
target giving the Bayesian model's `lphiA` more freedom to vary by season (e.g., an
IP5-specific unstructured year effect, analogous to the XYTE three-stage model's `eps_J`) rather
than treating either estimator's numbers as the presumptively correct one.

## 8. Caveats

- **IP5-only, stocked-only.** Unlike `GIEL_ClassicCJS.R`/`GIEL_ClassicCJS_PondGrouped.R`, this
  model excludes IP2, IP6, and all netting-tagged ("wild") fish. It is a narrow, targeted
  cross-check on one pond's die-off dynamics, not a general population survival model.
- **`gap` is a full annual time factor, not a die-off dummy** — the name is a holdover from the
  variable's original conception; it produces essentially the same estimates as an unconstrained
  `Phi(~time)` model (with a size offset), which is why the FY2020 and FY2021 estimates behave
  like the boundary/separation cells documented for other fully time-varying GIEL models
  (`GIEL_ClassicCJS_PondGrouped_Summary.md` §3).
- **FY2021 → FY2022 survival (0.888, SE = 0) is not a credible estimate** — see §3. Any reporting
  of this model's survival trajectory should flag this interval as unidentifiable rather than as
  "recovery."
- **The Bayesian comparison in §7 uses a simplified 12-month `dE` application** for the FY2020
  season; the more precise, nEvt-weighted version (documented in `AGENTS.md`) gives ≈0.041 instead
  of ≈0.021 for that season. Both are near the empirical floor and the distinction only matters
  for that one cell.
- **This model omits everything the Bayesian model includes beyond survival/detection**: no
  untagged/recruit pool, no netting-based abundance information, no maturation/stage-transition
  structure, and (for IP5 specifically) no netting-tagged fish at all.
- **The §7.1 stocked-only Bayesian re-run was done ad hoc in the R console** and is saved
  separately from the main `.RData` for this script, in
  `data/GIEL_IP5_StockedOnlyBayesianSensitivity.RData` (`annual_surv_stocked_wholeSeason`: 6,000
  posterior draws × 10 FY columns; `bpost_stocked`: posterior draws of `b_post[1,1]`;
  `stocked_summary`/`bpost_stocked_summary`: posterior mean + 95% CI convenience tables; `notes`:
  derivation string) — this file is **not** produced by any script and must be regenerated by hand
  (filtering `BWRobustDesign_MCMC_giel.RData`'s `rd_Fish` to `pond_new == 2 & origin == 1` and
  applying `b_post[1, 1]` across the stocking-FY season, per §7.1) if the underlying Bayesian fit
  is ever refit.
- **The remaining classic-CJS-vs-Bayesian disagreement (FY2018, FY2019, FY2021) is a model-structure
  finding, not a fish-subsetting artifact** (§7.1) — restricting the Bayesian posterior to
  stocked-only fish changed every season's estimate by ≤ 0.006, so any future reconciliation effort
  should target the survival-model architecture (e.g., letting `lphiA` vary more freely by cohort
  at IP5, or fitting the classic CJS's fully time-varying Phi structure inside the Bayesian model)
  rather than the fish population each model draws on.
- **The disagreement is not limited to the two crash-affected rows** (§7.2): excluding FY2020/FY2021,
  the remaining classic CJS estimates are all well-identified (no boundary values, tight SEs) yet
  still swing from 0.90 down to 0.056 (FY2018) and back to 0.94 (FY2019) — real, precisely-estimated
  year-to-year variability that the Bayesian model's single, modest `sigma_phi` will shrink toward
  the pond mean. **Neither estimator should be treated as unconditionally "more accurate"**: the
  classic CJS can track real shocks precisely but degenerates at the boundary in genuine crash
  intervals, while the Bayesian model avoids boundary degeneracy but likely underrepresents sharp,
  episodic drops (including a possible undocumented FY2018 event, not yet cross-checked against the
  raw known-alive series).

## 9. Key files

| File | Contents |
|---|---|
| `model/GIEL_ClassicCJS_RobustDesign_Annual.R` | Build + fit script (IP5 stocked-only data prep, 6-model AIC comparison, subprocess-isolated parametric bootstrap GOF, Bayesian comparison, empirical floor check) |
| `data/GIEL_ClassicCJS_RobustDesign_Annual.RData` | `pond_check` (stocked-fish TL completeness by pond), `giel_ip5` (fish-level data + encounter history), `dp_ip5`/`ddl_ip5` (processed data, design data), `model_specs_ip5` (6-model formula list), `fits_ip5` (all 6 fitted `crm` objects), `aic_ip5` (AIC ranking), `top_name_ip5`/`top_ip5` (AIC-best model, `gap_pdot`), `preds_ip5` (real-scale Phi/p estimates + CIs), `n_suppressed`, `deviance_obs`/`dev_boot`/`c_hat_mean`/`c_hat_median`/`boot_p_value`/`n_crashed_or_failed` (bootstrap GOF), `annual_surv_bayes` (Bayesian comparison, §7), `floor_2020` (empirical floor, §6) |
| `model/GIEL_ClassicCJS_RobustDesign_Annual_bootGOF.png` | Bootstrap deviance histogram with the observed deviance marked |
| `model/GIEL_ClassicCJS_RobustDesign.R` | Earlier, undocumented monthly robust-design attempt that motivated this script (see §0); not part of the analysis pipeline |
| `data/GIEL_IP5_StockedOnlyBayesianSensitivity.RData` | Stocked-only Bayesian re-run objects (§7.1): `annual_surv_stocked_wholeSeason` (6,000 × 10 posterior draws), `bpost_stocked` (posterior draws of `b_post[1,1]`), `stocked_summary`/`bpost_stocked_summary` (mean + 95% CI), `notes` (derivation string). Not produced by any script — see §8 for how to regenerate. |
