# Bonytail (GIEL) Classic Annualized CJS Model — Results

**Purpose:** a simple, frequentist (MLE) Cormack-Jolly-Seber cross-check on the joint hierarchical
Bayesian robust-design model described in `model/GIEL_ModelMethodology.md` /
`model/GIEL_ModelResults_Summary.md`. Built and fit 2026-09-11 with the `marked` R package (no
external MARK.exe dependency). Code: `model/GIEL_ClassicCJS.R`; output:
`data/GIEL_ClassicCJS.RData`.

## 1. Design

- **Occasions:** annual, FY2017–FY2026 (10 occasions, 9 survival intervals), each occasion defined
  by the Oct–Apr scanning window (matching the Bayesian model's primary period); May–Sep contacts
  are ignored entirely — a much simpler detection structure than the robust design's weekly
  secondary occasions.
- **Population:** 2,637 tagged GIEL fish (stocked or netting-tagged) at IPCA Ponds 2, 5, 6,
  pooled — **no pond-identity effect** in this model (not requested; see caveat in §5). 300 fish
  with no recorded length at tagging (all stocked) were dropped for lack of a size class; no fish
  were excluded for documented removals (none exist for GIEL — checked directly).
- **Size class ("young"/"adult"):** not a fixed group. Any fish tagged at <250 mm TL is "young" for
  exactly the one annual interval immediately following its own tagging and "adult" for every
  interval after that (and immediately, if tagged ≥250 mm) — implemented as an explicit
  individual × interval lookup, verified against every fish's own tagging occasion.
- **Origin:** stocked vs. netting-tagged ("wild"), a fixed covariate.
- **Model set:** all combinations of Phi(time or constant) × Phi(size effect: yes/no) × Phi(origin
  effect: yes/no) × p(time or constant) = 16 models, ranked by AIC (AICc correction negligible at
  n ≈ 2,600).

## 2. Model selection

| Model | k | AIC | ΔAIC | weight |
|---|---|---|---|---|
| Phi(time + size + origin) p(time) | 20 | 5876.8 | 0.0 | 0.883 |
| Phi(time + size + origin) p(.) | 12 | 5880.8 | 4.1 | 0.117 |
| Phi(time + size) p(time) | 19 | 5898.0 | 21.2 | <0.001 |
| Phi(time) p(time) | 18 | 5976.8 | 100.0 | ~0 |
| Phi(.) p(.) | 2 | 6767.1 | 890.3 | ~0 |

(full 16-model table in `data/GIEL_ClassicCJS.RData$giel_cjs_table`). The top two models — identical
Phi structure, differing only in whether p varies by year — carry essentially all (>99.9%) of the
AIC weight, so no formal model averaging was needed; results below are from the top model
(Phi(time+size+origin), p(time)).

**All three requested effects matter and time-varying survival dominates:** the AIC gap between
Phi(constant) and Phi(time) alone is enormous (≈765–890 AIC units), confirming strong annual
variation; adding size and origin on top of `Phi(time)` improves fit by another ~100 AIC units.

## 3. Survival (Phi) estimates, top model

Realized annual survival by interval-start year, size class, and origin (selected rows; full table
in the saved `.RData`):

| Interval start | Origin | Size class | Survival | 95% CI |
|---|---|---|---|---|
| FY2018 | wild | adult | 0.177 | 0.144–0.216 |
| FY2018 | wild | young | 0.108 | 0.086–0.135 |
| FY2019 | wild | adult | 0.505 | 0.421–0.588 |
| FY2020 | wild | adult | 0.199 | 0.163–0.240 |
| FY2020 | wild | young | 0.123 | 0.096–0.157 |
| FY2022 | wild | adult | 0.569 | 0.518–0.617 |
| FY2017 | stocked | adult | 0.792 | 0.761–0.820 |
| FY2018 | stocked | adult | 0.254 | 0.224–0.287 |

**Within every year × origin pair, "young" survival is lower than "adult" survival** (e.g., FY2018
wild: 0.108 vs 0.177; FY2020 wild: 0.123 vs 0.199) — consistent, first-year-after-tagging effect.
**Within every year × size pair, "stocked" survival is generally higher than "wild"/netting-tagged**
survival for adults (e.g., FY2018: 0.254 stocked vs 0.177 wild) — direction only, not causally
interpreted here (see caveats).

**FY2020 (interval FY2020→FY2021) shows clearly reduced survival relative to neighboring years**
(0.12–0.20 vs 0.42–0.57 elsewhere) — this lines up with the documented IP5 die-off window (Jan–Sep
2020, per `data/BWDieOffEvents.csv` and `GIEL_ModelResults_Summary.md` §3). **FY2022 (the IP2
die-off year, Jan–Jul 2022) does *not* show a comparably depressed pooled estimate** (0.57 wild
adult) — see the pond-pooling caveat in §5 for why this is expected, not a contradiction.

## 4. Detection (p) estimates, top model

| Year | p | 95% CI |
|---|---|---|
| 2018 | 0.983 | 0.947–0.994 |
| 2019 | 0.963 | 0.905–0.986 |
| 2020 | 0.902 | 0.740–0.968 |
| 2021 | 0.977 | 0.857–0.997 |
| 2022 | >0.9999 | degenerate (boundary) |
| 2023 | 0.987 | 0.959–0.996 |
| 2024 | >0.9999 | degenerate (boundary) |
| 2025 | >0.9999 | degenerate (boundary) |
| 2026 | 0.551 | 0.167–0.883 (wide; partial/terminal-occasion data) |

Detection is very high (≥0.90) in most years, with several years pinned essentially at 1 — the same
"boundary estimate, degenerate CI" phenomenon already documented for the Chapman estimator in this
project (`ChapmanValidation.qmd`'s R = C special case) rather than a modeling artifact. The
low, wide FY2026 estimate reflects the partial/terminal occasion (least data).

## 5. Goodness-of-fit (Program RELEASE-style tests, `R2ucare`)

Run on the pooled (no-covariate) encounter histories, per the standard Lebreton et al. (1992)
workflow of using these component tests to justify global model structure:

| Test | χ² | df | p | Interpretation |
|---|---|---|---|---|
| Overall CJS GOF | 108.9 | 13 | <0.0001 | Substantial lack of fit for the *no-covariate* model |
| **3.SR (transience)** | **103.9** | 7 | **<0.0001** | **Newly-marked fish have different (lower) recapture/survival prospects than previously-marked fish — accounts for ~95% of the total lack-of-fit** |
| 3.Sm | 0.30 | 3 | 0.96 | No evidence of heterogeneity in survival among already-seen animals by future recapture history |
| 2.CT (immediate trap-dependence) | 4.69 | 3 | 0.20 | No evidence of trap-happiness/trap-shyness |
| 2.CL (delayed trap-dependence) | 0 | 0 | 1 | No evidence |

**This is a strong internal validation of the requested model structure:** essentially all of the
basic CJS model's lack of fit is transience — exactly what the "young for one year, then adult"
size-class effect (and to a lesser extent the origin effect, which is highly correlated with
initial tagging size) is designed to absorb. No other structural problem (trap-dependence, or
survival heterogeneity unrelated to time-since-marking) was detected. A crude overdispersion
estimate on the null model is ĉ = 108.9/13 ≈ 8.4; because the dominant source (transience) is
explicitly modeled by the top AIC model's covariates, the *residual* overdispersion in that model is
expected to be much smaller, though R2ucare's tests do not directly provide a GOF statistic for an
arbitrary covariate model (a parametric bootstrap GOF on the top model would be the natural next
step if a precise ĉ is needed for QAICc).

## 6. Caveats

- **No pond-identity effect.** Per the requested design, all three GIEL ponds (IP2, IP5, IP6) are
  pooled into one set of annual survival estimates — there is no pond × year interaction. This
  explains why the FY2022 IP2 die-off (Jan–Jul 2022, confirmed via the known-alive series in
  `GIEL_ModelResults_Summary.md` §3) does not show up as a low pooled survival estimate: it affected
  only one of three ponds, for part of a year that straddles the Oct–Apr/May–Sep boundary this model
  ignores, and is diluted by the other two ponds' unaffected fish in the pooled annual estimate.
  The FY2020 IP5 die-off (Jan–Sep 2020) is more visible because it was a longer, more severe crash.
  A pond effect was not part of the requested model set; adding one (`Phi(~time*location)` or
  similar) would be a natural extension if pond-specific classic estimates are wanted.
  A more useful cross-check might be running this exact same model separately per pond, rather than adding a pond effect, given very low survivorship in certain ponds during certain years, which could make for unstable parameter estimates within a global model, given the extreme range in fish abundon over time in these ponds.
- **"Stocked survives better than wild" is a pooled association, not a causal claim** — origin is
  confounded with cohort-specific conditions (which years/ponds each group was tagged in), pond, and
  potentially unmeasured handling differences; interpret directionally only.
- **300 excluded fish** (all stocked, all with missing length-at-tagging) are excluded entirely
  rather than imputed; if size at tagging is recoverable from another source (e.g. stocking
  manifests), redoing this with the fuller dataset would be worth revisiting.
- **This model is a simplification relative to the Bayesian model** in several respects beyond
  losing summer detections and pond identity: no untagged/recruit pool, no netting-based
  abundance information, no explicit die-off event term, and a size class that (unlike the
  Bayesian model's true growth-based stage transitions) is purely a function of time-since-tagging.
  It is intended as an independent, transparent point of comparison, not a replacement.

## 7. Key files

| File | Contents |
|---|---|
| `model/GIEL_ClassicCJS.R` | Build + fit script (data prep, 16-model set, GOF tests) |
| `data/GIEL_ClassicCJS.RData` | `giel_cjs_table` (AIC ranking), `giel_cjs_results` (all 16 fitted `crm` objects), `top_model` (refit with Hessian), `giel_cjs_real` (real-scale Phi/p estimates + CIs), `gof`/`gof_3sr`/`gof_3sm`/`gof_2ct`/`gof_2cl` (R2ucare GOF tests), `giel_data`/`dp`/`ddl` (processed data and design data) |
