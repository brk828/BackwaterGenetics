# Bonytail Growing-Season Growth: Pond/Year Deviations and Die-off Cross-Reference

**Companion to:** `model/GIEL_GrowthSeasonThreshold.R` / `data/GIEL_GrowthSeasonThreshold.RData`
**Date:** 2026-09-12

> **Note:** the script was rescoped on 2026-09-12 to use capture/recapture pairs across the
> full bonytail size range (not just fish first tagged <250 mm TL) and refit as a Bayesian
> von Bertalanffy (Fabens) growth-increment model. Sections 1-7 below describe the earlier,
> smaller-fish-only linear-increment model and have not yet been rewritten for the rescoped
> analysis -- read their specific pond-year deviation numbers as describing that retired
> model. Sections 0-0c (immediately below) are current as of the rescope and hold the
> **final, adopted** P(TL >= 250 mm) predictions (§0c) -- BK decided (2026-09-12) to move
> forward with the current model's predictions rather than pursue further `Linf`/`K`
> restructuring; the identifiability caveat in §0b must be carried forward wherever §0c's
> numbers are cited.

## 0. Data Quality: Flagged Negative-Growth Records (for database correction)

While building growth-increment pairs from every consecutive pair of measured `StudyBWNFWG`
records per fish (all 3 bonytail ponds, all sizes), **11 pairs show TL2 < TL1** (the fish
appears to have shrunk). Bonytail do not shrink, so these are almost certainly measurement
or transcription errors in the underlying NFWG length records, not real biology. All 11 are
dropped from the growth model (`gs_data`, `gs_fit_data`) and are saved separately as
`NegativeGrowthPairs` in `data/GIEL_GrowthSeasonThreshold.RData` for database review/correction.

| PIT | Pond | Origin (TL1 record) | Date1 | TL1 (mm) | Date2 | TL2 (mm) | Increment (mm) | Days elapsed |
|---|---|---|---|---|---|---|---|---|
| 3B2398A924 | 5 | Stocking | 2022-03-09 | 250 | 2022-12-14 | 210 | -40 | 280 |
| 003D771ED9 | 6 | Capture | 2021-02-05 | 318 | 2021-04-01 | 315 | -3 | 55 |
| 003D795E74 | 6 | Capture | 2021-02-05 | 377 | 2021-04-01 | 375 | -2 | 55 |
| 003D772A7C | 6 | Capture | 2021-12-14 | 362 | 2022-03-09 | 360 | -2 | 85 |
| 003D772A94 | 6 | Capture | 2021-12-14 | 382 | 2022-03-10 | 378 | -4 | 86 |
| 003D772A9F | 6 | Capture | 2021-12-14 | 348 | 2022-03-09 | 345 | -3 | 85 |
| 003D772AA3 | 6 | Capture | 2021-12-14 | 340 | 2022-03-10 | 338 | -2 | 86 |
| **003D772AA4** | **6** | **Capture** | **2021-12-14** | **362** | **2022-03-09** | **258** | **-104** | **85** |
| 003D772AAA | 6 | Capture | 2021-12-14 | 364 | 2022-03-09 | 362 | -2 | 85 |
| 003D772ABB | 6 | Capture | 2021-12-14 | 388 | 2022-03-09 | 385 | -3 | 85 |
| 003D772A79 | 6 | Capture | 2023-12-19 | 420 | 2024-01-30 | 419 | -1 | 42 |

**Patterns worth noting for database review:**
- **9 of 11** involve Pond 6 captures measured in **Dec 2021 or Dec 2023**, re-measured
  ~6-12 weeks later (Feb-Mar) at a slightly *smaller* TL (1-4 mm smaller in 8 of these 9
  cases). A 1-4 mm discrepancy this consistent, this soon after the first measurement,
  looks like measurement/rounding noise (different netting crews or board precision)
  rather than a data-entry swap -- worth a lower-priority check.
  **`003D772AA4` is the exception and the highest-priority record to check**: same
  Dec 2021 → Mar 2022 window as the other Pond 6 pairs, but a 104 mm drop (362 → 258 mm)
  is far outside plausible measurement error and should be checked against the raw NFWG
  record for a transcription error (e.g., a digit swap, 358→258, or a mismatched PIT/length
  pairing).
- **`3B2398A924`** (Pond 5, stocked at 250 mm Mar 2022, recaptured at 210 mm Dec 2022, a
  full growing season later with a -40 mm apparent change) is also worth checking --
  unlike the Pond 6 cluster, this spans a full season, so a real fish this size should not
  measure smaller than at stocking regardless of measurement noise.
- Because only 1 of these 11 (`3B2398A924`) spans a growing season (`gs_days = 182`), the
  other 10 would have been excluded from the growth-rate fit anyway by the existing
  winter-only (`gs_days = 0`) filter -- dropping all 11 explicitly here is a data-quality
  correction, not a change to which rows previously entered the model, except for
  `3B2398A924`.

## 0b. Investigation: Linf posterior sits below the observed maximum size

The fitted model's global `Linf` (417.7 mm, SD 6.3) is below the largest observed
post-growth length in the data (TL2 = 462 mm) -- **every posterior draw** has
`Linf < max(TL2)`, which is not biologically sensible for a true asymptotic size. Three
alternative structures were tested to see whether this is fixable, none fit as well as the
original single-Linf, loose-prior model:

| Model | WAIC | Linf (mean, SD) |
|---|---|---|
| **Global Linf, loose prior N(635, 50) — current model** | **3309.3** | **417.7 (6.3)** |
| Per-pond Linf, loose prior N(635, 50) | 3311.7 | Pond 2: 627.2 (50.8, ≈ prior); Pond 5: 474.4 (27.7); Pond 6: 410.6 (6.6) |
| Global Linf, tight prior N(635, 15) | 3458.7 | 559.5 (20.8) |

**Diagnosis.** Only 3 of the 361 fitted rows have `TL1 > 400 mm` (all Pond 6, 2022 cohort,
still growing 17-42 mm per season at that size), against hundreds of smaller/mid-size fish
that pull toward a much smaller apparent asymptote. Letting each pond have its own `Linf`
does not resolve this: Pond 6 in isolation converges to an *even lower* `Linf` (411 mm,
tied WAIC) because the same small-fish-vs-few-large-fish tension exists within that single
pond. Forcing `Linf` toward the literature range (62-65 cm) with a tighter prior costs
~149 WAIC units and produces clear systematic residuals -- small fish (TL1 < 200 mm)
underpredicted by ~50 mm on average, mid-size fish (300-400 mm) overpredicted by ~17 mm --
a substantially worse overall fit to the bulk of the data.

This pattern matches the well-documented bias of the Fabens estimator toward low `Linf` /
high `K` when individual fish vary in their own growth trajectories (not all following one
population-wide curve) and few recaptures exist near the true asymptote -- exactly the
situation here (3 of 361 rows above 400 mm). **Recommendation:** do not read the fitted
`Linf`/`K` as literal biological growth-curve parameters. The retained model (single global
`Linf`, loose prior) is used only for interpolated increment predictions within the
well-populated size range (roughly 117-350 mm, where residuals are small); it should not be
used to extrapolate growth near or beyond 400 mm. Getting a reliable, unbiased `Linf`
estimate would require more recaptures of fish already above ~400-450 mm, which the current
data do not provide.

Alternative-fit objects (`data/GIEL_growth_pondLinf_fit.rds`, `data/GIEL_growth_tightLinf_fit.rds`)
are saved for reference but are not part of the main pipeline output.

## 0c. Final Results: P(TL >= 250 mm after one Apr-Sep growth season | TL1)

**Decision (2026-09-12, BK):** move forward reporting with the current model's predictions
(single global `Linf`, loose prior, `origin + pondyear` K structure -- the WAIC-best fit
from the 361-row, full-size-range dataset). No further restructuring of `Linf`/`K` is
planned; see §0b for why the alternatives tried did not improve on this fit.

| TL1 (mm) | Expected increment (mm) | P(TL2 >= 250 mm) |
|---|---|---|
| 150 | 99.9 | 0.491 |
| 175 | 90.5 | 0.574 |
| 200 | 81.2 | 0.676 |
| 225 | 71.9 | 0.800 |
| 249 | 63.0 | 0.917 |

> **Caveat to carry forward with these numbers wherever they are cited:** the model's
> `Linf` (417.7 mm) and implied `K` should **not** be read as literal biological growth-curve
> parameters (see §0b) -- they are a fitted compromise driven mostly by hundreds of
> small/mid-size (117-350 mm) recapture pairs, with almost no data (3 of 361 rows) above
> 400 mm to anchor a true asymptote. The table above is valid because every TL1 value in it
> (150-249 mm) sits well inside the well-populated, low-residual size range the model was
> checked against -- these are **interpolated** predictions, not extrapolations. Any future
> use of this table, or of `growth_season_threshold_probs` in `data/GIEL_GrowthSeasonThreshold.RData`,
> should repeat this caveat and should not use the model to predict growth for TL1 values
> approaching or beyond ~400 mm.

## 1. Background

The growing-season model (`model/GIEL_GrowthSeasonThreshold.R`) uses 29 recapture pairs
(fish first tagged <250 mm TL at IP2, IP5, or IP6, with a later physical recapture, after
dropping 12 pairs whose window fell entirely in Oct-Mar). Pond and year of first capture
are combined into a single partially-pooled `PondYear` random-slope term:

```
increment ~ Normal(rate + u[pondyear], sigma) * gs_days
u[pondyear] ~ Normal(0, sigma_py)
```

`sigma_py = 0.242 mm/day` (posterior mean) — a ~40% CV relative to the population rate of
0.561 mm/day — so pond-year does move the fitted rate meaningfully, and WAIC supports
keeping it (see `AGENTS.md`, GIEL Growth section). This note unpacks that variation.

**Important caveat on the `PondYear` label itself:** several recapture pairs have a long
elapsed window (up to 728 days = 4 growing seasons) and are labeled by the *year of first
capture*, not by which season(s) actually contributed the growth. The table below lists
the true `season_years` (the Apr-Sep season(s) with nonzero overlap) for every point, since
this matters for lining up with die-off timing.

## 2. Per-pond-year deviations, ranked

| PondYear | Posterior mean `u` | SD | n points |
|---|---|---|---|
| **P6_2017** | **-0.315** | 0.119 | 3 |
| P6_2018 | -0.096 | 0.112 | 5 |
| P2_2018 | -0.087 | 0.183 | 1 |
| P5_2017 | -0.026 | 0.127 | 2 |
| P6_2022 | -0.007 | 0.122 | 4 |
| P5_2018 | +0.011 | 0.121 | 5 |
| P6_2021 | +0.030 | 0.113 | 6 |
| **P6_2019** | **+0.239** | 0.207 | 1 |
| **P2_2021** | **+0.249** | 0.174 | 2 |

Units are mm/day of growing-season time added to (or subtracted from) the population rate
of 0.561 mm/day. Groups with n = 1-2 (P2_2018, P5_2017, P6_2019, P2_2021) carry wide SDs
and should be read as individual-fish signal, not cohort-level signal.

## 3. All 29 points, with implied simple rate (increment / gs_days)

This is the raw data behind the table above — no partial pooling, just each pair's own
increment divided by its own growing-season exposure, so outliers are visible directly.

| Pond | PondYear | Origin | Date1 | TL1 | Date2 | TL2 | gs_days | Season(s) | Increment | Implied rate (mm/d) |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | P2_2018 | Capture | 2018-12-27 | 210 | 2019-12-18 | 279 | 182 | 2019 | 69 | 0.379 |
| 2 | P2_2021 | Capture | 2021-03-09 | 140 | 2022-01-06 | 302 | 182 | 2021 | 162 | 0.890 |
| 2 | P2_2021 | Capture | 2021-03-10 | 132 | 2022-01-06 | 331 | 182 | 2021 | 199 | 1.093 |
| 5 | P5_2017 | Stocking | 2017-03-21 | 244 | 2018-03-28 | 370 | 182 | 2017 | 126 | 0.692 |
| 5 | P5_2017 | Stocking | 2017-03-21 | 248 | 2018-12-20 | 427 | 364 | 2017,2018 | 179 | 0.492 |
| 5 | P5_2018 | Capture | 2018-03-28 | 249 | 2018-12-20 | 367 | 182 | 2018 | 118 | 0.648 |
| 5 | P5_2018 | Capture | 2018-03-28 | 230 | 2018-12-20 | 351 | 182 | 2018 | 121 | 0.665 |
| 5 | P5_2018 | Capture | 2018-03-28 | 179 | 2018-12-20 | 325 | 182 | 2018 | 146 | 0.802 |
| 5 | P5_2018 | Capture | 2018-03-28 | 209 | 2019-12-05 | 386 | 364 | 2018,2019 | 177 | 0.486 |
| 5 | P5_2018 | Capture | 2018-12-20 | 242 | 2019-12-05 | 349 | 182 | 2019 | 107 | 0.588 |
| 6 | P6_2017 | Stocking | 2017-03-21 | 242 | 2017-11-29 | 340 | 182 | 2017 | 98 | 0.538 |
| 6 | P6_2017 | Stocking | 2017-03-21 | 245 | 2017-11-29 | 354 | 182 | 2017 | 109 | 0.599 |
| 6 | **P6_2017** | Stocking | 2017-03-21 | 240 | **2021-04-01** | 365 | **728** | **2017,2018,2019,2020,2021** | 125 | **0.172** |
| 6 | P6_2018 | Capture | 2018-12-19 | 172 | 2022-03-09 | 400 | 546 | 2019,2020,2021 | 228 | 0.418 |
| 6 | P6_2018 | Capture | 2018-12-19 | 203 | 2021-04-01 | 359 | 364 | 2019,2020,2021 | 156 | 0.429 |
| 6 | P6_2018 | Capture | 2018-12-19 | 241 | 2019-12-05 | 340 | 182 | 2019 | 99 | 0.544 |
| 6 | P6_2018 | Capture | 2018-12-19 | 186 | 2019-12-05 | 325 | 182 | 2019 | 139 | 0.764 |
| 6 | P6_2018 | Capture | 2018-12-19 | 223 | 2019-12-05 | 323 | 182 | 2019 | 100 | 0.549 |
| 6 | **P6_2019** | Capture | 2019-12-05 | 145 | 2021-04-01 | 346 | 182 | **2020,2021** | 201 | **1.104** |
| 6 | P6_2021 | Capture | 2021-02-05 | 117 | 2022-03-10 | 278 | 182 | 2021 | 161 | 0.885 |
| 6 | P6_2021 | Capture | 2021-04-01 | 249 | 2022-03-10 | 337 | 182 | 2021 | 88 | 0.484 |
| 6 | P6_2021 | Capture | 2021-12-14 | 165 | 2022-12-15 | 304 | 182 | 2022 | 139 | 0.764 |
| 6 | P6_2021 | Capture | 2021-12-14 | 162 | 2022-12-15 | 329 | 182 | 2022 | 167 | 0.918 |
| 6 | P6_2021 | Capture | 2021-12-14 | 237 | 2023-12-19 | 408 | 364 | 2022,2023 | 171 | 0.470 |
| 6 | P6_2021 | Capture | 2021-12-14 | 155 | 2023-12-19 | 360 | 364 | 2022,2023 | 205 | 0.563 |
| 6 | P6_2022 | Capture | 2022-03-09 | 231 | 2023-12-19 | 379 | 364 | 2022,2023 | 148 | 0.407 |
| 6 | P6_2022 | Capture | 2022-03-10 | 225 | 2022-12-15 | 343 | 182 | 2022 | 118 | 0.648 |
| 6 | P6_2022 | Capture | 2022-12-15 | 167 | 2023-12-19 | 325 | 182 | 2023 | 158 | 0.868 |
| 6 | P6_2022 | Capture | 2022-12-15 | 170 | 2023-12-19 | 305 | 182 | 2023 | 135 | 0.742 |

Population mean rate (all 29 points pooled, no grouping): 0.561 mm/day.

## 4. What is driving the two most extreme deviations

**P6_2017 (slow, -0.315, n = 3).** Two of the three points are unremarkable one-season
(182-day) 2017 stockings with rates 0.538 and 0.599 mm/day — close to or above the
population mean. **The negative deviation is driven almost entirely by the third fish**
(PIT `0015AC37D8`), whose single recapture spans 728 growing-season days — parts of *five*
calendar years (2017-2021) — averaging only 0.172 mm/day, less than a third of the
population rate. Because this one fish carries 4x the exposure of the other two combined,
it dominates the group's fitted rate under the linear-in-`gs_days` model. This is a single
data point's multi-year average, not three fish agreeing on a slow year.

**P6_2019 (fast, +0.239, n = 1).** A single fish (PIT `003D795E5E`) recaptured after
almost exactly the 2020 growing season (`season_years = "2020,2021"`, but 2021's
contribution is ~0 days — Date2 is 2021-04-01, the first day of that season) grew 201 mm in
182 growing-season days: **1.10 mm/day, the single fastest rate in the dataset**, roughly
double the population mean. With n = 1 this is pure individual variation until more data
exist; it happens to be the only Pond 6 point that isolates the 2020 season specifically.

**P2_2021 (fast, +0.249, n = 2, internally consistent).** Both points are 2021-season
recaptures at Pond 2, with rates 0.890 and 1.093 mm/day — unlike the two n=1 groups above,
these two independent fish agree with each other, so this is a real pond-year signal, not
a single outlier. Notably, **this is the growing season immediately before IP2's 2022
die-off** (see below) — worth flagging as a candidate for WQ follow-up: fast growth in
survivors the year before a crash is at least consistent with reduced competition from a
population already in decline, though this is speculative without density or WQ data for
that specific window.

## 5. Cross-reference with documented bonytail die-offs

`data/BWDieOffEvents.csv` and the project's Bayesian robust-design model record confirmed
near-total mortality events at **IP5 (2020)** and **IP2 (2022)** (see `AGENTS.md`, Bonytail
die-off section). Monthly known-alive counts (`TotalCountIPGIEL`) pin down the crash months
precisely — both fit the "summer heat, then early-autumn DO crash" seasonal pattern you
described:

**IP2 (Pond 2), 2021-2022:**

| Month | Known alive |
|---|---|
| 2022-06 | 54 |
| 2022-07 | 39 |
| **2022-08** | **0** |

**IP5 (Pond 5), 2019-2020:**

| Month | Known alive |
|---|---|
| 2020-07 | 52 |
| 2020-08 | 42 |
| **2020-09** | **0** |

Both crashes land in the Jul-Sep window — consistent with the hypothesized summer
temperature stress carrying into an early-autumn dissolved-oxygen crash from vegetation
die-off.

**A mechanical, but noteworthy, consequence for the growth model:** there is **no Pond 2
recapture pair covering the 2022 growing season, and no Pond 5 recapture pair covering the
2020 growing season** anywhere in the 41-fish dataset (the closest are P2_2021, ending
January 2022 just before the crash, and the P5_2017/P5_2018 multi-season points, which stop
at 2018-2019). This isn't additional evidence of the die-offs — it's the expected
consequence of them: by the time growth could have been measured across those seasons, the
ponds had almost no survivors left to recapture. It does mean the growth model is silent
about growth conditions during the crash seasons themselves at IP2 and IP5.

**Pond 6 has no documented die-off**, and its monthly known-alive counts from 2017-2023
show smooth, gradual decline with no month-to-month acute drop resembling the IP2/IP5
signature (largest single-month change in this period is ~10-15 fish out of 60-260, versus
IP2's 39→0 and IP5's 42→0 in one month). So the P6_2017 slow outlier and P6_2019 fast
outlier described above have no independent population-level corroboration at Pond 6 --
they read as individual variation, not evidence of an unflagged Pond 6 event, given
currently available data.

## 6. Hypotheses to check once water quality data are available

1. **IP2 Jul-Aug 2022 and IP5 Aug-Sep 2020** are the two calendar windows to prioritize for
   temperature and dissolved-oxygen records — both die-offs completed within a single
   month, in the Jul-Sep range consistent with the proposed summer-heat/autumn-DO-crash
   mechanism.
2. **P2_2021's above-average growth (0.89-1.09 mm/day, 2 fish agreeing)**, in the season
   immediately preceding IP2's crash, is worth checking against any WQ decline already
   underway in 2021 (sub-lethal stress can sometimes still permit growth, or reduced
   competition from an already-declining population could accelerate it) — currently only
   a hypothesis.
3. **P6_2017's single slow multi-season fish (0.172 mm/day, 2017-2021)** has no
   corroborating Pond 6 mortality signal, and with n = 1 spanning 4 seasons it cannot be
   localized to any specific year — if finer-grained WQ data become available, this fish's
   window is too broad to test against any single event.
4. **P6_2019's single fast fish (1.10 mm/day, 2020 season)** is not correlated with any
   known Pond 6 stress and, being the pond's only 2020-isolated data point, cannot be
   distinguished from ordinary individual growth variation.

## 7. Limitations

- Recapture pairs are an opportunistic subset of physically recaptured fish (29 of 936
  originally tagged <250 mm across all three ponds), not a random sample — capture effort
  itself may correlate with pond conditions in ways this analysis cannot detect.
- Several `PondYear` groups have only 1-2 points; their deviations are effectively
  individual-fish measurements, not cohort estimates, and should not be over-interpreted as
  pond-level or year-level effects.
- The linear-in-`gs_days` growth model, fit through the origin, does not allow rate to vary
  *within* a multi-season point (e.g., faster in year 1, slower in year 2) — a single
  average is all that's identifiable from two length measurements, which is why the
  728-day and 546-day points can't be attributed to a specific season within their span.
- No water quality data exist yet for any of these ponds; every WQ linkage above is a
  hypothesis to test, not a finding.
