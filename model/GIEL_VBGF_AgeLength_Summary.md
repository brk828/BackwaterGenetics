# Bonytail Age-Based von Bertalanffy Growth Model

**Script:** `model/GIEL_VBGF_AgeLength.R` | **Data:** `data/GIEL_VBGF_AgeLength.RData`
**Date:** 2026-09-14 (BK request); revised twice same day:
1. Anchor `L0` on a literature hatch-length value (see "L0 anchoring" section below).
2. **Shared-K, release-anchored revision** (see that section below) -- the current
   version. This drops the origin-specific growth rate entirely and anchors stocked
   fish's growth at their own release size/date rather than an assumed hatch date,
   based on a direct empirical test showing no credible growth-rate difference
   between stocked and wild-recruit bonytail once in the pond.

Two earlier fits are archived for comparison: the original `t0`-formulation fit
(`data/GIEL_VBGF_AgeLength_t0.RData`) and the `L0`-anchored, hatch-age,
origin-differentiated fit (`data/GIEL_VBGF_AgeLength_originDiff.RData`).

## Motivation

The existing bonytail growth model (`model/GIEL_GrowthSeasonThreshold.R`) fits a Fabens
growth-*increment* model to opportunistic recapture pairs, with no absolute age information.
That model's `Linf` (417.7 mm) sits below the largest observed length in the data (462 mm) --
an unresolved identifiability failure traced to having almost no recaptured fish near the true
size asymptote (see that script's header and `GIEL_GrowthSeasonThreshold_Summary.md` §0b).
BK proposed two previously-unused age data sources for bonytail that could fix this by
anchoring the growth curve with real ages instead of only paired size differences.

## Age data sources

1. **Known stocking-cohort year class.** `StudyBWNFWG$year_class` is populated for two GIEL
   stockings: 2017-03-21 (899 records, all 3 ponds, year class 2012) and 2023-04-12 (300
   records, Pond 2 only, year class 2018). Both cohorts were already **~5.2 years old at
   release** (long captive rearing before stocking is typical for bonytail hatchery
   production), at TL 234-299 mm for the 2012 cohort. A third stocking (2022-03-09, 300 fish,
   Pond 5) has **no recorded year class** -- flagged for database follow-up, not fixed here.
   `year_class` is already propagated in the source data to every later NFWG record
   (recaptures) for these PITs, so every measured record of a known-cohort fish is a valid
   age-length point, not just the stocking record itself.
2. **Assumed young-of-year (YOY) fall recruits.** Untagged fish first captured Oct-Dec with no
   `year_class` anywhere in their record history show real bimodal size distributions in
   several fiscal years. A 2-component Gaussian mixture (`mclust::Mclust`, `G = 1:2` selected
   by BIC) was fit per FY on log(TL); where `G = 2` was preferred, fish classified to the
   lower-mean cluster were assigned `year_class` = calendar year of capture (assumed spring
   birth that year).

   | FY | n candidates | Decision | n assigned YOY | Mean classification uncertainty |
   |---|---|---|---|---|
   | 2018 | 9 | too few fish, excluded | -- | -- |
   | 2019 | 101 | included (marginal) | 70 | 0.077 |
   | 2020 | 334 | included | 79 | 0.085 |
   | 2021 | 0 | no fall captures (die-off aftermath) | -- | -- |
   | 2022 | 54 | included | 27 | 0.056 |
   | 2023 | 94 | included | 16 | 0.011 |
   | 2024 | 356 | included | 305 | 0.0004 |

   FY2019 has a narrower gap between cluster edges than the other years (low-cluster max 228 mm
   vs high-cluster min 233 mm) and should be read with slightly less confidence; all other
   included years have low classification uncertainty.

## Combined dataset

The two sources were merged into one year-class map (no PIT overlap by construction, 1,199
known + 497 assumed = 1,696 fish) and used to compute age at **every** NFWG record for those
fish: `age_years = (collection_date - April 1 of year_class) / 365.25`. April 1 is a nominal
spawning-date assumption, chosen for consistency with this project's existing Apr-Sep
"growing season" boundary convention used elsewhere (joint robust-design / YCB3S models) --
it is not a literature-sourced bonytail spawning date.

Result: **1,601 age-length rows** across **1,476 distinct fish** (1,365 fish contribute one
row, 111 contribute 2-3 repeat measurements, treated as independent rows -- the same
simplification used in the increment model).

| Source | n rows | n fish | Age range (yr) |
|---|---|---|---|
| Assumed YOY fall capture | 525 | 497 | 0.66 - 3.94 |
| Known stocking cohort | 1,076 | 979 | 4.97 - 10.7 |

**Important structural caveat:** the two age sources occupy almost entirely disjoint age
windows (max recruit age 3.94 yr vs min known-cohort age 4.97 yr) -- see the model results
below for why this matters.

## L0 anchoring (revision, 2026-09-14)

The initial draft of this model used a standard `t0` parameterization because it assumed no
literature-sourced bonytail hatch-length value existed to anchor `L0` directly (unlike the
gray triggerfish study; `documents/Chamberlin2025_GrayTriggerfishVBGM_Summary.md`). BK
pointed out this was incorrect: Hamman (1982) reports a mean bonytail hatch length of
**6.8 mm TL** (stated range 6.5-7.5 mm) --

> Hamman, Roger L. "Induced Spawning and Culture of Bonytail Chub." *The Progressive
> Fish-Culturist*, vol. 44, no. 4, 1982, pp. 201-203.
> doi:10.1577/1548-8659(1982)44[201:ISACOB]2.0.CO;2.

The model below was revised to use the standard `L0`/size-at-hatch VBGF form, with `L0`
given an informative prior centered on this value (`Normal(6.8, 0.25)` mm; the SD is chosen
so the stated 6.5-7.5 mm range spans roughly a 95% interval), replacing the earlier
essentially-uninformed `t0 ~ Normal(0, 1)` prior. All results in this document are from the
revised, `L0`-anchored fit; the superseded `t0` results are noted briefly for comparison
where relevant. The archived `t0` fit lives in `data/GIEL_VBGF_AgeLength_t0.RData`.

## L0-anchored, origin-differentiated results (superseded)

With the `L0`/hatch-anchor form above (age since assumed hatch date for every fish, origin-
specific `logK`), WAIC clearly favored the full (origin + pondyear) model over pondyear-only,
origin-only, and fully pooled (14,561.0 vs 14,715.9 / 15,382.7 / 16,095.6 -- a large, clean
separation). The final fit (barker-sampler block, 3 chains, 100k iter/30k burn-in/thin 10,
R-hat <= 1.06) gave `Linf` = 486.2 mm (95% CI 463.5-509.4), `logK[StockedCohort]` implying
K = 0.221/yr (0.178-0.260), `logK[AssumedRecruit]` implying K = 0.540/yr (0.443-0.622),
`L0` = 6.79 mm (posterior ~ prior, as expected), residual sigma = 22.7 mm. This was a
substantial change from the earlier `t0` version (`Linf` 703 mm, `t0` = -0.69 yr) because
anchoring the age-0 intercept at a real value requires faster early growth to reach observed
sizes over the 0.66-10.7 yr age range, pulling `Linf` down -- see the archived
`data/GIEL_VBGF_AgeLength_originDiff.RData` and prior git history of this file for full detail
on that comparison.

**This origin-differentiated result was itself superseded the same day** by the
shared-K, release-anchored revision below, after BK raised a structural concern with it (see
next section).

## Shared-K, release-anchored revision (current)

### Motivation: hatchery growth vs. pond growth are conflated in the hatch-anchored model

The L0-anchored fit above found a large, "significant" origin effect: stocked-cohort fish
appeared to grow at less than half the rate of wild recruits (K = 0.221 vs 0.540/yr). BK
raised a structural concern: bonytail stocking cohorts spend **~5.2 years in captive/hatchery
rearing before release**, and hatchery growth conditions are not representative of pond
conditions. If stocked and wild fish actually grow at the **same rate once both are in the
pond**, then modeling stocked-fish growth as one continuous curve from an assumed hatch date
would blend the (different, slower) hatchery-phase growth with pond-phase growth, producing a
spurious "origin effect" that is really a hatchery-vs-pond artifact -- not a true difference in
how bonytail grow once released.

### Direct empirical test (natural experiment)

Almost every pond-year cell in the companion growth-increment dataset
(`GIEL_GrowthSeasonThreshold.R`'s `gs_fit_data`) is single-origin (either all newly-stocked or
all previously-established fish), so origin and pond-year/cohort are almost completely
confounded in general. **Exactly one cell has both:** IP5 (Pond 5), 2022-23, where the 2022
stocking cohort (all released at TL = 250 mm, recaptured Mar-Dec 2022) overlaps in size with
36 previously-established fish measured Dec 2022 (TL1 210-376 mm) and recaptured through
Nov 2023. Restricting both groups to one-growing-season (182-day) recapture pairs and
regressing increment on `TL1 + Origin` (n = 143):

| Term | Estimate | SE | p |
|---|---|---|---|
| TL1 | -0.526 mm/mm | 0.134 | < 0.001 |
| Origin (Stocking vs. Capture) | -18.3 mm | 11.3 | 0.108 |

The origin effect is small and not statistically significant -- much smaller than the ~2.4x K
ratio implied by the hatch-anchored age model. This supports the hypothesis that most of the
apparent origin effect elsewhere is a pond-year/cohort confound, not a real pond-growth-rate
difference. (Caveat noted at the time: only 6 of the 36 established fish fell near the
stocked cohort's exact 250 mm starting size, and the two groups' growing seasons are one
calendar year apart -- 2022 vs. 2023 -- so this is suggestive, not conclusive, on its own.)

### Model revision

Rather than continuing to assume one growth regime from hatch to now for stocked fish, the
model now anchors each stocked fish's growth at **its own release size and date** (a
Fabens-style reset), and drops the origin-specific growth rate entirely in favor of a single
shared `K` (still varying by pond-year/cohort via `v[pondyear]`):

```
Recruit-origin row i:  mu_i = L0 + (Linf - L0) * (1 - exp(-K_i * age_i))
                        (age_i = years since assumed hatch date, unchanged; L0 anchored on
                         Hamman 1982 as before)
Stocked-origin row i:  mu_i = TLrel_i + (Linf - TLrel_i) * (1 - exp(-K_i * t_i))
                        (t_i = years since THIS FISH's own release-anchor date; TLrel_i =
                         that fish's release-anchor TL -- both known data, not modeled)
K_i = exp(logK + v[pondyear_i])          <- single shared logK (primary model)
    = exp(logK[origin_i] + v[pondyear_i]) <- origin-specific logK (sensitivity check only)
TL_i ~ Normal(mu_i, sigma)
```

This directly implements "adjust the effective age of stocked fish for growth purposes": a
stocked fish's growth clock resets at release rather than counting hatchery years, and the
same curve applies to everyone once they are in the pond. It also removes the model's least
defensible assumption (a constant, unobserved hatchery growth rate).

**Release-anchor caveat:** the true stocking-event record has a measured TL for the
2017-03-21 cohort (year class 2012, 899 fish; TL 234-299 mm at t=0, an exact release anchor)
but **not** for the 2023-04-12 cohort (year class 2018, 300 fish -- no `total_length` recorded
at that stocking event). For that cohort the anchor falls back to each fish's own first later
record with a measured TL (in practice Dec 2023, ~244 days post-release, TL 297-406 mm,
n = 80 of 300 fish) -- growth between true release and that first measurement is unobserved
and not attributed to anyone. Practically, this means all 80 of those fish carry **zero**
weight for estimating K (an anchor row is a tautology: `mu = TLrel_i` regardless of K at
t = 0), rather than misleadingly informing a blended rate as they effectively did before. All
real growth information on stocked-origin fish now comes from the 2012 cohort's repeat
measurements: only **89 of 979** release-anchored fish (` n_informative_stocked` in the saved
data) have any post-release recapture at all (81 with one extra measurement, 8 with two).
This is a much smaller *informative* sample than the raw row count suggests -- appropriately
so, since the earlier richness was informing a blended hatchery+pond rate, not a pond-only one.

Priors unchanged: `Linf ~ Normal(635, 50)`, `logK ~ Normal(log(0.15), 0.5)`,
`L0 ~ Normal(6.8, 0.25)` mm (Hamman 1982). Same barker-sampler block over `{Linf, logK[...]}`
as before (requires `buildDerivs = TRUE`); this model converges even better than the
hatch-anchored version (R-hat <= 1.02 vs <= 1.06), likely because the release anchor is a hard
data constraint rather than an estimated parameter interacting with the rest of the curve.

### Results: origin effect vanishes on the full dataset, not just the one IP5 cell

Three models compared by WAIC: the primary shared-K model, a fully-pooled (no pond-year)
variant, and a sensitivity check that re-adds an origin-specific `logK` on top of
release-anchoring (to test whether origin still matters once anchoring is fixed, using the
*entire* 1,601-row dataset rather than just the 143-fish IP5 natural experiment):

| Model | WAIC |
|---|---|
| **shared K + pondyear (primary)** | **13,646.35** |
| origin + pondyear (sensitivity check) | 13,646.87 |
| shared K, pooled (no pondyear) | 14,771.76 |

**The origin sensitivity check is a statistical tie with the shared-K model (Delta WAIC =
0.52)** -- a dramatic reversal from the hatch-anchored version's ~800-unit gap favoring
origin. Pond-year/cohort variation still matters a great deal (Delta WAIC = 1,125 vs. pooled).
In the origin-check model itself, implied K is nearly identical between origins:

| Origin | K (yr^-1) | 95% CI |
|---|---|---|
| StockedCohort | 0.718 | 0.47 - 1.01 |
| AssumedRecruit | 0.768 | 0.63 - 0.91 |

P(K[Stocked] < K[Recruit]) = 0.63 -- essentially a coin flip, not the P = 1.0 (non-overlapping
CIs) found under hatch-anchoring. **This confirms BK's hypothesis on the full dataset:** once
growth is anchored at release rather than hatch, stocked and wild-recruit bonytail show no
credible difference in pond growth rate.

**Primary (shared-K) model posterior** (barker-sampler block, 3 chains, 100k iter/30k
burn-in/thin 10, R-hat <= 1.02 for every parameter):

| Parameter | Mean | SD | 95% CI |
|---|---|---|---|
| `Linf` (mm) | 406.1 | 5.3 | 395.9 - 416.7 |
| `logK` | -0.233 | 0.077 | -0.386 - -0.086 |
| `L0` (mm) | 6.80 | 0.25 | 6.32 - 7.29 |
| `sigma` (mm) | 17.0 | 0.3 | 16.4 - 17.6 |
| `sigma_v` | 0.267 | 0.064 | 0.175 - 0.419 |

Implied shared `K` = 0.792/yr (0.68-0.92). Two further improvements over the hatch-anchored,
origin-differentiated fit: **residual sigma dropped from 22.7 mm to 17.0 mm** (a
better-fitting model with fewer parameters -- not just a simpler one), and `Linf`'s CI is much
narrower (5.3 mm SD vs. 11.8 mm) because release-anchoring uses hard data constraints instead
of an additional estimated parameter per origin.

### New caveat introduced by this revision: Linf is now below the observed maximum length

**`Linf`'s posterior is entirely below the observed maximum TL in the dataset (442 mm) --
the single highest posterior draw across all three chains is 427.5 mm, i.e. 0% posterior
support for `Linf` >= 442 mm.** This is the same identifiability failure already flagged for
the increment (Fabens) model (`GIEL_GrowthSeasonThreshold_Summary.md` section 0b) and
explicitly avoided by the previous (hatch-anchored) age-based fit, whose `Linf` sat safely
above the observed max. It has re-appeared here despite the much better origin/WAIC story.
**Likely cause:** real growth information on large/old fish is thin under release-anchoring --
only 89 of 979 stocked fish have any post-release repeat measurement, and the recruit data
(which drives most of the shared `K` estimate) tops out at 408 mm. Four individual fish (3
stocked, 1 recruit; TL 408-442 mm) are underpredicted by 20-43 mm at posterior-mean
parameters. **Read this model's `Linf` as a fitted compromise for interpolation within the
well-populated size range, not a literal asymptote** -- exactly the same caution already
standing for the increment model and (to a lesser degree) for the GAM-based growth-season
threshold work. The *shape* conclusion (no origin effect once anchored at release) is not
undermined by this -- it comes from the WAIC/K comparison, which is well-identified -- but the
specific `Linf`/`K` values should not be over-interpreted as final growth-curve parameters.

### Critical caveat: age/origin confound is reduced, not eliminated

The original hatch-anchored version had essentially zero age overlap between origins (recruit
ages 0.66-3.94 yr vs. known-cohort hatch-ages 4.97-10.7 yr). Release-anchoring substantially
reduces this: time since release spans 0-9.4 yr (2012 cohort) and 0-3.4 yr (2018 cohort), both
overlapping the recruit age range (0.66-3.94 yr). However, the *informative* (repeat-measured)
stocked-fish data are still sparse (89 fish) and concentrated in the 2012 cohort, so this is an
improvement in principle more than a fully resolved confound in practice -- the IP5 natural
experiment above remains the cleanest direct test because it holds pond and (approximately)
calendar time fixed as well as size.

## Caveats to carry forward

1. **April 1 birth-date assumption** for recruits (unchanged) is a convention borrowed from
   this project's Apr-Sep growing-season boundary, not a literature-sourced bonytail spawning
   date. No longer used for stocked fish (release-anchored instead).
2. **YOY classification is a statistical assumption, not a validated age determination** --
   fish are assigned age-0 status based on falling in a fall size distribution's lower mode,
   not by any direct ageing method (otoliths, etc.). FY2018 is excluded (too few fish) and
   FY2019 is marginal (see table above).
3. **Release-anchor gap for the 2018 year-class cohort** (see "Release-anchor caveat" above)
   -- no measured TL at the true 2023-04-12 stocking event, so those 300 fish's growth from
   release to first measurement (~244 days) is unobserved; none of them are currently
   informative for K.
4. **Repeat measurements per fish are treated as independent rows**, the same simplification
   used in the increment model; a full individual-random-effects growth model would be more
   rigorous but was not attempted.
5. The missing `year_class` for the 2022-03-09 Pond 5 stocking (300 fish) means that cohort
   contributes no age-length data at all under this approach -- worth a database fix if this
   model is extended.
6. **`L0` is anchored on a single point estimate from one hatchery-culture paper** (Hamman
   1982, n and variance not reported in the summary read for this analysis) rather than a
   wild-population estimate for these specific bonytail broodstock; the SD (0.25 mm) reflects
   the paper's stated 6.5-7.5 mm range, not a formally reported measurement uncertainty.
7. **`Linf` is now below the observed maximum length** (see dedicated caveat above) -- treat
   `Linf`/`K` as fitted-compromise values for interpolation, not literal growth-curve
   parameters, matching the standing caveat already in place for the increment model.
8. **The IP5 natural-experiment test (n=143, one pond-year) is thin evidence on its own**;
   the full-dataset origin-check WAIC tie is the stronger piece of evidence for "no origin
   effect," but both point the same direction.
