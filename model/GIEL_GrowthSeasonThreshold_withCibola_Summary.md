# Bonytail Growth Model: Backwater + Cibola Combined Fit (2026-09-14)

**Status: alternative model, not yet adopted as the sole production model.**
Produced by `model/GIEL_GrowthSeasonThreshold_withCibola.R` (output
`data/GIEL_GrowthSeasonThreshold_withCibola.RData`), a fork of
`model/GIEL_GrowthSeasonThreshold.R`. The original backwater-only model/file are
unchanged and remain the primary reference until BK decides otherwise.

## Motivation

Both existing GIEL growth models (`GIEL_GrowthSeasonThreshold.R`,
`GIEL_VBGF_AgeLength.R`) found `Linf` collapses to ~406-418 mm, below their own
observed maximum recaptured lengths, from a lack of fish recaptured above ~400 mm.
A prior historical-only check at Cibola High Levee Pond (a different but
environmentally similar off-channel GIEL pond; see
`model/GIEL_Cibola_LinfSensitivity.R`) found real, if 20-30-year-old, growth pairs
from larger fish (TL1 up to 504 mm) giving an independently-fitted `Linf` of 662 mm
(95% CI 581-749) -- well above the backwater estimate. This analysis asks two
follow-up questions: (1) is the underlying growth rate (Fabens `K`) actually
different between the two locations, in the size range where both have data, and
(2) if not clearly different, does pooling Cibola's large-fish data into the
production model help anchor `Linf`.

## Step 1: overlap-range site-comparison test

Restricted both datasets to the TL1 range where both have real data (189-442 mm,
the intersection of backwater's 117-442 mm and Cibola's 189-504 mm): 349 backwater
pairs + 45 Cibola pairs (pre-2005), n = 394. Fit a Fabens model with shared `Linf`
and either a single pooled `K` or a site-specific `K` (Backwater vs. Cibola), each
with a site-year random effect (pond-year for backwater, Cibola-year for Cibola),
same priors as all prior GIEL VBGF work (`Linf ~ Normal(635, sd 50)`,
`logK ~ Normal(log(0.15), sd 0.5)`).

**Result: genuinely borderline, not a clean answer.**

| Metric | Value |
|---|---|
| WAIC, pooled K | 3667.1 |
| WAIC, site-specific K | 3668.9 (Delta = 1.9, favors pooled but below the conventional ~2-unit threshold for a meaningful improvement) |
| K[Backwater] | 0.896/yr (95% CI 0.651-1.151) |
| K[Cibola] | 0.630/yr (95% CI 0.390-0.878) |
| logK[Cibola] - logK[Backwater], 95% CI | -0.805 to +0.007 (barely fails to exclude zero) |
| P(K[Cibola] < K[Backwater]) | 0.974 |

WAIC (this project's standing model-selection convention) narrowly favors pooling,
but there is a ~97% one-sided posterior probability that Cibola fish grew somewhat
slower in this size range historically -- plausibly real (different era, climate,
stocking source, or density) but underpowered to confirm with only 45 Cibola pairs
in the overlap window.

**BK decision (2026-09-14):** treat as "not significant" per the project's standard
WAIC-based term-selection convention, and proceed to pool the two sites under a
single shared `K` (varying only by cohort/origin, not by site) in the combined
production model. This is a judgment call given the borderline evidence, not an
unambiguous finding, and should be revisited if a larger Cibola-overlap sample
becomes available.

## Step 2: combined production model

Pools all 361 backwater growth-increment rows (full size range, same data as
`GIEL_GrowthSeasonThreshold.R`) with all 52 Cibola pre-2005 pairs (TL1 189-504 mm)
into one dataset, n = 413, TL1 117-504 mm, TL2 up to 510 mm. Same 4-way WAIC
term-selection (origin x site-year-group) as the original model; Cibola rows all
carry `Origin = "Capture"` (matching backwater's existing "Capture" level for
wild/netting-tagged fish) and their own per-year grouping factor.

**WAIC again selects the full origin + group model** (3844.4 vs. 3849.7 group-only,
3879.8 pooled, 3880.1 origin-only) -- the same structure selected for backwater
alone, now with Cibola's year groups folded into the same random-effect
distribution. Converges cleanly (max R-hat 1.03, 0 of 22 monitored nodes > 1.1).

| Parameter | Combined model | Backwater-only model (original) |
|---|---|---|
| Linf | 435.5 mm (95% CI 421.5-451.5) | 417.4 mm |
| K[Capture-origin] | 1.00/yr | (see original summary) |
| K[Stocking-origin] | 0.77/yr | (see original summary) |
| Posterior P(Linf < max observed size) | 1.0 (max obs = 510 mm) | 1.0 (max obs = 462 mm) |

**Linf moved up by about 18 mm (417 -> 436 mm) but the identifiability problem is
NOT resolved** -- every posterior draw for `Linf` is still below the combined
dataset's own observed maximum (510 mm). This is expected, not a bug: `Linf` is
shared across all groups and both sites (per the Step 1 decision not to give Cibola
its own rate, and there was never a case for a site-specific `Linf`), and the 361
backwater rows still numerically dominate the 52 Cibola rows in fitting a single
shared curve. Cibola's own separate historical fit (662 mm `Linf`, fit without any
backwater data at all) has no direct pull on this combined model's `Linf` beyond its
52 rows' proportional contribution to the pooled likelihood.

## Derived predictions: P(TL >= 250 mm after one growth season | TL1)

| TL1 (mm) | Combined model | Backwater-only model (original) |
|---|---|---|
| 150 | 0.643 | 0.491 |
| 175 | 0.770 | 0.574 |
| 200 | 0.879 | 0.676 |
| 225 | 0.951 | 0.800 |
| 249 | 0.985 | 0.917 |

All predictions increase under the combined model, driven by both the modestly
larger `Linf` and a higher fitted Capture-origin `K` once Cibola's large,
faster-in-absolute-mm-but-slower-in-rate fish join the pool (the Fabens increment at
a given size depends on both parameters jointly, not `Linf` alone).

## Comparison figure

`model/GIEL_GrowthModel_CombinedVsOriginal.png` overlays posterior-median VBGF
curves (95% ribbons) for the backwater-only and combined models from a common
starting size (TL0 = 100 mm). The combined curve sits modestly above the
backwater-only curve throughout, with its flatten-out point (`Linf`) about 18 mm
higher -- a visible but modest shift, not a resolution of the underlying
data-sparsity problem.

## Interpretation and recommendation

Adding Cibola's historical large-fish data measurably moves `Linf` in the expected
(upward) direction, consistent with the standalone Cibola-only sensitivity check,
but the shift is small relative to the gap between the current estimate (436 mm) and
either the biologically-plausible literature range (~635 mm, this project's prior
center) or Cibola's own separately-fit value (662 mm). A well-anchored `Linf` would
require either many more large-fish recaptures (from either location) or a tighter
external prior (e.g., literature/otolith-based) than currently used. Until such data
exist, **derived quantities in the size range actually covered by data (TL1 up to
~350-400 mm) are usable from either model, but `Linf`/`K` themselves should not be
over-interpreted as literal growth-curve parameters from either fit.**

## Standing caveats (repeat wherever this combined model is cited)

1. Cibola data are 20-30 years old (pre-2005), a different pond/location, under
   different management, water source, and possibly different stocking provenance
   than the current IP2/IP5/IP6 backwater ponds.
2. The site-comparison test (Step 1) is a borderline, underpowered result (45
   Cibola pairs in the overlap range) -- pooling `K` across sites was a judgment
   call given WAIC's narrow preference, not an unambiguous finding.
3. `Linf` is still identifiability-limited even after adding Cibola (posterior
   P(Linf < max observed) = 1) -- treat 435.5 mm as a modest improvement over
   417.4 mm, not a resolved estimate.
4. The Apr-Sep growing-season window is assumed, not verified, for Cibola.
5. Every Cibola row is `Origin = "Capture"` -- there is no Cibola stocking-origin
   growth data to compare against backwater's stocked cohorts.

## Output

`data/GIEL_GrowthSeasonThreshold_withCibola.RData`: `bw_prepped`, `cibola_prepped`,
`combined_data` (413-row fitting dataset), `overlap_data_saved` (394-row
site-comparison dataset), `fit_pooledK_overlap`, `fit_siteK_overlap`, `siteTest`
(summary list), the four combined-model WAIC candidates (`comb_full`,
`comb_noorigin`, `comb_nogroup`, `comb_neither`), `comb_waic_tab`, `comb_winner`,
`winning_key_comb`, `winner_use_group_comb`,
`growth_season_threshold_probs_combined`, plus prior constants and `max_obs_comb`.
