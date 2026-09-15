# Cibola High Levee Pond Linf Sensitivity Check (Bonytail)

**Status: a historical, standalone sensitivity check only.** This analysis is separate
from, and does not modify or feed into, the backwater bonytail growth models
(`model/GIEL_GrowthSeasonThreshold.R`, `model/GIEL_VBGF_AgeLength.R`). It was run
2026-09-14 at BK's request after determining that Cibola High Levee Pond (LID 10222,
Cibola NWR, Reach 4) has no *recent* large-fish growth-recapture pairs (see
conversation log): the only physical GIEL recaptures there are from 1993-2005 (large
fish) plus 3 small-fish pairs from 2015-2016, and no GIEL were physically recaptured
at Cibola at all after 2017.

## Purpose

Both current backwater growth models independently found that `Linf` collapses to a
value (406-418 mm) below the observed maximum post-growth length in their own data,
and diagnosed this as a Fabens-estimator bias caused by having almost no recaptured
fish above ~400 mm TL. Cibola's pre-2005 records include real recaptures of much
larger fish (TL1 up to 504 mm, TL2 up to 510 mm) and offer an independent, if old and
site-different, check on whether that low `Linf` is plausible.

## Data

- Source: `data/NFWGAnalysis.RData` (`NFWGAnalysis`), `location_id == 10222`.
- 71 consecutive-measured-record growth pairs at Cibola overall (Date2 spanning
  1994-04-07 to 2016-11-15); restricted here to the **pre-2005 subset** (Date2 <
  2005-01-01), 67 of 71 pairs, for internal consistency (one pond, one management
  era) with the case for including Cibola at all (large-fish coverage).
- Growing-season (Apr 1-Sep 30) exposure computed with the same calendar-overlap
  logic as the backwater models -- **applied here as an unverified assumption**,
  since Cibola is a different, more southerly location.
- 13 of 67 pairs are negative-growth (TL2 < TL1, almost certainly measurement error,
  same pattern already documented for the backwater data) and are dropped; 2 more are
  winter-only (`gs_days = 0`) and dropped. **52 rows used for fitting**, TL1 range
  189-504 mm, TL2 range up to 510 mm.

## Model

Single-population Fabens/VBGF increment model (NIMBLE), no origin or pond-year
grouping (one pond, one era, insufficient fish for grouping):

```
increment = TL2 - TL1 ~ Normal(mu, sigma)
mu[i] = (Linf - TL1[i]) * (1 - exp(-K * t_years[i]))
```

Priors matched to the backwater models for direct comparability: `Linf ~ Normal(635,
sd 50)`, `logK ~ Normal(log(0.15), sd 0.5)`.

## Results

**Primary fit** (`Linf` prior sd = 50, matching the backwater models): converges
cleanly, R-hat <= 1.001 for all parameters (unlike the backwater fits, which show a
curved `Linf`-`logK` ridge). `Linf = 662 mm (95% CI 581-749)`, `K = 0.316/yr (95% CI
0.237-0.420)`. **Posterior P(Linf < 510 mm, the dataset's own observed max) = 0** --
every posterior draw supports an asymptote above the observed data, the opposite
identifiability failure from the backwater fits (which had zero posterior support for
`Linf` above their own observed max).

**Prior-sensitivity fit** (`Linf` prior sd = 300, otherwise identical): `Linf` moves
further, to 834 mm (sd 155) -- so the 662 mm point estimate is itself meaningfully
prior-influenced and should not be read as precise. What is robust across both prior
widths: `Linf` sits well above 510 mm and well above the backwater models' ~406-418 mm
in every posterior draw at both prior widths.

## Interpretation

This supports the qualitative conclusion that the backwater growth models' low `Linf`
is likely an artifact of insufficient large-fish coverage in their own data, not a
true reflection of species/system asymptotic size -- Cibola's real (if old)
large-fish recaptures are consistent with a substantially larger asymptote. It does
**not** provide a precise alternative `Linf` value (the two priors tested give 662 mm
and 834 mm), and it does not resolve the backwater models' own identifiability
problem, since it is a different pond and a different (much older) era of data.

## Standing caveats (repeat in any report or document citing this check)

1. Data are 20-30 years old (1993-2005), a different pond/location, under different
   management, water source, and possibly different stocking provenance than the
   current IP2/IP5/IP6 backwater ponds.
2. The Apr-Sep growing-season window is assumed for Cibola, not verified.
3. Single-pond, single-era, ungrouped model -- no origin or cohort effects estimated.
4. `Linf`/`K` here are themselves prior-sensitive (662 mm at prior sd 50 vs. 834 mm at
   prior sd 300) -- treat as a qualitative direction-of-bias check, not a replacement
   point estimate.
5. Does not modify `GIEL_GrowthSeasonThreshold.R`, `GIEL_VBGF_AgeLength.R`, or any of
   their saved `.RData`/summary files.

## Comparison figure

`model/GIEL_Cibola_LinfComparison.png` overlays posterior-median VBGF growth curves
(with 95% posterior ribbons) for the backwater model (`GIEL_GrowthSeasonThreshold.R`
winning fit, `origin + pondyear`, marginalized over the pond-year random effect and
using the "Capture"-origin `logK`) and the Cibola historical primary fit (`Linf` prior
sd = 50), both projected from a common reference size (TL0 = 100 mm at t = 0) on the
shared growing-season-year time axis both models use. The two curves track closely
through about TL 350-400 mm (the size range both datasets cover reasonably well) and
then diverge sharply: the backwater curve flattens at its ~417 mm `Linf` while the
Cibola curve continues rising toward ~661 mm, consistent with the backwater model
having no data to constrain growth beyond ~400-460 mm. Jittered rug tick marks along
the left edge show each dataset's observed `TL1` (starting size) values -- backwater
(n = 361, `gs_fit_data$TL1`, blue) spans roughly 117-442 mm and is densely populated
below ~350 mm; Cibola (n = 52, `cibola_fit_data$TL1`, orange) spans 189-504 mm and is
the only one of the two with meaningful support above ~400 mm, visually confirming
that the divergence in the fitted curves begins exactly where the backwater data run
out and only the Cibola data continue. Not saved to the pipeline's `output/` folder
(this is a one-off diagnostic, not a report figure); regenerate by loading
`data/GIEL_GrowthSeasonThreshold.RData` and `data/GIEL_Cibola_LinfSensitivity.RData`
together (see script comments / conversation log 2026-09-14 for the exact plotting code).

## Output

`data/GIEL_Cibola_LinfSensitivity.RData`: `cibola_pairs` (all 71 pairs, any date),
`cibola_hist` (pre-2005, pre-negative/winter filtering), `cibola_fit_data` (52 rows
used for fitting), `CibolaNegativeGrowthPairs` (13 flagged rows), primary fit
(`cibola_hist_fit`, `rhat_cib`, `samps_cib`, `cib_summ`) and wide-prior sensitivity fit
(`cibola_hist_fit_wide`, `rhat_cib_wide`, `samps_cib_wide`, `cib_summ_wide`), plus
prior constants and `max_obs`.
