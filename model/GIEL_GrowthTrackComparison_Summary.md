# Bonytail Growth-Model Track Comparison (2026-09-14)

**Script:** `model/GIEL_GrowthTrackComparison.R` | **Data:** `data/GIEL_GrowthTrackComparison.RData`

Step 2 of the growth-model sensitivity plan feeding the joint robust-design model's size-tied
maturation hazard (see `AGENTS.md` "GIEL Size-Tied Maturation Hazard" and the plan file
`.posit/assistant/plans/2026-09-14-1444-giel-growth-model-tracks-feeding-robust-design-
maturation-hazard.md`). Compares the three growth-model tracks selected for this pass (Track
A2, age-based + Cibola, is out of scope):

| Track | Model | Data | Source |
|---|---|---|---|
| **F1** | Fabens increment | IPCA-only | `GIEL_GrowthSeasonThreshold.R` |
| **F2** | Fabens increment | IPCA + Cibola | `GIEL_GrowthSeasonThreshold_withCibola.R` |
| **A1** | Age-based VBGF (shared-K, release-anchored) | IPCA-only | `GIEL_VBGF_AgeLength.R` |

## P(TL >= 250mm after one Apr-Sep growth season | TL1)

| TL1 (mm) | F1 | F2 | A1 | Spread (max - min) |
|---|---|---|---|---|
| 150 | 0.491 | 0.643 | 0.268 | 0.375 |
| 175 | 0.574 | 0.770 | 0.521 | 0.249 |
| 200 | 0.676 | 0.879 | 0.788 | 0.203 |
| 225 | 0.800 | 0.951 | 0.951 | 0.151 |
| 249 | 0.917 | 0.985 | 0.995 | 0.078 |

**The three tracks disagree substantially at small starting sizes and converge as TL1
approaches 250mm** (mechanically expected -- less growth is needed to cross the threshold).
At TL1 = 150mm the tracks span 0.27 (A1) to 0.64 (F2), a 0.375 spread; by TL1 = 249mm they
agree within 0.08. **A1 (age-based) gives the lowest probability at small sizes** (0.268 at
150mm, notably below both Fabens tracks) but is essentially indistinguishable from F2 by
200-249mm; **F2 (Fabens + Cibola) gives the highest probability at every TL1**, consistent
with its larger fitted `Linf` (435.5mm vs. F1's 417.4mm) and higher Capture-origin `K`.

This is the divergence flagged ad hoc during the 2026-09-14 planning conversation (there
estimated informally as "0.27 age-model vs. 0.49 increment-model" at TL1=150mm, matching F1
here) -- now computed for all three tracks on identical grids and saved as a reusable table.

## Posterior growth curves (common start TL0 = 100mm)

`model/GIEL_GrowthTrackComparison.png` overlays posterior-median VBGF curves (95% ribbons)
for all three tracks from a common TL0 = 100mm, with dashed lines at each track's posterior
median `Linf`. All three curves track closely through roughly 200-300mm, then diverge as each
flattens toward its own fitted `Linf` (F1 417mm, F2 436mm, A1 406mm) -- none of these values
should be read as a literal biological asymptote (see caveat below); the divergence in the
curves past ~300mm reflects each track's own data-sparsity-driven compromise, not a
resolved difference in true growth rate.

`model/GIEL_GrowthTrackComparison_ProbCurves.png` shows the three P(TL>=250mm) curves from
the table above as a single overlay plot.

## Standing caveat carried from all three source models

**All three tracks' `Linf` sits at or below that track's own observed maximum TL** -- this is
the same Fabens-type identifiability failure (too few large-fish recaptures to pin down the
true asymptote) documented independently for F1/F2 (`GIEL_GrowthSeasonThreshold_Summary.md`
§0b, `GIEL_GrowthSeasonThreshold_withCibola_Summary.md`) and A1
(`GIEL_VBGF_AgeLength_Summary.md`). None of the three `Linf`/`K` posteriors should be
over-interpreted as literal growth-curve parameters; all three P(TL>=250mm) curves are best
read as track-specific, data-driven compromises valid for interpolation in the well-populated
size range (roughly 100-350mm), not as competing "true" biological estimates. The spread
between tracks reported above is therefore a **lower bound on total growth-model
uncertainty** -- it reflects disagreement between three specific model structures/datasets,
not full parameter uncertainty within any one of them (though each track's own probability
table already incorporates its own posterior uncertainty in the underlying simulation).

## Additional caveat specific to Track A1

Track A1's residual `sigma` (17.0mm) was estimated from absolute age-length residuals
spanning ages up to ~10 years, not from short (182-day) one-season increments directly;
applying it to a one-season projection likely overstates uncertainty in A1's curve relative
to a residual sd fit on short-interval increments alone (see
`GIEL_VBGF_AgeLength_Summary.md`).

## Purpose / next step

This comparison quantifies how much choice of growth-model track could move the maturation
hazard used downstream in the joint robust-design model (Steps 3-6 of the sensitivity plan):
build a `psi_lookup` hazard table from each track's posterior (Step 3), wire it into
`model/BWRobustDesign_data.R`/`BWRobustDesign_defs.R` (Step 4), and refit the GIEL-only
robust-design model once per track (Step 5) to see whether this growth-model disagreement
propagates into materially different survival/maturation/recruitment conclusions, or is
absorbed without much effect (Step 6).

## Output

`data/GIEL_GrowthTrackComparison.RData`: `GrowthTrackComparisonTable`, `GrowthTrackSpread`,
`curve_f1`, `curve_f2`, `curve_a1`, `curves_all`, `linf_lines`.
