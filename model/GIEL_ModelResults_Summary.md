# Bonytail (GIEL) Population Model — Results Summary

**Purpose of this file:** self-contained briefing on the Bonytail (Bonytail Chub, species code
GIEL) population monitoring model for the BackwaterGenetics project, for an AI assistant that
needs context without re-reading the full project memory file. Written 2026-09-10.

## 1. Background

Bonytail (GIEL) are PIT-tagged and monitored in three artificial ponds at the Imperial Ponds
Conservation Area (IPCA), part of the USFWS Imperial National Wildlife Refuge near Yuma, AZ:

| Pond | LID | Surface area |
|---|---|---|
| IPCA (Pond 2) — "IP2" | 1044 | 12.5 acres |
| IPCA (Pond 5) — "IP5" | 1047 | 21.3 acres |
| IPCA (Pond 6) — "IP6" | 8.8 acres (3.6 ha) |

(The companion species, Razorback Sucker/XYTE, occupies Yuma Cove Backwater and IPCA Ponds 1, 3,
4, and has its own three-stage density-dependent model — not covered here.)

Fish are PIT-tagged at stocking or at netting events, and re-detected via passive PIT antennas
plus periodic netting. The data span FY2017 (first stocking after a 2014–2015 pond renovation)
through FY2026 (partial).

## 2. Model structure decision: 2-stage, not 3-stage

A parallel three-stage density-dependent model was built for Razorback Sucker at Yuma Cove
(splitting juveniles into small/large substages, with recruit density driving growth). Before
extending that architecture to Bonytail, a spot-check (2026-09-10) tested whether it was
warranted:

- **Size vs. survival supports only one break, at 250 mm** (the existing stage threshold used in
  the standard 2-stage model), not two. Netting-tagged juveniles re-sighted ≥365 days later:
  250–349 mm and ≥350 mm both survive much better than <250 mm, but do **not** differ from each
  other (log-odds difference 0.25 ± 0.20, z = 1.2, p = 0.23). No evidence for a second size break.
- **No smooth density-survival gradient to fit.** Known-alive adult density swings by an order of
  magnitude year to year, but as discrete crash/rebound events (see §3), not a continuum.
- **Recruit density too sparse for a growth-density submodel.** Untagged juveniles at fall netting
  are 0 in most pond-seasons; one extreme outlier (IP2 FY2024, 24.1/acre) is an order of magnitude
  above everything else.

**Decision:** reuse the existing 2-stage joint robust-design model (shared with Razorback Sucker,
stage split at 250 mm for GIEL) rather than building a Bonytail-specific three-stage variant.

## 3. Confirmed die-off events

Marsh et al. (2024) documents two Bonytail die-offs attributed to well pump failure / human error.
The PIT-scan known-alive time series (continuously and heavily scanned throughout both windows, so
these are real deaths, not detection gaps) pins down tighter dates than the source's calendar-year
citations:

| Pond | Window | Known-alive trajectory |
|---|---|---|
| IP5 | 2020-01 to 2020-09 (FY2020) | 126 → 75 → 59 → 56 → 55 → 55 → 52 → 42 → 0 (from a post-tagging peak of 145 in Dec 2019); stayed at 0–1 fish through FY2021 |
| IP2 | 2022-01 to 2022-07 (FY2022) | 121 → 85 → 76 → 61 → 55 → 54 → 39; scanning dropped to 0 h in Aug–Sep 2022 right as the population hit 0, so the window is truncated at the last well-scanned month rather than extended further |

These windows were encoded as a per-pond, per-season fixed "event" indicator (`Evt`/`nEvt`) and a
survival-reducing offset parameter `dE[pond]` was added to the existing joint robust-design NIMBLE
model (ported from a mechanism first built for the XYTE three-stage model). Recorded in
`data/BWDieOffEvents.csv`.

## 4. GIEL-only model fit

**What was fit:** the existing joint 2-stage robust-design model (`model/BWRobustDesign_defs.R`),
restricted to the three GIEL ponds only (`RUN_PONDS <- c(3, 6, 7)` → model indices IP2=1, IP5=2,
IP6=3). Monthly time step, Oct–Apr primary detection occasions, states Juvenile/Adult/Dead (split
at 250 mm), 3 MCMC chains, 30,000 iterations / 10,000 burn-in / thin 10, ~60 minutes wall time.
Output: `data/BWRobustDesign_MCMC_giel.RData`.

**Convergence:** clean. Of 673 monitored parameters, 42 have R-hat > 1.1, but all of these are
either a known floating-point artifact on the known-alive-constrained `N_tag`/`Tal` nodes (not a
real mixing problem) or recruitment nuisance parameters (`RecJ`, `r`, `U`) in the last 2–3 seasons.
The parameters that matter for inference — `lphiA` (adult survival), `dE` (die-off offsets),
detection, and hyperparameters — all converge well (R-hat ≤ 1.1).

**Key results:**

- **Die-off offsets:** `dE[IP2] = -1.74 (95% CI -3.15, -0.55)`, `dE[IP5] = -1.45 (-2.07, -0.87)`,
  both clearly negative on the logit-survival scale. `dE[IP6]` (no recorded event at this pond) is
  `-2.42 (-6.76, -0.09)` — essentially prior-driven, as expected.
- **Event-adjusted realized annual adult survival vs. direct empirical floors** (fraction of adults
  known alive at one October re-seen after the next October):
  - IP5 FY2020: model 0.041 (0.017–0.077) vs. floor 0/39 (0%) — **excellent match**, the
    near-whole-season (9 of 12 months) event window let the model recover the true near-total loss.
  - IP2 FY2022: model 0.089 (0.044–0.147) vs. floor 0/14 (0%) — **still biased high**. The
    event window here covers only 7 of 12 months; hierarchical shrinkage pulls the fitted rate back
    toward the cross-pond/season average, an under-correction pattern already seen for Razorback
    Sucker's analogous IP1 FY2022 crash in the companion three-stage model. The common factor
    appears to be *partial-season* event windows, not species or pond differences.
  - **IP2 is also chronically low-survival in seasons with no recorded event** (empirical floor
    ~0.05–0.09 in FY2018, FY2019, FY2021, FY2025) — the single documented 2022 crash does not
    explain this pattern; IP2 may have additional unflagged mortality episodes worth investigating
    separately from Marsh et al.'s two-event record.
  - Pond IP6 (no die-off) shows season-to-season adult survival tracking its empirical floor almost
    exactly across every season (e.g., FY2018 model 0.50 vs. floor 0.51; FY2025 model 0.56 vs. floor
    0.57) — the base random-effects structure fits well when nothing anomalous is happening.
- **Pond-level vs. season-level variation:** within-pond season-to-season variation in adult
  survival (SD ≈ 0.69 on the logit scale) is larger than the variation across pond means (SD ≈
  0.45). This favors the current model's structure — pond × season random draws around a shared
  species-level hyper-mean (`mu_phi = 2.71`, `sigma_phi = 0.85` on the logit scale) — over a design
  with fixed pond intercepts and a shared covariate slope.
- **Stocked vs. wild-hatched survival split — not currently testable.** Sample sizes for
  origin = "tagged at netting" (our proxy for recruited/wild fish) are substantial at all three
  ponds (IP2: 857 vs. 599 stocked; IP5: 159 vs. 600 stocked; IP6: 423 vs. 300 stocked), so a
  sample-size limitation is not the blocker. The model architecture is: origin only affects a
  temporary post-release survival offset (`b_post`, active for a few months after tagging), with no
  origin term on steady-state established-adult survival (`lphiA`) at all. Answering this question
  would require a model extension, not just re-reading this fit.
- **Data sufficiency check:** of 27 pond-season cells with any known-alive fish at their October
  census, 2 are the flagged die-off events, leaving 25 non-event cells — of which only 19 have ≥15
  known-alive fish (a workable sample for a floor comparison). Even restricted to these 19 "clean"
  cells, empirical annual survival swings from ~0.05 to ~0.92 with no apparent smooth pattern. This
  is a thin, noisy base for fitting even a single covariate slope across the three ponds.

## 5. Comparison against the externally-proposed BONY integrated population model (IPM)

A separate document (Sound Science LLC, "BONY Proposed IPM," draft 10 July 2026) proposes a
richer model: 4 life stages with early wild sub-stages, environmental covariates (dissolved
oxygen, temperature, water level) on survival and reproduction, release/habitat covariates on
post-stocking survival, and origin-based survival splits by pond. See
`model/BONY_ProposedIPM_Critique.md` for the full structural critique. Headline conclusions,
now updated with the fit results above:

- The proposed 3–4 wild sub-stage structure is not supported by our size-vs-survival data (§2);
  our simpler 250-mm, 2-stage split matches what re-sighting data show.
- No evidence supports a density-dependent survival or growth covariate for GIEL (§2); the model's
  biggest survival swings are step-function die-off events, not a gradient an environmental
  covariate would smoothly track.
- The two known die-off events are well-represented as explicit, dated event terms (our `dE`
  mechanism) rather than requiring a continuous environmental covariate — and the proposed
  document's own text attributes both known crashes to operational pump failures, not
  environmental drivers.
- A companion analysis of a 2013–2016 Lake Havasu open-water bonytail telemetry study
  (`documents/C64_2016_LakeHavasu_Summary.md`) found that even purpose-built, telemetry-grade
  survival data could not identify a bare (covariate-free) survival curve for stocked bonytail in
  open water. Our own pond data are sparser for the contrasts the proposed IPM's covariates would
  need (§4's "19 informative cells, swinging widely" finding) — reinforcing skepticism that the
  proposed model's fuller covariate stack (environmental, origin-by-pond, release covariates) is
  identifiable from currently available GIEL data.
- Reproduction (F) via genetic parentage is not currently implemented in this project's pipeline
  for any GIEL pond (parentage data currently only exist for XYTE at Yuma Cove and IP Pond 1); a
  University of New Mexico lab reportedly holds some bonytail parentage data independent of our
  Dowling-lab panel, but its coverage of IP2/IP5/IP6 specifically is unconfirmed.

## 6. Key files

| File | Contents |
|---|---|
| `data/BWDieOffEvents.csv` | Hand-entered die-off event windows (all ponds/species), used to build the `Evt`/`nEvt`/`dE` mechanism |
| `data/BWRobustDesign_data.RData` | NIMBLE inputs for the joint 2-stage robust-design model (all 7 ponds; rebuild before any refit) |
| `data/BWRobustDesign_MCMC_giel.RData` | Posterior samples/summary from the GIEL-only (IP2/IP5/IP6) fit described here. Only 9 objects saved (`rd_code`, `rd_Fish`, `rd_Hist`, `rd_Netting`, `rd_ponds`, `rd_samples`, `rd_settings`, `rd_summary`, `rhat_all`) — load into a fresh environment to avoid confusion with stale objects of the same name from other model-building sessions. `rd_summary` covers only 26 core hyperparameters; full posteriors (673 nodes) are in `rd_samples`. |
| `model/BWRobustDesign_traceplots_giel.pdf` | Trace plots for the GIEL-only fit |
| `model/BONY_ProposedIPM_Critique.md` | Full structural critique of the proposed IPM, including a §8 with detailed evidence for each open question, underlying this summary |
| `documents/BONY_Proposed_IPM_Summary.md` | Factual (non-critical) summary of the proposed IPM document itself |
| `documents/C64_2016_LakeHavasu_Summary.md` | Summary of the Lake Havasu open-water bonytail telemetry study, cited as precedent against complex covariate models for bonytail survival |

## 7. Open follow-ups

- Investigate whether IP2 has additional unflagged mortality episodes beyond the documented 2022
  event (its non-event-season survival is chronically low).
- Consider a pond-and-season-specific `dE` (rather than one shared scalar per pond) if a closer
  match to partial-season crash floors (as at IP2 FY2022) is wanted.
- To test stocked-vs-wild survival, add an origin effect to established-adult survival (`lphiA`)
  and refit — not answerable from the current model.
- Given the thin, noisy non-event data (19 informative pond-seasons), any future environmental or
  origin covariate should be justified against this data-sufficiency finding before being added.
