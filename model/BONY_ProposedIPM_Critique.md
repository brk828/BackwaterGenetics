# Proposed BONY IPM (Sound Science LLC, 9-2-26) — Structure Outline and Critique

**Purpose.** This document outlines the structure of the externally-proposed bonytail integrated
population model (`documents/BONY Proposed IPM - Model Description 9-2-26.docx`; see
`documents/BONY_Proposed_IPM_Summary.md` for a plain factual summary) and records a critique to
revisit once the GIEL-only joint robust-design fit (job `C1FED42C` launched 2026-09-10, see
`AGENTS.md`) returns results. The goal is not to re-derive the proposed model but to note (a) where
its structure assumes things our data don't support or that belong to razorback sucker (XYTE) rather
than bonytail (GIEL), and (b) where our fitted models already speak to a structural choice the
proposal has to make on paper.

---

## 1. Structure outline (as proposed)

### 1.1 Scope split
- **Off-channel (OCH) state-space IPM** — CHLP, IPCA (Ponds 2/5/6 only), Yuma Cove, ephemeral
  ponds "where data support." Full abundance + survival + recruitment model.
- **Open-water CJS** — mainstem releases, survival-only, weekly time-since-release scale, no
  persistent population assumed.

Our project's scope is OCH-only (YCB + all 6 IPCA ponds); we have no open-water release data, so
the CJS component is entirely out of scope for comparison here.

### 1.2 Stages and rates
Four annual stages: newly-stocked adult → established adult, with early wild stages (S1–2, S2–4)
and reproduction F closing the loop only in isolated ponds. Rates: φ³ (post-stocking, ~year 1),
φ⁴ (established-adult annual), F (reproduction), plus detection p and observation error σ.

### 1.3 Graphical-model / joint-likelihood structure
Three sub-models share parameters:
- Reproduction (parentage genetics) → F
- OCH abundance (state-space, PIT-scan + net counts) → N³, N⁴, φ³, φ⁴, σ
- Encounter (CJS on individual PIT histories) → φ³, φ⁴, p

φ³/φ⁴ are shared across the abundance and encounter likelihoods (the actual "integration"). Total
abundance is recovered by scaling tagged abundance (capture-recapture) with the netted untagged
fraction, assuming equal catchability of tagged/untagged fish; a closed-population fallback is
proposed if that assumption fails. Close-kin mark-recapture from the same genetic panel is floated
as a future extension, not implemented.

### 1.4 Site coverage / inclusion rule
Coverage (stocking, tagging, detections, netting effort, antenna deployment), not raw fish counts,
decides which ponds are modeled. Only IP2, IP5, IP6, and CHLP qualify. IP1, IP3, and Yuma Cove are
explicitly dropped for the BONY analysis (recapture ≤ 5.3%); pre-2016 Pond 2 and post-2016 CHLP
records are also censored.

### 1.5 Reproduction (F) estimation
F estimated from parentage assignment of genotyped larvae/YOY against a near-complete adult
genotype panel — participation, reproductive skew, and relative per-adult output are estimated
directly from parentage; absolute F magnitude is pinned by trammel-net recruit counts. Supporting
evidence cited (participation 89–99%, collapsing to 27–33% in a poor year) is from Osborne et al.
(2018), a different ("Nine Mile") pond system, not CHLP/IPCA.

### 1.6 Covariate structure
- `logit φ⁴`, `log F` ~ environmental covariates (dissolved oxygen, temperature, water level) +
  pond intercept + year random effect (shared slopes across ponds, pond-specific intercepts).
- `logit φ³` ~ release/habitat covariates (reach, size, method, predator load, timing, site).
- Origin (stocked vs. wild-hatched) is carried as a survival split, allowed to vary by pond.
- Sex is explicitly **not** a survival covariate (unreliable field assignment; confounded with
  release size), despite citing a sex-biased-survival result from elsewhere in the literature.
- Site-level (whole-pond) covariates are flagged as unidentifiable with only a handful of ponds
  and are reported descriptively, not fitted.
- Year random effects absorb the 2020/2022 IPCA die-offs as one-off operational shocks (well pump
  failure), explicitly *not* attributed to the environmental covariates.

### 1.7 Modeling progression
Staged build: open-water CJS → OCH CJS with pond fixed effects → covariate slopes on survival →
full state-space abundance + recruitment integration, applied only as far as each site's data
support.

---

## 2. Critique: relevance to bonytail vs. razorback sucker

Our 2026-09-10 GIEL structure assessment (`AGENTS.md`, "Bonytail (GIEL) Model Structure
Assessment") is the direct empirical counterpoint to several structural choices below.

- **Three vs. two size stages.** The proposed IPM's stage diagram (S1–2, S2–4, S3–4, S4–4) reads
  like the XYTE three-stage architecture we built and then abandoned for cross-pond use (v12,
  2026-09-10) — multiple juvenile sub-stages driven by size at maturation. Our GIEL spot-check
  found the data support only **one break at 250 mm** (existing `StageTL`); 250–349 mm and ≥350 mm
  netting-tagged juveniles do not differ in re-sighting (diff 0.25 ± 0.20, z = 1.2, p = 0.23). A
  3-stage or continuous-growth structure for GIEL would not be identifiable from re-sighting data,
  independent of any covariate model layered on top.
- **Density-dependent growth/survival covariate.** The XYTE three-stage model's central finding
  (`a1 = -1.9`, strong density effect on recruit size) does not transfer: GIEL recruit density is
  near-zero in most pond-seasons with one extreme outlier (IP2 FY2024, 24.1/acre), and known-alive
  adult density swings are step-function crash/rebound events, not a smooth gradient. A
  `cA · density` or growth-density term the proposed IPM's environmental-covariate slot might
  invite (dissolved oxygen and temperature as proxies for stocking density effects) has no leverage
  independent of the already-known die-off events for GIEL — this is a RASU-shaped intuition
  (YCB's dense juvenile cohorts) that likely does not carry over.
- **Sex-biased survival citation.** The proposed document cites "pond experiments have documented
  lower survival for stocked males than females (Marsh et al. 2024)" as motivation for *excluding*
  sex from the BONY model but flagging the concern. That finding is our own YCB **razorback
  sucker** result (24M vs. 80F surviving the 2013 stocking cohort;
  `documents/2025Report_YCB_PopulationMonitoring_Summary.md`) — Marsh et al. (2024) is a
  multi-species synthesis and the cited sex-survival contrast is not established for bonytail in
  our data. Worth flagging in review: is there an independent BONY-specific source for this claim,
  or is it XYTE evidence imported by association? Either way our GIEL model has no sex-survival
  signal to report against it.
- **Established-adult survival split by origin (stocked vs. wild-hatched), allowed to vary by
  pond.** This is a reasonable ask in principle but is far better supported for XYTE at YCB, where
  the majority of the current population is naturally recruited and origin contrasts are
  data-rich. For GIEL, recruit density is 0 in most pond-seasons (per the same assessment above),
  so wild-hatched adults will be rare-to-absent in most pond-years; a per-pond stocked-vs-wild
  contrast may only be estimable for whichever pond(s) had a productive spawning year, if any.
- **Site list includes IP1, IP3, and Yuma Cove for a BONY model.** Per our `Study Locations` table,
  Ponds 1, 3, 4 are XYTE-stocked and Ponds 2, 5, 6 are GIEL-stocked; Yuma Cove is XYTE-only. It's
  unclear why a bonytail-specific document discusses dropping IP1/IP3/Yuma Cove for insufficient
  BONY recapture — either there was minor incidental historical BONY presence at those sites not
  reflected in our data, or this passage was carried over from a joint/RASU-parallel document
  without being re-scoped to BONY. Worth confirming with the authors before treating the "5.3%
  recapture" figure as informative for our three GIEL ponds.
- **CHLP inclusion.** Cibola High Levee Pond is outside `BackwaterGenetics`'s scope entirely (not
  in `Study Locations`); any comparison of our GIEL fit against the proposed model can only speak
  to IP2/5/6, not the CHLP component of their analysis.

## 3. Critique: data we don't have on hand

- **Larvae/YOY genetic parentage for GIEL — correction (BK, 2026-09-10).** The original draft of
  this critique stated that parentage data exist only for YCB and IP1 XYTE per `Analysis Rules and
  Caveats`, and treated the proposed F sub-model as a capability we entirely lack for bonytail.
  That's too strong: **a University of New Mexico genetics lab holds some bonytail parentage data**
  independent of the Dowling-lab panel this project's `Analysis Rules and Caveats` describes, so
  the proposal is not out of line to assume some parentage-based F estimation is possible in
  principle. What's still an open question, not a settled gap: whether that UNM dataset (a) covers
  IP2/5/6 specifically (vs. other refuge ponds/CHLP), (b) spans enough years and enough of each
  pond's adult pool to satisfy the "near-complete genotyped adult pool" requirement the proposal's
  F-estimation section leans on, and (c) is in a form (individual larvae/YOY assignments, not just
  summary statistics) that a state-space IPM could actually consume. `BackwaterGenetics`'s own
  pipeline still has no F term or GIEL parentage join, so *for this project's model* the capability
  gap is real; whether it's a project-scope gap or a true data gap for the LCR MSCP program as a
  whole depends on the UNM dataset's coverage, which should be checked before assuming F is
  estimable for all three ponds.
- **Environmental covariates (dissolved oxygen, temperature, water level).** We do not currently
  collect or store any of these for IPCA ponds. Task asks us to *assume* they become available and
  evaluate their likely relevance:
  - The proposed document's own text undercuts their importance for the events that matter most:
    it explicitly attributes the 2020 (IP5) and 2022 (IP2) die-offs to **well pump failure**, i.e.
    an operational/mechanical cause, and states these are "not environmental" and are absorbed by
    a year random effect rather than a covariate. If the known catastrophic mortality events are
    by the proposal's own account not environmentally driven, DO/temperature/water-level
    covariates would only be informative for the *remaining, presumably smaller* variation in
    survival between crashes — and our GIEL assessment above found no evidence of a smooth
    density/environment gradient in survival outside the documented crash events (adult density
    varies by crash/rebound, not a continuum). It's plausible environmental covariates would
    explain little incremental variance beyond what our simple pond × season random effect (or the
    proposed model's year random effect) already absorbs.
  - A concrete test once real DO/temp/water-level series exist: check whether they predict
    survival *within* non-crash seasons, or whether they mostly just co-vary with (i.e., serve as
    an indirect flag for) the known pump-failure windows already captured by an explicit event
    indicator (our `dE`/`Evt` mechanism, ported to the joint 2-stage model 2026-09-10). If the
    latter, an explicit event term is more transparent and identifiable with only a few ponds/years
    than a continuous covariate whose slope is pooled across ponds by necessity (per the proposal's
    own few-ponds caveinta about τ being weakly identified).
  - Practical corollary: if/when this covariate data is pursued, prioritize water level and DO
    around the two known crash windows (IP5 Jan–Sep 2020, IP2 Jan–Jul 2022) specifically to
    confirm or refute the "operational, not environmental" attribution, rather than assuming it.

## 4. What our fitted models already say about the right structure

- **Two-stage (not three-stage) size structure for GIEL** — supported directly (see §2). Our
  `model/BWRobustDesign_*` joint model already uses `StageTL = 250` mm for GIEL; no change
  indicated.
- **No density-dependent survival/growth covariate for GIEL** — supported directly (see §2); we
  are not planning to add one.
- **Discrete die-off events, not a smooth environmental gradient, drive the largest survival
  swings** — our `dE`/`Evt`/`nEvt` mechanism (ported from the XYTE three-stage model to the joint
  2-stage model on 2026-09-10, `data/BWDieOffEvents.csv`) already implements the proposal's own
  "one-off shock absorbed by a year effect" idea, but as an explicit, dated event term rather than
  a generic year random effect — arguably a closer match to the proposal's own stated causal
  account (operational pump failures) than a continuous environmental-covariate slope would be.
- **Origin-based (stocked vs. wild) survival split** may be uninformative for most GIEL pond-years
  given how sparse recruitment is (see §2) — this is an empirical question the GIEL fit itself can
  partly answer: once results are in, check whether any pond-season has enough non-stocked
  (`origin` code) adults to even attempt a stocked-vs-wild contrast.
- **Reproduction (F) is not currently estimable for GIEL at all**, by parentage or otherwise — this
  is the sharpest capability gap, not a structural disagreement, and should be stated plainly in
  any comparison rather than treated as a modeling choice.

## 5. Is there enough data to support this level of complexity at all?

Separate from whether individual structural choices are right for bonytail (§2) is a more basic
question: even granting every assumption, do bonytail stocking/release/capture/scanning data —
either program-wide or within our three GIEL ponds — carry enough information to identify the full
stack of stages, covariates, and sub-models the proposal lays out? Two lines of evidence bear on
this directly.

- **The open-water φ³ + covariate submodel has already been tried on far richer open-water data
  and failed.** `documents/C64_2016_LakeHavasu_Summary.md` (Humphrey, Kesner, and Marsh 2016)
  reports a dedicated 3-year acoustic/radio telemetry study (85 tagged fish, 6 iterations) plus PIT
  scanning at four Reach 3 sites, purpose-built to inform exactly this kind of post-stocking
  survival/covariate analysis. The authors explicitly attempted to fit continuous PIT-scanning
  data to CJS and Barker mark-recapture models — a **survival-only** model, with no covariates
  yet layered on — and reported that a "complex and biologically realistic model with time varying
  survival and re-sight rates resulted in highly variable parameter estimates for survival and
  re-sight rates (nearly 0 to 1)" even among structural variants with acceptable likelihoods. If a
  telemetry-grade dataset purpose-designed for this question could not identify a bare CJS/Barker
  survival curve, it's a real question whether the routine stocking/PIT-scan records the proposed
  IPM would draw on for φ³ (lower resolution, not telemetry-grade, not purpose-designed) can
  support φ³ **plus** a multi-covariate linear predictor (reach, size, method, predator load,
  timing, site) on top of it. The proposal should show this submodel converges and recovers
  sensible, non-degenerate covariate slopes before treating its outputs as informative — the
  Havasu precedent suggests that's not guaranteed.
- **Our own GIEL data is sparser than Havasu's, not richer, for the contrasts the model needs.**
  The 2026-09-10 GIEL structure assessment found recruit density near-zero in most pond-seasons
  (one order-of-magnitude outlier aside) and adult-density variation dominated by two discrete
  crash events rather than a gradient (§2). That's the same "not enough contrast to fit a slope"
  problem Havasu's telemetry data ran into for open-water survival, showing up again in the pond
  data for growth/density covariates. Before adding *any* covariate (environmental or otherwise)
  to the GIEL model, check whether the pond-season count (30 for our 4-XYTE-pond build; fewer for
  3 GIEL ponds over 9–10 seasons) and the number of "non-event" seasons with real variation in the
  covariate of interest are enough to estimate a slope distinguishable from noise — the proposal's
  own "few ponds ⇒ don't pool intercepts" argument (§1.6) applies with at least as much force to
  "few informative pond-seasons ⇒ don't add slopes."
- **Practical implication.** The staged build the proposal describes (§1.7) is the right instinct —
  climb only as far as the data support — but the Havasu study is a documented case of that climb
  stalling at the very first rung (bare survival) for open water. Before committing to the full
  stack, it would be worth explicitly fitting the same staged progression to the current GIEL
  IP2/5/6 + CHLP data and reporting where it stalls, rather than presenting the complete DAG (§1.3)
  as the target architecture from the outset.

## 6. Critique: nonnative predator load as a φ³-only covariate

The proposal's release/habitat covariate set for post-stocking survival φ³ includes "nonnative-
predator load" (§1.6, Table 2 in `documents/BONY_Proposed_IPM_Summary.md`). Two things about this
are worth flagging given what open-water data already show:

- **The effect size predator load would need to explain is close to a ceiling, not a gradient.**
  Per `documents/C64_2016_LakeHavasu_Summary.md`, apparent mortality/loss in open water reached
  55–100% across all four Reach 3 sites, mostly within a month of release, "no matter the size" of
  the fish (per your note) and regardless of which specific predator assemblage was present
  (avian predation dominant at some sites, fish predation documented at others). When an outcome is
  this close to complete and this fast everywhere it's been measured, there may not be enough
  surviving contrast *within* open water to fit a meaningful predator-load slope — the data mostly
  distinguish "open water" from "predator-free pond," a near-binary contrast the model's two-limb
  scope split (OCH vs. open-water CJS, §1.1) already encodes structurally. A continuous predator-
  load covariate on φ³ is asking the data to resolve gradations within a regime that, empirically,
  looks close to uniformly lethal.
- **Restricting the covariate to φ³ only may not be the limiting assumption you'd expect.** The
  proposal's own scope split already puts all "predator load" variation into the open-water
  component (ponds are treated as "effectively predator-free," §1.6/§4), so φ⁴ never sees any
  predator-load contrast to begin with — the choice to exclude it from φ⁴'s linear predictor
  is largely moot rather than a meaningful modeling restriction, since φ⁴ is essentially never
  observed under nonzero predator load in this design. The substantive question is whether predator
  load varies enough *within* the open-water φ³ data (across reaches/sites/years) to be
  distinguishable from a simple site or reach fixed effect — which is closer to what the covariate
  is actually testing in practice than "does predator load matter," a question this project's own
  data already answer via the OCH-vs-open-water contrast itself (backwaters: bonytail persist for
  years; open water: bonytail are gone in weeks).
- **What the joint/stepwise structure actually buys here is unclear.** Because predator load is
  (by the model's own design) confined to the one submodel (open-water CJS) that does not share
  detection or abundance parameters with the OCH state-space side, a predator-load effect on φ³
  cannot inform, and is not informed by, anything estimated in the pond-side model. It's a
  standalone regression sitting inside a larger joint-likelihood framework, buying none of the
  cross-stream sharing that justifies an IPM's added complexity elsewhere (e.g., φ³/φ⁴ shared
  between the pond abundance and encounter sub-models, §1.3). If the open-water CJS component
  can't even resolve a bare survival curve (§5), adding a predator-load slope to it is unlikely to
  produce a defensible, generalizable coefficient — and given the near-ceiling mortality Havasu
  documents, it's questionable whether the resulting number would represent anything more than
  "open water is not survivable," which is already known without fitting a slope.

## 7. Open questions to revisit once the GIEL fit (job `C1FED42C`) returns

1. Does per-pond adult survival (`lphiA`-equivalent) show pond-specific baselines consistent with
   the proposed model's pond-fixed-effect + shared-covariate-slope design, or does it look more
   like ours (pond × season random draws around a species-level hyper-mean, no covariates)?
2. Do the confirmed die-off windows (IP5 2020-01–2020-09, IP2 2022-01–2022-07) fully explain the
   crash-and-recover pattern in known-alive counts, or is there residual unexplained variance that
   an environmental covariate (if ever collected) might legitimately explain?
3. How much of GIEL's adult survival variation (if any) is left after removing the two known
   crashes — is there enough remaining signal to argue for *any* covariate structure (environmental
   or otherwise), or is a simple pond/season random effect sufficient, as our exploratory
   assessment already suggests?
4. Does IP1FY22-style pond×season deviation logic (attempted for XYTE IP1, largely unsuccessful —
   see v12 result in `AGENTS.md`) generalize any better for GIEL's two confirmed, better-isolated
   crash windows (which, unlike IP1's, were pinned to single, better-covered ponds/seasons)?
5. Is stocked-vs-wild survival separable in any GIEL pond-season, and if so, does it change the
   pond-season baseline appreciably, or does the sparse wild sample size make it noise?
6. How many pond-seasons of "non-event" GIEL data actually exist (i.e., excluding the two confirmed
   crash windows), and is that enough to argue for or against *any* future covariate on GIEL
   survival — including a hypothetical predator-load or environmental term — given the Havasu
   precedent that even telemetry-grade open-water data failed to identify a bare survival curve
   (§5)?

## 8. Answers from the completed GIEL-only fit (`data/BWRobustDesign_MCMC_giel.RData`, 2026-09-10)

Fit: 3 ponds (IP2, IP5, IP6; model indices 1/2/3), `NSP = 2` (species 2 slot prior-only, ignore),
3 chains, 30k iter / 10k burn / thin 10, 60 min wall time. `dE[1]` (IP2) `-1.74 (-3.15, -0.55)`,
`dE[2]` (IP5) `-1.45 (-2.07, -0.87)`, `dE[3]` (IP6, no event, prior-only) `-2.42 (-6.76, -0.09)` — all
converge (Rhat ≤ 1.1). 42 of 673 monitored nodes have Rhat > 1.1, all either the known-alive-pin
floating-point artifact on `N_tag`/`Tal` or genuine recruitment nuisance parameters (`RecJ`, `r`,
`U` in the last 2–3 seasons) — `lphiA`, `dE`, detection, and hyperparameters all converge cleanly.

1. **Pond baselines vs. random draws:** looks like our structure, not the proposal's. Per-pond mean
   `lphiA` (logit, averaged over season) is IP2 2.20, IP5 2.89, IP6 3.04 — a real but modest
   ~0.8-logit spread — while within-pond season-to-season SD (0.69) exceeds the across-pond SD
   (0.45). Season variation dominates pond identity; `mu_phi = 2.71`, `sigma_phi = 0.85` (hyper-mean
   ± hierarchical draw) describes the data better than fixed pond intercepts + a shared covariate
   slope.
2. **Do the confirmed die-off windows fully explain the crash pattern? Partially, and unevenly.**
   Event-adjusted realized annual survival: IP5 FY2020 (9/12 event months) 0.041 (0.017–0.077)
   vs. the direct empirical floor of 0/39 known-alive adults re-seen after Oct 2020 — an excellent
   match. IP2 FY2022 (7/12 event months) 0.089 (0.044–0.147) vs. floor 0/14 — the model still sits
   well above a near-total-loss floor, the same under-correction pattern seen for XYTE IP1 FY2022
   (v10–v12 in `AGENTS.md`). More notably, **IP2's *non-event* seasons are also chronically low**
   (empirical floor 0.051 FY2018, 0.094 FY2019, 0.080 FY2021, 0.064 FY2025 — all comparable to or
   lower than the "crash" year), which the single flagged 2022 window does not address at all. This
   suggests IP2 may have additional unflagged mortality episodes, or a generally higher baseline
   mortality regime, that the current `dE`/`Evt` inventory (built from Marsh et al. 2024's two
   documented pump-failure events) doesn't capture — worth a closer look at IP2's known-alive series
   independent of this fit.
3. **Residual variation after removing the two crashes:** substantial. Even outside IP5-FY2020 and
   IP2-FY2022, floor-based annual survival ranges from ~0.05 (IP2, several seasons) to ~0.92 (IP5
   FY2019). The model already absorbs this via a free `lphiA[j,y]` per pond-season (a saturated
   random effect), so nothing is "unexplained" in a fitting sense — but a real covariate would need
   to reproduce this same swing pattern, and IP2's repeated low-survival seasons don't track the one
   documented event, arguing against a smooth environmental gradient and for either (a) a
   season/pond random effect being sufficient as a description, or (b) an as-yet-undocumented
   operational issue at IP2 specifically, not a data gap the current covariate wishlist would fill.
4. **Does the IP1FY22-style deviation logic generalize better here? Mixed, and the pattern is
   informative.** IP5's crash (event window covering 9 of 12 months that season) was matched well by
   `dE`; IP2's crash (7 of 12 months) showed the same under-correction as XYTE IP1's 7-month window.
   The common thread across both under-corrected cases (IP1/XYTE and IP2/GIEL) is a *partial-season*
   event window, not a species or pond-isolation difference — a single scalar `dE` applied to a
   window that doesn't span the (near-)whole season lets hierarchical shrinkage pull the discounted
   rate back toward the cross-pond/season mean. IP5's near-whole-season window avoided this. If a
   future refit wants a closer floor match for partial-season events, the fix is architectural (e.g.,
   a season-and-pond-specific `dE`, not a shared scalar), not a species-specific one.
5. **Stocked-vs-wild survival split — not testable from this fit as built.** Sample sizes for
   origin = 2 ("tagged at netting", our proxy for recruited/wild-hatched fish) are nontrivial at all
   three ponds — IP2 857 (vs. 599 stocked), IP5 159 (vs. 600 stocked), IP6 423 (vs. 300 stocked) — so
   there's no sample-size barrier to attempting a split. But **the joint 2-stage model does not
   currently have an origin term on established-adult survival (`lphiA`) at all**; `origin` only
   enters through the temporary post-release offset `b_post`, which applies for a few months after
   tagging, not to steady-state survival. This is a genuine architecture gap (also true of the XYTE
   side) rather than a data-sufficiency question this fit can resolve — answering it would require
   adding an origin effect to `lphiA` and refitting, not re-reading this fit's output.
6. **Informative non-event pond-seasons: 19 of 25.** 27 pond-season cells have any known-alive fish
   at their Oct census (FY2017, the stocking year, has none by construction); 2 are the flagged
   events, leaving 25 non-event cells, of which 19 have ≥15 known-alive fish (a reasonable floor
   sample). Combined with the finding in #3 that even these 19 "clean" cells swing from ~0.05 to
   ~0.92 survival with no apparent smooth pattern, this is a small and noisy base for fitting even
   one covariate slope shared across 3 ponds — reinforcing the Havasu-precedent argument in §5 that
   the data likely cannot support the proposed IPM's fuller covariate stack, and arguably not even a
   single environmental slope on GIEL adult survival at present.
