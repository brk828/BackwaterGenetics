# Bonytail (GIEL) Population Model — Methodology

**Purpose of this file:** companion to `model/GIEL_ModelResults_Summary.md`, aimed at an outside
Bayesian modeler who wants to understand (not just the results but) exactly what model was fit,
why the priors take the values they do, and how NIMBLE was configured to estimate it. Written
2026-09-11. Model code lives in `model/BWRobustDesign_defs.R` (shared with the Razorback Sucker
fits — the GIEL-only run is that same model code restricted to the 3 Bonytail ponds); the GIEL-only
fit itself was launched from `model/BWRobustDesign_NIMBLE.R` with `RUN_PONDS <- c(3, 6, 7)`.

## 1. Modeling framework

This is a **hierarchical Bayesian robust-design mark-recapture model**, in the Cormack-Jolly-Seber
(CJS) family, extended with:

- **Two latent life stages** (Juvenile/Adult, split at 250 mm total length) plus Dead, rather than
  the usual single "alive" state.
- **A robust-design detection layer**: within each monthly "primary" period, weekly ("secondary")
  scanner-coverage occasions are used to estimate detection probability as a per-week Bernoulli
  trial, summed to a Binomial per month, rather than treating "detected in month t" as a single
  binary outcome. This lets detection probability adapt to how much scanning effort (hours) was
  actually deployed.
- **An untagged ("wild"/recruit) sub-population** per pond, tracked as a deterministic-mean annual
  process (not individually), that is linked to the tagged population only through: (a) sharing the
  same survival/maturation parameters, and (b) a periodic netting event that samples both tagged and
  untagged fish and reports how many of the netted fish were already tagged — this is the
  information source that lets the model estimate untagged abundance at all.
- **A fixed, dated "die-off event" covariate** representing two documented Bonytail mass-mortality
  events (well pump failures), rather than a smooth environmental covariate.

Inference is fully joint: tagged-fish encounter histories, untagged-pool dynamics, and netting
counts are combined into one NIMBLE model and one MCMC run, so that all parameter uncertainty
(e.g., in survival or detection) automatically propagates into every derived quantity (abundance,
recruitment, event effects).

**Species/pond scope of the GIEL-only fit.** The general model code in `BWRobustDesign_defs.R`
supports an arbitrary subset of ponds and both species (XYTE and GIEL) jointly (used for the
7-pond fit). For the GIEL-only fit, `make_rd_inputs(c(3, 6, 7))` (pond indices for IP2/IP5/IP6)
remaps species indices to `1:length(unique(species present))` — since only GIEL is present, `NSP
= 1` and every species-indexed prior below is a single Bonytail-wide value, not shared with or
influenced by any Razorback Sucker data.

## 2. Data structure and time indexing

- **Study window:** FY2017 (first stocking after the 2014-2015 renovation) through FY2026
  (partial), 10 seasons.
- **Primary periods:** calendar months October-April each season (`T_PRIM` = 7 x 10 = 70 total
  primary periods). May-September is not a detection occasion — PIT scanning still occurs then,
  but the corresponding fish contacts are used only as *known-alive constraints* (see below), not
  as Binomial detection trials, because the original robust-design was built around the Oct-Apr
  netting/scanning protocol. `gap[t]` is 1 month within a season and 6 months from April to the
  following October.
- **Secondary (within-month) occasions:** each month is divided into up to 4 weeks (days 1-7,
  8-14, 15-21, 22-end). `K[i,t]` = number of weeks the scanner was actually running in that fish's
  pond that month (from `StudyBWEffort`, with a data-quality screen described in §7); `y[i,t]` =
  number of those weeks with >= 1 contact of that specific fish. This is what makes detection a
  **Binomial(K, p)** trial rather than a single Bernoulli per month.
- **Fish encounter codes (`stobs[i,t]`):** 0 = no informative observation that month; 1 = a
  direct observation as a juvenile (netting measurement < 250 mm); 2 = a direct observation as an
  adult (>= 250 mm); 3 = "known alive, stage unknown" — assigned to every month between a fish's
  first record and its last-ever PIT contact (`MaxScanDate`), including summer months, even when no
  stage measurement exists for that month. Code 3 is what lets May-September contacts inform the
  model without being treated as formal detection trials.
- **Removals:** harvest or inter-pond transfer ends a fish's encounter history *with the fish
  known alive* at that time (biomass removed from the pond, not a mortality event) — handled as a
  special terminal likelihood term, not as a missing/censored observation.
- **Netting events:** `data/BWNettingEvents.csv`-style records (harvested/returned/mortality counts
  by stage) supply, for each pond x month x stage, the total number of fish handled (`ev_n`) and how
  many of those were already-tagged fish (`ev_m`) — the single likelihood term that identifies
  untagged abundance (§4).

## 3. Latent-state (survival/maturation/detection) submodel

For each pond *j* and month *t*, a fish can be Juvenile, Adult, or Dead. NIMBLE code
(`rdPond` inside `BWRobustDesign_defs.R`) marginalizes each individual fish's unobserved
month-by-month state sequence via a **scaled forward-backward algorithm** (the HMM/CJS
standard), rather than sampling discrete latent states with data augmentation — this is far more
efficient in NIMBLE (no discrete-parameter MCMC, smooth log-likelihood surface, and it directly
returns the smoothed expected number of tagged fish alive by stage, `N_tag`, needed later for the
netting likelihood).

**Survival.** Monthly logit-scale survival:

```
logit(phiJ[j,y]) = lphiA[j,y] + dJ            # juvenile monthly survival
logit(phiA[j,y]) = lphiA[j,y]                 # adult monthly survival
```

with a **pond x season random effect** on adult survival, `lphiA[j,y] ~ Normal(mu_phi, sigma_phi)`,
shared hierarchically across all pond-seasons for the species. `dJ` is a single (not
pond/season-varying) offset representing the juvenile-vs-adult survival gap. A **post-release
penalty** `b_post[origin]` (origin = stocked vs. netting-tagged; 0 for established fish) is added
to `lphiA`/`lphiJ` for the first several months after tagging (`post[i,t] = 1` flag in the data),
representing handling/acclimation mortality distinct from steady-state survival. Over a multi-month
gap (the April-October transition), the per-month rate is simply raised to the `gap` power — an
assumption of a constant monthly rate within the gap.

**Known die-off events.** A fixed, **pond-specific, survival-reducing-only offset** `dE[j]` is
added to both `lphiJ` and `lphiA` in any primary month flagged by `Evt[j,t] = 1` (built from the
hand-curated `data/BWDieOffEvents.csv` date windows — see §6 for why this is a fixed dated
indicator rather than an estimated changepoint or a continuous environmental covariate). Ponds
with no recorded event (e.g. IP6) still carry a `dE[j]` parameter, but it is never multiplied by a
nonzero `Evt` value, so it is estimated entirely from its prior (a deliberate check built into the
results report — see `GIEL_ModelResults_Summary.md` §4, `dE[IP6]` sits at its prior mean with a
very wide CI).

**Maturation.** A single logit-scale monthly maturation rate `lpsi` (Juvenile -> Adult) applies
pond- and season-uniformly, again raised to `gap` when needed via `1 - (1-psi)^gap`.

**Detection.** Weekly per-fish detection probability:

```
logit(pJ[j,t]) = lp[j,1] + g_eff * eff[j,t] + d_moy[moy(t)]
logit(pA[j,t]) = lp[j,2] + g_eff * eff[j,t] + d_moy[moy(t)]
```

where `lp[j,stage]` is a pond x stage intercept, `eff[j,t]` is standardized log monthly scan hours
(so `g_eff` is a single shared effort-elasticity slope), and `d_moy` is a month-of-season fixed
effect (Oct = 0 reference level) capturing seasonal detectability patterns (e.g., spawning-season
aggregation). `y[i,t] ~ Binomial(K[i,t], p_stage[j,t])` for months with any scanning; when
`stobs[i,t]` codes a direct stage observation or "known alive" that month, the corresponding state
probability is fixed to 1 and the alternatives to 0 inside the forward pass, so it constrains the
hidden-state trajectory but is not treated as a second, redundant Binomial draw.

## 4. Untagged ("wild") pool submodel

Because untagged fish are not individually trackable, their abundance is represented as a
**deterministic-mean annual difference equation** per pond and stage, driven by the *same* survival
and maturation parameters estimated from the tagged fish (so the untagged pool borrows strength
from — and is constrained to be consistent with — the tagged-fish survival estimates, rather than
having independent survival):

```
Atot[j,y]      = N_tag[j, April_t, adult] + U[j, y, adult]     # total adults present in April
r[j,y]         = exp(lr[j,y]),  lr[j,y] ~ Normal(mu_r, sigma_r)  # recruits produced per adult
RecJ[j,y]      = r[j,y] * Atot[j,y]                              # new juveniles entering the pool
U[j,y+1,J]     = U_after_removal[j,y,J] * phiUJ[j,y] * (1-psiY) + RecJ[j,y]
U[j,y+1,A]     = U_after_removal[j,y,A] * phiUA[j,y] + matured juveniles + nu
```

`phiUJ`/`phiUA` apply the same monthly survival rates as the tagged model (raised to 12 months,
split between event and non-event months using `nEvt[j,y]`, the count of calendar months in that
season flagged as die-off months). `nu` is a small per-season "unexplained adult appearances" term
(catches any adults that enter a pond's untagged pool other than by aging in from the juvenile
stage — e.g., unreported supplemental stocking or model misspecification slack) with an
Exponential(1) prior, deliberately weakly informative and non-negative. `U[j,1,stage] = (1-Uk[j]) *
U0[j,stage]` where `Uk[j] = 1` for every IP pond (all three GIEL ponds were rotenone-renovated and
restocked from empty in 2017, so their untagged pool truly starts at zero by construction — no
`U0` parameter is estimated for them; `Uk = 0` only applied to Yuma Cove Backwater in the combined
7-pond fit, which was never emptied). `dconstraint()` nodes enforce `U[j,y,stage] >=
remU[j,y,stage]` (known removals of untagged fish can never exceed the modeled pool), and initial
values are iteratively repaired (`repair_inits()`) until this holds everywhere, since NIMBLE
requires a finite starting log-density.

## 5. Netting join (identifies untagged abundance)

The single likelihood link between tagged and untagged fish is the netting event: of `ev_n[e]`
fish handled in event *e* (pond, month, stage), `ev_m[e]` were already PIT-tagged. Given the
smoothed expected tagged count `Tal[e] = N_tag[pond,t,stage]` from the forward-backward pass and
the available untagged pool `Uav[e]` (after subtracting any untagged fish already removed earlier
that season), the model assumes fish are handled in proportion to their availability:

```
pm[e] = Tal[e] / (Tal[e] + Uav[e])         # (with small additive constants for numerical stability)
ev_m[e] ~ Binomial(ev_n[e], pm[e])
```

This is the only place `U` (untagged abundance) touches the data directly, so untagged-abundance
precision is limited by how many netting events exist and how large `ev_n` is at each one — noted
explicitly as a "recruitment is a nuisance parameter, not a precisely estimated one" caveat in the
results summary.

## 6. Priors and their rationale

All priors are set on the same NIMBLE model code used for both species (XYTE and GIEL); the table
below gives the values as they appear in `rd_code` and the reasoning behind each choice. None of
these were tuned against the GIEL data post hoc — they were chosen (or carried over unchanged) from
the earlier XYTE robust-design development, based on domain knowledge and standard weakly
informative practice, and then just re-used because a GIEL-only refit with `NSP = 1` uses the same
code path.

| Parameter | Prior | Scale / meaning | Rationale |
|---|---|---|---|
| `mu_phi` | Normal(3, sd = 1.5) | logit mean monthly adult survival | `logit(3) ≈ 0.95` monthly survival, i.e. a plausible pond-fish annual survival around exp(12*log(0.95)) ≈ 0.54 at the prior mean — centered near what is biologically expected for adult fish in a predator-free managed pond, but with sd = 1.5 logit units (a very wide range, roughly 0.5-0.999 monthly survival within 2 SD) so the data, not the prior, drive the pond-season estimates. |
| `sigma_phi` | Half-Normal(0, sd = 1), truncated >= 0 | pond x season SD of `lphiA` around `mu_phi` | Weakly informative half-normal on a variance-type parameter, standard practice to avoid the Jeffreys-type impropriety of a flat or inverse-gamma prior at small sample sizes; allows anywhere from near-homogeneous pond-seasons to substantial (multiple logit-units) heterogeneity. |
| `dJ` | Normal(-0.5, sd = 1) | logit offset, juvenile vs adult monthly survival | Centered slightly negative (juveniles usually survive somewhat worse than established adults) but wide enough (sd = 1) to let the data determine the actual juvenile survival deficit, including potentially near-zero. |
| `b_post[origin]` | Normal(-1, sd = 1), truncated <= 0 (established-origin fixed at 0) | logit offset, temporary post-release survival penalty | Truncated to be survival-reducing only (post-release mortality is a real, directional phenomenon — a positive offset would be biologically implausible and unidentifiable from the wide prior otherwise) with a moderately negative center reflecting known post-stocking mortality risk documented elsewhere in this project (see the Lake Havasu telemetry precedent, `documents/C64_2016_LakeHavasu_Summary.md`). |
| `lpsi` | Normal(-2.4, sd = 1) | logit monthly juvenile->adult maturation rate | `logit(-2.4) ≈ 0.083` monthly, i.e. roughly 1 year expected time-to-maturation at the prior mean — a biologically reasonable starting point for a species reaching ~250 mm within 1-2 years of stocking, but wide enough to be revised substantially by the netting-based direct stage observations. |
| `g_eff` | Normal(0, sd = 1) | detection elasticity to standardized log scan-hours | Centered at zero (no assumed direction) since more scanning *should* increase detection but the standardized covariate's sign/scale was unknown a priori; sd = 1 is wide relative to the realistic range of the standardized covariate. |
| `d_moy[2:NMOY]` | Normal(0, sd = 1) each (month 1 = Oct fixed at 0 as reference) | month-of-season fixed effects on detection | Standard reference-cell weakly informative prior for a categorical fixed effect; Oct is the arbitrary reference level (first primary month of the biological year). |
| `mu_r` | Normal(-1.2, sd = 1.5) | mean log recruits (juvenile-stage entrants) per adult | `exp(-1.2) ≈ 0.3` juveniles recruited per adult per season at the prior mean — a weak, order-of-magnitude prior on a term that is only loosely identified by the netting join, deliberately wide (sd = 1.5 on the log scale spans roughly two orders of magnitude) since the underlying data are known to be thin (see `GIEL_ModelResults_Summary.md` and the recruitment-nuisance caveat). |
| `sigma_r` | Half-Normal(0, sd = 1), truncated >= 0 | pond-season SD of log recruitment rate | Same rationale as `sigma_phi` — weakly informative, non-negative. |
| `nu` | Exponential(rate = 1) | mean unexplained adult appearances per pond-season | Non-negative by construction (Exponential support), small prior mean (1) reflecting that this term should normally be near zero and only picks up genuine model-untracked adult inflows (e.g., an unrecorded stocking); an Exponential rather than Normal/Half-Normal was chosen specifically because a hard floor at exactly 0 with a long right tail matches the expectation "usually nothing, occasionally something." |
| `lp[j,stage]` (pond x stage detection intercept) | Normal(1, sd = 1.5) | logit baseline weekly detection probability | `logit(1) ≈ 0.73`, reflecting that PIT scanning at these sites is generally quite effective when running, but sd = 1.5 is wide enough to let poorly-scanned ponds/stages (e.g., small juveniles) fall much lower. |
| `dE[j]` | Half-Normal(0, sd = 3), truncated <= 0 | fixed pond-specific die-off logit-survival offset, applied only in flagged event months | Deliberately very wide and one-directional: this parameter only has data to inform it in ponds/seasons with a recorded event (`Evt = 1`); a wide, survival-reducing-only prior lets the magnitude of a real crash be estimated essentially unconstrained by data outside the event window, while for ponds with **no** recorded event (IP6 here) the posterior is expected to — and does — collapse back to (a slightly negative, since truncated) version of the prior itself, which is used explicitly as an internal validity check in the results summary. |
| `U0[j,stage]` | Uniform(0, 5000) | initial (FY2017) untagged pool size, IP ponds only | A wide, essentially uninformative prior appropriate for ponds where `Uk[j] = 1` forces `U0` to not actually enter the likelihood at all (all three GIEL ponds were emptied by rotenone renovation before 2017 stocking) — retained in the model code for generality (it *is* used for Yuma Cove Backwater in the full 7-pond fit) but is a structural no-op for the GIEL-only run. |

**Species-indexed vs. pond-indexed parameters.** `mu_phi`, `sigma_phi`, `dJ`, `b_post`, `lpsi`,
`g_eff`, `d_moy`, `mu_r`, `sigma_r`, `nu` are indexed by species (here, a single GIEL-wide value
each, since `NSP = 1` for the GIEL-only subset — see §1). `lphiA[j,y]`, `lp[j,stage]`, `dE[j]`,
`U0[j,stage]`, `lr[j,y]` are indexed by pond (and season, where applicable), drawn hierarchically
around the species-level parameters where a hierarchical structure exists (`lphiA`, `lr`) or given
their own independent, wide priors where it does not (`lp`, `dE`, `U0` — chosen deliberately, since
detection quality and die-off exposure are pond-specific engineering/event facts, not expected to
share statistical strength the way survival biology does across ponds of the same species).

## 7. Data screening that feeds the likelihood

Before any fitting, `model/BWRobustDesign_data.R` applies a **scanner-week screen**: an Oct-Apr
week with effort deployed but zero contacts of *any* fish is dropped from `K` (treated as
unavailable, not as "scanned with zero detections") if either (a) its hours are below 25% of that
pond's median weekly hours (a probable partial/failed deployment), or (b) at least 20 tagged fish
are known alive that month and the implied probability of observing zero contacts under the
pond-season's typical detection rate is below 1e-3 (a probable scanner malfunction). This avoids
biasing detection-probability estimates downward from data-quality artifacts rather than genuine
low detectability; dropped weeks are retained in `DroppedWeeks` for auditing. This is a data-prep
step, not part of the NIMBLE model itself, but is essential context for interpreting `K`, `lp`, and
`g_eff`.

## 8. MCMC implementation in NIMBLE

**Chains and length.** 3 parallel chains (one R worker process per chain via `parallel::parLapply`
onto a `makeCluster()`), each run for 30,000 iterations with the first 10,000 discarded as burn-in
and the remainder thinned by 10 — 2,000 retained posterior draws per chain, 6,000 total. Distinct
fixed seeds per chain (`4127`, `8351`, `2960` for the general driver; the GIEL-only run used the
same three-seed scheme). Wall time for the GIEL-only (3-pond, 1-species) fit was ~60 minutes;
compare to ~3 hours for the full 7-pond joint fit.

**Marginalization strategy.** Rather than including each fish's monthly true state as a discrete
latent parameter (which NIMBLE can sample but which mixes poorly and scales badly with the number
of fish x months), the entire per-fish, per-pond likelihood is computed by the custom
`rdPond()` nimbleFunction: a scaled forward-backward recursion over each unique encounter history
(after collapsing identical histories across fish into one weighted row via `make_rd_inputs()`, to
avoid redundant computation), returning (a) the total log-likelihood for that pond's data and (b)
the smoothed posterior-mean number of tagged fish alive by stage at each month (`N_tag`), needed by
the netting-join likelihood. The scalar log-likelihood is injected into the NIMBLE graph via a
"zeros trick" custom distribution (`dLLtrick`): a dummy data node `zero[j] ~ dLLtrick(ll)` whose
log-density is defined to just equal the supplied `ll`, a standard technique for incorporating an
arbitrary externally-computed likelihood into a NIMBLE/BUGS-style model graph.

**Sampler configuration.** NIMBLE's default configuration would assign a univariate adaptive
random-walk Metropolis sampler to nearly every continuous stochastic node. Several groups of
parameters are highly correlated in this model (survival-offset parameters that only ever appear
summed together in the likelihood; detection parameters that trade off similarly) and mix far
better as `RW_block` (adaptive multivariate random-walk) samplers instead of one-at-a-time
univariate updates:

- Per species: `{dJ, b_post[,1], b_post[,2], lpsi}` blocked together (all jointly determine the
  realized juvenile survival, which is only identified as a sum/combination in the likelihood).
- Per species: `{g_eff, d_moy[2:NMOY]}` blocked together (detection covariates that trade off
  against each other and the pond-level `lp` intercepts).
- Per pond: `{lphiA[j, 1:NS]}` (all season-specific survival deviations for that pond) blocked
  together.
- Per pond: `{lp[j, 1], lp[j, 2]}` (juvenile/adult detection intercepts) blocked together.

All other nodes (species hyperparameters `mu_phi`/`sigma_phi`/`mu_r`/`sigma_r`/`nu`, pond-level
`dE`, `U0`, `lr`) keep NIMBLE's default univariate adaptive samplers, since they are less directly
correlated with each other in the likelihood. Block-sampler adaptation interval is set to 200
iterations (`control = list(adaptInterval = 200)`), NIMBLE's mechanism for periodically re-tuning
the proposal covariance/step size during the early part of the chain.

**Initialization.** `make_inits()` draws chain-specific starting values from tight, biologically
plausible distributions around the prior means (e.g., `mu_phi ~ Normal(3, 0.3)`, well inside the
Normal(3, 1.5) prior), rather than from the full prior spread, purely to give the sampler a
numerically stable starting point; recruitment (`mu_r`) is deliberately started high (`log(3)` per
adult, well above the prior mean of `-1.2`/`exp(-1.2) ≈ 0.3`) so the deterministic untagged-pool
skeleton starts comfortably above its `U >= remU` floor, then the posterior naturally pulls it back
down as the chain runs. `repair_inits()` is a small deterministic loop (not a sampler) that,
starting from `make_inits()`'s draw, iteratively raises whichever pond-season's untagged pool
component (`U0`, `lr`, or `nu`) is causing a `dconstraint()` violation, re-checking the model's
overall log-density after each nudge, until a finite starting log-probability is achieved (capped
at 200 attempts, after which it errors out rather than silently proceeding with an invalid start) —
this is necessary because NIMBLE (like any MCMC engine) cannot begin sampling from a
zero-probability point in parameter space, and the joint constraint structure here makes a
biologically-motivated random initial draw occasionally infeasible without such a repair step.

**Monitored nodes.** All species- and pond-level parameters listed in §6, plus derived quantities
`lphiA`, `r`, `U` (untagged pool), `RecJ` (recruits), `N_tag` (tagged fish alive by stage),
`Tal` (tagged count at each netting event), and `Atot` (total spawning-season adult abundance) —
i.e., both the primary parameters and the abundance/survival quantities users actually want to
report on are full posterior objects, not post hoc point-estimate plug-ins.

## 9. Convergence diagnostics used

- **Gelman-Rubin R-hat** (`coda::gelman.diag(..., multivariate = FALSE)`) computed across all 3
  chains for every monitored node (up to several hundred, depending on how many pond-seasons and
  netting events are in scope). A node is flagged if R-hat > 1.1 (the summary report additionally
  distinguishes R-hat > 1.05 as a stricter secondary threshold).
- **Effective sample size** (implicitly via `MCMCvis::MCMCsummary()`, which reports `n.eff`
  alongside posterior mean/SD/quantiles for the core hyperparameter set).
- **Trace plots** for every core hyperparameter, saved to a PDF (`MCMCvis::MCMCtrace()`), for
  visual inspection of chain mixing and stationarity beyond the summary statistics.
- **A known, deliberately-tolerated R-hat artifact:** many `N_tag`/`Tal` nodes register as R-hat >
  1.1 purely because of floating-point noise on quantities that are essentially pinned to a fixed
  value by the known-alive (`stobs = 3`) data constraint in a given month — i.e., all three chains
  agree to ~1e-6 to 1e-8, and `gelman.diag()`'s ratio-of-variances statistic becomes numerically
  unstable (not scientifically meaningful) at that level of within-chain agreement. This was
  verified directly for the analogous XYTE three-stage model by checking within-chain SDs on
  flagged nodes before applying the same reasoning here; it should **not** be mistaken for a real
  mixing failure, and any newly-flagged node should be checked this way before being treated as a
  convergence problem.
- **Substantive, non-artifactual mixing concerns are reported separately from the artifact above**
  — for the GIEL-only fit, per `GIEL_ModelResults_Summary.md` §4, the parameters that matter for
  inference (`lphiA`, `dE`, detection, hyperparameters) all converged cleanly at R-hat <= 1.1, and
  the ~42 flagged nodes were confirmed to be either the artifact above or late-season recruitment
  nuisance parameters (`RecJ`, `r`, `U` in the final 2-3 seasons, where data are thinnest) rather
  than a modeling problem.

## 10. Practical reproducibility notes

- **Rebuild-before-refit discipline.** `model/BWRobustDesign_NIMBLE.R` loads
  `data/BWRobustDesign_data.RData` but does not itself rebuild it; a `BUILD_INFO`-based mtime guard
  at the top of the script stops the run (or warns, for older pre-guard data files) if the input
  script/data are newer than the saved NIMBLE-ready data, to prevent a previously-encountered
  stale-data bug (see `AGENTS.md`, "Known Code Issues") from recurring silently.
- **Pond subsetting is a first-class option**, not a hack: `RUN_PONDS` is set before sourcing the
  driver script (e.g., `c(3, 6, 7)` for GIEL-only), and `make_rd_inputs()` re-derives `NSP`
  (species count actually present), `pond_new` indices, and all constants/data arrays fresh for
  that subset — so a GIEL-only fit genuinely never sees XYTE data or shares hyperparameters with a
  XYTE fit, rather than merely filtering a joint posterior after the fact.
- **All settings, seeds, and run time are saved alongside the posterior** (`rd_settings` in the
  output `.RData`), so a given results file is self-documenting for exactly how it was produced.

## 11. Where this differs from the Razorback Sucker (XYTE) model

For contrast (see `model/YCB3S_data.R`/`YCB3S_defs.R` and the "Three-Stage Density-Dependent Model"
notes in `AGENTS.md` for full detail): the XYTE model at Yuma Cove Backwater uses **three** alive
stages (Small/Large/Adult, splitting the juvenile stage further) with a **density-dependent
recruit-size submodel** (`logit(pLg) = a0 + a1 * recruit_density`) and season-specific juvenile
survival deviations (`eps_J`), because XYTE recruitment is large, episodic, and shows a genuine
size-density gradient. Both the size-vs-survival check and the density-vs-survival check reported
in `GIEL_ModelResults_Summary.md` §2 concluded Bonytail data do **not** support that added
complexity — hence the decision to reuse the simpler, already-existing 2-stage joint model
described in this document rather than port the three-stage architecture to GIEL. The `dE`
die-off-event mechanism described in §3 above was originally built for the XYTE three-stage model
and ported into this 2-stage model afterward; the underlying mathematical device (a fixed,
dated, pond-specific logit-survival offset activated by a data-derived indicator) is identical in
both.
