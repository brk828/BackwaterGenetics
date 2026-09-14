# Summary: Bayesian VBGM for Gray Triggerfish (Chamberlin et al. 2025)

**Citation:** Chamberlin, D.W., Siders, Z.A., Potts, J.C., Rogers, W.D., Taylor, M.A., and
Patterson, W.F. III. 2025. Bayesian estimation of von Bertalanffy growth parameters for gray
triggerfish, *Balistes capriscus*, incorporating multiple readers and ageing structures. Can.
J. Fish. Aquat. Sci. 82: 1-12. doi:10.1139/cjfas-2024-0315.

## What the paper does

A sex-specific Bayesian hierarchical von Bertalanffy growth model (VBGM) for Gulf of Mexico
gray triggerfish, fit in Stan, comparing age estimates from three ageing protocols/structures
(whole otolith opaque-zone counts, a new dorsal-spine protocol, and the historical/old
dorsal-spine protocol) with up to 3 independent readers per structure. Core innovations:

1. **Separates observation error (ageing/reader imprecision) from process error (true
   growth variability).** Repeated age reads per fish/protocol are modeled as draws from a
   location-scale Student-t distribution centered on a latent "true age," with a
   reader-imprecision parameter (`sigma_ageing`) distinct from the growth-model's own residual
   SD (`sigma_growth`). This lets `sigma_growth` reflect only real biological variability, with
   ageing noise quantified and propagated separately rather than inflating growth-process
   uncertainty. Including multiple readers *increases* total reported uncertainty (relative to
   a single-reader model) because it makes previously-hidden observation error explicit,
   rather than reducing it.
2. **Size-at-birth VBGM parameterization** (`Li = Linf - (Linf - L0)*exp(-k*ti)`), fit
   simultaneously across protocols with a *shared* `L0` (birth size) but independent
   `Linf`, `k`, `sigma` per protocol; sex enters as an additive effect on the log scale
   (`log(Linf,s) = log(Linf) + betaL,s`, etc.) with a shared intercept across sexes -- this
   handles individuals of unknown sex without discarding them (they just contribute to the
   shared intercept, not either sex's additive term).
3. **Priors:** weakly informative normal priors on log-transformed `Linf`, `k`, `sigma`,
   `sigma_ageing`, centered on a quick maximum-likelihood fit with a 50% CV for prior width.
   The *one* strongly informative prior is on `L0` (birth size), pinned to an independently
   measured literature value (2.2 mm +/- 0.04 mm SD from reared larvae) -- i.e., they anchor
   the well-known end of the growth curve tightly and let the rest float. A sensitivity run
   with much more diffuse ("hyper weakly informative") priors on `Linf`/`k` changed final
   estimates only in the last significant digit, confirming the model is well identified
   given their data.
4. **Reports differences via joint posterior comparisons** (MAP-based p-values), not by
   comparing `Linf`/`k` point estimates independently, because `Linf` and `k` are highly
   correlated in the VBGM and sample size / presence-or-absence of small or large fish
   biases naive comparisons (citing Gwinn et al. 2010; Wilson et al. 2015).

## Key findings

- Otolith- and new-spine-protocol-derived ages agreed well; the historical old-spine
  protocol systematically underestimated age for fish >= 5 years, which inflated the
  fitted Brody growth coefficient `k` (old-spine `k` ~0.13-0.18 /yr higher than the
  validated protocols) -- age underestimation biases `k` upward without changing `Linf`
  much, because ageing error compresses the time axis, not the observed size range.
  Overestimating `k` overestimates natural-mortality proxies and stock productivity.
- Gray triggerfish show significant sexually dimorphic growth: males reach a larger
  asymptotic size than females (`Linf,M` ~502 mm vs `Linf,F` ~447 mm from the validated
  new-spine protocol), though the sexes track together until age 2-3.
- Studies lacking small fish near birth size produced strongly negative `t0` and
  implausibly large implied `L0` (10-16x the true measured value) with dampened `k`
  estimates -- a specific, cited example of size-selective-sampling bias in VBGM fits
  (Gwinn et al. 2010).

## Relevance to the GIEL growth-season threshold model (our project)

We lack age data entirely (no otoliths/spines for live PIT-tagged fish), so the paper's
central multi-reader ageing-error machinery does not transfer directly. What *does* transfer:

1. **Our `Linf`/`K` non-identifiability is a generic, well-documented VBGF problem, not a
   coding issue.** The paper explicitly notes `Linf` and `k` are highly correlated in this
   model family, and that lack of large fish near the asymptote biases `Linf` low / `k`
   high (Gwinn et al. 2010; Wilson et al. 2015) -- this is exactly our finding (`Linf`
   posterior at 417.7 mm sitting below the observed max `TL2` = 462 mm, with only 3 of 361
   rows having `TL1` > 400 mm). Their fix -- anchoring the *well-sampled* end of the curve
   with a strong, literature-grounded prior (`L0`) rather than forcing a prior on the
   poorly-sampled end (`Linf`) -- matches what we already tried and found *worse* (tight
   `Linf` prior cost ~149 WAIC units with clear residual bias). This is a useful citation to
   add to our own caveat language: the issue is fish-size coverage near the asymptote, not
   the absence of age data per se.
2. **Report derived predictions, not raw `Linf`/`K`, when parameters are correlated/poorly
   identified.** The paper recommends comparing size-at-age (or, in our case,
   probability-of-reaching-a-threshold) directly rather than comparing `Linf`/`k` point
   estimates, since those aren't independently interpretable when correlated and
   sample-coverage-dependent. This validates our existing approach of reporting
   `P(TL2 >= 250mm | TL1)` from `growth_season_threshold_probs` as the primary deliverable,
   with the `Linf`/`K` identifiability caveat attached, rather than treating `Linf`/`K`
   as literal biological parameters.
3. **Consider separating measurement/observation error from process error.** We already
   flagged 11 `TL2 < TL1` "negative growth" pairs as likely measurement/transcription noise
   (1-4 mm discrepancies for most, one 104 mm outlier) and excluded them rather than modeling
   them. The paper's approach suggests an alternative worth considering if more repeat/paired
   measurements become available: model a small explicit measurement-noise term separately
   from the pond-year growth-rate variance (`sigma_py`), rather than only hand-excluding
   outliers. With our current sample size this is likely not separately identifiable, but the
   ~1-4 mm noise floor from the flagged records could inform a fixed, literature/data-anchored
   measurement-SD if that extension is pursued later.
4. **Sex-specific growth is common in many fish species and worth checking as a covariate.**
   The paper's central finding is a strong, `Linf`-level sex effect. If sex is recorded for
   PIT-tagged GIEL fish (`StudyBWAnalysis$sex`) with adequate sample size, adding sex as an
   additive covariate (paralleling how `origin` is already handled in the GIEL growth model)
   would be a natural, low-cost extension to test for -- currently sex is not in the growth
   model at all.
5. **Priors-from-a-quick-MLE-fit + informative-prior-only-on-the-well-known-end** is a
   generally reusable pattern. Our own `Linf ~ Normal(635, 50)` prior is centered on the
   requested/literature TL range but (per points 1-2 above) sits on the poorly-sampled end of
   the curve, which is why it doesn't stabilize identifiability the way the paper's `L0` prior
   does for their model -- confirming that no simple prior tweak on `Linf` alone will fix our
   identifiability issue without more large-fish recaptures.

**Bottom line:** the paper does not offer a fix for our core data limitation (too few
recaptures near the size asymptote), but it independently confirms that this specific failure
mode is expected and well documented for the VBGF, provides citable support (Gwinn et al.
2010; Wilson et al. 2015) for the caveat language already attached to our `P(TL>=250mm)`
predictions, and surfaces one concrete, currently-untested extension worth checking: a sex
covariate on GIEL growth.
