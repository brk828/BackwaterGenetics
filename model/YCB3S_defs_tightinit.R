# Diagnostic variant of model/YCB3S_defs.R (2026-09-23).
#
# The YOY-classified refit (job BAD9373E, data/YCB3S_MCMC_m12.RData) failed to
# converge: chains 1 and 2 landed in a basin matching the pre-YOY reference
# fit (dJ[1] ~ -2.6, lpsiSL[1]/lpsiLA[1] ~ -2.5, sigma_phi[1] ~ 0.56-0.57),
# while chain 3 found a distinct alternate mode with near-instant S->L->A
# maturation (lpsiSL[1]/lpsiLA[1] flipped positive, ~+1.6 to +2.3),
# compensated by a much more negative juvenile post-release offset
# (b_post[1,2]), a weaker dJ[1], and larger sigma_phi[1]/sigma_J[1]. The
# divergence concentrated in eps_J[1, 8:9] (FY2024/FY2025 -- the two seasons
# whose recruit counts were trimmed by the carryover/YOY reclassification).
#
# This variant overrides make_3s_inits() so all 3 chains start close to the
# converged (chain 1/2) basin instead of the wide generic priors, to test
# whether the alternate mode is an avoidable mixing failure (a bad-inits
# artifact) rather than a real, comparable second posterior mass. Only
# species 1 (XYTE/YCB) is anchored; species 2 (GIEL) is prior-only in this
# YCB-only build and is left on the original generic inits. See AGENTS.md,
# "Young-of-Year (YOY) vs. Carryover Recruit Classification" section.

source("model/YCB3S_defs.R")

.make_3s_inits_orig <- make_3s_inits

make_3s_inits <- function(inp, seed) {
  inits <- .make_3s_inits_orig(inp, seed)
  cst <- inp$constants; NP <- cst$NPOND; NS <- cst$NS

  # Distinct sub-stream so the 3 chains still start with some spread around
  # the anchor point, rather than all 3 landing on identical values.
  set.seed(seed + 90000L)

  inits$mu_phi[1]    <- rnorm(1, 4.60, 0.08)
  inits$sigma_phi[1] <- runif(1, 0.48, 0.62)
  inits$dJ[1]        <- rnorm(1, -2.60, 0.10)
  inits$b_post[1, 1] <- rnorm(1, -3.03, 0.10)
  inits$b_post[1, 2] <- rnorm(1, -1.25, 0.10)
  inits$sigma_J[1]   <- runif(1, 0.90, 1.15)
  inits$lpsiSL[1]    <- rnorm(1, -2.55, 0.10)
  inits$lpsiLA[1]    <- rnorm(1, -2.53, 0.10)

  # eps_raw (pre-centering) seeded near the chain 1/2 eps_J trajectory by
  # season (FY2017-FY2026); centering inside the model (eps_J = eps_raw -
  # mean(eps_raw[j, ])) means only the relative pattern across seasons
  # matters, not the absolute level.
  epsJ_ref <- c(-1.40, 0.60, 0.55, 1.40, -1.50, 0.43, -0.65, 0.14, 0.25, 0.19)
  inits$eps_raw[1, ] <- epsJ_ref[seq_len(NS)] + rnorm(NS, 0, 0.15)

  # Regenerate lphiA consistent with the anchored mu_phi (same formula as the
  # original function, just re-evaluated after the override above).
  inits$lphiA <- matrix(rnorm(NP * NS, inits$mu_phi[cst$sp], 0.2), NP, NS)

  inits
}
