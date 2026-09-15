## GIEL_MaturationHazard_Build.R
##
## Step 3 of the growth-model sensitivity plan: convert each growth-model
## track's posterior draws into a size- and time-specific bonytail (GIEL)
## maturation hazard surface (psi_lookup), for use by the joint robust-design
## model's size-tied J->A maturation mechanism (see AGENTS.md "GIEL Size-Tied
## Maturation Hazard" and the plan file .posit/assistant/plans/2026-09-14-
## 1444-giel-growth-model-tracks-feeding-robust-design-maturation-hazard.md).
##
## Tracks: "F1" (Fabens, IPCA-only), "F2" (Fabens, IPCA+Cibola), "A1"
## (age-based VBGF, IPCA-only). Track A2 is out of scope.
##
## TL bins: 20mm bins spanning 100-249mm (8 bins), confirmed default from the
## plan. delta_t grid: 1-24 months post-entry. For each bin's midpoint TL0
## and each delta_t, project forward via the Fabens-renewal identity
## TL(t) = Linf - (Linf - TL0) * exp(-K_i * t_years), drawing a fresh
## individual/pond-year K_i = exp(logK + Normal(0, sigma_v)) per posterior
## draw, plus observation noise ~ Normal(0, sigma). Compute
## F(bin, delta_t) = P(TL(delta_t) >= 250mm), then convert to a discrete
## monthly hazard:
##   psi_lookup[bin, delta_t] = (F(delta_t) - F(delta_t-1)) / (1 - F(delta_t-1))
## with F(0) = 0. Beyond delta_t = 24, maturation is treated as certain
## (psi_lookup -> ~1), handled downstream in the NIMBLE code, not stored here.
##
## Run once per track by setting TRACK before sourcing, e.g.:
##   TRACK <- "F1"; source("model/GIEL_MaturationHazard_Build.R")
## Writes data/GIEL_MaturationHazard_<TRACK>.RData.

library(tidyverse)

if (!exists("TRACK")) TRACK <- "F1"
stopifnot(TRACK %in% c("F1", "F2", "A1"))
cat("Building maturation hazard table for track:", TRACK, "\n")

## ---- 1. Load posterior draws for the chosen track -------------------------

if (TRACK == "F1") {
  e <- new.env(); load("data/GIEL_GrowthSeasonThreshold.RData", envir = e)
  samps <- do.call(rbind, e$gs_winner$samples)
  Linf_draws <- samps[, "Linf"]
  logK_draws <- samps[, "logK[1]"]
  sigma_draws <- samps[, "sigma"]
  sigma_v_draws <- if (e$winner_use_group) samps[, "sigma_v"] else rep(0, length(Linf_draws))
  track_desc <- "Fabens increment model, IPCA-only (GIEL_GrowthSeasonThreshold.R)"
} else if (TRACK == "F2") {
  e <- new.env(); load("data/GIEL_GrowthSeasonThreshold_withCibola.RData", envir = e)
  samps <- do.call(rbind, e$comb_winner$samples)
  Linf_draws <- samps[, "Linf"]
  logK_draws <- samps[, "logK[1]"]
  sigma_draws <- samps[, "sigma"]
  sigma_v_draws <- if (e$winner_use_group_comb) samps[, "sigma_v"] else rep(0, length(Linf_draws))
  track_desc <- "Fabens increment model, IPCA + Cibola (GIEL_GrowthSeasonThreshold_withCibola.R)"
} else if (TRACK == "A1") {
  e <- new.env(); load("data/GIEL_VBGF_AgeLength.RData", envir = e)
  samps <- do.call(rbind, e$al_final$samples)
  Linf_draws <- samps[, "Linf"]
  logK_draws <- samps[, "logK[1]"]
  sigma_draws <- samps[, "sigma"]
  sigma_v_draws <- samps[, "sigma_v"]
  track_desc <- "Age-based VBGF (shared-K, release-anchored), IPCA-only (GIEL_VBGF_AgeLength.R)"
}

n_draws <- length(Linf_draws)
cat("Loaded", n_draws, "posterior draws for track", TRACK, "--", track_desc, "\n")

## ---- 2. TL bins and delta_t grid -------------------------------------------

TL_BIN_EDGES <- c(100, 120, 140, 160, 180, 200, 220, 240, 249)
N_TLBIN <- length(TL_BIN_EDGES) - 1
TL_bin_mid <- (TL_BIN_EDGES[-length(TL_BIN_EDGES)] + TL_BIN_EDGES[-1]) / 2
TL_bin_label <- paste0("[", TL_BIN_EDGES[-length(TL_BIN_EDGES)], ",", TL_BIN_EDGES[-1], ")")

DELTA_T_MAX <- 24  # months post-entry
delta_t_grid <- 1:DELTA_T_MAX

cat("\nTL bins (n =", N_TLBIN, "):\n")
print(data.frame(bin = seq_len(N_TLBIN), label = TL_bin_label, midpoint = TL_bin_mid))

## ---- 3. Simulate F(bin, delta_t) = P(TL(delta_t) >= 250mm) ----------------

set.seed(5566)
## For each TL bin, draw ONE K_i and ONE eps_i per simulated "individual"
## (n_draws of them), shared across the whole delta_t grid (common random
## numbers). This is required for F(bin, delta_t) to be a proper monotone
## non-decreasing cumulative probability in delta_t: for fixed K_i > 0, an
## individual's projected TL_t = TL0 + (Linf-TL0)*(1-exp(-K_i*t)) + eps_i is
## monotonically increasing in t, so the crossing indicator (and therefore
## its average across individuals) cannot decrease as delta_t grows. Drawing
## fresh noise independently at every (bin, delta_t) cell -- as an earlier
## version of this script did -- introduces pure Monte Carlo noise that can
## make F appear to dip from one delta_t to the next, which then produces a
## negative discrete hazard after clipping and breaks the psi_lookup <->
## F_mat reconstruction identity used as a sanity check below.
F_mat <- matrix(NA_real_, nrow = N_TLBIN, ncol = DELTA_T_MAX)
for (b in seq_len(N_TLBIN)) {
  TL0 <- TL_bin_mid[b]
  v_new <- rnorm(n_draws, 0, sigma_v_draws)
  K_i <- exp(logK_draws + v_new)
  eps <- rnorm(n_draws, 0, sigma_draws)
  for (dt in delta_t_grid) {
    t_years <- dt / 12
    TL_t <- TL0 + (Linf_draws - TL0) * (1 - exp(-K_i * t_years)) + eps
    F_mat[b, dt] <- mean(TL_t >= 250)
  }
}
dimnames(F_mat) <- list(TL_bin_label, paste0("m", delta_t_grid))

cat("\nCumulative F(bin, delta_t) at delta_t = 1, 6, 12, 24 months:\n")
print(round(F_mat[, c(1, 6, 12, 24)], 3))

## ---- 4. Convert to discrete monthly hazard psi_lookup ---------------------
## psi_lookup[b, dt] = (F(dt) - F(dt-1)) / (1 - F(dt-1)), F(0) = 0

psi_lookup <- matrix(NA_real_, nrow = N_TLBIN, ncol = DELTA_T_MAX)
for (b in seq_len(N_TLBIN)) {
  F_prev <- 0
  for (dt in delta_t_grid) {
    F_now <- F_mat[b, dt]
    denom <- 1 - F_prev
    psi_lookup[b, dt] <- if (denom > 1e-9) max(0, min(1, (F_now - F_prev) / denom)) else 1
    F_prev <- F_now
  }
}
dimnames(psi_lookup) <- list(TL_bin_label, paste0("m", delta_t_grid))

## Clip away from exact 0/1 before the logit-offset (lpsi_adj) mechanism in
## the NIMBLE model touches this table -- logit(0)/logit(1) are +-Inf, which
## algebraically still resolves correctly once lpsi_adj (finite) is added
## and ilogit() applied (ilogit(-Inf) = 0, ilogit(Inf) = 1), but clipping
## avoids relying on IEEE Inf arithmetic inside compiled NIMBLE code.
PSI_CLIP <- 1e-6
psi_lookup <- pmin(pmax(psi_lookup, PSI_CLIP), 1 - PSI_CLIP)

cat("\npsi_lookup (monthly hazard) at delta_t = 1, 6, 12, 24 months:\n")
print(round(psi_lookup[, c(1, 6, 12, 24)], 3))

## Sanity check: reconstructing F from psi_lookup should reproduce F_mat
F_check <- matrix(NA_real_, nrow = N_TLBIN, ncol = DELTA_T_MAX)
for (b in seq_len(N_TLBIN)) {
  surv <- 1
  for (dt in delta_t_grid) {
    surv <- surv * (1 - psi_lookup[b, dt])
    F_check[b, dt] <- 1 - surv
  }
}
max_recon_err <- max(abs(F_check - F_mat))
cat("\nMax reconstruction error (F_check vs F_mat):", signif(max_recon_err, 6), "\n")
## Tolerance relaxed from 1e-8 to allow for the PSI_CLIP applied above (which
## deliberately perturbs psi_lookup away from the raw F_mat-implied values at
## the extremes); a clean run before clipping showed <2e-16 error, so any
## discrepancy above ~PSI_CLIP magnitude here would indicate a real bug.
stopifnot(max_recon_err < 10 * PSI_CLIP)

## ---- 5. Save ---------------------------------------------------------------

OUT_FILE <- paste0("data/GIEL_MaturationHazard_", TRACK, ".RData")
BUILD_INFO <- list(track = TRACK, track_desc = track_desc, built = Sys.time(),
                    TL_BIN_EDGES = TL_BIN_EDGES, DELTA_T_MAX = DELTA_T_MAX)

save(psi_lookup, F_mat, TL_BIN_EDGES, TL_bin_mid, TL_bin_label,
     delta_t_grid, DELTA_T_MAX, N_TLBIN, TRACK, track_desc, BUILD_INFO,
     file = OUT_FILE)

cat("\nSaved:", OUT_FILE, "\n")
