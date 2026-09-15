## GIEL_Cibola_LinfSensitivity.R
##
## Historical sensitivity check on bonytail (GIEL) Linf, using pre-2005
## growth-recapture pairs from Cibola High Levee Pond (LID 10222, Cibola
## NWR, Reach 4) -- a DIFFERENT location from the three backwater ponds
## (IP2/IP5/IP6) used in GIEL_GrowthSeasonThreshold.R / GIEL_VBGF_AgeLength.R,
## included here ONLY because the backwater growth-increment data have almost
## no fish above ~400 mm TL and every fitted Linf collapses to a value below
## the observed max TL2 (see AGENTS.md "Linf-below-observed-max-size
## investigation" and the parallel finding in the age-based model). Cibola
## has real growth pairs from larger fish, but only from 1993-2005 (~20-30
## years old) plus 3 small-fish pairs from 2015-2016; no GIEL were physically
## recaptured there after 2017, so there is no way to build a RECENT
## large-fish comparison from this location -- see conversation log
## 2026-09-14. This script therefore fits an explicitly HISTORICAL,
## standalone sensitivity model -- it is NOT merged with or used to update
## the backwater growth models, and its Linf/K should not be read as
## describing current backwater bonytail growth. It exists purely to ask:
## "if we had had access to a size range extending well above 400 mm, would
## the implied asymptote have looked like the ~406-418 mm the backwater-only
## data currently imply?"
##
## Data: data/NFWGAnalysis.RData, NFWGAnalysis table, location_id == 10222
## (Cibola High Levee Pond), species == "GIEL", event %in% c("stocking",
## "capture"). Every consecutive pair of measured records per fish (same
## logic as GIEL_GrowthSeasonThreshold.R) gives 71 growth-increment rows,
## Date2 spanning 1994-04-07 to 2016-11-15. This script restricts to the
## PRE-2005 subset (Date2 < 2005-01-01, 67 of 71 rows) -- the era with real
## coverage of large fish (TL1 up to 504 mm, TL2 up to 510 mm) -- rather than
## the full history, to keep the check internally consistent (one pond, one
## management era) and because the only 3 post-2005 pairs are recent (2015-
## 2016), small fish (TL1 172-289 mm), and would not change the large-fish
## picture anyway.
##
## Growing-season exposure: same Apr 1 - Sep 30 window and calendar-overlap
## logic as the backwater growth models, applied here as a simplifying
## ASSUMPTION -- Cibola is a different, more southerly location and may not
## share the same growing-season timing as the IP ponds; not verified here.
## 13 of 67 pre-2005 pairs have TL2 < TL1 (negative growth, almost certainly
## measurement/transcription error, consistent with the pattern already
## documented for the backwater data) and are dropped; 2 more are winter-only
## (gs_days = 0) and dropped. 52 rows remain for fitting.
##
## Model: single-population Fabens/VBGF increment model, NO origin or
## pond-year grouping (unlike the backwater models) -- this dataset is a
## single pond/era with too few fish to support any grouping structure, and
## the point of this check is the single implied Linf/K, not sub-group
## contrasts:
##   increment = TL2 - TL1 ~ Normal(mu, sigma)
##   mu[i] = (Linf - TL1[i]) * (1 - exp(-K * t_years[i]))
## Priors matched to the backwater models for direct comparability:
## Linf ~ Normal(635, sd 50), logK ~ Normal(log(0.15), sd 0.5).
##
## RESULT (2026-09-14): this fit converges cleanly (R-hat <= 1.001, unlike
## the backwater fits' Linf/K ridge) and gives Linf = 662 mm (95% CI 581-749)
## and K = 0.316/yr (0.237-0.420) -- well ABOVE the observed max size in this
## dataset (510 mm; posterior P(Linf < 510) = 0), the opposite problem from
## the backwater fits (which never had any posterior support for Linf above
## their own observed max). This directly supports the hypothesis that the
## backwater models' low Linf (~406-418 mm) is an artifact of insufficient
## large-fish coverage, not a true reflection of a smaller adult size in this
## species/system.
##
## PRIOR-SENSITIVITY CAVEAT: a second fit with a much wider Linf prior
## (Normal(635, sd 300), otherwise identical) moves the estimate further, to
## Linf = 834 mm (sd 155) -- so the 662 mm point estimate is still
## meaningfully prior-influenced and should not be read as a precise value.
## What is robust across both prior widths is that Linf sits well above the
## backwater estimate and above this dataset's own observed max size in
## every posterior draw at both prior widths -- i.e., the direction and
## qualitative conclusion (backwater Linf is likely biased low) is robust;
## the exact magnitude of the "true" Linf is not pinned down by this data
## alone.
##
## STANDING CAVEATS -- repeat these in any report or document that cites this
## sensitivity check:
##  1. Data are 20-30 years old (1993-2005), from a different pond/location
##     under different management, water source, and possibly stocking
##     provenance than the current IP2/IP5/IP6 backwater ponds.
##  2. Growing-season window (Apr-Sep) is assumed, not verified, for Cibola.
##  3. This is a SINGLE-POND, single-era, ungrouped model -- no origin or
##     cohort effects are estimated, unlike the backwater growth models.
##  4. Linf/K here are themselves prior-sensitive (see above) and should be
##     read as a qualitative check ("is the backwater Linf plausible relative
##     to a location with real large-fish data?"), not as a replacement point
##     estimate for Linf.
##  5. This script does NOT feed back into GIEL_GrowthSeasonThreshold.R or
##     GIEL_VBGF_AgeLength.R -- no backwater-model files are modified here.

library(tidyverse)
library(lubridate)
library(nimble)
library(coda)

load("data/NFWGAnalysis.RData")   # NFWGAnalysis, LocationTable, LocationPolygons

CIBOLA_LID <- 10222   # Cibola High Levee Pond (see LocationTable)

## ---- 1. Build Cibola GIEL growth pairs (all consecutive measured records) -

cibola_recs <- NFWGAnalysis |>
  filter(
    location_id == CIBOLA_LID, species == "GIEL",
    event %in% c("stocking", "capture"),
    !is.na(total_length), !is.na(collection_date), !is.na(PITIndex)
  ) |>
  select(PITIndex, event, collection_date, total_length) |>
  arrange(PITIndex, collection_date)

cat("Qualifying physical NFWG records (any size), GIEL, Cibola High Levee Pond:",
    nrow(cibola_recs), "\n")

cibola_pairs <- cibola_recs |>
  group_by(PITIndex) |>
  filter(n() >= 2) |>
  arrange(collection_date, .by_group = TRUE) |>
  mutate(
    Date1 = collection_date, TL1 = total_length, Origin = event,
    Date2 = lead(collection_date), TL2 = lead(total_length)
  ) |>
  ungroup() |>
  filter(!is.na(Date2)) |>
  transmute(PITIndex, Origin, Date1, TL1, Date2, TL2,
            days_elapsed = as.numeric(Date2 - Date1)) |>
  arrange(Date1)

cat("Consecutive-pair growth-increment rows (all dates):", nrow(cibola_pairs), "\n")
cat("Date2 range:", paste(as.character(range(cibola_pairs$Date2)), collapse = " to "), "\n")

## ---- 2. Restrict to the pre-2005 historical, large-fish-coverage era ------

CIBOLA_CUTOFF <- as.Date("2005-01-01")
cibola_hist <- cibola_pairs |> filter(Date2 < CIBOLA_CUTOFF)

cat("\nPre-2005 pairs:", nrow(cibola_hist), "of", nrow(cibola_pairs), "\n")
cat("TL1 range, pre-2005 pairs:", paste(range(cibola_hist$TL1), collapse = "-"), "mm\n")
cat("TL2 range, pre-2005 pairs:", paste(range(cibola_hist$TL2), collapse = "-"), "mm\n")

## ---- 3. Growing-season (Apr 1 - Sep 30) exposure, negative-growth flag ----

growing_season_days <- function(d1, d2) {
  mapply(function(a, b) {
    yrs <- year(a):year(b)
    total <- 0
    for (yr in yrs) {
      gs_start <- ymd(sprintf("%d-04-01", yr))
      gs_end <- ymd(sprintf("%d-09-30", yr))
      os <- max(a, gs_start); oe <- min(b, gs_end)
      if (oe >= os) total <- total + as.numeric(oe - os)
    }
    total
  }, d1, d2)
}

cibola_hist <- cibola_hist |>
  mutate(increment = TL2 - TL1, gs_days = growing_season_days(Date1, Date2))

CibolaNegativeGrowthPairs <- cibola_hist |>
  filter(increment < 0) |>
  select(PITIndex, Origin, Date1, TL1, Date2, TL2, increment, days_elapsed, gs_days) |>
  arrange(Date1)

cat("\nNegative-growth pairs flagged and excluded (TL2 < TL1):",
    nrow(CibolaNegativeGrowthPairs), "\n")

cibola_hist <- cibola_hist |> filter(increment >= 0)

cibola_fit_data <- cibola_hist |> filter(gs_days > 0)

cat("Winter-only pairs excluded (gs_days = 0):",
    nrow(cibola_hist) - nrow(cibola_fit_data), "\n")
cat("Rows used for the historical Linf sensitivity fit:", nrow(cibola_fit_data), "\n")
cat("Max observed TL (TL1 or TL2) in fit data:",
    max(c(cibola_fit_data$TL1, cibola_fit_data$TL2)), "mm\n")

## ---- 4. Single-population Fabens/VBGF model (NIMBLE) ----------------------
##
## Priors matched to GIEL_GrowthSeasonThreshold.R / GIEL_VBGF_AgeLength.R for
## direct comparability: Linf ~ Normal(635, sd 50), logK ~ Normal(log(0.15),
## sd 0.5). No origin or pond-year grouping (single pond, single era, too few
## fish to support grouping).

LINF_PRIOR_MEAN <- 635
LINF_PRIOR_SD   <- 50
LOGK_PRIOR_MEAN <- log(0.15)
LOGK_PRIOR_SD   <- 0.5

fit_cibola_vbgf <- function(linf_sd, seed_set, niter = 60000, nburnin = 20000, thin = 10) {
  N <- nrow(cibola_fit_data)
  t_years <- cibola_fit_data$gs_days / 365
  TL1 <- cibola_fit_data$TL1
  increment <- cibola_fit_data$increment

  code <- nimbleCode({
    for (i in 1:N) {
      mu[i] <- (Linf - TL1[i]) * (1 - exp(-K * t_years[i]))
      increment[i] ~ dnorm(mu[i], sd = sigma)
    }
    logK ~ dnorm(logK_mean, sd = logK_sd)
    K <- exp(logK)
    Linf ~ dnorm(linf_mean, sd = linf_sd)
    sigma ~ dunif(0, 200)
  })
  constants <- list(N = N, t_years = t_years, TL1 = TL1,
                     linf_mean = LINF_PRIOR_MEAN, linf_sd = linf_sd,
                     logK_mean = LOGK_PRIOR_MEAN, logK_sd = LOGK_PRIOR_SD)
  inits <- list(logK = LOGK_PRIOR_MEAN, Linf = LINF_PRIOR_MEAN, sigma = 50)

  m <- nimbleModel(code, constants = constants, data = list(increment = increment),
                    inits = inits, calculate = TRUE)
  conf <- configureMCMC(m, monitors = c("Linf", "logK", "K", "sigma"))
  mcmc <- buildMCMC(conf)
  cm <- compileNimble(m)
  cmcmc <- compileNimble(mcmc, project = m)
  runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin,
          nchains = 3, samplesAsCodaMCMC = TRUE, setSeed = seed_set)
}

## Primary fit: prior matched to the backwater models (sd = 50)
cibola_hist_fit <- fit_cibola_vbgf(linf_sd = LINF_PRIOR_SD, seed_set = c(7331, 4162, 9048))
rhat_cib <- gelman.diag(cibola_hist_fit, multivariate = FALSE)$psrf[, 1]
samps_cib <- do.call(rbind, cibola_hist_fit)

cat("\n--- Primary fit (Linf ~ Normal(635, sd 50)) ---\n")
cat("R-hat:\n"); print(round(rhat_cib, 3))
cib_summ <- t(apply(samps_cib, 2, function(x)
  c(mean = mean(x), sd = sd(x), q2.5 = quantile(x, 0.025), q97.5 = quantile(x, 0.975))))
print(round(cib_summ, 3))

max_obs <- max(c(cibola_fit_data$TL1, cibola_fit_data$TL2))
cat("\nMax observed TL (TL1 or TL2):", max_obs, "mm\n")
cat("Posterior P(Linf < max observed size):",
    round(mean(samps_cib[, "Linf"] < max_obs), 4), "\n")

## Prior-sensitivity fit: much wider Linf prior (sd = 300)
cibola_hist_fit_wide <- fit_cibola_vbgf(linf_sd = 300, seed_set = c(111, 222, 333))
rhat_cib_wide <- gelman.diag(cibola_hist_fit_wide, multivariate = FALSE)$psrf[, 1]
samps_cib_wide <- do.call(rbind, cibola_hist_fit_wide)

cat("\n--- Sensitivity fit (Linf ~ Normal(635, sd 300)) ---\n")
cat("R-hat:\n"); print(round(rhat_cib_wide, 3))
cib_summ_wide <- t(apply(samps_cib_wide, 2, function(x)
  c(mean = mean(x), sd = sd(x), q2.5 = quantile(x, 0.025), q97.5 = quantile(x, 0.975))))
print(round(cib_summ_wide, 3))

cat("\nCONCLUSION: primary fit gives Linf =", round(mean(samps_cib[, "Linf"]), 1),
    "mm (95% CI", round(cib_summ["Linf", "q2.5.2.5%"], 1), "-",
    round(cib_summ["Linf", "q97.5.97.5%"], 1),
    "mm), well above both the observed max size in this dataset (", max_obs,
    "mm) and the backwater models' fitted Linf (~406-418 mm). The wide-prior",
    "sensitivity fit moves further still (Linf =", round(mean(samps_cib_wide[, "Linf"]), 1),
    "mm) rather than settling back down -- direction of the finding (backwater",
    "Linf likely biased low by lack of large-fish data) is robust to prior width;",
    "the precise magnitude of Linf is not. See header comments for full caveats.\n")

## ---- 5. Save results --------------------------------------------------------

save(
  cibola_pairs, cibola_hist, cibola_fit_data, CibolaNegativeGrowthPairs,
  CIBOLA_LID, CIBOLA_CUTOFF,
  cibola_hist_fit, rhat_cib, samps_cib, cib_summ,
  cibola_hist_fit_wide, rhat_cib_wide, samps_cib_wide, cib_summ_wide,
  LINF_PRIOR_MEAN, LINF_PRIOR_SD, LOGK_PRIOR_MEAN, LOGK_PRIOR_SD, max_obs,
  file = "data/GIEL_Cibola_LinfSensitivity.RData"
)
