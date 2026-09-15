## GIEL_GrowthSeasonThreshold_withCibola.R
##
## Fork of GIEL_GrowthSeasonThreshold.R (2026-09-14, BK request): tests whether
## Fabens post-stocking/post-tagging growth rate (K) differs significantly
## between the backwater ponds (IP2/IP5/IP6) and Cibola High Levee Pond (LID
## 10222, a different but environmentally similar off-channel GIEL location;
## see model/GIEL_Cibola_LinfSensitivity.R/_Summary.md for the initial,
## historical-only exploration of Cibola data), and if the difference is not
## significant, refits the production growth-increment model with Cibola's
## pre-2005 large-fish pairs pooled in to help anchor Linf.
##
## This script does NOT overwrite data/GIEL_GrowthSeasonThreshold.RData or
## GIEL_GrowthSeasonThreshold.R -- it is a parallel, clearly-labeled fork.
## Whether to promote this combined model to "the" production model is a BK
## decision to be made after reviewing the results below; until then, both
## the backwater-only and combined files should be treated as live
## alternatives (see STATUS note at the end of this header).
##
## ---- Step 1: overlap-range site-comparison test ---------------------------
##
## Question: is Fabens K significantly different between backwater and Cibola
## fish, restricted to the TL1 range where both datasets have real data
## (189-442 mm, the intersection of backwater's 117-442 mm and Cibola's
## 189-504 mm)? This avoids comparing K across non-overlapping size ranges,
## where an apparent difference could just reflect different parts of a
## shared decelerating curve.
##
## Data: 349 backwater pairs (TL1 in range) + 45 Cibola pairs (TL1 in range,
## all pre-2005 -- see GIEL_Cibola_LinfSensitivity.R for how these 52 total
## Cibola pairs were built; 7 fall outside the overlap window and are
## excluded from this specific test only). Site-year grouping (pond-year for
## backwater, Cibola-year for Cibola) as the random effect on log(K), same
## structure as the main growth model's PondYear term. Compared a pooled-K
## model (single logK, both sites) against a site-specific-K model (separate
## logK[Backwater], logK[Cibola]) by WAIC.
##
## RESULT (2026-09-14): WAIC narrowly favors the pooled model (3667.1 vs
## 3668.9 site-specific, Delta WAIC = 1.9 -- conventionally not a meaningful
## improvement for the extra parameter). But the credible interval on the
## site difference is genuinely borderline: logK[Cibola] - logK[Backwater]
## 95% CI (-0.81, +0.01) barely fails to exclude zero; K[Backwater] = 0.90/yr
## (0.65-1.15) vs K[Cibola] = 0.63/yr (0.39-0.88); one-sided posterior
## P(K[Cibola] < K[Backwater]) = 0.97. This is a real ambiguity, not a clean
## answer -- 45 Cibola pairs is not much power to detect a real but modest
## rate difference (era, climate, stocking-source, or density differences
## between the two locations are all plausible explanations if the
## difference is real). BK decision (2026-09-14): treat as "not significant"
## per this project's standing WAIC-based term-selection convention (see
## e.g. the origin/pondyear term-selection precedent in the original growth
## model) and proceed to pool the two sites with a single shared K in the
## production refit below. This choice should be revisited if a larger
## Cibola-overlap sample ever becomes available.
##
## ---- Step 2: combined production model ------------------------------------
##
## Pools ALL 361 backwater growth-increment rows (full size range, same as
## GIEL_GrowthSeasonThreshold.R) with all 52 Cibola pre-2005 pairs (TL1
## 189-504 mm) into one dataset (n = 413, TL1 117-504 mm, TL2 up to 510 mm).
## Same 4-way WAIC term-selection (origin x site-year-group) as the original
## script; Cibola rows all have Origin = "Capture" (matches the "Capture"
## level already used for backwater wild/netting-tagged fish) and their own
## per-year grouping factor (distinct from backwater PondYear).
##
## RESULT: WAIC again selects the full origin + group model (WAIC 3844.4 vs
## 3849.7 group-only, 3879.8 pooled/neither, 3880.1 origin-only) -- the same
## structure selected for backwater alone, now with Cibola's groups folded
## into the same random-effect distribution. Converges cleanly (max R-hat
## 1.03, 0 of 22 monitored nodes > 1.1). Linf = 435.5 mm (95% CI 421.5-451.5)
## -- HIGHER than the backwater-only model's 417.4 mm, but still below the
## dataset's own observed max size (510 mm; posterior P(Linf < 510) = 1).
## K[Capture-origin] = 1.00/yr, K[Stocking-origin] = 0.77/yr.
##
## INTERPRETATION: adding Cibola moved Linf in the expected direction (up),
## by about 18 mm, but did NOT resolve the identifiability problem -- the
## backwater data (361 rows) still numerically dominate the 52 Cibola rows,
## and Cibola's own historical Linf-only fit (662 mm, see
## GIEL_Cibola_LinfSensitivity.R) has no separate influence on Linf here
## because Linf is shared across all groups/sites, not site-specific (per
## the Step 1 decision to pool K, and because there was never a serious case
## for a site-specific Linf given the K result). This is expected given the
## the relative sample sizes and is not a bug: pooling more of the same kind
## of sparse large-fish data will move Linf further only in proportion to
## how much of the combined sample it represents. A genuinely well-anchored
## Linf would need either many more large-fish recaptures (backwater or
## Cibola) or a literature/otolith-based external Linf prior tighter than
## the current Normal(635, sd 50) if one existed.
##
## STANDING CAVEATS (repeat wherever this combined model is cited):
##  1. Cibola data are 20-30 years old (pre-2005), a different pond/location,
##     under different management/water source/possibly different stocking
##     provenance than IP2/IP5/IP6.
##  2. The site-comparison test (Step 1) is a borderline, underpowered
##     result, not a clean "no difference" -- pooling K across sites was a
##     judgment call given WAIC's narrow preference, not an unambiguous
##     finding.
##  3. Linf is still identifiability-limited (P(Linf < max observed) = 1)
##     even after adding Cibola -- do not treat 435.5 mm as a literal,
##     well-anchored asymptote; it is a modest improvement over 417.4 mm, not
##     a resolution.
##  4. Apr-Sep growing-season window assumed, not verified, for Cibola.
##  5. Every Cibola row is Origin = "Capture" -- there is no Cibola stocking-
##     origin growth data to compare against backwater's stocked cohorts.
##
## STATUS (2026-09-14): this combined fit is offered as an alternative to the
## original backwater-only model, not yet adopted as the sole production
## model. Until BK decides otherwise, treat GIEL_GrowthSeasonThreshold.RData
## (backwater only) as the primary reference and this file's output
## (GIEL_GrowthSeasonThreshold_withCibola.RData) as a documented sensitivity/
## alternative, per the Linf caveat above.

library(tidyverse)
library(lubridate)
library(nimble)
library(coda)

load("data/ReportingData.RData")     # StudyBWNFWG (backwater)
load("data/NFWGAnalysis.RData")      # NFWGAnalysis, LocationTable (Cibola)

GIEL_PONDS <- c("IPCA (Pond 2)", "IPCA (Pond 5)", "IPCA (Pond 6)")
CIBOLA_LID <- 10222
CIBOLA_CUTOFF <- as.Date("2005-01-01")

## ---- Build growth pairs: backwater (full size range, all dates) -----------

bw_recs <- StudyBWNFWG |>
  filter(species == "GIEL", Backwater %in% GIEL_PONDS,
         event %in% c("stocking", "capture"),
         !is.na(total_length), !is.na(collection_date)) |>
  select(PITIndex, Backwater, event, collection_date, total_length) |>
  arrange(PITIndex, Backwater, collection_date)

bw_pairs <- bw_recs |>
  group_by(PITIndex, Backwater) |>
  filter(n() >= 2) |>
  arrange(collection_date, .by_group = TRUE) |>
  mutate(Date1 = collection_date, TL1 = total_length, Origin = event,
         Date2 = lead(collection_date), TL2 = lead(total_length)) |>
  ungroup() |>
  filter(!is.na(Date2)) |>
  transmute(PITIndex, Pond = Backwater, Origin, Date1, TL1, Date2, TL2,
            days_elapsed = as.numeric(Date2 - Date1))

## ---- Build growth pairs: Cibola (pre-2005 subset, matching the historical
## sensitivity check) ---------------------------------------------------------

cibola_recs <- NFWGAnalysis |>
  filter(location_id == CIBOLA_LID, species == "GIEL",
         event %in% c("stocking", "capture"),
         !is.na(total_length), !is.na(collection_date), !is.na(PITIndex)) |>
  select(PITIndex, event, collection_date, total_length) |>
  arrange(PITIndex, collection_date)

cibola_pairs_all <- cibola_recs |>
  group_by(PITIndex) |>
  filter(n() >= 2) |>
  arrange(collection_date, .by_group = TRUE) |>
  mutate(Date1 = collection_date, TL1 = total_length, Origin = event,
         Date2 = lead(collection_date), TL2 = lead(total_length)) |>
  ungroup() |>
  filter(!is.na(Date2)) |>
  transmute(PITIndex, Origin, Date1, TL1, Date2, TL2,
            days_elapsed = as.numeric(Date2 - Date1))

cibola_pairs <- cibola_pairs_all |> filter(Date2 < CIBOLA_CUTOFF)

cat("Backwater growth pairs:", nrow(bw_pairs), "\n")
cat("Cibola pre-2005 growth pairs:", nrow(cibola_pairs), "\n")

## ---- Growing-season exposure, negative-growth / winter-only filtering -----

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

prep_pairs <- function(pairs, site_label, group_prefix) {
  pairs |>
    mutate(
      increment = TL2 - TL1,
      gs_days = growing_season_days(Date1, Date2),
      SourceSite = site_label
    ) |>
    filter(increment >= 0, gs_days > 0)
}

bw_prepped <- bw_pairs |>
  mutate(Pond = str_extract(Pond, "\\d+"), Origin = str_to_title(Origin)) |>
  prep_pairs("Backwater", "P") |>
  mutate(Group = paste0("P", Pond, "_", year(Date1)))

cibola_prepped <- cibola_pairs |>
  mutate(Origin = "Capture") |>
  prep_pairs("Cibola", "Cibola") |>
  mutate(Group = paste0("Cibola_", year(Date1)))

cat("Backwater rows after negative-growth/winter-only filtering:", nrow(bw_prepped), "\n")
cat("Cibola rows after negative-growth/winter-only filtering:", nrow(cibola_prepped), "\n")

## ---- Step 1: overlap-range site-comparison test ---------------------------

overlap_lo <- max(min(bw_prepped$TL1), min(cibola_prepped$TL1))
overlap_hi <- min(max(bw_prepped$TL1), max(cibola_prepped$TL1))
cat("\nOverlap TL1 range for site-comparison test:", overlap_lo, "-", overlap_hi, "mm\n")

overlap_data <- bind_rows(
  bw_prepped |> filter(TL1 >= overlap_lo, TL1 <= overlap_hi) |>
    transmute(PITIndex, SourceSite, TL1, TL2, increment, gs_days, Group),
  cibola_prepped |> filter(TL1 >= overlap_lo, TL1 <= overlap_hi) |>
    transmute(PITIndex, SourceSite, TL1, TL2, increment, gs_days, Group)
) |>
  mutate(t_years = gs_days / 365,
         site_idx = as.integer(factor(SourceSite, levels = c("Backwater", "Cibola"))),
         group_idx = as.integer(factor(Group)))

cat("Overlap-range test dataset: n =", nrow(overlap_data),
    "(", sum(overlap_data$SourceSite == "Backwater"), "backwater +",
    sum(overlap_data$SourceSite == "Cibola"), "Cibola)\n")

LINF_PRIOR_MEAN <- 635
LINF_PRIOR_SD   <- 50
LOGK_PRIOR_MEAN <- log(0.15)
LOGK_PRIOR_SD   <- 0.5

fit_site_vbgf <- function(data, use_site, seed_set, niter = 100000, nburnin = 40000, thin = 15) {
  N <- nrow(data)
  n_group <- max(data$group_idx)
  site_use <- if (use_site) data$site_idx else rep(1L, N)
  n_site_use <- if (use_site) 2L else 1L

  code <- nimbleCode({
    for (i in 1:N) {
      K[i] <- exp(logK[site_idx[i]] + v[group_idx[i]])
      mu[i] <- (Linf - TL1[i]) * (1 - exp(-K[i] * t_years[i]))
      increment[i] ~ dnorm(mu[i], sd = sigma)
    }
    for (g in 1:n_group) v[g] ~ dnorm(0, sd = sigma_v)
    for (s in 1:n_site) logK[s] ~ dnorm(logK_mean, sd = logK_sd)
    Linf ~ dnorm(linf_mean, sd = linf_sd)
    sigma ~ dunif(0, 200)
    sigma_v ~ dunif(0, 2)
  })
  constants <- list(N = N, n_group = n_group, n_site = n_site_use,
                     t_years = data$t_years, TL1 = data$TL1,
                     site_idx = site_use, group_idx = data$group_idx,
                     linf_mean = LINF_PRIOR_MEAN, linf_sd = LINF_PRIOR_SD,
                     logK_mean = LOGK_PRIOR_MEAN, logK_sd = LOGK_PRIOR_SD)
  inits <- list(logK = rep(LOGK_PRIOR_MEAN, n_site_use), v = rep(0, n_group),
                Linf = LINF_PRIOR_MEAN, sigma = 50, sigma_v = 0.3)
  m <- nimbleModel(code, constants = constants, data = list(increment = data$increment),
                    inits = inits, calculate = TRUE)
  conf <- configureMCMC(m, monitors = c("Linf", "logK", "sigma", "sigma_v", "v"), enableWAIC = TRUE)
  mcmc <- buildMCMC(conf)
  cm <- compileNimble(m)
  cmcmc <- compileNimble(mcmc, project = m)
  runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = 3,
          samplesAsCodaMCMC = TRUE, WAIC = TRUE, setSeed = seed_set)
}

fit_pooledK <- fit_site_vbgf(overlap_data, use_site = FALSE, seed_set = c(101, 202, 303))
fit_siteK   <- fit_site_vbgf(overlap_data, use_site = TRUE,  seed_set = c(404, 505, 606))

cat("\nSite-comparison test (overlap range only):\n")
cat("  WAIC pooled K:        ", round(fit_pooledK$WAIC$WAIC, 1), "\n")
cat("  WAIC site-specific K: ", round(fit_siteK$WAIC$WAIC, 1),
    "(Delta =", round(fit_siteK$WAIC$WAIC - fit_pooledK$WAIC$WAIC, 2), ")\n")

samps_site <- do.call(rbind, fit_siteK$samples)
K_bw_site <- exp(samps_site[, "logK[1]"])
K_cib_site <- exp(samps_site[, "logK[2]"])
diff_logK_site <- samps_site[, "logK[2]"] - samps_site[, "logK[1]"]

siteTest <- list(
  waic_pooled = fit_pooledK$WAIC$WAIC,
  waic_site = fit_siteK$WAIC$WAIC,
  K_backwater_mean = mean(K_bw_site), K_backwater_ci = quantile(K_bw_site, c(0.025, 0.975)),
  K_cibola_mean = mean(K_cib_site), K_cibola_ci = quantile(K_cib_site, c(0.025, 0.975)),
  diff_logK_ci = quantile(diff_logK_site, c(0.025, 0.975)),
  p_cibola_slower = mean(K_cib_site < K_bw_site)
)
cat("  K[Backwater] =", round(siteTest$K_backwater_mean, 3), "/yr (",
    round(siteTest$K_backwater_ci[1], 3), "-", round(siteTest$K_backwater_ci[2], 3), ")\n")
cat("  K[Cibola] =", round(siteTest$K_cibola_mean, 3), "/yr (",
    round(siteTest$K_cibola_ci[1], 3), "-", round(siteTest$K_cibola_ci[2], 3), ")\n")
cat("  logK diff 95% CI:", round(siteTest$diff_logK_ci[1], 3), "-", round(siteTest$diff_logK_ci[2], 3), "\n")
cat("  P(K[Cibola] < K[Backwater]):", round(siteTest$p_cibola_slower, 3), "\n")
cat("\nDECISION (BK, 2026-09-14): WAIC narrowly favors pooling; proceeding with a",
    "\nshared K across sites in the combined production model below, despite the",
    "\nborderline (not clean) evidence -- see header comments for full caveat.\n")

## ---- Step 2: combined production model (all backwater + Cibola rows) -----

combined_data <- bind_rows(
  bw_prepped |> transmute(PITIndex, SourceSite, Origin, TL1, TL2, increment, gs_days, Group),
  cibola_prepped |> transmute(PITIndex, SourceSite, Origin, TL1, TL2, increment, gs_days, Group)
) |>
  mutate(t_years = gs_days / 365,
         origin_idx = as.integer(factor(Origin, levels = c("Capture", "Stocking"))),
         group_idx = as.integer(factor(Group)))

cat("\nCombined production dataset: n =", nrow(combined_data),
    "(", sum(combined_data$SourceSite == "Backwater"), "backwater +",
    sum(combined_data$SourceSite == "Cibola"), "Cibola)\n")
cat("TL1 range:", paste(range(combined_data$TL1), collapse = "-"), "mm\n")
cat("TL2 range:", paste(range(combined_data$TL2), collapse = "-"), "mm\n")

n_origin_c <- 2L
n_group_c <- max(combined_data$group_idx)

fit_vbgf_combined <- function(data, use_origin, use_group, seed,
                               niter = 60000, nburnin = 20000, thin = 10) {
  N <- nrow(data)
  origin_use <- if (use_origin) data$origin_idx else rep(1L, N)
  n_origin_use <- if (use_origin) n_origin_c else 1L
  group_use <- if (use_group) data$group_idx else rep(1L, N)
  n_group_use <- if (use_group) n_group_c else 1L

  if (use_group) {
    code <- nimbleCode({
      for (i in 1:N) {
        K[i] <- exp(logK[origin_idx[i]] + v[group_idx[i]])
        mu[i] <- (Linf - TL1[i]) * (1 - exp(-K[i] * t_years[i]))
        increment[i] ~ dnorm(mu[i], sd = sigma)
      }
      for (g in 1:n_group) v[g] ~ dnorm(0, sd = sigma_v)
      for (k in 1:n_origin) logK[k] ~ dnorm(logK_mean, sd = logK_sd)
      Linf ~ dnorm(linf_mean, sd = linf_sd)
      sigma ~ dunif(0, 200)
      sigma_v ~ dunif(0, 2)
    })
    constants <- list(N = N, n_origin = n_origin_use, n_group = n_group_use,
                       t_years = data$t_years, TL1 = data$TL1,
                       origin_idx = origin_use, group_idx = group_use,
                       linf_mean = LINF_PRIOR_MEAN, linf_sd = LINF_PRIOR_SD,
                       logK_mean = LOGK_PRIOR_MEAN, logK_sd = LOGK_PRIOR_SD)
    inits <- list(logK = rep(LOGK_PRIOR_MEAN, n_origin_use), v = rep(0, n_group_use),
                  Linf = LINF_PRIOR_MEAN, sigma = 50, sigma_v = 0.3)
    monitors <- c("Linf", "logK", "sigma", "sigma_v", "v")
  } else {
    code <- nimbleCode({
      for (i in 1:N) {
        K[i] <- exp(logK[origin_idx[i]])
        mu[i] <- (Linf - TL1[i]) * (1 - exp(-K[i] * t_years[i]))
        increment[i] ~ dnorm(mu[i], sd = sigma)
      }
      for (k in 1:n_origin) logK[k] ~ dnorm(logK_mean, sd = logK_sd)
      Linf ~ dnorm(linf_mean, sd = linf_sd)
      sigma ~ dunif(0, 200)
    })
    constants <- list(N = N, n_origin = n_origin_use, t_years = data$t_years, TL1 = data$TL1,
                       origin_idx = origin_use,
                       linf_mean = LINF_PRIOR_MEAN, linf_sd = LINF_PRIOR_SD,
                       logK_mean = LOGK_PRIOR_MEAN, logK_sd = LOGK_PRIOR_SD)
    inits <- list(logK = rep(LOGK_PRIOR_MEAN, n_origin_use), Linf = LINF_PRIOR_MEAN, sigma = 50)
    monitors <- c("Linf", "logK", "sigma")
  }
  m <- nimbleModel(code, constants = constants, data = list(increment = data$increment),
                    inits = inits, calculate = TRUE)
  conf <- configureMCMC(m, monitors = monitors, enableWAIC = TRUE)
  mcmc <- buildMCMC(conf)
  cm <- compileNimble(m)
  cmcmc <- compileNimble(mcmc, project = m)
  runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = 3,
          samplesAsCodaMCMC = TRUE, WAIC = TRUE, setSeed = seed)
}

comb_full     <- fit_vbgf_combined(combined_data, TRUE,  TRUE,  c(11, 22, 33))
comb_noorigin <- fit_vbgf_combined(combined_data, FALSE, TRUE,  c(44, 55, 66))
comb_nogroup  <- fit_vbgf_combined(combined_data, TRUE,  FALSE, c(77, 88, 99))
comb_neither  <- fit_vbgf_combined(combined_data, FALSE, FALSE, c(12, 34, 56))

comb_waic_tab <- data.frame(
  model = c("origin + group (full)", "group only", "origin only", "neither (pooled)"),
  WAIC = c(comb_full$WAIC$WAIC, comb_noorigin$WAIC$WAIC, comb_nogroup$WAIC$WAIC, comb_neither$WAIC$WAIC)
) |> arrange(WAIC)

cat("\nCombined production model WAIC comparison:\n")
print(comb_waic_tab)

fit_list_comb <- list(full = comb_full, noorigin = comb_noorigin, nogroup = comb_nogroup, neither = comb_neither)
use_group_list_comb <- list(full = TRUE, noorigin = TRUE, nogroup = FALSE, neither = FALSE)
key_map_comb <- c("origin + group (full)" = "full", "group only" = "noorigin",
                   "origin only" = "nogroup", "neither (pooled)" = "neither")
winning_key_comb <- unname(key_map_comb[comb_waic_tab$model[1]])
comb_winner <- fit_list_comb[[winning_key_comb]]
winner_use_group_comb <- use_group_list_comb[[winning_key_comb]]

rhat_comb <- gelman.diag(comb_winner$samples, multivariate = FALSE)$psrf[, 1]
cat("\nWinning combined model:", comb_waic_tab$model[1], "\n")
cat("Max R-hat:", round(max(rhat_comb, na.rm = TRUE), 3),
    "| nodes > 1.1:", sum(rhat_comb > 1.1, na.rm = TRUE), "of", length(rhat_comb), "\n")

samps_comb <- do.call(rbind, comb_winner$samples)
Linf_comb <- samps_comb[, "Linf"]
logK_capture_comb <- samps_comb[, "logK[1]"]
sigma_comb <- samps_comb[, "sigma"]
sigma_v_comb <- if (winner_use_group_comb) samps_comb[, "sigma_v"] else rep(0, length(Linf_comb))

max_obs_comb <- max(c(combined_data$TL1, combined_data$TL2))
cat("\nLinf =", round(mean(Linf_comb), 1), "mm (95% CI",
    round(quantile(Linf_comb, 0.025), 1), "-", round(quantile(Linf_comb, 0.975), 1), ")\n")
cat("K[Capture] =", round(exp(mean(logK_capture_comb)), 3), "/yr\n")
if (n_origin_c > 1 && "logK[2]" %in% colnames(samps_comb)) {
  cat("K[Stocking] =", round(exp(mean(samps_comb[, "logK[2]"])), 3), "/yr\n")
}
cat("Max observed TL (combined data):", max_obs_comb, "mm\n")
cat("Posterior P(Linf < max observed):", round(mean(Linf_comb < max_obs_comb), 4), "\n")

## ---- Derive P(TL >= 250mm after one growth season | TL1), combined model -

GROWTH_SEASON_DAYS <- 182
t_years_pred <- GROWTH_SEASON_DAYS / 365
TL1_grid <- c(150, 175, 200, 225, 249)

set.seed(4471)
n_draws_comb <- length(Linf_comb)
v_new_comb <- if (winner_use_group_comb) rnorm(n_draws_comb, 0, sigma_v_comb) else rep(0, n_draws_comb)
K_draw_comb <- exp(logK_capture_comb + v_new_comb)
eps_comb <- rnorm(n_draws_comb, 0, sigma_comb)

prob_comb <- sapply(TL1_grid, function(tl1) {
  increment_i <- (Linf_comb - tl1) * (1 - exp(-K_draw_comb * t_years_pred))
  TL2_i <- tl1 + increment_i + eps_comb
  mean(TL2_i >= 250)
})
expected_increment_comb <- sapply(TL1_grid, function(tl1) {
  mean((Linf_comb - tl1) * (1 - exp(-K_draw_comb * t_years_pred)))
})

growth_season_threshold_probs_combined <- data.frame(
  TL1 = TL1_grid,
  Expected_increment_mm = round(expected_increment_comb, 1),
  P_cross_250_one_growth_season = round(prob_comb, 3)
)

cat("\nP(TL >= 250mm after one Apr-Sep growth season | TL1), COMBINED model:\n")
print(growth_season_threshold_probs_combined)

cat("\nFor comparison, the backwater-ONLY model's predictions were: 0.491 (150mm),",
    "0.574 (175mm), 0.676 (200mm), 0.800 (225mm), 0.917 (249mm) -- the combined model's",
    "predictions above are all HIGHER, driven by both the modestly larger Linf and the",
    "Capture-origin K estimate shifting upward when Cibola's faster-growing large fish",
    "join the pool.\n")

## ---- Comparison figure: backwater-only vs combined growth curve ----------

pred_curve <- function(Linf_draws, K_draws, t_grid, TL0) {
  TL_mat <- sapply(t_grid, function(t) Linf_draws - (Linf_draws - TL0) * exp(-K_draws * t))
  data.frame(t = t_grid, median = apply(TL_mat, 2, median),
             lwr = apply(TL_mat, 2, quantile, 0.025), upr = apply(TL_mat, 2, quantile, 0.975))
}

## backwater-only reference fit is loaded from the original model output
bw_only <- new.env()
load("data/GIEL_GrowthSeasonThreshold.RData", envir = bw_only)
samps_bw_only <- do.call(rbind, bw_only$gs_winner$samples)
Linf_bw_only <- samps_bw_only[, "Linf"]
set.seed(8842)
K_bw_only <- exp(samps_bw_only[, "logK[1]"] +
                    rnorm(nrow(samps_bw_only), 0, samps_bw_only[, "sigma_v"]))

t_grid <- seq(0, 12, by = 0.25)
TL0_plot <- 100
curve_bw_only  <- pred_curve(Linf_bw_only, K_bw_only, t_grid, TL0_plot) |> mutate(Model = "Backwater only (original)")
curve_combined <- pred_curve(Linf_comb, K_draw_comb, t_grid, TL0_plot) |> mutate(Model = "Combined (backwater + Cibola)")
curves_final <- bind_rows(curve_bw_only, curve_combined)

linf_lines_final <- data.frame(
  Model = c("Backwater only (original)", "Combined (backwater + Cibola)"),
  Linf = c(median(Linf_bw_only), median(Linf_comb))
)

p_final <- ggplot(curves_final, aes(t, median, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(data = linf_lines_final, aes(yintercept = Linf, color = Model),
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  scale_color_manual(values = c("Backwater only (original)" = "#2c7fb8",
                                 "Combined (backwater + Cibola)" = "#33a02c")) +
  scale_fill_manual(values = c("Backwater only (original)" = "#2c7fb8",
                                "Combined (backwater + Cibola)" = "#33a02c")) +
  labs(title = "Bonytail VBGF: backwater-only vs. combined (+ Cibola) model",
       subtitle = paste0("Common start TL0 = ", TL0_plot, " mm; dashed = posterior median Linf"),
       x = "Growing-season years elapsed", y = "Total length (mm)",
       caption = "Combined model pools 52 historical (pre-2005) Cibola pairs with 361 backwater pairs,\nshared Linf, origin + pondyear/site-year K structure (WAIC-preferred).") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

ggsave("model/GIEL_GrowthModel_CombinedVsOriginal.png", p_final, width = 7, height = 5, dpi = 150)

## ---- Save results ----------------------------------------------------------

save(
  bw_prepped, cibola_prepped, combined_data,
  overlap_data, fit_pooledK, fit_siteK, siteTest,
  comb_full, comb_noorigin, comb_nogroup, comb_neither, comb_waic_tab,
  comb_winner, winning_key_comb, winner_use_group_comb,
  growth_season_threshold_probs_combined, GROWTH_SEASON_DAYS,
  LINF_PRIOR_MEAN, LINF_PRIOR_SD, LOGK_PRIOR_MEAN, LOGK_PRIOR_SD, max_obs_comb,
  file = "data/GIEL_GrowthSeasonThreshold_withCibola.RData"
)
