## GIEL_GrowthSeasonThreshold.R
##
## Side analysis, general bonytail (GIEL) growth model (2026-09-12 rescope,
## BK request). Originally scoped to only fish first tagged <250 mm TL, to
## directly target the probability of crossing the 250 mm size-class
## threshold. Rescoped per BK: fit a general growth model using capture and
## recapture records across the FULL size spectrum (not just small fish),
## so the growth curve is identified from real curvature (small fish growing
## fast, large fish decelerating) rather than fit to a narrow size window
## where Linf and K cannot be told apart. The 250 mm threshold probability
## is now a downstream prediction from this general model (section 4),
## computed after first locking down the general growth curve.
##
## Data source: data/ReportingData.RData (StudyBWNFWG). Uses only physical
## NFWG capture/stocking records with a measured total_length -- PIT scanner
## contacts do not record length, so this analysis is limited to fish with a
## later physical recapture, not the full known-alive population. Every
## CONSECUTIVE pair of measured records for a given fish x backwater is used
## as a growth-increment observation (not just the first tag -> first
## recapture pair), to make full use of fish with more than one remeasurement
## and to span as much of the size range as possible. This means a small
## number of fish contribute more than one row (314 distinct fish -> 362
## rows); rows from the same fish are not independent, a simplification
## noted as a caveat below.
##
## Growing-season exposure: growth is assumed to occur only during Apr 1 -
## Sep 30 each year (matching the joint robust-design / YCB3S models'
## Apr->Oct between-year gap), computed as the exact calendar overlap of
## [Date1, Date2] with that window in every year spanned. Pairs whose window
## falls entirely within Oct-Mar (gs_days = 0) are dropped -- overwinter
## growth is minimal and these pairs carry no information about the growing-
## season rate. In this dataset gs_days happens to take only 0 or exact
## multiples of 182 (one Apr-Sep season) -- fish are captured/recaptured at
## discrete netting/stocking events that land close to season boundaries, so
## there are no partial-season pairs to worry about; the gs_days > 0 filter
## alone satisfies BK's "must include at least one full growing season, not
## just winter" requirement without a separate minimum-exposure threshold.
##
## Model: Bayesian von Bertalanffy growth-increment (Fabens 1965) model fit
## in NIMBLE:
##   increment = TL2 - TL1 ~ Normal(mu, sigma)
##   mu[i] = (Linf - TL1[i]) * (1 - exp(-K[i] * t_years[i]))
##   K[i]  = exp(logK[origin[i]] + v[pondyear[i]])
##   v[pondyear] ~ Normal(0, sigma_v)          # log-K pond-year deviation
## t_years = gs_days / 365. Fit through the origin (mu = 0 when t_years = 0)
## is automatic under the Fabens form. Priors: Linf ~ Normal(635 mm, sd 50)
## -- weakly informative around the requested 62-65 cm TL asymptotic size;
## logK[k] ~ Normal(log(0.15), sd 0.5) -- weakly informative, centered near
## the requested 0.1-0.2/yr range (median K = 0.15/yr, ~90% prior interval
## roughly 0.06-0.4/yr). Pond and year are combined into a single PondYear
## grouping factor acting multiplicatively on K (log scale) as a first step
## (per BK: lock down the general Linf/K first, then assess whether pond
## and/or year matter, and whether they can be separated -- the expanded
## data now have a real, if unbalanced, pond x year cross-tab: Pond 2 is
## almost entirely 2017, Pond 5 is 2017 and 2022 (its two stocking cohorts),
## Pond 6 spans 2017-2022, so full pond x year separation is left for a
## follow-up once this general model is validated). Stepwise term reduction
## (origin, PondYear) is by WAIC across the same four nested structures used
## previously (full, origin only, pondyear only, neither).
##
## TL1 now spans 117-442 mm (vs 117-249 mm in the original <250 mm-only
## scope) -- still short of the ~635 mm Linf prior, but with real curvature
## contrast between small and mid-size fish. See section 3 output / AGENTS.md
## for whether this resolves the earlier Linf/K identifiability problem.
##
## Linf-below-observed-max-size investigation (2026-09-12, see AGENTS.md for
## full writeup): the winning model's Linf (418 mm) sits below several
## observed post-growth lengths (max TL2 = 462 mm) -- every posterior draw
## has Linf < max(TL2). Diagnosed as a data-sparsity / Fabens-bias problem,
## not fixed by restructuring: (1) only 3 of 361 rows have TL1 > 400 mm (all
## Pond 6, 2022), and those 3 still grow 17-42 mm/season, pulling against
## the hundreds of small/mid-size rows that want a much smaller apparent
## asymptote; (2) giving each pond its own Linf (tested, not kept in the
## pipeline) does not help -- Pond 6 alone converges to an even lower Linf
## (411 mm) with a tied WAIC, because the same small-vs-large-fish tension
## exists within that one pond; (3) forcing Linf toward the literature range
## with a tight prior (sd 15 instead of 50, tested, not kept) costs ~149
## WAIC and produces large systematic residuals (small fish underpredicted
## by ~50 mm, mid-size fish overpredicted by ~17 mm) -- a much worse overall
## fit. This matches the well-documented Fabens estimator bias toward low
## Linf / high K when individual fish vary in their own growth trajectories
## and few recaptures exist near the true asymptote. Conclusion: the fitted
## Linf/K here should not be read as literal biological growth-curve
## parameters; the model is retained (loose prior, single global Linf) only
## for its interpolated increment predictions within the well-populated size
## range (roughly 117-350 mm, where residuals are small), not for
## extrapolation near or beyond 400+ mm.
##
## STATUS (2026-09-12, BK decision): moving forward reporting this model's
## P(TL >= 250 mm | TL1) predictions as final (section 4 / growth_season_
## threshold_probs) -- no further Linf/K restructuring is planned. Any
## report or document that cites these predictions must repeat the Linf/K
## identifiability caveat above (see also model/GIEL_GrowthSeasonThreshold_
## Summary.md section 0c, which holds the adopted results table plus the
## required caveat language).
##
## Caveat: recapture pairs are an opportunistic subset of physically
## recaptured fish, not a random resample of all tagged/stocked GIEL at
## these 3 ponds -- if capture probability correlates with fish size or
## growth rate, these estimates may be biased. Consecutive-pair rows from
## the same fish are treated as independent observations (a common
## simplification for Fabens-style fitting with limited data; a full
## individual-random-effects growth model, e.g. Laslett/GROTAG-style, would
## be more rigorous but is not attempted here).
##
## Negative-growth pairs (TL2 < TL1) are dropped before fitting -- these are
## almost certainly measurement/transcription errors (fish do not shrink),
## not real biology, and would otherwise bias sigma upward and distort the
## Fabens fit. Every dropped pair is written to
## model/GIEL_GrowthSeasonThreshold_Summary.md as a flagged-record list for
## database correction (see NegativeGrowthPairs below).

library(tidyverse)
library(lubridate)
library(nimble)
library(coda)

load("data/ReportingData.RData")

GIEL_PONDS <- c("IPCA (Pond 2)", "IPCA (Pond 5)", "IPCA (Pond 6)")

## ---- 1. Build recapture pairs (all sizes, all consecutive remeasurements) -

giel_recs <- StudyBWNFWG |>
  filter(
    species == "GIEL", Backwater %in% GIEL_PONDS,
    event %in% c("stocking", "capture"),
    !is.na(total_length), !is.na(collection_date)
  ) |>
  select(PITIndex, Backwater, event, collection_date, total_length) |>
  arrange(PITIndex, Backwater, collection_date)

cat("Qualifying physical NFWG records (any size), GIEL, IP2/IP5/IP6:", nrow(giel_recs), "\n")

# every consecutive pair of measured records per fish x backwater (not just
# first tag -> first recapture) -- spans the full size range and uses fish
# with >1 remeasurement more fully
recap_pairs <- giel_recs |>
  group_by(PITIndex, Backwater) |>
  filter(n() >= 2) |>
  arrange(collection_date, .by_group = TRUE) |>
  mutate(
    Date1 = collection_date, TL1 = total_length, Origin = event,
    Date2 = lead(collection_date), TL2 = lead(total_length)
  ) |>
  ungroup() |>
  filter(!is.na(Date2)) |>
  transmute(
    PITIndex, Pond = Backwater, Origin,
    Date1, TL1, Date2, TL2,
    days_elapsed = as.numeric(Date2 - Date1)
  ) |>
  arrange(Pond, Date1)

cat("Fish x backwater with >=2 measured records:",
    giel_recs |> count(PITIndex, Backwater) |> filter(n >= 2) |> nrow(), "\n")
cat("Consecutive-pair growth-increment rows (all sizes, any gap length):", nrow(recap_pairs), "\n")
cat("TL1 range across all pairs:", paste(range(recap_pairs$TL1), collapse = "-"), "mm\n")

## ---- 2. Growing-season (Apr 1 - Sep 30) exposure ------------------------

growing_season_days <- function(d1, d2) {
  mapply(function(a, b) {
    yrs <- year(a):year(b)
    total <- 0
    for (yr in yrs) {
      gs_start <- ymd(sprintf("%d-04-01", yr))
      gs_end <- ymd(sprintf("%d-09-30", yr))
      os <- max(a, gs_start)
      oe <- min(b, gs_end)
      if (oe >= os) total <- total + as.numeric(oe - os)
    }
    total
  }, d1, d2)
}

gs_data <- recap_pairs |>
  mutate(
    Pond = str_extract(Pond, "\\d+"),
    Origin = str_to_title(Origin),
    increment = TL2 - TL1,
    gs_days = growing_season_days(Date1, Date2),
    PondYear = paste0("P", Pond, "_", year(Date1))
  )

## Flag and drop negative-growth pairs (TL2 < TL1) -- almost certainly
## measurement/transcription error, not real shrinkage. Kept as a separate
## object for database review rather than silently discarded.
NegativeGrowthPairs <- gs_data |>
  filter(increment < 0) |>
  select(PITIndex, Pond, Origin, Date1, TL1, Date2, TL2, increment, days_elapsed, gs_days) |>
  arrange(Pond, Date1)

cat("Negative-growth pairs flagged and excluded (TL2 < TL1):", nrow(NegativeGrowthPairs), "\n")
if (nrow(NegativeGrowthPairs) > 0) print(NegativeGrowthPairs)

gs_data <- gs_data |> filter(increment >= 0)

gs_fit_data <- gs_data |>
  filter(gs_days > 0) |>                       # drop winter-only pairs
  mutate(
    origin_idx = as.integer(factor(Origin, levels = c("Capture", "Stocking"))),
    py_idx = as.integer(factor(PondYear))
  )

cat("Winter-only pairs excluded (gs_days = 0):", nrow(gs_data) - nrow(gs_fit_data), "\n")
cat("Rows used for growing-season model:", nrow(gs_fit_data), "\n")

n_origin <- length(unique(gs_fit_data$origin_idx))
n_py <- max(gs_fit_data$py_idx)
py_levels <- levels(factor(gs_fit_data$PondYear))

## ---- 3. Bayesian von Bertalanffy (Fabens) increment model ----------------
##
## Priors (weakly informative, per BK request 2026-09-12):
##   Linf ~ Normal(635 mm, sd 50)         (62-65 cm TL asymptotic size)
##   logK[k] ~ Normal(log(0.15), sd 0.5)  (median K = 0.15/yr; ~90% interval
##                                          roughly 0.06-0.4/yr)

LINF_PRIOR_MEAN <- 635   # mm; midpoint of requested 62-65 cm
LINF_PRIOR_SD   <- 50
LOGK_PRIOR_MEAN <- log(0.15)  # per year; midpoint of requested 0.1-0.2/yr
LOGK_PRIOR_SD   <- 0.5

fit_vbgf_nimble <- function(t_years, TL1, origin_idx, py_idx, increment,
                             use_origin = TRUE, use_group = TRUE,
                             n_origin, n_py, seed = 4127,
                             niter = 60000, nburnin = 20000, thin = 10) {
  N <- length(increment)
  origin_use <- if (use_origin) origin_idx else rep(1L, N)
  n_origin_use <- if (use_origin) n_origin else 1L
  py_use <- if (use_group) py_idx else rep(1L, N)
  n_py_use <- if (use_group) n_py else 1L

  if (use_group) {
    code <- nimbleCode({
      for (i in 1:N) {
        K[i] <- exp(logK[origin_idx[i]] + v[py_idx[i]])
        mu[i] <- (Linf - TL1[i]) * (1 - exp(-K[i] * t_years[i]))
        increment[i] ~ dnorm(mu[i], sd = sigma)
      }
      for (j in 1:n_py) {
        v[j] ~ dnorm(0, sd = sigma_v)
      }
      for (k in 1:n_origin) {
        logK[k] ~ dnorm(logK_mean, sd = logK_sd)
      }
      Linf ~ dnorm(linf_mean, sd = linf_sd)
      sigma ~ dunif(0, 200)
      sigma_v ~ dunif(0, 2)
    })
    constants <- list(N = N, n_origin = n_origin_use, n_py = n_py_use,
                       t_years = t_years, TL1 = TL1,
                       origin_idx = origin_use, py_idx = py_use,
                       linf_mean = LINF_PRIOR_MEAN, linf_sd = LINF_PRIOR_SD,
                       logK_mean = LOGK_PRIOR_MEAN, logK_sd = LOGK_PRIOR_SD)
    inits <- list(logK = rep(LOGK_PRIOR_MEAN, n_origin_use), v = rep(0, n_py_use),
                  Linf = LINF_PRIOR_MEAN, sigma = 50, sigma_v = 0.3)
    monitors <- c("Linf", "logK", "sigma", "sigma_v", "v")
  } else {
    code <- nimbleCode({
      for (i in 1:N) {
        K[i] <- exp(logK[origin_idx[i]])
        mu[i] <- (Linf - TL1[i]) * (1 - exp(-K[i] * t_years[i]))
        increment[i] ~ dnorm(mu[i], sd = sigma)
      }
      for (k in 1:n_origin) {
        logK[k] ~ dnorm(logK_mean, sd = logK_sd)
      }
      Linf ~ dnorm(linf_mean, sd = linf_sd)
      sigma ~ dunif(0, 200)
    })
    constants <- list(N = N, n_origin = n_origin_use, t_years = t_years, TL1 = TL1,
                       origin_idx = origin_use,
                       linf_mean = LINF_PRIOR_MEAN, linf_sd = LINF_PRIOR_SD,
                       logK_mean = LOGK_PRIOR_MEAN, logK_sd = LOGK_PRIOR_SD)
    inits <- list(logK = rep(LOGK_PRIOR_MEAN, n_origin_use), Linf = LINF_PRIOR_MEAN, sigma = 50)
    monitors <- c("Linf", "logK", "sigma")
  }

  m <- nimbleModel(code, constants = constants, data = list(increment = increment),
                    inits = inits, calculate = TRUE)
  conf <- configureMCMC(m, monitors = monitors, enableWAIC = TRUE)
  mcmc <- buildMCMC(conf)
  cm <- compileNimble(m)
  cmcmc <- compileNimble(mcmc, project = m)
  runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin,
          nchains = 2, samplesAsCodaMCMC = TRUE, WAIC = TRUE,
          setSeed = c(seed, seed + 1))
}

rhat <- function(fit) gelman.diag(fit$samples, multivariate = FALSE)$psrf[, 1]
summ <- function(fit, pars) round(summary(fit$samples)$statistics[pars, c("Mean", "SD")], 3)

t_years_gs <- gs_fit_data$gs_days / 365
TL1_gs <- gs_fit_data$TL1

gs_full     <- fit_vbgf_nimble(t_years_gs, TL1_gs, gs_fit_data$origin_idx, gs_fit_data$py_idx, gs_fit_data$increment,
                                use_origin = TRUE,  use_group = TRUE,
                                n_origin = n_origin, n_py = n_py, seed = 5511)
gs_noorigin <- fit_vbgf_nimble(t_years_gs, TL1_gs, gs_fit_data$origin_idx, gs_fit_data$py_idx, gs_fit_data$increment,
                                use_origin = FALSE, use_group = TRUE,
                                n_origin = n_origin, n_py = n_py, seed = 8127)
gs_nogroup  <- fit_vbgf_nimble(t_years_gs, TL1_gs, gs_fit_data$origin_idx, gs_fit_data$py_idx, gs_fit_data$increment,
                                use_origin = TRUE,  use_group = FALSE,
                                n_origin = n_origin, n_py = n_py, seed = 2298)
gs_neither  <- fit_vbgf_nimble(t_years_gs, TL1_gs, gs_fit_data$origin_idx, gs_fit_data$py_idx, gs_fit_data$increment,
                                use_origin = FALSE, use_group = FALSE,
                                n_origin = n_origin, n_py = n_py, seed = 7743)

gs_waic_tab <- data.frame(
  model = c("origin + pondyear (full)", "pondyear only", "origin only", "neither (fully pooled)"),
  WAIC = c(gs_full$WAIC$WAIC, gs_noorigin$WAIC$WAIC, gs_nogroup$WAIC$WAIC, gs_neither$WAIC$WAIC)
) |> arrange(WAIC)

print(gs_waic_tab)

## Select winning model automatically by lowest WAIC
fit_list <- list(full = gs_full, noorigin = gs_noorigin, nogroup = gs_nogroup, neither = gs_neither)
use_group_list <- list(full = TRUE, noorigin = TRUE, nogroup = FALSE, neither = FALSE)
winning_name <- gs_waic_tab$model[1]
key_map <- c(
  "origin + pondyear (full)" = "full", "pondyear only" = "noorigin",
  "origin only" = "nogroup", "neither (fully pooled)" = "neither"
)
winning_key <- unname(key_map[winning_name])
gs_winner <- fit_list[[winning_key]]
winner_use_group <- use_group_list[[winning_key]]

cat("\nWinning model (lowest WAIC):", winning_name, "\n")
cat("\nR-hat, winning model:\n")
print(round(rhat(gs_winner), 3))
cat("\nParameter summary (Linf, logK[1], sigma", if (winner_use_group) ", sigma_v" else "", "):\n", sep = "")
pars_report <- c("Linf", "logK[1]", "sigma")
if (winner_use_group) pars_report <- c(pars_report, "sigma_v")
par_summ_win <- summ(gs_winner, pars_report)
print(par_summ_win)
cat("Implied K[1] (exp(mean logK[1])):",
    round(exp(par_summ_win["logK[1]", "Mean"]), 3), "/yr\n")

if (winner_use_group) {
  v_summ <- summ(gs_winner, paste0("v[", seq_len(n_py), "]"))
  rownames(v_summ) <- py_levels
  cat("\nPond-year log-K deviations (v):\n")
  print(v_summ)
}

## ---- 4. Derive P(TL >= 250mm after one growth season | TL1) -------------

samps_gs <- do.call(rbind, gs_winner$samples)
Linf_gs <- samps_gs[, "Linf"]
logK_gs <- samps_gs[, "logK[1]"]
sigma_gs <- samps_gs[, "sigma"]
sigma_v_gs <- if (winner_use_group) samps_gs[, "sigma_v"] else rep(0, length(Linf_gs))

GROWTH_SEASON_DAYS <- 182  # Apr 1 - Sep 30
t_years_pred <- GROWTH_SEASON_DAYS / 365

set.seed(9931)
n_draws <- length(Linf_gs)
v_new <- if (winner_use_group) rnorm(n_draws, 0, sigma_v_gs) else rep(0, n_draws)
K_draw <- exp(logK_gs + v_new)
eps <- rnorm(n_draws, 0, sigma_gs)

TL1_grid <- c(150, 175, 200, 225, 249)
prob_gs <- sapply(TL1_grid, function(tl1) {
  increment_i <- (Linf_gs - tl1) * (1 - exp(-K_draw * t_years_pred))
  TL2_i <- tl1 + increment_i + eps
  mean(TL2_i >= 250)
})
expected_increment <- sapply(TL1_grid, function(tl1) {
  mean((Linf_gs - tl1) * (1 - exp(-K_draw * t_years_pred)))
})

growth_season_threshold_probs <- data.frame(
  TL1 = TL1_grid,
  Expected_increment_mm = round(expected_increment, 1),
  P_cross_250_one_growth_season = round(prob_gs, 3)
)

cat("\nP(TL >= 250mm after one Apr-Sep growth season | TL1), von Bertalanffy model:\n")
print(growth_season_threshold_probs)

cat("\nPosterior: Linf =", round(mean(Linf_gs), 1), "mm (sd", round(sd(Linf_gs), 1),
    "), K =", round(mean(K_draw), 3), "/yr (sd", round(sd(K_draw), 3), ")\n")

cat("\nCAVEAT (adopted 2026-09-12, see header comments / AGENTS.md): Linf/K above are a",
    "\nfitted compromise, not literal growth-curve parameters (only 3 of", nrow(gs_fit_data),
    "\nrows have TL1 > 400 mm). The P(TL>=250mm) predictions above are valid because TL1 =",
    "\n150-249 mm is well inside the well-populated, low-residual size range this model was",
    "\nchecked against -- repeat this caveat in any report that cites them.\n")

## ---- 5. Save results ------------------------------------------------------

save(
  recap_pairs, gs_data, gs_fit_data, NegativeGrowthPairs,
  gs_full, gs_noorigin, gs_nogroup, gs_neither, gs_waic_tab,
  gs_winner, winning_name, winner_use_group,
  growth_season_threshold_probs, GROWTH_SEASON_DAYS,
  LINF_PRIOR_MEAN, LINF_PRIOR_SD, LOGK_PRIOR_MEAN, LOGK_PRIOR_SD,
  file = "data/GIEL_GrowthSeasonThreshold.RData"
)
