## GIEL_VBGF_AgeLength.R
##
## Bonytail (GIEL) age-based von Bertalanffy growth model (2026-09-14, BK
## request). Builds a real age-length dataset -- as opposed to the size-only
## growth-increment model in GIEL_GrowthSeasonThreshold.R -- from two age
## sources, then refits the VBGF on absolute age instead of paired growth
## increments.
##
## ============================================================================
## REVISION 2026-09-14 (later same day): SHARED-K, RELEASE-ANCHORED MODEL
## ============================================================================
## BK hypothesized that stocked (hatchery-reared) and wild-recruit bonytail
## grow at the SAME rate once both are in the pond, and that the large
## origin-specific K difference in the original age-based fit (K[Stocked]
## 0.221/yr vs K[Recruit] 0.540/yr) was an artifact of blending ~5.2 years of
## slow captive/hatchery growth with subsequent pond growth into one
## continuous hatch-anchored curve for stocked fish, rather than a real
## difference in pond growth rate.
##
## Direct empirical test (see conversation log / GIEL_VBGF_AgeLength_Summary.md
## section "Natural-experiment test"): IP5 (Pond 5) in 2022-23 is the ONLY
## pond-year cell in the companion growth-increment dataset
## (GIEL_GrowthSeasonThreshold.R's gs_fit_data) where both a stocked cohort
## and previously-established (recruit-origin) fish have one-growing-season
## recapture pairs. Controlling for starting TL, the origin effect on
## increment was small and not significant (-18.3 mm/season, SE 11.3,
## p = 0.11, n = 143) -- much smaller than the ~2.4x K ratio implied by the
## original hatch-anchored age model. This supports the hypothesis: once in
## the pond, growth looks the same regardless of origin, and the origin
## effect elsewhere is confounded with pond-year (almost every other
## pond-year cell is single-origin).
##
## FIX: stop assuming hatch-to-now is one continuous growth regime for
## stocked fish. Anchor stocked-fish records at their own RELEASE size/date
## (a Fabens-style reset), not at an assumed hatch length/date, and share a
## single K (with the existing pond-year/cohort random effect v[j]) across
## origins -- i.e., DROP the origin-specific logK term entirely from the
## primary model:
##   Recruit-origin row i:  mu_i = L0 + (Linf - L0) * (1 - exp(-K_i * age_i))
##     (age_i = years since assumed hatch date, unchanged from before; L0
##      anchored on Hamman 1982 as before)
##   Stocked-origin row i:  mu_i = TLrel_i + (Linf - TLrel_i) * (1 - exp(-K_i * t_i))
##     (t_i = years since that fish's own release-anchor date; TLrel_i = that
##      fish's release-anchor TL, both DATA, not modeled)
##   K_i = exp(logK + v[pondyear_i])  -- single shared logK, no origin index
##
## Release anchor caveat: the true stocking record has a recorded TL for the
## 2017-03-21 cohort (year_class 2012, 899 fish; TL 234-299 mm at t=0, exact
## release anchor) but NOT for the 2023-04-12 cohort (year_class 2018, 300
## fish -- no total_length recorded at the stocking event itself). For that
## cohort the release anchor falls back to each fish's first later record
## with a measured TL (in practice, Dec 2023, ~244 days after release, TL
## 297-406 mm, n = 80 of 300 fish) -- i.e., growth between true release and
## that first measurement is not observed and cannot be attributed. This
## degrades but does not eliminate the design: all 80 of these fish
## contribute only their own anchor row (no repeat measurements exist for
## this recent cohort yet), so they add no K information at all under this
## scheme (an anchor row is a tautology: mu = TLrel_i regardless of K, since
## at t = 0 the exponential term is 0) -- they simply carry no weight, rather
## than misleadingly informing a blended hatchery+pond rate as before. All
## real K information for stocked-origin fish now comes from the 899-fish
## 2012 cohort's 97 repeat-measurement rows (89 of 899 fish were recaptured
## at least once after release: 81 fish x 1 extra measurement + 8 fish x 2).
## This is a much smaller informative sample than the previous hatch-anchored
## approach appeared to use (n = 1076 total rows) -- but that apparent
## richness was the problem, not a feature: it was informative about the
## hatch-to-now AVERAGE rate, not the pond rate this model now targets.
##
## Side benefit: releasing the hatch-age assumption for stocked fish also
## resolves most of the original age/origin confound (recruit ages 0.66-3.94
## yr had NO overlap with stocked-cohort hatch-ages 4.97-10.7 yr). Time since
## release for the 2012 cohort now spans 0-9.4 yr and for the 2018 cohort
## 0-3.4 yr, both of which overlap the recruit age range.
##
## RESULT: see model/GIEL_VBGF_AgeLength_Summary.md for the refit posterior,
## a sensitivity check re-adding an origin term on top of release-anchoring
## (to confirm the origin effect on K is now negligible on the full
## dataset, not just the one IP5 cell), and WAIC comparison. The prior
## origin-differentiated, hatch-anchored fit is archived at
## data/GIEL_VBGF_AgeLength_originDiff.RData for comparison.
##
## Age sources:
##   1. Known stocking-cohort year class. StudyBWNFWG$year_class is populated
##      for two GIEL stockings: 2017-03-21 (899 records, all 3 ponds, year
##      class 2012) and 2023-04-12 (300 records, Pond 2 only, year class
##      2018). Both cohorts were ~5.2 years old at release (long captive
##      rearing before stocking, typical of bonytail hatchery production).
##      The third stocking (2022-03-09, 300 fish, Pond 5) has NO recorded
##      year class -- a gap flagged for database follow-up, not fixed here.
##      year_class is already propagated in the source data to every later
##      NFWG record (stocking + recapture) for these PITs.
##   2. Assumed young-of-year (YOY) fall recruits. Untagged fish first
##      captured Oct-Dec (no year_class anywhere in their record history)
##      show real bimodal size distributions in several fiscal years,
##      consistent with a distinguishable small (age-0) cohort vs older fish.
##      A 2-component Gaussian mixture (mclust::Mclust, G = 1:2 by BIC) is
##      fit per FY on log(TL) among these first-capture untagged fall fish;
##      when G = 2 is selected, fish classified to the lower-mean cluster are
##      assigned year_class = calendar year of capture (assumed spring birth
##      that same year). FY2018 has too few fish (n = 9) for a reliable
##      mixture fit and is excluded (no YOY assignment that year); FY2021 has
##      no fall captures at all (consistent with the documented IP2 2022
##      die-off's preceding population crash). FY2019's mixture fit is
##      marginal (mean classification uncertainty 0.102, vs <0.09 for all
##      other fitted years) and is included but should be read with slightly
##      less confidence than FY2020/2022/2023/2024.
##
## Both sources are combined into one year_class map (no PIT overlaps
## between them by construction) and used to compute age at every NFWG
## record (stocking + capture) for every fish carrying either a real or an
## assumed year class: age_years = (collection_date - April 1 of year_class)
## / 365.25. April 1 is used as a nominal spawning-date convention for
## consistency with this project's existing Apr-Sep "growing season"
## boundary (used throughout the joint robust-design / YCB3S models); it is
## an assumption, not a literature-sourced bonytail spawning date.
##
## Resulting dataset: 1,601 age-length rows (525 from assumed-YOY recruits,
## age 0.66-3.94 yr; 1,076 from known stocking cohorts). 1,365 fish
## contribute exactly one row, 111 contribute 2-3 (repeat measurements are
## treated as independent rows, the same simplification used in
## GIEL_GrowthSeasonThreshold.R).
##
## SUPERSEDED CAVEAT (hatch-anchored version only, see revision note above):
## age and origin were almost completely confounded by construction
## (assumed-recruit ages all < 4 yr, known-cohort hatch-ages all > 4.97 yr).
## The release-anchored revision uses time-since-release (not time-since-
## hatch) for stocked fish, which mostly resolves this (see revision note).
##
## Priors: Linf ~ Normal(635, 50); logK ~ Normal(log(0.15), 0.5) (single
## shared logK in the primary model; logK[origin] in the sensitivity
## variant); L0 ~ Normal(6.8, 0.25) mm, anchored on Hamman (1982) mean
## bonytail hatch length (6.8 mm; stated range 6.5-7.5 mm) -- BK flagged
## this literature source (Hamman, R.L. 1982. "Induced Spawning and Culture
## of Bonytail Chub." The Progressive Fish-Culturist 44(4):201-203,
## doi:10.1577/1548-8659(1982)44[201:ISACOB]2.0.CO;2) after an earlier draft
## of this script incorrectly assumed no such value existed; this mirrors
## the gray triggerfish paper's L0-anchoring approach (see
## documents/Chamberlin2025_GrayTriggerfishVBGM_Summary.md).
##
## Sampler note (carried over from the hatch-anchored version): Linf and
## logK are correlated in the posterior via a curved ridge (both enter
## through exp() and a product) that a linear-Gaussian RW_block handles
## poorly (R-hat 1.3-1.4 after 80k iter) but a gradient-informed **barker**
## sampler resolves well (R-hat <= 1.06); requires `buildDerivs = TRUE`.
##
## RESULT (t0 version, superseded): WAIC favored the full (origin +
## pondyear) model (14506.9); Linf = 703 mm (95% CI 630-783); t0 = -0.69 yr;
## K[StockedCohort] = 0.115/yr, K[AssumedRecruit] = 0.20/yr.
##
## RESULT (L0-anchored, hatch-age, origin-differentiated -- superseded by
## the release-anchored revision below): WAIC favored the full (origin +
## pondyear) model (14561); Linf = 486 mm (95% CI 464-509); K[StockedCohort]
## = 0.221/yr (0.178-0.260), K[AssumedRecruit] = 0.540/yr (0.443-0.622).
## Archived at data/GIEL_VBGF_AgeLength_originDiff.RData.
##
## RESULT (release-anchored, shared K -- current, 2026-09-14): CONFIRMS the
## hypothesis on the full dataset. WAIC is a statistical tie between the
## primary shared-K model (13646.35) and a sensitivity check re-adding an
## origin-specific logK (13646.87, Delta = 0.52 -- noise) -- a dramatic
## reversal from the hatch-anchored version's ~800-unit WAIC gap favoring
## origin. In the origin-check model, implied K is nearly identical between
## origins: StockedCohort 0.718/yr (0.47-1.01), AssumedRecruit 0.768/yr
## (0.63-0.91), heavily overlapping, P(K[Stocked] < K[Recruit]) = 0.63 (~coin
## flip). The shared-K primary model gives Linf = 406.1 mm (95% CI
## 395.9-416.7, much tighter than the hatch-anchored version's 464-509) and
## K = 0.792/yr (0.68-0.92); residual sigma dropped from 22.7 mm to 17.0 mm
## (a better-fitting model with FEWER parameters). CAVEAT: this Linf is
## BELOW the observed max TL (442 mm) -- the maximum posterior draw across
## all chains is 427.5 mm, i.e. 0% posterior support for Linf >= 442 -- the
## same identifiability failure seen in the increment (Fabens) model
## (GIEL_GrowthSeasonThreshold_Summary.md section 0b), now reappearing here
## despite the much better origin/WAIC story. Likely cause: real growth
## information for large/old fish is thin under release-anchoring (only 89
## of 979 stocked fish have any post-release repeat measurement; recruit
## data tops out at 408 mm). Read Linf here as a fitted compromise, not a
## literal asymptote, exactly as already cautioned for the increment model.
## See model/GIEL_VBGF_AgeLength_Summary.md for the full writeup.

library(tidyverse)
library(lubridate)
library(mclust)
library(nimble)
library(coda)

load("data/ReportingData.RData")

GIEL_PONDS <- c("IPCA (Pond 2)", "IPCA (Pond 5)", "IPCA (Pond 6)")

## ---- 1. Known stocking-cohort year classes --------------------------------

giel_nfwg <- StudyBWNFWG |> filter(species == "GIEL", Backwater %in% GIEL_PONDS)

known_yc <- giel_nfwg |>
  filter(!is.na(year_class)) |>
  distinct(PITIndex, year_class) |>
  group_by(PITIndex) |>
  summarise(year_class = as.integer(first(year_class)), .groups = "drop")

cat("Known-year-class PITs:", nrow(known_yc), "\n")
print(known_yc |> dplyr::count(year_class))

## ---- 2. Assumed age-0 (YOY) fall recruits via per-FY mixture model --------

recruit_candidates <- giel_nfwg |>
  filter(FirstRecord == "yes", event == "capture", month(collection_date) %in% 10:12,
         !(PITIndex %in% known_yc$PITIndex), !is.na(total_length)) |>
  mutate(cal_year = year(collection_date), fy = cal_year + 1)

cat("\nCandidate untagged fall first-captures (no known year class):", nrow(recruit_candidates), "\n")

mix_groups <- recruit_candidates |> group_by(fy) |> group_split()
yoy_rows <- list()
mixture_log <- list()
for (grp in mix_groups) {
  fy_val <- grp$fy[1]
  x <- log(grp$total_length)
  if (length(x) < 15) {
    mixture_log[[as.character(fy_val)]] <- data.frame(fy = fy_val, n = length(x), decision = "too few fish (excluded)")
    next
  }
  m <- Mclust(x, G = 1:2, verbose = FALSE)
  if (m$G == 1) {
    mixture_log[[as.character(fy_val)]] <- data.frame(fy = fy_val, n = length(x), decision = "unimodal (excluded)")
    next
  }
  low_idx <- which(m$parameters$mean == min(m$parameters$mean))
  is_low <- m$classification == low_idx
  yoy_rows[[as.character(fy_val)]] <- grp[is_low, ] |>
    mutate(assigned_year_class = cal_year, mix_unc = m$uncertainty[is_low])
  mixture_log[[as.character(fy_val)]] <- data.frame(
    fy = fy_val, n = length(x), decision = "included (2-component mixture)",
    n_assigned_yoy = sum(is_low), mean_uncertainty = round(mean(m$uncertainty[is_low]), 3)
  )
}

MixtureClassificationLog <- bind_rows(mixture_log)
cat("\nPer-FY mixture classification decisions:\n")
print(MixtureClassificationLog)

AssumedYOY <- bind_rows(yoy_rows) |> as_tibble() |>
  select(PITIndex, Backwater, collection_date, total_length, cal_year, mix_unc)

cat("\nTotal fish assigned assumed age-0 (YOY) status:", nrow(AssumedYOY), "\n")

assumed_yc_map <- AssumedYOY |> distinct(PITIndex, year_class = cal_year)

## ---- 3. Combine into one year-class map and build age-length dataset -----

YearClassMap <- bind_rows(
  known_yc |> mutate(source = "known_stocking_cohort"),
  assumed_yc_map |> mutate(source = "assumed_YOY_fall_capture")
)

stopifnot(length(intersect(known_yc$PITIndex, assumed_yc_map$PITIndex)) == 0)

AgeLength <- giel_nfwg |>
  select(-year_class) |>
  inner_join(YearClassMap, by = "PITIndex") |>
  filter(!is.na(total_length)) |>
  mutate(
    birth_date = as.Date(paste0(year_class, "-04-01")),
    age_years = as.numeric(collection_date - birth_date) / 365.25
  )

stopifnot(all(AgeLength$age_years >= 0))

cat("\nTotal age-length observations:", nrow(AgeLength), "\n")
cat("Age range by source:\n")
print(AgeLength |> group_by(source) |> summarise(n = dplyr::n(), min_age = min(age_years), max_age = max(age_years)))

## ---- 3b. Release anchors for stocked-cohort fish --------------------------
##
## Each stocked-cohort fish's growth is now anchored at its own release
## record (TL + date), not at an assumed hatch date/length. The true
## stocking-event record has a measured TL for the 2012 year-class cohort
## (2017-03-21, 899 fish); it does NOT for the 2018 year-class cohort
## (2023-04-12, 300 fish -- no total_length recorded at that stocking
## event), so for that cohort the anchor falls back to each fish's own
## first later record with a measured TL (see revision note at top of file).
## Either way the anchor is simply "this fish's own earliest measured
## record" -- computed uniformly below.

release_anchor <- AgeLength |>
  filter(source == "known_stocking_cohort") |>
  group_by(PITIndex) |>
  slice_min(collection_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(PITIndex, release_date = collection_date, release_TL = total_length)

cat("\nRelease anchors built for", nrow(release_anchor), "of", nrow(known_yc), "known-cohort fish\n")
cat("Release anchor TL/date by year_class:\n")
print(AgeLength |> filter(source == "known_stocking_cohort") |>
        distinct(PITIndex, year_class) |>
        inner_join(release_anchor, by = "PITIndex") |>
        group_by(year_class) |> summarise(n = dplyr::n(), min_TL = min(release_TL), max_TL = max(release_TL),
                                           min_date = min(release_date), max_date = max(release_date)))

VBGF_data <- AgeLength |>
  left_join(release_anchor, by = "PITIndex") |>
  mutate(
    Pond = str_extract(Backwater, "\\d+"),
    AgeSource = if_else(source == "known_stocking_cohort", "StockedCohort", "AssumedRecruit"),
    PondYear = paste0("P", Pond, "_", year_class),
    origin_idx = as.integer(factor(AgeSource, levels = c("StockedCohort", "AssumedRecruit"))),
    py_idx = as.integer(factor(PondYear)),
    is_recruit = if_else(AgeSource == "AssumedRecruit", 1L, 0L),
    t_since_release = if_else(AgeSource == "StockedCohort",
                               as.numeric(collection_date - release_date) / 365.25, NA_real_),
    ## t = growth-clock time: years since hatch (recruit) or years since
    ## release (stocked); TL_release_const = release TL (stocked) or 0
    ## (recruit, unused placeholder -- zeroed out by is_recruit in the model)
    t = if_else(AgeSource == "StockedCohort", t_since_release, age_years),
    TL_release_const = if_else(AgeSource == "StockedCohort", release_TL, 0)
  )

stopifnot(all(VBGF_data$t >= -1e-9, na.rm = TRUE))
stopifnot(all(!is.na(VBGF_data$t)))

n_origin <- length(unique(VBGF_data$origin_idx))
n_py <- max(VBGF_data$py_idx)
py_levels <- levels(factor(VBGF_data$PondYear))

n_informative_stocked <- VBGF_data |> filter(AgeSource == "StockedCohort") |>
  dplyr::count(PITIndex) |> filter(n > 1) |> nrow()
cat("\nStocked-cohort fish with >=1 post-release repeat measurement (K-informative):",
    n_informative_stocked, "of", sum(VBGF_data$AgeSource == "StockedCohort" & !duplicated(VBGF_data$PITIndex)), "\n")

## ---- 4. Bayesian von Bertalanffy growth model, release-anchored, shared K -
##
## mu_i = anchor_i + (Linf - anchor_i) * (1 - exp(-K_i * t_i))
##   anchor_i = L0 (recruit) or TL_release_i (stocked, known data)
##   K_i = exp(logK + v[pondyear_i])  [primary model: single shared logK]
##       = exp(logK[origin_i] + v[pondyear_i])  [sensitivity variant only]
## Priors: Linf ~ Normal(635, 50); logK ~ Normal(log(0.15), 0.5);
## L0 ~ Normal(6.8, 0.25) mm (Hamman 1982, as before).

LINF_PRIOR_MEAN <- 635; LINF_PRIOR_SD <- 50
LOGK_PRIOR_MEAN <- log(0.15); LOGK_PRIOR_SD <- 0.5
L0_PRIOR_MEAN <- 6.8; L0_PRIOR_SD <- 0.25

fit_vbgf_release_nimble <- function(t, TL, is_recruit, TL_release, py_idx, origin_idx,
                                     use_origin = FALSE, use_group = TRUE,
                                     n_origin, n_py, seed = 4127,
                                     niter = 60000, nburnin = 20000, thin = 10, nchains = 2,
                                     block = FALSE) {
  N <- length(TL)
  origin_use <- if (use_origin) origin_idx else rep(1L, N)
  n_origin_use <- if (use_origin) n_origin else 1L
  py_use <- if (use_group) py_idx else rep(1L, N)
  n_py_use <- if (use_group) n_py else 1L

  code <- nimbleCode({
    for (i in 1:N) {
      anchor[i] <- is_recruit[i] * L0 + (1 - is_recruit[i]) * TL_release[i]
      K[i] <- exp(logK[origin_idx[i]] + v[py_idx[i]])
      mu[i] <- anchor[i] + (Linf - anchor[i]) * (1 - exp(-K[i] * t[i]))
      TL[i] ~ dnorm(mu[i], sd = sigma)
    }
    for (j in 1:n_py) v[j] ~ dnorm(0, sd = sigma_v)
    for (k in 1:n_origin) logK[k] ~ dnorm(logK_mean, sd = logK_sd)
    Linf ~ dnorm(linf_mean, sd = linf_sd)
    L0 ~ dnorm(l0_mean, sd = l0_sd)
    sigma ~ dunif(0, 200)
    sigma_v ~ dunif(0, 2)
  })
  constants <- list(N = N, n_origin = n_origin_use, n_py = n_py_use,
                     t = t, is_recruit = is_recruit, TL_release = TL_release,
                     origin_idx = origin_use, py_idx = py_use,
                     linf_mean = LINF_PRIOR_MEAN, linf_sd = LINF_PRIOR_SD,
                     logK_mean = LOGK_PRIOR_MEAN, logK_sd = LOGK_PRIOR_SD,
                     l0_mean = L0_PRIOR_MEAN, l0_sd = L0_PRIOR_SD)
  inits <- list(logK = rep(LOGK_PRIOR_MEAN, n_origin_use), v = rep(0, n_py_use),
                Linf = LINF_PRIOR_MEAN, L0 = L0_PRIOR_MEAN, sigma = 50, sigma_v = 0.3)
  monitors <- c("Linf", "logK", "L0", "sigma", "sigma_v", "v")

  m <- nimbleModel(code, constants = constants, data = list(TL = TL),
                    inits = inits, calculate = TRUE, buildDerivs = block)
  conf <- configureMCMC(m, monitors = monitors, enableWAIC = TRUE, print = FALSE)

  if (block) {
    # Same Linf/logK curved-ridge issue as the hatch-anchored version (see
    # revision note); barker sampler resolves it (requires buildDerivs = TRUE).
    conf$removeSamplers(c("Linf", paste0("logK[1:", n_origin_use, "]")))
    conf$addSampler(target = c("Linf", paste0("logK[", 1:n_origin_use, "]")),
                     type = "barker", control = list(scale = 1))
  }

  mcmc <- buildMCMC(conf)
  cm <- compileNimble(m)
  cmcmc <- compileNimble(mcmc, project = m)
  seeds <- seed + seq_len(nchains) - 1
  runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin,
          nchains = nchains, samplesAsCodaMCMC = TRUE, WAIC = TRUE, setSeed = seeds)
}

rhat <- function(fit) gelman.diag(fit$samples, multivariate = FALSE)$psrf[, 1]

t_v <- VBGF_data$t
TL_v <- VBGF_data$total_length
is_recruit_v <- VBGF_data$is_recruit
TL_release_v <- VBGF_data$TL_release_const
origin_v <- VBGF_data$origin_idx
py_v <- VBGF_data$py_idx

## Three-way comparison: primary (shared K + pondyear), pooled (shared K, no
## pondyear), and a sensitivity check re-adding an origin-specific logK on
## top of release-anchoring (to test whether origin still matters once
## anchoring is fixed, on the full dataset rather than just the one shared
## IP5 pond-year cell used for the initial hypothesis test).
al_shared        <- fit_vbgf_release_nimble(t_v, TL_v, is_recruit_v, TL_release_v, py_v, origin_v,
                                             use_origin = FALSE, use_group = TRUE, n_origin = n_origin, n_py = n_py, seed = 6213)
al_shared_pooled <- fit_vbgf_release_nimble(t_v, TL_v, is_recruit_v, TL_release_v, py_v, origin_v,
                                             use_origin = FALSE, use_group = FALSE, n_origin = n_origin, n_py = n_py, seed = 8721)
al_origin_check  <- fit_vbgf_release_nimble(t_v, TL_v, is_recruit_v, TL_release_v, py_v, origin_v,
                                             use_origin = TRUE, use_group = TRUE, n_origin = n_origin, n_py = n_py, seed = 3391)

al_waic_tab <- data.frame(
  model = c("shared K + pondyear (primary)", "shared K, pooled (no pondyear)", "origin + pondyear (sensitivity check)"),
  WAIC = c(al_shared$WAIC$WAIC, al_shared_pooled$WAIC$WAIC, al_origin_check$WAIC$WAIC)
) |> arrange(WAIC)
print(al_waic_tab)

## Refit the primary (shared K + pondyear) model with the barker-sampler
## block for a trustworthy final posterior.
al_final <- fit_vbgf_release_nimble(t_v, TL_v, is_recruit_v, TL_release_v, py_v, origin_v,
                                     use_origin = FALSE, use_group = TRUE, n_origin = n_origin, n_py = n_py, seed = 2247,
                                     niter = 100000, nburnin = 30000, thin = 10, nchains = 3, block = TRUE)

cat("\nFinal (barker-sampler, shared-K) fit R-hat:\n")
print(round(rhat(al_final), 3))
cat("\nFinal fit WAIC:", al_final$WAIC$WAIC, "\n")

samps_final <- do.call(rbind, al_final$samples)
pars <- c("Linf", "logK[1]", "L0", "sigma", "sigma_v")
al_param_summary <- data.frame(
  param = pars,
  mean = round(colMeans(samps_final[, pars]), 4),
  sd = round(apply(samps_final[, pars], 2, sd), 4),
  q2.5 = round(apply(samps_final[, pars], 2, quantile, 0.025), 4),
  q97.5 = round(apply(samps_final[, pars], 2, quantile, 0.975), 4)
)
print(al_param_summary)

al_v_summary <- data.frame(
  PondYear = py_levels,
  mean_v = round(colMeans(samps_final[, paste0("v[", 1:n_py, "]")]), 3),
  sd_v = round(apply(samps_final[, paste0("v[", 1:n_py, "]")], 2, sd), 3)
)
print(al_v_summary)

## Also refit the origin-check sensitivity model with the barker block, so
## its logK[1]/logK[2] posteriors are trustworthy for the writeup (not just
## used for WAIC above).
al_origin_check_final <- fit_vbgf_release_nimble(t_v, TL_v, is_recruit_v, TL_release_v, py_v, origin_v,
                                                   use_origin = TRUE, use_group = TRUE, n_origin = n_origin, n_py = n_py, seed = 9142,
                                                   niter = 100000, nburnin = 30000, thin = 10, nchains = 3, block = TRUE)
cat("\nOrigin-check sensitivity fit R-hat:\n")
print(round(rhat(al_origin_check_final), 3))
samps_origin_check <- do.call(rbind, al_origin_check_final$samples)
origin_check_pars <- c("Linf", "logK[1]", "logK[2]", "L0", "sigma", "sigma_v")
al_origin_check_summary <- data.frame(
  param = origin_check_pars,
  mean = round(colMeans(samps_origin_check[, origin_check_pars]), 4),
  sd = round(apply(samps_origin_check[, origin_check_pars], 2, sd), 4),
  q2.5 = round(apply(samps_origin_check[, origin_check_pars], 2, quantile, 0.025), 4),
  q97.5 = round(apply(samps_origin_check[, origin_check_pars], 2, quantile, 0.975), 4)
)
print(al_origin_check_summary)
cat("\nP(logK[StockedCohort] < logK[AssumedRecruit]):",
    mean(samps_origin_check[, "logK[1]"] < samps_origin_check[, "logK[2]"]), "\n")

## ---- 5. Save results -------------------------------------------------------

save(
  known_yc, recruit_candidates, MixtureClassificationLog, AssumedYOY,
  YearClassMap, AgeLength, release_anchor, VBGF_data, n_informative_stocked,
  al_shared, al_shared_pooled, al_origin_check, al_waic_tab, al_final,
  al_param_summary, al_v_summary,
  al_origin_check_final, al_origin_check_summary,
  LINF_PRIOR_MEAN, LINF_PRIOR_SD, LOGK_PRIOR_MEAN, LOGK_PRIOR_SD, L0_PRIOR_MEAN, L0_PRIOR_SD,
  file = "data/GIEL_VBGF_AgeLength.RData"
)
