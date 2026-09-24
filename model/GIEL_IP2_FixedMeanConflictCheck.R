## GIEL_IP2_FixedMeanConflictCheck.R
##
## Conflict check for the BONY (GIEL) fixed-mean age-0 (YOY) vs. holdover
## size-at-netting priors proposed for IPCA Ponds 2, 5, and 6 (name retained
## from the original IP2-only request; scope is now all three GIEL ponds).
##
## Background: for BONY at IPCA, BK confirmed (1) there is not enough data to
## support a density-dependent growth/survival submodel like the one built
## for RASU at YCB (model/YCBSizeDensityEDA.qmd), and (2) the May 1 spawn-date
## anchor (Jonez and Sumner 1954, via BK) is BONY/IPCA-specific. Given that,
## the plan is to classify untagged fall-netted (Oct-Dec) GIEL recruits into
## age-0 (YOY) vs. holdover using FIXED size-at-age means/SDs computed from
## the existing shared-K, release-anchored VBGF posterior
## (data/GIEL_VBGF_AgeLength.RData, object `al_final`), anchored at May 1 for
## age-0 and May 1 of the PRIOR year for a one-year-old holdover, rather than
## letting an unconstrained 2-component mixture (as used for RASU/YCB in
## RecruitYOYClassification.R) fit both means freely. BK expects little
## conflict between the growth-model prediction and observed netted sizes,
## since IPCA is data-poor and the growth model already fits reasonably well.
##
## This script tests that expectation with two cheap, per-event diagnostics
## (chosen over a full parametric-bootstrap goodness-of-fit test as the first
## screening pass -- see the 2026-09-23 conversation log in AGENTS.md):
##   1. Per-component z-check: hard-assign each fish to whichever fixed
##      component (age-0 or holdover) has higher density, then compare the
##      assigned group's observed mean TL to the fixed predicted mean via
##      z = (obs_mean - pred_mean) / (pred_sd / sqrt(n_assigned)). Flag
##      |z| > 2.
##   2. Extreme-residual fish: flag individual fish whose distance from BOTH
##      fixed components (min(|z0|, |z1|)) exceeds 3 -- these don't fit
##      either predicted cluster and may indicate an unmodeled age class
##      (e.g. a multi-year-old holdover) or a data problem.
##
## IMPORTANT ASYMMETRY: the "Holdover" fixed mean/SD is calibrated for
## EXACTLY one year old (age ~1.6 yr at a typical Dec netting). But the
## untagged-recruit pool used here (fish with no known year class, first
## captured at a fall netting) can include truly older fish that simply
## evaded capture in prior years -- age-2, age-3, etc. Those fish are
## expected to run LARGER than the age-1 anchor (VBGF growth is decelerating
## but still positive at these ages), so a positive z-test flag on a
## "Holdover" group is consistent with ordinary multi-age mixing and is NOT
## by itself evidence of a growth-model conflict. A NEGATIVE flag (observed
## smaller than the age-1 prediction), or any flag on the "Age-0" group in
## either direction, is the more informative signal of real conflict (wrong
## spawn-date anchor, a bad growth-model year, or a data issue). The output
## table below annotates this distinction explicitly (`interpretation`
## column) rather than leaving every flagged row looking equally concerning.
##
## Data sources: data/ReportingData.RData (StudyBWNFWG, for the untagged
## fall-netted recruit candidates -- same candidate definition used in
## model/GIEL_VBGF_AgeLength.R's step 2) and data/GIEL_VBGF_AgeLength.RData
## (al_final, for the posterior draws). The latter file is a full workspace
## dump from an interactive session (see AGENTS.md caution) -- loaded into a
## throwaway environment here so it cannot clobber anything in this script.
##
## Output: data/GIEL_IP2_FixedMeanConflictCheck.RData (EventFixedPriors,
## FishClassified, EventZCheck, ExtremeResidualFish) and
## model/GIEL_IP2_FixedMeanConflictCheck.png (per-event histograms with the
## fixed component densities overlaid).

source("LabFunctions.R")
packages(dplyr)
packages(lubridate)
packages(ggplot2)
packages(coda)

load("data/ReportingData.RData")

## Posterior draws for the shared-K, release-anchored VBGF fit, loaded into a
## scratch environment so the full-workspace-dump RData can't overwrite
## anything above (see file header).
al_env <- new.env()
load("data/GIEL_VBGF_AgeLength.RData", envir = al_env)
al_final <- al_env$al_final
rm(al_env)

samps <- as.matrix(al_final$samples)
stopifnot(all(c("Linf", "logK[1]", "L0", "sigma", "sigma_v") %in% colnames(samps)))
Linf_s   <- samps[, "Linf"]
logK_s   <- samps[, "logK[1]"]
L0_s     <- samps[, "L0"]
sigma_s  <- samps[, "sigma"]
sigmav_s <- samps[, "sigma_v"]
N_DRAW <- length(Linf_s)
cat("Loaded", N_DRAW, "posterior draws from al_final (Linf, logK[1], L0, sigma, sigma_v)\n")

EXTREME_Z <- 3    # per-fish: flag if min(|z0|, |z1|) exceeds this
EVENT_Z   <- 2    # per-event: flag if the assigned-group z-test exceeds this
SEED <- 9142

# ---------------------------------------------------------------------------
# 1. Untagged fall-netted (Oct-Dec) GIEL recruit candidates, IP2/IP5/IP6 --
#    same definition as model/GIEL_VBGF_AgeLength.R step 2 (fish with no
#    known stocking year class, first physical capture, Oct-Dec, measured TL)
# ---------------------------------------------------------------------------
GIEL_PONDS <- c("IPCA (Pond 2)", "IPCA (Pond 5)", "IPCA (Pond 6)")

giel_nfwg <- StudyBWNFWG |> filter(species == "GIEL", location %in% GIEL_PONDS)

known_yc <- giel_nfwg |>
  filter(!is.na(year_class)) |>
  distinct(PITIndex, year_class)

RecruitCandidates <- giel_nfwg |>
  filter(FirstRecord == "yes", event == "capture", month(collection_date) %in% 10:12,
         !(PITIndex %in% known_yc$PITIndex), !is.na(total_length)) |>
  mutate(cal_year = year(collection_date), ym = format(collection_date, "%Y-%m"))

cat("\nUntagged fall-netted recruit candidates (IP2/IP5/IP6):", nrow(RecruitCandidates), "\n")

# ---------------------------------------------------------------------------
# 2. Fixed size-at-age predictions per netting event, anchored to May 1
#    (age-0) / May 1 of the prior year (one-year-old holdover)
# ---------------------------------------------------------------------------
may1_this <- function(d) as.Date(paste0(year(d), "-05-01"))
may1_prev <- function(d) as.Date(paste0(year(d) - 1, "-05-01"))

set.seed(SEED)
predict_fixed <- function(age_years) {
  # Posterior-predictive mean/SD for a NEW pond-year (fresh v ~ N(0, sigma_v))
  # plus residual observation noise -- same recipe used for the original
  # IP2 age-0/holdover prior table.
  v_draw <- rnorm(N_DRAW, 0, sigmav_s)
  K <- exp(logK_s + v_draw)
  mu <- L0_s + (Linf_s - L0_s) * (1 - exp(-K * age_years))
  TL_draw <- rnorm(N_DRAW, mu, sigma_s)
  c(mean = mean(TL_draw), sd = sd(TL_draw))
}

EventFixedPriors <- RecruitCandidates |>
  group_by(location, ym) |>
  summarise(n = dplyr::n(), rep_date = median(collection_date), .groups = "drop") |>
  arrange(location, ym) |>
  mutate(age0_age = as.numeric(rep_date - may1_this(rep_date)) / 365.25,
         age1_age = as.numeric(rep_date - may1_prev(rep_date)) / 365.25)

fixed_preds <- t(sapply(seq_len(nrow(EventFixedPriors)), function(i) {
  c(predict_fixed(EventFixedPriors$age0_age[i]), predict_fixed(EventFixedPriors$age1_age[i]))
}))
colnames(fixed_preds) <- c("mean0", "sd0", "mean1", "sd1")
EventFixedPriors <- bind_cols(EventFixedPriors, as.data.frame(fixed_preds))

cat("\nFixed age-0 / holdover(age-1) size predictions by netting event:\n")
print(EventFixedPriors, n = Inf)

# ---------------------------------------------------------------------------
# 3. Per-fish classification against the fixed means, and the two
#    diagnostics: extreme-residual fish, and per-event/per-assigned-class
#    z-check
# ---------------------------------------------------------------------------
FishClassified <- RecruitCandidates |>
  left_join(EventFixedPriors |> select(location, ym, mean0, sd0, mean1, sd1),
            by = c("location", "ym")) |>
  mutate(
    z0 = (total_length - mean0) / sd0,
    z1 = (total_length - mean1) / sd1,
    assigned = if_else(abs(z0) <= abs(z1), "Age-0 (YOY)", "Holdover"),
    min_abs_z = pmin(abs(z0), abs(z1)),
    extreme = min_abs_z > EXTREME_Z
  )

ExtremeResidualFish <- FishClassified |>
  filter(extreme) |>
  select(location, ym, PITIndex, collection_date, total_length, z0, z1, min_abs_z) |>
  arrange(location, ym, desc(min_abs_z))

cat("\n--- Diagnostic 1: extreme-residual fish (min(|z0|,|z1|) >", EXTREME_Z, ") ---\n")
cat("Count:", nrow(ExtremeResidualFish), "of", nrow(FishClassified), "classified fish\n")
if (nrow(ExtremeResidualFish) > 0) print(ExtremeResidualFish, n = Inf)

EventZCheck <- FishClassified |>
  group_by(location, ym, assigned) |>
  summarise(n_assigned = dplyr::n(),
            obs_mean = mean(total_length), obs_sd = sd(total_length),
            pred_mean = if (first(assigned) == "Age-0 (YOY)") first(mean0) else first(mean1),
            pred_sd   = if (first(assigned) == "Age-0 (YOY)") first(sd0)   else first(sd1),
            .groups = "drop") |>
  mutate(
    z_test = (obs_mean - pred_mean) / (pred_sd / sqrt(n_assigned)),
    flag_z = !is.na(z_test) & abs(z_test) > EVENT_Z,
    direction = if_else(is.na(z_test), NA_character_, if_else(z_test > 0, "high", "low")),
    interpretation = case_when(
      !flag_z ~ "consistent with fixed prediction",
      assigned == "Holdover" & direction == "high" ~
        "expected: consistent with a multi-age holdover pool (fish older than the age-1 anchor), not a growth-model conflict",
      assigned == "Holdover" & direction == "low" ~
        "POSSIBLE CONFLICT: holdover fish smaller than the age-1 prediction",
      assigned == "Age-0 (YOY)" ~
        paste0("POSSIBLE CONFLICT: age-0 fish running ", direction,
               "er than predicted -- check spawn-date anchor / growth rate for this event"),
      TRUE ~ NA_character_
    )
  ) |>
  arrange(location, ym, assigned)

cat("\n--- Diagnostic 2: per-event z-check (assigned-group mean vs. fixed prediction) ---\n")
print(EventZCheck, n = Inf)

cat("\nEvents/classes flagged as a POSSIBLE (not just 'expected multi-age') conflict:\n")
print(EventZCheck |> filter(flag_z, grepl("POSSIBLE", interpretation)), n = Inf)

# ---------------------------------------------------------------------------
# 4. Figure: per-event histograms with the two fixed component densities
# ---------------------------------------------------------------------------
pond_short <- function(x) sub("IPCA \\((Pond \\d+)\\)", "\\1", x)

EventFixedPriors <- EventFixedPriors |>
  mutate(panel_label = paste0(pond_short(location), " ", ym, " (n=", n, ")"))
FishClassified <- FishClassified |>
  left_join(EventFixedPriors |> select(location, ym, panel_label), by = c("location", "ym"))

dens_grid <- EventFixedPriors |>
  rowwise() |>
  reframe(
    panel_label = panel_label,
    TL = seq(50, 500, by = 2),
    dens0 = dnorm(TL, mean0, sd0),
    dens1 = dnorm(TL, mean1, sd1)
  ) |>
  tidyr::pivot_longer(c(dens0, dens1), names_to = "component", values_to = "density") |>
  mutate(component = if_else(component == "dens0", "Age-0 (YOY)", "Holdover (age-1)"))

p_conflict <- ggplot(FishClassified, aes(x = total_length)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 15, fill = "grey70", color = "white") +
  geom_line(data = dens_grid, aes(x = TL, y = density, color = component), linewidth = 0.8) +
  facet_wrap(~panel_label, scales = "free_y", ncol = 4) +
  labs(x = "Total length (mm)", y = "Density",
       color = NULL,
       title = "GIEL fall-netted recruits vs. fixed May-1-anchored size-at-age predictions",
       subtitle = "IPCA Ponds 2, 5, 6 -- untagged first captures, Oct-Dec") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))
print(p_conflict)
ggsave("model/GIEL_IP2_FixedMeanConflictCheck.png", p_conflict, width = 13, height = 9, dpi = 150)

# ---------------------------------------------------------------------------
# 5. Save
# ---------------------------------------------------------------------------
save(EventFixedPriors, FishClassified, EventZCheck, ExtremeResidualFish,
     file = "data/GIEL_IP2_FixedMeanConflictCheck.RData")
cat("\nSaved data/GIEL_IP2_FixedMeanConflictCheck.RData\n")
cat("Saved model/GIEL_IP2_FixedMeanConflictCheck.png\n")
