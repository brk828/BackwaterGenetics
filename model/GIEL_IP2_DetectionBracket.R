# Model-free (bracketed known-alive) detection comparison: juvenile vs. adult
# PIT-scanner detection at IPCA Pond 2 (GIEL).
#
# Purpose: test whether the joint robust-design model's finding of higher
# apparent juvenile-than-adult survival at Pond 2 (dJ[GIEL] > 0, driven
# largely by Pond 2 per AGENTS.md's "IP2 adult-survival investigation") is an
# artifact of lower juvenile detection probability, rather than a real
# survival difference. This analysis does NOT fit a survival model at all --
# it estimates detection probability directly from fish whose alive status in
# a given month is not in question, so any month without a contact is
# unambiguously a missed detection, not a death.
#
# Bracket definition: a fish-month is included only if (a) the fish had
# already entered the pond (stocked/tagged) by the end of that month, and
# (b) the fish was contacted again strictly AFTER that month (MaxScanDate >
# month_end). Condition (b) guarantees the fish was alive during the month in
# question (closed population, no interim mortality/removal evidence) and
# also excludes the trivial terminal month, whose "detected" status is
# guaranteed by definition of MaxScanDate and would otherwise inflate
# detection estimates.
#
# Stage assignment: last-observation-carried-forward from physical NFWG total
# length measurements (StudyBWNFWG), classified at the StageTL threshold used
# throughout the project (GIEL: TL < 250 mm = juvenile, >= 250 mm = adult).
# This uses only observed lengths, not model-imputed stage.

source("LabFunctions.R")
packages(dplyr)
packages(lubridate)
packages(tidyr)
packages(ggplot2)
packages(lme4)

load("data/ReportingData.RData")

LOC <- "IPCA (Pond 2)"
STAGE_TL <- 250

# ---------------------------------------------------------------------------
# 1. Fish at Pond 2 (GIEL), entry date + last-contact date
# ---------------------------------------------------------------------------
Pond2Fish <- StudyBWAnalysis |>
  filter(location == LOC, species == "GIEL") |>
  transmute(PITIndex, first_date = as.Date(first_date), event,
            entry_TL = total_length, MaxScanDate = as.Date(MaxScanDate))

cat("Pond 2 GIEL fish:", nrow(Pond2Fish), "\n")

# ---------------------------------------------------------------------------
# 2. TL observation history (physical NFWG records only -- length is never
#    recorded by the PIT scanner) for last-observation-carried-forward
#    stage assignment
# ---------------------------------------------------------------------------
TLObsP2 <- StudyBWNFWG |>
  filter(location == LOC, species == "GIEL", !is.na(total_length)) |>
  transmute(PITIndex, obs_date = as.Date(collection_date), TL = total_length) |>
  distinct() |>
  arrange(PITIndex, obs_date)

# ---------------------------------------------------------------------------
# 3. Monthly grid per fish, from entry month through the month before
#    MaxScanDate's month (the strict bracket -- see header)
# ---------------------------------------------------------------------------
month_floor <- function(d) floor_date(d, "month")

Pond2Fish <- Pond2Fish |>
  mutate(entry_month = month_floor(first_date),
         last_bracket_month = month_floor(MaxScanDate %m-% months(1)))

cat("Fish with no possible bracketed month (last contact within 1 month of entry):",
    sum(Pond2Fish$last_bracket_month < Pond2Fish$entry_month), "of", nrow(Pond2Fish), "\n")

FishMonths <- Pond2Fish |>
  filter(last_bracket_month >= entry_month) |>
  rowwise() |>
  reframe(PITIndex = PITIndex,
          month_start = seq(entry_month, last_bracket_month, by = "month")) |>
  mutate(month_end = ceiling_date(month_start, "month") - 1)

cat("Bracketed known-alive fish-months (pre-stage-assignment):", nrow(FishMonths), "\n")
cat("Fish contributing at least 1 bracketed month:", n_distinct(FishMonths$PITIndex),
    "of", nrow(Pond2Fish), "\n")

# ---------------------------------------------------------------------------
# 4. Last-observation-carried-forward stage per fish-month
# ---------------------------------------------------------------------------
StageByMonth <- FishMonths |>
  left_join(TLObsP2, by = "PITIndex", relationship = "many-to-many") |>
  filter(is.na(obs_date) | obs_date <= month_end) |>
  group_by(PITIndex, month_start) |>
  slice_max(obs_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(stage = if_else(TL < STAGE_TL, "Juvenile", "Adult")) |>
  select(PITIndex, month_start, month_end, TL_carried = TL, stage)

FishMonths <- FishMonths |>
  left_join(StageByMonth, by = c("PITIndex", "month_start", "month_end")) |>
  filter(!is.na(stage))   # drop the rare fish with no TL record at/before entry

cat("Fish-months with a defined stage:", nrow(FishMonths), "\n")
cat("Stage-transition fish-months (juvenile at entry, later reclassified adult):\n")
print(FishMonths |> group_by(PITIndex) |> summarise(n_stage = n_distinct(stage)) |>
        dplyr::count(n_stage))

# ---------------------------------------------------------------------------
# 5. Detection: any PIT contact at Pond 2 that fish that calendar month
# ---------------------------------------------------------------------------
ContactsP2 <- StudyBWContacts |>
  filter(Backwater == LOC, !is.na(PITIndex)) |>
  transmute(PITIndex, month_start = floor_date(Date, "month")) |>
  distinct()

FishMonths <- FishMonths |>
  left_join(ContactsP2 |> mutate(detected = 1L), by = c("PITIndex", "month_start")) |>
  mutate(detected = coalesce(detected, 0L))

# ---------------------------------------------------------------------------
# 6. Monthly scan effort at Pond 2 (prorate deployments to calendar days,
#    same method used in IPGIELScannerDiagnostics.qmd / BWRobustDesign_data.R)
# ---------------------------------------------------------------------------
EffortDaily <- StudyBWEffort |>
  filter(Location == LOC, ScanTimeHrs > 0, !is.na(Deploy), !is.na(Retrieve),
         Retrieve > Deploy) |>
  transmute(EID, dep = as.Date(Deploy), ret = as.Date(Retrieve),
            hrs_per_day = ScanTimeHrs / pmax(1, as.numeric(ret - dep) + 1)) |>
  rowwise() |>
  reframe(day = seq(dep, ret, by = "day"), hrs = hrs_per_day) |>
  mutate(month_start = floor_date(day, "month")) |>
  group_by(month_start) |>
  summarise(effort_hrs = sum(hrs), .groups = "drop")

FishMonths <- FishMonths |>
  left_join(EffortDaily, by = "month_start") |>
  mutate(effort_hrs = coalesce(effort_hrs, 0))

# ---------------------------------------------------------------------------
# 7. Monthly summary: proportion of bracketed known-alive fish detected, by
#    stage
# ---------------------------------------------------------------------------
MonthlyByStage <- FishMonths |>
  group_by(month_start, stage) |>
  summarise(n_alive = n(), n_detected = sum(detected),
            prop_detected = n_detected / n_alive,
            effort_hrs = first(effort_hrs), .groups = "drop")

cat("\nOverall proportion detected by stage (all bracketed fish-months pooled):\n")
print(FishMonths |> group_by(stage) |>
        summarise(n_fish_months = n(), n_fish = n_distinct(PITIndex),
                  n_detected = sum(detected),
                  prop_detected = mean(detected), .groups = "drop"))

cat("\nMedian of monthly proportion-detected, by stage (months with >=5 known-alive fish):\n")
print(MonthlyByStage |> filter(n_alive >= 5) |> group_by(stage) |>
        summarise(n_months = n(), median_prop = median(prop_detected),
                  IQR_lo = quantile(prop_detected, .25),
                  IQR_hi = quantile(prop_detected, .75), .groups = "drop"))

# ---------------------------------------------------------------------------
# 8. Formal comparison: mixed-effects logistic regression on the fish-month
#    level, adjusting for scan effort and clustering by fish (repeated
#    monthly observations per fish) and by month (shared scanner conditions)
# ---------------------------------------------------------------------------
FishMonths <- FishMonths |>
  mutate(log_effort = log1p(effort_hrs),
         stage = factor(stage, levels = c("Adult", "Juvenile")))

fit_detect <- glmer(detected ~ stage + log_effort + (1 | PITIndex) + (1 | month_start),
                     data = FishMonths, family = binomial)

cat("\n--- Mixed-effects logistic model: detected ~ stage + log_effort + (1|PITIndex) + (1|month) ---\n")
print(summary(fit_detect))

or_ci <- exp(confint(fit_detect, parm = "stageJuvenile", method = "Wald"))
cat("\nOdds ratio (Juvenile vs Adult detection), adjusted for effort:\n")
print(c(OR = exp(fixef(fit_detect)["stageJuvenile"]), or_ci))

# ---------------------------------------------------------------------------
# 9. Plots
# ---------------------------------------------------------------------------
p_monthly <- MonthlyByStage |>
  filter(n_alive >= 5) |>
  ggplot(aes(x = month_start, y = prop_detected, color = stage)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(size = n_alive), alpha = 0.7) +
  labs(x = NULL, y = "Proportion of bracketed known-alive fish detected",
       color = "Stage", size = "Fish alive",
       title = "IPCA Pond 2 (GIEL): monthly detection by stage",
       subtitle = "Bracketed known-alive fish-months only; months with < 5 fish alive omitted") +
  theme_minimal()

p_box <- FishMonths |>
  ggplot(aes(x = stage, y = detected)) +
  stat_summary(fun = mean, geom = "col", aes(fill = stage), alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(x = NULL, y = "Fraction of fish-months detected",
       title = "IPCA Pond 2 (GIEL): overall detection rate by stage",
       subtitle = "Bracketed known-alive fish-months; mean \u00b1 SE") +
  theme_minimal() + theme(legend.position = "none")

print(p_monthly)
print(p_box)

save(Pond2Fish, TLObsP2, FishMonths, MonthlyByStage, fit_detect,
     file = "data/GIEL_IP2_DetectionBracket.RData")

cat("\nSaved data/GIEL_IP2_DetectionBracket.RData\n")
