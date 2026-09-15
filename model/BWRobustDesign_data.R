# Data preparation for the joint stage-structured robust-design model of
# backwater Razorback Sucker (XYTE) and Bonytail (GIEL) populations.
# B. Kesner / Posit Assistant, September 2026
#
# Builds all inputs for model/BWRobustDesign_NIMBLE.R from
# data/ReportingData.RData (DataWrangling.R output) and data/BWNettingEvents.csv
# (hand-maintained untagged netting counts), and saves
# data/BWRobustDesign_data.RData.
#
# Time structure: seasons = fiscal years FY2017..FY2026; primary periods =
# calendar months Oct-Apr (7 per season, 70 total); secondary occasions =
# weeks within a month (days 1-7, 8-14, 15-21, 22-end).

source("LabFunctions.R")
packages(dplyr)
packages(lubridate)
packages(tidyr)
packages(stringr)
packages(readr)

load("data/ReportingData.RData")

if (!exists("DIEOFF_FILE")) DIEOFF_FILE <- "data/BWDieOffEvents.csv"

# GIEL size-tied maturation hazard track (2026-09-14; see AGENTS.md "GIEL
# Size-Tied Maturation Hazard" and the plan file .posit/assistant/plans/
# 2026-09-14-1444-giel-growth-model-tracks-feeding-robust-design-maturation-
# hazard.md). One of "F1" (Fabens, IPCA-only), "F2" (Fabens, IPCA+Cibola),
# "A1" (age-based VBGF, IPCA-only) -- built by model/GIEL_MaturationHazard_Build.R.
if (!exists("MATURATION_TRACK")) MATURATION_TRACK <- "F1"
stopifnot(MATURATION_TRACK %in% c("F1", "F2", "A1"))
MATURATION_HAZARD_FILE <- paste0("data/GIEL_MaturationHazard_", MATURATION_TRACK, ".RData")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
START_DATE  <- as.Date("2016-10-01")   # first primary period (FY2017)
FIRST_FY    <- 2017
N_SEASONS   <- 10                      # FY2017..FY2026
N_MOY       <- 7                       # Oct..Apr
T_PRIM      <- N_SEASONS * N_MOY       # 70 primary periods
N_SEC       <- 4                       # weekly secondaries per month
POST_MONTHS <- 6                       # post-release survival offset window
StageTL     <- c(XYTE = 350, GIEL = 250)   # juvenile < threshold <= adult

PondTable <- tibble(
  pond      = 1:7,
  location  = c("Yuma Cove backwater", paste0("IPCA (Pond ", 1:6, ")")),
  LID       = c(592L, 1043:1048),
  Backwater = c("YCB", paste0("IP", 1:6)),
  species   = c("XYTE", "XYTE", "GIEL", "XYTE", "XYTE", "GIEL", "GIEL"),
  sp        = c(1L, 1L, 2L, 1L, 1L, 2L, 2L)
)

# ---------------------------------------------------------------------------
# Primary period table
# ---------------------------------------------------------------------------
Primary <- tibble(t = 1:T_PRIM) |>
  mutate(
    season      = (t - 1) %/% N_MOY + 1,
    moy         = (t - 1) %% N_MOY + 1,
    FY          = FIRST_FY + season - 1,
    month       = ifelse(moy <= 3, moy + 9, moy - 3),
    year        = ifelse(month > 9, FY - 1, FY),
    month_start = make_date(year, month, 1),
    month_end   = ceiling_date(month_start, "month") - 1,
    gap         = ifelse(moy == N_MOY, 6L, 1L)   # months to next primary
  )

fy_of_date <- function(d) ifelse(month(d) > 9, year(d) + 1L, year(d))
season_of_date <- function(d) fy_of_date(d) - FIRST_FY + 1L

# Primary period at or before a date: Oct-Apr -> that month; May-Sep -> Apr of
# the same FY. Dates before START_DATE or beyond the last primary return NA.
prim_floor <- function(d) {
  fy  <- fy_of_date(d)
  m   <- month(d)
  moy <- ifelse(m > 9, m - 9, m + 3)
  t   <- (fy - FIRST_FY) * N_MOY + pmin(moy, N_MOY)
  t[d < START_DATE | t > T_PRIM] <- NA_integer_
  as.integer(t)
}

week_of_month <- function(d) pmin(N_SEC, (day(d) - 1L) %/% 7L + 1L)

# ---------------------------------------------------------------------------
# Removals (harvest / transfer / recovered mortality), per fish x backwater
# ---------------------------------------------------------------------------
Removals <- bind_rows(
  StudyBWTransfersAnalysis |>
    mutate(Backwater = recode(Backwater, "Imperial NWR Pond 3" = "IPCA (Pond 3)")) |>
    filter(!PITIndex %in% SuspectTransfers$PITIndex) |>
    transmute(location = Backwater, PITIndex, rem_date = as.Date(TransferDate),
              rem_type = ifelse(str_detect(TransferLocation, "^IPCA"), "transfer", "harvest"),
              rem_TL = coalesce(CaptureTL, TransferTL)),
  StudyBWNFWG |>
    filter(event == "capture", disposition == "retained",
           !PITIndex %in% SuspectTransfers$PITIndex) |>
    transmute(location, PITIndex, rem_date = as.Date(collection_date),
              rem_type = "harvest", rem_TL = total_length),
  StudyBWNFWG |>
    filter(event == "capture", disposition == "recovered") |>
    transmute(location, PITIndex, rem_date = as.Date(collection_date),
              rem_type = "mortality", rem_TL = total_length)
) |>
  group_by(location, PITIndex) |>
  slice_min(rem_date, n = 1, with_ties = FALSE) |>
  ungroup()

# ---------------------------------------------------------------------------
# Fish table: one row per PIT x backwater, FY2017 onward
# ---------------------------------------------------------------------------
FirstContactInWindow <- StudyBWContacts |>
  filter(!is.na(PITIndex), Date >= START_DATE) |>
  group_by(PITIndex, location = Backwater) |>
  summarise(first_contact = min(Date), .groups = "drop")

TLObs <- StudyBWNFWG |>
  filter(!is.na(total_length)) |>
  transmute(location, PITIndex, obs_date = as.Date(collection_date), TL = total_length)

Fish <- StudyBWAnalysis |>
  inner_join(PondTable |> select(location, pond, pond_species = species, sp), by = "location") |>
  mutate(first_date = as.Date(first_date))

SpeciesMismatch <- Fish |> filter(species != pond_species)
if (nrow(SpeciesMismatch) > 0) {
  warning(nrow(SpeciesMismatch), " fish excluded: species does not match pond species ",
          "(see SpeciesMismatch)")
  Fish <- Fish |> filter(species == pond_species)
}

Fish <- Fish |>
  left_join(FirstContactInWindow, by = c("PITIndex", "location")) |>
  left_join(Removals, by = c("PITIndex", "location")) |>
  mutate(
    origin = case_when(
      first_date < START_DATE ~ 3L,          # established before window
      event == "stocking"     ~ 1L,          # stocked
      TRUE                    ~ 2L           # tagged at netting
    ),
    entry_date = if_else(origin == 3L, first_contact, first_date)
  ) |>
  # established fish never contacted in the window carry no information
  filter(!is.na(entry_date)) |>
  # fish tagged and removed in the same event carry no information as
  # individuals; they are counted as untagged removals in the netting table
  mutate(same_event_removal = !is.na(rem_date) & origin == 2L &
           rem_date <= first_date + 3)

SameEventRemovals <- Fish |> filter(same_event_removal)
Fish <- Fish |> filter(!same_event_removal)

# Stage at entry: TL at tagging for new fish; for established fish the latest
# TL within 12 months before entry, otherwise adult
PriorTL <- Fish |>
  select(PITIndex, location, entry_date) |>
  inner_join(TLObs, by = c("PITIndex", "location"), relationship = "many-to-many") |>
  filter(obs_date <= entry_date) |>
  group_by(PITIndex, location) |>
  slice_max(obs_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(PITIndex, location, obs_date, TL)

StageAtEntry <- Fish |>
  select(PITIndex, location, entry_date, species, origin, total_length) |>
  left_join(PriorTL, by = c("PITIndex", "location")) |>
  mutate(
    thr = StageTL[species],
    recent_obs = origin == 3L & !is.na(obs_date) & obs_date >= entry_date %m-% months(12),
    TL_use = case_when(origin == 3L & recent_obs ~ TL,
                       origin == 3L ~ NA_real_,
                       TRUE ~ total_length),
    stage0 = case_when(
      origin == 3L & !recent_obs ~ 2L,
      !is.na(TL_use) & TL_use <  thr ~ 1L,
      !is.na(TL_use) & TL_use >= thr ~ 2L,
      origin == 1L ~ 2L,                     # stocked, TL missing
      TRUE         ~ 1L                      # netting-tagged, TL missing
    ),
    stage0_imputed = is.na(TL_use) & origin != 3L
  ) |>
  select(PITIndex, location, stage0, stage0_imputed, TL_entry = TL_use)

Fish <- Fish |>
  left_join(StageAtEntry, by = c("PITIndex", "location")) |>
  mutate(
    entry_t = prim_floor(entry_date),
    rem_t   = prim_floor(rem_date),
    removed = !is.na(rem_t),
    end_t   = if_else(removed, rem_t, T_PRIM),
    rem_type = if_else(removed, rem_type, NA_character_)
  ) |>
  filter(!is.na(entry_t), end_t >= entry_t) |>
  arrange(pond, entry_t, PITIndex) |>
  mutate(fish = row_number())

# ---------------------------------------------------------------------------
# GIEL size-tied maturation hazard: TL bin at entry (2026-09-14). Only
# GIEL juveniles (stage0 == 1) with a real (non-imputed) TL_entry get the
# size-tied hazard; everyone else (XYTE, adults-at-entry, and GIEL juveniles
# with imputed/missing TL) keeps the existing scalar lpsi[sp] pathway --
# hasTL gates this in RD_constants/rd_code below. TL bin edges match the
# 20mm, 100-249mm bins used by model/GIEL_MaturationHazard_Build.R for every
# track, so a fixed table here is safe regardless of which track's hazard
# file is loaded.
GIEL_TLBIN_EDGES <- c(100, 120, 140, 160, 180, 200, 220, 240, 249)
N_TLBIN_GIEL <- length(GIEL_TLBIN_EDGES) - 1

Fish <- Fish |>
  mutate(
    hasTL = species == "GIEL" & stage0 == 1L & !stage0_imputed & !is.na(TL_entry),
    TL_entry_clip = pmin(pmax(TL_entry, GIEL_TLBIN_EDGES[1]), GIEL_TLBIN_EDGES[N_TLBIN_GIEL + 1] - 1e-6),
    TL_entry_bin = if_else(hasTL,
                            findInterval(TL_entry_clip, GIEL_TLBIN_EDGES, all.inside = TRUE),
                            1L)
  ) |>
  select(-TL_entry_clip)

cat("\nGIEL size-tied maturation hazard: fish with known TL at entry (hasTL) by pond x bin:\n")
print(Fish |> filter(hasTL) |>
        mutate(bin_label = paste0("[", GIEL_TLBIN_EDGES[TL_entry_bin], ",",
                                   GIEL_TLBIN_EDGES[TL_entry_bin + 1], ")")) |>
        dplyr::count(pond, bin_label))
cat("GIEL juveniles (stage0==1) WITHOUT usable TL (scalar-lpsi fallback):",
    sum(Fish$species == "GIEL" & Fish$stage0 == 1L & !Fish$hasTL), "of",
    sum(Fish$species == "GIEL" & Fish$stage0 == 1L), "\n")

N_FISH <- nrow(Fish)

# ---------------------------------------------------------------------------
# Effort: weekly availability and monthly scan hours per pond
# ---------------------------------------------------------------------------
WeekGrid <- Primary |>
  select(t, month_start, month_end) |>
  crossing(k = 1:N_SEC) |>
  mutate(wk_start = month_start + (k - 1) * 7,
         wk_end   = if_else(k == N_SEC, month_end, month_start + k * 7 - 1))

EffortRecords <- StudyBWEffort |>
  inner_join(PondTable |> select(location, pond), by = c("Location" = "location")) |>
  filter(ScanTimeHrs > 0, !is.na(Deploy), !is.na(Retrieve), Retrieve > Deploy) |>
  transmute(pond, EID,
            dep = as.Date(Deploy), ret = as.Date(Retrieve),
            hrs_per_day = ScanTimeHrs / pmax(1, as.numeric(ret - dep) + 1)) |>
  filter(ret >= START_DATE)

# Overlap of each effort record with each pond-week
EffortWeekly <- EffortRecords |>
  cross_join(WeekGrid) |>
  mutate(ov_days = as.numeric(pmin(ret, wk_end) - pmax(dep, wk_start)) + 1) |>
  filter(ov_days > 0) |>
  group_by(pond, t, k) |>
  summarise(eff_hrs = sum(ov_days * hrs_per_day), .groups = "drop")

# Weeks with contacts in a pond imply a functioning scanner regardless of
# effort record quality
ContactWeeks <- StudyBWContacts |>
  filter(Date >= START_DATE, month(Date) %in% c(10:12, 1:4)) |>
  inner_join(PondTable |> select(location, pond), by = c("Backwater" = "location")) |>
  mutate(t = prim_floor(Date), k = week_of_month(Date)) |>
  filter(!is.na(t)) |>
  distinct(pond, t, k)

eff_wk  <- array(0L, dim = c(7, T_PRIM, N_SEC))
eff_hrs <- matrix(0, 7, T_PRIM)
for (r in seq_len(nrow(EffortWeekly))) {
  eff_wk[EffortWeekly$pond[r], EffortWeekly$t[r], EffortWeekly$k[r]] <- 1L
  eff_hrs[EffortWeekly$pond[r], EffortWeekly$t[r]] <-
    eff_hrs[EffortWeekly$pond[r], EffortWeekly$t[r]] + EffortWeekly$eff_hrs[r]
}
for (r in seq_len(nrow(ContactWeeks))) {
  eff_wk[ContactWeeks$pond[r], ContactWeeks$t[r], ContactWeeks$k[r]] <- 1L
}

# ---------------------------------------------------------------------------
# Scanner-week screening (2026-09-08; identical to model/YCB3S_data.R).
# Effort records can call a week "scanned" when the antenna produced no
# information. Two rules drop such weeks from availability (eff_wk -> K_mat)
# and from the monthly hours covariate; weeks with any contact are never
# dropped (y <= K must hold).
#  (a) low effort: hours < MIN_WK_FRAC x pond median weekly hours and no
#      contact of any fish (partial deployments).
#  (b) dead scanner: no contact of any fish although >= MIN_ALIVE tagged fish
#      were known alive and the pond-season weekly per-fish detection rate
#      (median over weeks with contacts) makes P(zero fish) < P0_MAX
#      (e.g. YCB Apr 2020, IP1 Oct-Nov 2021, IP5 Mar 2022 wk 1, IP3 Nov 2024 wk 1).
# ---------------------------------------------------------------------------
MIN_WK_FRAC <- 0.25
MIN_ALIVE   <- 20L
P0_MAX      <- 1e-3

ContactFishWeeks <- StudyBWContacts |>
  filter(!is.na(PITIndex), Date >= START_DATE, month(Date) %in% c(10:12, 1:4)) |>
  inner_join(PondTable |> select(location, pond), by = c("Backwater" = "location")) |>
  mutate(t = prim_floor(Date), k = week_of_month(Date)) |>
  filter(!is.na(t)) |>
  distinct(pond, t, k, PITIndex) |>
  dplyr::count(pond, t, k, name = "n_fish")

KnownAliveMonth <- Fish |>
  select(pond, entry_t, end_t, MaxScanDate) |>
  cross_join(Primary |> select(t, month_start)) |>
  filter(entry_t <= t, end_t >= t, !is.na(MaxScanDate), MaxScanDate >= month_start) |>
  dplyr::count(pond, t, name = "known_alive")

WeekScreen <- expand_grid(pond = 1:7, t = 1:T_PRIM, k = 1:N_SEC) |>
  mutate(avail = eff_wk[cbind(pond, t, k)]) |>
  filter(avail == 1L) |>
  left_join(EffortWeekly, by = c("pond", "t", "k")) |>
  left_join(ContactFishWeeks, by = c("pond", "t", "k")) |>
  left_join(KnownAliveMonth, by = c("pond", "t")) |>
  left_join(Primary |> select(t, season), by = "t") |>
  mutate(across(c(eff_hrs, n_fish, known_alive), ~ replace_na(.x, 0))) |>
  group_by(pond) |>
  mutate(med_hrs = median(eff_hrs[eff_hrs > 0])) |>
  group_by(pond, season) |>
  mutate(p_season = median((n_fish / known_alive)[n_fish > 0 & known_alive >= 5])) |>
  ungroup() |>
  mutate(
    P0 = (1 - p_season)^known_alive,
    low_effort = n_fish == 0 & eff_hrs < MIN_WK_FRAC * med_hrs,
    dead       = n_fish == 0 & known_alive >= MIN_ALIVE & !is.na(P0) & P0 < P0_MAX,
    drop       = low_effort | dead
  )

DroppedWeeks <- WeekScreen |>
  filter(drop) |>
  left_join(PondTable |> select(pond, Backwater), by = "pond") |>
  left_join(Primary |> select(t, month_start), by = "t") |>
  select(Backwater, pond, t, month_start, k, eff_hrs, known_alive, p_season, P0,
         low_effort, dead)
for (r in seq_len(nrow(DroppedWeeks))) {
  eff_wk[DroppedWeeks$pond[r], DroppedWeeks$t[r], DroppedWeeks$k[r]] <- 0L
  eff_hrs[DroppedWeeks$pond[r], DroppedWeeks$t[r]] <-
    eff_hrs[DroppedWeeks$pond[r], DroppedWeeks$t[r]] - DroppedWeeks$eff_hrs[r]
}
eff_hrs[eff_hrs < 1e-6] <- 0
cat("Scanner-week screening: dropped", nrow(DroppedWeeks), "weeks (",
    sum(DroppedWeeks$low_effort), "low effort,", sum(DroppedWeeks$dead & !DroppedWeeks$low_effort),
    "dead scanner ); pond-months left with no available week:",
    DroppedWeeks |> distinct(pond, t) |>
      filter(mapply(function(j, tt) all(eff_wk[j, tt, ] == 0L), pond, t)) |> nrow(), "\n")

# Detection covariate: standardized log hours. Months with contacts but no
# effort record get the pond's median non-zero hours to avoid a spurious 0.
for (j in 1:7) {
  has_contact <- apply(eff_wk[j, , , drop = FALSE], 2, function(m) any(m > 0))
  fill <- median(eff_hrs[j, eff_hrs[j, ] > 0])
  eff_hrs[j, has_contact & eff_hrs[j, ] == 0] <- fill
}
log_eff <- log(eff_hrs + 1)
eff_std <- (log_eff - mean(log_eff[eff_hrs > 0])) / sd(log_eff[eff_hrs > 0])
eff_std[eff_hrs == 0] <- 0

# ---------------------------------------------------------------------------
# Detection arrays y[i, t] (weeks detected) and K[i, t] (weeks available)
# ---------------------------------------------------------------------------
Detections <- StudyBWContacts |>
  filter(!is.na(PITIndex), Date >= START_DATE, month(Date) %in% c(10:12, 1:4)) |>
  inner_join(Fish |> select(PITIndex, location, fish, entry_date, rem_date),
             by = c("PITIndex", "Backwater" = "location")) |>
  filter(Date >= entry_date, is.na(rem_date) | Date <= rem_date) |>
  mutate(t = prim_floor(Date), k = week_of_month(Date)) |>
  filter(!is.na(t)) |>
  distinct(fish, t, k)

y_mat <- matrix(0L, N_FISH, T_PRIM)
det_counts <- Detections |> dplyr::count(fish, t)
y_mat[cbind(det_counts$fish, det_counts$t)] <- det_counts$n

# K: weeks available to fish i in month t (scanner running, fish at large)
wk_start_mat <- matrix(as.numeric(WeekGrid$wk_start), T_PRIM, N_SEC, byrow = TRUE)
wk_end_mat   <- matrix(as.numeric(WeekGrid$wk_end),   T_PRIM, N_SEC, byrow = TRUE)
K_mat <- matrix(0L, N_FISH, T_PRIM)
for (i in seq_len(N_FISH)) {
  j  <- Fish$pond[i]
  ts <- Fish$entry_t[i]:Fish$end_t[i]
  ent <- as.numeric(Fish$entry_date[i])
  rem <- if (Fish$removed[i]) as.numeric(Fish$rem_date[i]) else Inf
  avail <- eff_wk[j, ts, , drop = FALSE][1, , , drop = TRUE]
  avail <- matrix(avail, length(ts), N_SEC)
  avail <- avail == 1L & wk_end_mat[ts, , drop = FALSE] >= ent &
    wk_start_mat[ts, , drop = FALSE] <= rem
  K_mat[i, ts] <- rowSums(avail)
}
stopifnot(all(y_mat <= K_mat))

last_det <- apply(y_mat, 1, function(v) if (any(v > 0)) max(which(v > 0)) else 0L)
Fish$last_det_t <- as.integer(last_det)

# ---------------------------------------------------------------------------
# Stage observations at recapture (0 none, 1 juvenile, 2 adult)
# ---------------------------------------------------------------------------
StageObs <- StudyBWNFWG |>
  filter(event == "capture", FirstRecord == "no", !is.na(total_length),
         collection_date >= START_DATE) |>
  inner_join(Fish |> select(PITIndex, location, fish, species, entry_date, end_t),
             by = c("PITIndex", "location", "species")) |>
  mutate(obs_date = as.Date(collection_date), t = prim_floor(obs_date),
         stage = ifelse(total_length >= StageTL[species], 2L, 1L)) |>
  filter(!is.na(t), obs_date > entry_date, t <= end_t) |>
  group_by(fish, t) |>
  summarise(stage = max(stage), .groups = "drop")

stobs_mat <- matrix(0L, N_FISH, T_PRIM)
stobs_mat[cbind(StageObs$fish, StageObs$t)] <- StageObs$stage
# A fish cannot revert to juvenile: drop juvenile observations after an adult one
for (i in unique(StageObs$fish)) {
  ts_a <- which(stobs_mat[i, ] == 2L)
  if (length(ts_a) > 0) {
    later_j <- which(stobs_mat[i, ] == 1L & seq_len(T_PRIM) > min(ts_a))
    stobs_mat[i, later_j] <- 0L
  }
}
# Entry stage adult with later juvenile observation is a threshold artefact;
# trust the entry TL and drop those observations
adult_entry <- which(Fish$stage0 == 2L)
stobs_mat[adult_entry, ][stobs_mat[adult_entry, ] == 1L] <- 0L

# Known-alive constraint (code 3 = alive, stage unknown): any own-backwater
# contact, including May-Sep contacts that are not used as detections, proves
# the fish was alive through the primary period at or before that date
Fish <- Fish |>
  mutate(last_alive_t = if_else(MaxScanDate > Primary$month_end[T_PRIM], T_PRIM,
                                prim_floor(MaxScanDate)),
         last_alive_t = pmin(coalesce(last_alive_t, entry_t), end_t))
for (i in seq_len(N_FISH)) {
  if (Fish$last_alive_t[i] < Fish$entry_t[i]) next
  rng <- Fish$entry_t[i]:Fish$last_alive_t[i]
  rng <- rng[stobs_mat[i, rng] == 0L]
  stobs_mat[i, rng] <- 3L
}

# ---------------------------------------------------------------------------
# Post-release offset indicator for the transition starting at t
# ---------------------------------------------------------------------------
post_mat <- matrix(0L, N_FISH, T_PRIM)
for (i in seq_len(N_FISH)) {
  if (Fish$origin[i] == 3L) next
  cutoff <- Fish$entry_date[i] %m+% months(POST_MONTHS)
  ts <- Fish$entry_t[i]:Fish$end_t[i]
  post_mat[i, ts[Primary$month_start[ts] < cutoff]] <- 1L
}

# ---------------------------------------------------------------------------
# Netting events (pond x calendar month), by stage
# ---------------------------------------------------------------------------
NettingCSV <- read_csv("data/BWNettingEvents.csv", show_col_types = FALSE) |>
  mutate(EventMonth = ym(EventMonth), stage = ifelse(Stage == "A", 2L, 1L)) |>
  inner_join(PondTable |> select(location, pond), by = c("Backwater" = "location")) |>
  group_by(pond, EventMonth, stage) |>
  summarise(harvU_csv = sum(HarvestedUntagged) + sum(MortsUntagged),
            retU_pool = sum(ReturnedUntagged * ReturnedInPool),
            retU_excl = sum(ReturnedUntagged * !ReturnedInPool),
            .groups = "drop")

CaptureRecords <- StudyBWNFWG |>
  filter(event == "capture", collection_date >= START_DATE) |>
  inner_join(PondTable |> select(location, pond, pond_species = species), by = "location") |>
  filter(species == pond_species) |>
  mutate(cap_date = as.Date(collection_date),
         EventMonth = floor_date(cap_date, "month"),
         stage_cap = case_when(!is.na(total_length) ~ ifelse(total_length >= StageTL[species], 2L, 1L),
                               FirstRecord == "yes" ~ 1L,
                               TRUE ~ 2L),
         prev_tagged = FirstRecord == "no",
         retained = disposition == "retained") |>
  group_by(pond, EventMonth, PITIndex) |>
  summarise(stage_cap = max(stage_cap), prev_tagged = any(prev_tagged),
            retained = any(retained), cap_date = min(cap_date),
            TL = suppressWarnings(max(total_length, na.rm = TRUE)), .groups = "drop") |>
  mutate(TL = ifelse(is.finite(TL), TL, NA))

NettingNFWG <- CaptureRecords |>
  group_by(pond, EventMonth, stage = stage_cap) |>
  summarise(m = sum(prev_tagged),
            newtag = sum(!prev_tagged & !retained),
            harv_newtag = sum(!prev_tagged & retained),
            remT = sum(prev_tagged & retained),
            event_date = min(cap_date),
            meanTL_new = mean(TL[!prev_tagged], na.rm = TRUE),
            .groups = "drop")

Netting <- NettingNFWG |>
  full_join(NettingCSV, by = c("pond", "EventMonth", "stage")) |>
  mutate(across(c(m, newtag, harv_newtag, remT, harvU_csv, retU_pool, retU_excl),
                ~ replace_na(as.integer(.x), 0L))) |>
  group_by(pond, EventMonth) |>
  mutate(event_date = if (all(is.na(event_date))) EventMonth[1] else min(event_date, na.rm = TRUE)) |>
  ungroup() |>
  mutate(
    harvU = harvU_csv + harv_newtag,              # untagged fish removed
    n     = m + newtag + harvU + retU_pool,       # handled fish in the pool
    remU  = newtag + harvU,                       # leave the untagged pool
    t_e   = prim_floor(event_date),
    season = (t_e - 1L) %/% N_MOY + 1L
  ) |>
  filter(!is.na(t_e)) |>
  arrange(pond, event_date, stage) |>
  group_by(pond, season, stage) |>
  mutate(remU_before = cumsum(remU) - remU) |>   # untagged removed earlier in season
  ungroup() |>
  mutate(e = row_number())

N_EVENTS <- nrow(Netting)

# Untagged removals per pond x season x stage (for the U process)
remU_arr <- array(0L, dim = c(7, N_SEASONS, 2))
for (r in seq_len(N_EVENTS)) {
  remU_arr[Netting$pond[r], Netting$season[r], Netting$stage[r]] <-
    remU_arr[Netting$pond[r], Netting$season[r], Netting$stage[r]] + Netting$remU[r]
}

# Tagged removals with stage at removal (reporting)
RemovedTagged <- Fish |>
  filter(removed) |>
  mutate(season = (rem_t - 1L) %/% N_MOY + 1L,
         EventMonth = floor_date(rem_date, "month")) |>
  left_join(CaptureRecords |> select(pond, PITIndex, EventMonth, stage_cap),
            by = c("pond", "PITIndex", "EventMonth")) |>
  mutate(rem_stage = coalesce(stage_cap, stage0)) |>
  select(fish, PITIndex, pond, species, season, rem_date, rem_type, rem_TL, rem_stage)

# ---------------------------------------------------------------------------
# Known die-off events (pond x calendar month indicator), ported from
# model/YCB3S_data.R (2026-09-10). Hand-maintained table of independently
# documented mass-mortality events (not inferred from the detection data).
# Applied as a fixed monthly logit-survival offset dE[pond] in rd_code/rdPond,
# shared across juvenile and adult. Under this model's m7 calendar,
# prim_floor() collapses May-Sep dates to April, so a summer event only
# partially covered by the Oct-Apr primaries flags the whole April primary
# (same approximation used in the XYTE three-stage model before its m12
# calendar was added); nEvt below is computed from calendar months directly
# so the untagged-pool annual survival discount is exact regardless.
# ---------------------------------------------------------------------------
DieOffsRaw <- read_csv(DIEOFF_FILE, show_col_types = FALSE) |>
  mutate(StartMonth = ym(StartMonth), EndMonth = ym(EndMonth)) |>
  inner_join(PondTable |> select(location, pond, pond_species = species), by = c("Backwater" = "location")) |>
  filter(Species == pond_species)

Evt  <- matrix(0L, 7, T_PRIM)
nEvt <- matrix(0L, 7, N_SEASONS)
if (nrow(DieOffsRaw) > 0) {
  for (r in seq_len(nrow(DieOffsRaw))) {
    mos <- seq(DieOffsRaw$StartMonth[r], DieOffsRaw$EndMonth[r], by = "month")
    ts  <- prim_floor(mos)
    Evt[DieOffsRaw$pond[r], ts[!is.na(ts)]] <- 1L
    ys <- season_of_date(mos)
    for (yy in unique(ys)) nEvt[DieOffsRaw$pond[r], yy] <- nEvt[DieOffsRaw$pond[r], yy] + sum(ys == yy)
  }
}

# ---------------------------------------------------------------------------
# Length-weight relationships for biomass (log W = a + b log TL)
# ---------------------------------------------------------------------------
LengthWeight <- StudyBWNFWG |>
  filter(!is.na(weight), !is.na(total_length), weight > 0, total_length > 50) |>
  group_by(species) |>
  summarise(fit = list(lm(log(weight) ~ log(total_length))), n = n(), .groups = "drop") |>
  mutate(a = sapply(fit, function(f) unname(coef(f)[1])),
         b = sapply(fit, function(f) unname(coef(f)[2])),
         sigma = sapply(fit, function(f) summary(f)$sigma)) |>
  select(species, n, a, b, sigma)

# ---------------------------------------------------------------------------
# Stocking table (reporting)
# ---------------------------------------------------------------------------
Stocking <- Fish |>
  filter(origin == 1L) |>
  group_by(pond, species, EventMonth = floor_date(first_date, "month")) |>
  summarise(n = n(), meanTL = mean(total_length, na.rm = TRUE),
            n_juv = sum(stage0 == 1L), n_adult = sum(stage0 == 2L), .groups = "drop")

# ---------------------------------------------------------------------------
# GIEL size-tied maturation hazard lookup table (psi_lookup), from whichever
# track model/GIEL_MaturationHazard_Build.R was last run for (MATURATION_TRACK
# above). 20mm TL bins x 24 months-since-entry, clipped away from exact 0/1.
# ---------------------------------------------------------------------------
MaturationHazard <- new.env()
load(MATURATION_HAZARD_FILE, envir = MaturationHazard)
stopifnot(identical(MaturationHazard$TL_BIN_EDGES, GIEL_TLBIN_EDGES))
psi_lookup   <- MaturationHazard$psi_lookup       # N_TLBIN_GIEL x DELTA_T_MAX
DELTA_T_MAX  <- MaturationHazard$DELTA_T_MAX
cat("\nGIEL maturation hazard table loaded from", MATURATION_HAZARD_FILE,
    "(track", MaturationHazard$TRACK, "):", nrow(psi_lookup), "TL bins x",
    ncol(psi_lookup), "months\n")

# ---------------------------------------------------------------------------
# NIMBLE constants / data (fish contiguous within pond)
# ---------------------------------------------------------------------------
pond_off <- c(0L, cumsum(tabulate(Fish$pond, nbins = 7)))
stopifnot(all(diff(Fish$pond) >= 0))

RD_constants <- list(
  NPOND = 7L, NSP = 2L, T = T_PRIM, NS = N_SEASONS, NMOY = N_MOY, NFISH = N_FISH,
  NEV = N_EVENTS,
  sp = PondTable$sp, season = Primary$season, moy = Primary$moy, gap = Primary$gap,
  tOct = Primary$t[Primary$moy == 1],
  pond_off = pond_off,
  entry_t = Fish$entry_t, end_t = Fish$end_t, removed = as.integer(Fish$removed),
  stage0 = Fish$stage0, origin = Fish$origin,
  count_from = Fish$entry_t + as.integer(Fish$origin == 2L),
  eff = eff_std, Evt = Evt, nEvt = nEvt,
  ev_pond = Netting$pond, ev_t = Netting$t_e, ev_season = Netting$season,
  ev_stage = Netting$stage, ev_n = Netting$n, ev_remU_before = Netting$remU_before,
  remU = remU_arr,
  U_known = ifelse(PondTable$pond == 1, 0L, 1L),   # IP ponds start with no untagged fish
  tlbin = Fish$TL_entry_bin, hasTL = as.integer(Fish$hasTL),
  psi_lookup = psi_lookup, N_TLBIN = N_TLBIN_GIEL, DELTA_T_MAX = DELTA_T_MAX
)

RD_data <- list(
  y = y_mat, K = K_mat, stobs = stobs_mat, post = post_mat,
  ev_m = Netting$m
)

BUILD_INFO_INPUTS <- c("model/BWRobustDesign_data.R", "data/ReportingData.RData",
                       "data/BWNettingEvents.csv", DIEOFF_FILE, MATURATION_HAZARD_FILE)
BUILD_INFO <- list(
  built = Sys.time(),
  inputs = BUILD_INFO_INPUTS,
  input_mtime = file.info(BUILD_INFO_INPUTS)$mtime,
  DIEOFF_FILE = DIEOFF_FILE,
  MATURATION_TRACK = MATURATION_TRACK,
  MATURATION_HAZARD_FILE = MATURATION_HAZARD_FILE
)

save(Fish, Primary, PondTable, Netting, NettingCSV, Stocking, Removals, RemovedTagged,
     SameEventRemovals, SpeciesMismatch, LengthWeight, DroppedWeeks, DieOffsRaw, Evt, nEvt,
     y_mat, K_mat, stobs_mat, post_mat, eff_wk, eff_hrs, eff_std,
     psi_lookup, GIEL_TLBIN_EDGES, N_TLBIN_GIEL, DELTA_T_MAX, MATURATION_TRACK,
     RD_constants, RD_data, BUILD_INFO,
     START_DATE, FIRST_FY, N_SEASONS, N_MOY, T_PRIM, N_SEC, POST_MONTHS, StageTL,
     prim_floor, fy_of_date, season_of_date,
     file = "data/BWRobustDesign_data.RData")

# ---------------------------------------------------------------------------
# Reconciliation printout
# ---------------------------------------------------------------------------
cat("\nFish by pond / origin (1 stocked, 2 netting-tagged, 3 established) / stage at entry:\n")
print(Fish |> dplyr::count(pond, origin, stage0) |>
        pivot_wider(names_from = c(origin, stage0), values_from = n, values_fill = 0,
                    names_glue = "orig{origin}_stage{stage0}"))
cat("\nRemoved tagged fish by pond and type:\n")
print(Fish |> filter(removed) |> dplyr::count(pond, rem_type))
cat("\nSame-event tagged-and-removed fish (counted as untagged harvest):",
    nrow(SameEventRemovals), "\n")
cat("\nNetting events:\n")
print(Netting |>
        select(pond, EventMonth, stage, t_e, season, m, newtag, harvU, retU_pool, retU_excl, n, remU) |>
        as.data.frame())
cat("\nDetection summary: fish-months with K>0:", sum(K_mat > 0),
    "; with y>0:", sum(y_mat > 0), "; stage observations:", sum(stobs_mat %in% 1:2),
    "; known-alive months (code 3):", sum(stobs_mat == 3L), "\n")
cat("\nLength-weight fits:\n"); print(LengthWeight)
cat("\nGIEL maturation hazard track:", MATURATION_TRACK, "(", MaturationHazard$track_desc, ")\n")
cat("\nDie-off events:", nrow(DieOffsRaw), "record(s); event pond-months (Evt):", sum(Evt),
    "; event calendar-months by pond x season (nEvt):\n")
print(DieOffsRaw |> select(Backwater, Species, StartMonth, EndMonth, Cause) |> as.data.frame())
if (any(nEvt > 0)) print(which(nEvt > 0, arr.ind = TRUE) |> as.data.frame() |>
      setNames(c("pond", "season")) |> mutate(Backwater = PondTable$Backwater[pond], FY = FIRST_FY + season - 1,
                                              n_months = nEvt[cbind(pond, season)]) |>
      select(Backwater, FY, n_months) |> as.data.frame())
