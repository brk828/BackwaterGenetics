# Data preparation for the three-stage, density-dependent robust-design model
# (Phase 1 of .posit/assistant/plans/2026-09-04-1335-ycb-three-stage-...md).
# B. Kesner / Posit Assistant, September 2026
#
# Forked from model/BWRobustDesign_data.R. Differences:
#   * three alive stages per species (S = small juvenile, L = large juvenile,
#     A = adult) using StageBreaks; stobs codes 1/2/3 = observed S/L/A,
#     4 = known alive (stage unknown)
#   * MODEL_PONDS selects which backwaters are built (YCB only for Phases 1-3;
#     Phase 4 default c(1, 2, 4, 5) = all XYTE ponds; set before sourcing to
#     override, e.g. MODEL_PONDS <- 1L to rebuild the YCB-only inputs)
#   * per pond x season covariates: observed recruit index (untagged juveniles
#     handled at the Oct-Dec netting, per acre), recruit size split (n0, nL)
#     for the density-dependent growth sub-model, mean weight by stage
#   * pond surface areas from data/BWPondAreas.csv (placeholder 1 if missing)
#   * known die-off events (data/BWDieOffEvents.csv) -> Evt[pond, t] indicator
#     and nEvt[pond, season] calendar-month count, feeding a fixed dE[pond]
#     logit-survival offset in rd3_code (2026-09-09)
#   * mtime guard so the saved inputs cannot silently go stale
#
# Calendar (Phase 3b): CALENDAR = "m7" (default) uses the Oct-Apr robust design
# (7 primaries per season, 6-month Apr->Oct gap; May-Sep contacts are only
# known-alive constraints collapsed to April). CALENDAR = "m12" makes every
# month a primary period (12 per season, gap == 1) so May-Sep scanning becomes
# detection occasions. Set CALENDAR (and optionally OUT_FILE) before sourcing.
#
# Output: data/YCB3S_data[<pondtag>][_m12].RData, where <pondtag> is empty for
# the YCB-only build (so the Phase 1-3 files keep their names), "_xyte" for
# c(1, 2, 4, 5), and "_p<ponds>" otherwise.

source("LabFunctions.R")
packages(dplyr)
packages(lubridate)
packages(tidyr)
packages(stringr)
packages(readr)
source("model/RecruitYOYClassification.R")

load("data/ReportingData.RData")

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
if (!exists("MODEL_PONDS")) MODEL_PONDS <- c(1L, 2L, 4L, 5L)   # 1 = YCB, 2/4/5 = IP1/IP3/IP4
MODEL_PONDS <- sort(as.integer(MODEL_PONDS))
if (!exists("CALENDAR")) CALENDAR <- "m7"
stopifnot(CALENDAR %in% c("m7", "m12"))
# Die-off event windows: default is the working table (v11, 2026-09-09: IP1
# Jul-Sep, IP3/IP4 Jun-Sep). Set DIEOFF_FILE <- "data/BWDieOffEvents_ip3MaySep.csv"
# before sourcing to build the IP3 May-Sep sensitivity variant instead.
if (!exists("DIEOFF_FILE")) DIEOFF_FILE <- "data/BWDieOffEvents.csv"
POND_TAG <- if (identical(MODEL_PONDS, 1L)) "" else
  if (identical(MODEL_PONDS, c(1L, 2L, 4L, 5L))) "_xyte" else
  paste0("_p", paste(MODEL_PONDS, collapse = ""))
if (!exists("OUT_FILE"))
  OUT_FILE <- paste0("data/YCB3S_data", POND_TAG, if (CALENDAR == "m12") "_m12", ".RData")

START_DATE  <- as.Date("2016-10-01")
FIRST_FY    <- 2017
N_SEASONS   <- 10
N_MOY       <- if (CALENDAR == "m12") 12L else 7L    # primaries per season (moy 1 = Oct)
DET_MONTHS  <- if (CALENDAR == "m12") 1:12 else c(10:12, 1:4)  # calendar months that are detection occasions
T_PRIM      <- N_SEASONS * N_MOY
N_SEC       <- 4
POST_MONTHS <- 6
N_STAGE     <- 3L
# S < break1 <= L < break2 <= A. GIEL breaks are provisional (Phase 5).
StageBreaks <- list(XYTE = c(300, 350), GIEL = c(200, 250))
RECRUIT_MONTHS <- 10:12                 # fall netting window for the recruit index
CSV_SMALL_MAX  <- 250                   # CSV untagged juveniles with MeanTL below this count as stage S in n0

PondTable <- tibble(
  pond      = 1:7,
  location  = c("Yuma Cove backwater", paste0("IPCA (Pond ", 1:6, ")")),
  LID       = c(592L, 1043:1048),
  Backwater = c("YCB", paste0("IP", 1:6)),
  species   = c("XYTE", "XYTE", "GIEL", "XYTE", "XYTE", "GIEL", "GIEL"),
  sp        = c(1L, 1L, 2L, 1L, 1L, 2L, 2L)
) |>
  filter(pond %in% MODEL_PONDS) |>
  mutate(pond_all = pond, pond = row_number())
N_POND <- nrow(PondTable)

stage_of <- function(TL, species) {
  b <- do.call(rbind, StageBreaks[species])
  as.integer(1L + (TL >= b[, 1]) + (TL >= b[, 2]))
}

# ---------------------------------------------------------------------------
# Primary period table and date helpers
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
    # months from this primary to the next: 1 within the season; the last
    # primary of a season jumps to October (6 for m7, 1 for m12)
    gap         = ifelse(moy == N_MOY, 13L - N_MOY, 1L),
    summer      = month %in% 5:9
  )

fy_of_date <- function(d) ifelse(month(d) > 9, year(d) + 1L, year(d))
season_of_date <- function(d) fy_of_date(d) - FIRST_FY + 1L

# Primary period at or before a date. Under m7, May-Sep dates collapse to April
# (pmin); under m12 pmin is a no-op and every date maps to its own month.
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
# Pond areas
# ---------------------------------------------------------------------------
Areas <- read_csv("data/BWPondAreas.csv", show_col_types = FALSE) |>
  select(location = Location, Acres)
PondTable <- PondTable |> left_join(Areas, by = "location")
if (any(is.na(PondTable$Acres))) {
  warning("Surface area missing for: ",
          paste(PondTable$Backwater[is.na(PondTable$Acres)], collapse = ", "),
          ". Using Area = 1 for ALL modeled ponds so densities are per pond, not per acre.")
  PondTable$Area <- 1
} else {
  PondTable$Area <- PondTable$Acres
}

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
  warning(nrow(SpeciesMismatch), " fish excluded: species does not match pond species")
  Fish <- Fish |> filter(species == pond_species)
}

Fish <- Fish |>
  left_join(FirstContactInWindow, by = c("PITIndex", "location")) |>
  left_join(Removals, by = c("PITIndex", "location")) |>
  mutate(
    origin = case_when(
      first_date < START_DATE ~ 3L,
      event == "stocking"     ~ 1L,
      TRUE                    ~ 2L
    ),
    entry_date = if_else(origin == 3L, first_contact, first_date)
  ) |>
  filter(!is.na(entry_date)) |>
  mutate(same_event_removal = !is.na(rem_date) & origin == 2L & rem_date <= first_date + 3)

SameEventRemovals <- Fish |> filter(same_event_removal)
Fish <- Fish |> filter(!same_event_removal)

# Stage at entry (3 classes)
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
    recent_obs = origin == 3L & !is.na(obs_date) & obs_date >= entry_date %m-% months(12),
    TL_use = case_when(origin == 3L & recent_obs ~ TL,
                       origin == 3L ~ NA_real_,
                       TRUE ~ total_length),
    stage0 = case_when(
      origin == 3L & !recent_obs ~ 3L,           # established, no recent TL: adult
      !is.na(TL_use)             ~ stage_of(TL_use, species),
      origin == 1L               ~ 3L,           # stocked, TL missing
      TRUE                       ~ 1L            # netting-tagged, TL missing
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

EffortWeekly <- EffortRecords |>
  cross_join(WeekGrid) |>
  mutate(ov_days = as.numeric(pmin(ret, wk_end) - pmax(dep, wk_start)) + 1) |>
  filter(ov_days > 0) |>
  group_by(pond, t, k) |>
  summarise(eff_hrs = sum(ov_days * hrs_per_day), .groups = "drop")

ContactWeeks <- StudyBWContacts |>
  filter(Date >= START_DATE, month(Date) %in% DET_MONTHS) |>
  inner_join(PondTable |> select(location, pond), by = c("Backwater" = "location")) |>
  mutate(t = prim_floor(Date), k = week_of_month(Date)) |>
  filter(!is.na(t)) |>
  distinct(pond, t, k)

eff_wk  <- array(0L, dim = c(N_POND, T_PRIM, N_SEC))
eff_hrs <- matrix(0, N_POND, T_PRIM)
for (r in seq_len(nrow(EffortWeekly))) {
  eff_wk[EffortWeekly$pond[r], EffortWeekly$t[r], EffortWeekly$k[r]] <- 1L
  eff_hrs[EffortWeekly$pond[r], EffortWeekly$t[r]] <-
    eff_hrs[EffortWeekly$pond[r], EffortWeekly$t[r]] + EffortWeekly$eff_hrs[r]
}
for (r in seq_len(nrow(ContactWeeks))) {
  eff_wk[ContactWeeks$pond[r], ContactWeeks$t[r], ContactWeeks$k[r]] <- 1L
}

# ---------------------------------------------------------------------------
# Scanner-week screening (2026-09-08). Effort records can call a week
# "scanned" when the antenna produced no information. Two rules drop such
# weeks from availability (eff_wk -> K_mat) and from the monthly hours
# covariate; weeks with any contact are never dropped (y <= K must hold).
#  (a) low effort: hours < MIN_WK_FRAC x pond median weekly hours and no
#      contact of any fish (partial deployments, e.g. YCB Apr 2020: 77.7 h
#      over the month, 205 fish at large, zero contacts).
#  (b) dead scanner: no contact of any fish although >= MIN_ALIVE tagged fish
#      were known alive and the pond-season weekly per-fish detection rate
#      (median over weeks with contacts) makes P(zero fish) < P0_MAX
#      (e.g. IP1 Oct-Nov 2021, IP5 Mar 2022 wk 1, IP3 Nov 2024 wk 1).
# ---------------------------------------------------------------------------
MIN_WK_FRAC <- 0.25
MIN_ALIVE   <- 20L
P0_MAX      <- 1e-3

ContactFishWeeks <- StudyBWContacts |>
  filter(!is.na(PITIndex), Date >= START_DATE, month(Date) %in% DET_MONTHS) |>
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

WeekScreen <- expand_grid(pond = 1:N_POND, t = 1:T_PRIM, k = 1:N_SEC) |>
  mutate(avail = eff_wk[cbind(pond, t, k)]) |>
  filter(avail == 1L) |>
  left_join(EffortWeekly, by = c("pond", "t", "k")) |>
  left_join(ContactFishWeeks, by = c("pond", "t", "k")) |>
  left_join(KnownAliveMonth, by = c("pond", "t")) |>
  left_join(Primary |> select(t, season, summer), by = "t") |>
  mutate(across(c(eff_hrs, n_fish, known_alive), ~ replace_na(.x, 0))) |>
  group_by(pond) |>
  mutate(med_hrs = median(eff_hrs[eff_hrs > 0])) |>
  # summer and Oct-Apr detection rates differ, so the reference rate for the
  # dead-scanner test is taken within pond x season x summer (no-op under m7)
  group_by(pond, season, summer) |>
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
  select(Backwater, pond, t, month_start, summer, k, eff_hrs, known_alive, p_season, P0,
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

for (j in 1:N_POND) {
  has_contact <- apply(eff_wk[j, , , drop = FALSE], 2, function(m) any(m > 0))
  fill <- median(eff_hrs[j, eff_hrs[j, ] > 0])
  eff_hrs[j, has_contact & eff_hrs[j, ] == 0] <- fill
}
log_eff <- log(eff_hrs + 1)
eff_std <- (log_eff - mean(log_eff[eff_hrs > 0])) / sd(log_eff[eff_hrs > 0])
eff_std[eff_hrs == 0] <- 0
# Partial-week deployments (e.g. YCB Apr 2021, 22 h) are extreme leverage
# points on the effort slope; floor the covariate
EFF_FLOOR <- -3
eff_std <- pmax(eff_std, EFF_FLOOR)

# ---------------------------------------------------------------------------
# Detection arrays y[i, t] (weeks detected) and K[i, t] (weeks available)
# ---------------------------------------------------------------------------
Detections <- StudyBWContacts |>
  filter(!is.na(PITIndex), Date >= START_DATE, month(Date) %in% DET_MONTHS) |>
  inner_join(Fish |> select(PITIndex, location, fish, entry_date, rem_date),
             by = c("PITIndex", "Backwater" = "location")) |>
  filter(Date >= entry_date, is.na(rem_date) | Date <= rem_date) |>
  mutate(t = prim_floor(Date), k = week_of_month(Date)) |>
  filter(!is.na(t)) |>
  distinct(fish, t, k)

y_mat <- matrix(0L, N_FISH, T_PRIM)
det_counts <- Detections |> dplyr::count(fish, t)
y_mat[cbind(det_counts$fish, det_counts$t)] <- det_counts$n

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

Fish$last_det_t <- as.integer(apply(y_mat, 1, function(v) if (any(v > 0)) max(which(v > 0)) else 0L))

# ---------------------------------------------------------------------------
# Stage observations at recapture (0 none, 1 S, 2 L, 3 A, 4 known alive)
# ---------------------------------------------------------------------------
StageObs <- StudyBWNFWG |>
  filter(event == "capture", FirstRecord == "no", !is.na(total_length),
         collection_date >= START_DATE) |>
  inner_join(Fish |> select(PITIndex, location, fish, species, entry_date, end_t),
             by = c("PITIndex", "location", "species")) |>
  mutate(obs_date = as.Date(collection_date), t = prim_floor(obs_date),
         stage = stage_of(total_length, species)) |>
  filter(!is.na(t), obs_date > entry_date, t <= end_t) |>
  group_by(fish, t) |>
  summarise(stage = max(stage), .groups = "drop")

stobs_mat <- matrix(0L, N_FISH, T_PRIM)
stobs_mat[cbind(StageObs$fish, StageObs$t)] <- StageObs$stage
# Stages are monotone: drop any observation lower than an earlier one or lower
# than the entry stage (threshold measurement artefacts)
for (i in unique(StageObs$fish)) {
  v <- stobs_mat[i, ]
  run_max <- Fish$stage0[i]
  for (t in Fish$entry_t[i]:Fish$end_t[i]) {
    if (v[t] > 0L) {
      if (v[t] < run_max) v[t] <- 0L else run_max <- v[t]
    }
  }
  stobs_mat[i, ] <- v
}
N_STAGE_OBS <- sum(stobs_mat %in% 1:3)

# Known-alive constraint (code 4)
Fish <- Fish |>
  mutate(last_alive_t = if_else(MaxScanDate > Primary$month_end[T_PRIM], T_PRIM,
                                prim_floor(MaxScanDate)),
         last_alive_t = pmin(coalesce(last_alive_t, entry_t), end_t))
for (i in seq_len(N_FISH)) {
  if (Fish$last_alive_t[i] < Fish$entry_t[i]) next
  rng <- Fish$entry_t[i]:Fish$last_alive_t[i]
  rng <- rng[stobs_mat[i, rng] == 0L]
  stobs_mat[i, rng] <- 4L
}

# ---------------------------------------------------------------------------
# Post-release offset indicator
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
NettingCSVraw <- read_csv("data/BWNettingEvents.csv", show_col_types = FALSE) |>
  mutate(EventMonth = ym(EventMonth)) |>
  inner_join(PondTable |> select(location, pond, pond_species = species), by = c("Backwater" = "location")) |>
  filter(Species == pond_species)

# CSV rows are staged J/A with a mean TL; assign J rows to S or L by MeanTL
# (S when MeanTL is below the first break or missing)
NettingCSV <- NettingCSVraw |>
  mutate(stage = case_when(
    Stage == "A" ~ 3L,
    !is.na(MeanTL) & MeanTL >= sapply(Species, function(s) StageBreaks[[s]][1]) ~ 2L,
    TRUE ~ 1L)) |>
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
         stage_cap = case_when(!is.na(total_length) ~ stage_of(total_length, species),
                               FirstRecord == "yes" ~ 1L,
                               TRUE ~ 3L),
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
    harvU = harvU_csv + harv_newtag,
    n     = m + newtag + harvU + retU_pool,
    remU  = newtag + harvU,
    t_e   = prim_floor(event_date),
    season = (t_e - 1L) %/% N_MOY + 1L
  ) |>
  filter(!is.na(t_e)) |>
  arrange(pond, event_date, stage) |>
  group_by(pond, season, stage) |>
  mutate(remU_before = cumsum(remU) - remU) |>
  ungroup() |>
  mutate(e = row_number())

N_EVENTS <- nrow(Netting)

remU_arr <- array(0L, dim = c(N_POND, N_SEASONS, N_STAGE))
for (r in seq_len(N_EVENTS)) {
  remU_arr[Netting$pond[r], Netting$season[r], Netting$stage[r]] <-
    remU_arr[Netting$pond[r], Netting$season[r], Netting$stage[r]] + Netting$remU[r]
}

RemovedTagged <- Fish |>
  filter(removed) |>
  mutate(season = (rem_t - 1L) %/% N_MOY + 1L,
         EventMonth = floor_date(rem_date, "month")) |>
  left_join(CaptureRecords |> select(pond, PITIndex, EventMonth, stage_cap),
            by = c("pond", "PITIndex", "EventMonth")) |>
  mutate(rem_stage = coalesce(stage_cap, stage0)) |>
  select(fish, PITIndex, pond, species, season, rem_date, rem_type, rem_TL, rem_stage)

# ---------------------------------------------------------------------------
# Recruit index and recruit size split (pond x season)
# ---------------------------------------------------------------------------
# Recruit index: untagged juvenile fish (S or L, i.e. below the adult break)
# handled at the Oct-Dec netting of season y, from NFWG first captures plus the
# CSV untagged counts (harvested, returned whether or not in the pool, morts).
# This is the observed density covariate; it is per acre once areas are known.
MIN_YOY_N <- 15L   # cohorts smaller than this fall back to a default (see below / RecruitYOYClassification.R)

# Pond x season flag: does the bulk-harvested/returned untagged YOY record
# for this season (data/BWNettingEvents.csv, Stage == "J") DOMINATE over the
# individually-tagged untagged cohort handled at the same fall netting (i.e.
# n_csv_bulk > n_tagged_raw)? A large dominant bulk record means the true
# recruit cohort was already comprehensively handled that way (that's
# exactly why it wasn't individually tagged -- too many fish, no time to tag
# each one), so individually tagged fish from the same event are
# presumptively established/carryover fish that finally got caught, not this
# season's YOY. A magnitude comparison (not just "any bulk record exists") is
# required: FY2017's CSV row is just 7 incidental morts alongside a real
# 113-fish individually-tagged recruit cohort and must NOT trigger the
# override, whereas FY2020's 2,877-fish bulk harvest (mean TL 152mm) against
# only 52 individually-tagged fish (all TL 330-585mm, clearly established
# adults) must. This flips the classifier default below from "assume YOY" to
# "assume not YOY" for dominant-bulk pond-seasons, AND (further below) caps
# any classified-YOY cluster to genuinely small TL in those same
# pond-seasons regardless of sample size (2026-09-24 fix, Part B; see
# AGENTS.md).
NettingCSV_bulk <- NettingCSVraw |>
  filter(Stage == "J", month(EventMonth) %in% RECRUIT_MONTHS) |>
  mutate(season = season_of_date(EventMonth),
         n_csv_evt = HarvestedUntagged + ReturnedUntagged + MortsUntagged) |>
  group_by(pond, season) |>
  summarise(n_csv_bulk = sum(n_csv_evt), .groups = "drop")

TaggedCohortN <- CaptureRecords |>
  filter(!prev_tagged, month(cap_date) %in% RECRUIT_MONTHS) |>
  mutate(season = season_of_date(cap_date)) |>
  group_by(pond, season) |>
  summarise(n_tagged_raw = n(), .groups = "drop")

BulkYOYSeasons <- NettingCSV_bulk |>
  left_join(TaggedCohortN, by = c("pond", "season")) |>
  mutate(n_tagged_raw = replace_na(n_tagged_raw, 0L)) |>
  filter(n_csv_bulk > n_tagged_raw) |>
  distinct(pond, season) |>
  mutate(has_bulk_csv = TRUE)

# A netting cohort's untagged first-captures are not all this season's
# young-of-year (YOY) recruits -- some are older untagged fish missed by a
# prior netting, or long-established fish first captured now. Classify each
# pond x season cohort with a Gaussian mixture on TL (mirrors the GIEL bimodal
# reporting table in PopulationMonitoring.qmd, generalized in
# model/RecruitYOYClassification.R) and keep only the YOY-classified fish
# before computing the recruit density / size-split covariates below.
#
# Classification runs on the FULL untagged first-capture cohort (any TL), not
# pre-filtered to sub-adult TL -- the stage_cap < 3L restriction is applied
# AFTER classification instead. Pre-filtering before classifying truncates
# away the very evidence (large, clearly-established fish) that lets the
# mixture fit recognize a cohort as non-YOY; a truncated remainder can also
# fall below MIN_YOY_N and default to "all YOY" even when the untruncated
# cohort is obviously bimodal and entirely non-YOY (2026-09-24 fix, Part A;
# see AGENTS.md for the FY2020 and FY2018 cases this corrected).
RecruitYOYClassified <- CaptureRecords |>
  filter(!prev_tagged, month(cap_date) %in% RECRUIT_MONTHS) |>
  mutate(season = season_of_date(cap_date)) |>
  left_join(BulkYOYSeasons, by = c("pond", "season")) |>
  mutate(yoy_default = !coalesce(has_bulk_csv, FALSE)) |>
  classify_yoy_df(TL, group_vars = c("pond", "season"), min_n = MIN_YOY_N,
                   default_col = "yoy_default") |>
  mutate(
    # The small-n default above only covers cohorts too small to fit a
    # mixture at all; it does NOT prevent the mixture from finding an
    # internal "smaller" cluster that is still not YOY-sized (FY2020: n = 52
    # >= MIN_YOY_N, mclust fits two established-fish modes at 343mm and
    # 472mm, and the smaller one is labeled is_yoy by construction even
    # though neither is close to that season's true 152mm bulk-harvested
    # recruit mean). When a bulk YOY record exists for this pond-season, only
    # trust a classified-YOY cluster if its own mean TL is genuinely small
    # (< CSV_SMALL_MAX, the same threshold already used to treat bulk-CSV
    # counts as stage S) -- otherwise it is just the smaller of two
    # established-fish modes and should not be counted as a recruit.
    is_yoy = is_yoy & !(coalesce(has_bulk_csv, FALSE) & mu1 >= CSV_SMALL_MAX)
  )

RecruitNFWG <- RecruitYOYClassified |>
  filter(is_yoy, stage_cap < 3L) |>
  group_by(pond, season) |>
  summarise(n_tagged_juv = n(),
            n0 = sum(!is.na(TL)),                       # with a TL, eligible for the size split
            nL = sum(!is.na(TL) & stage_cap == 2L),
            .groups = "drop")

RecruitCSV <- NettingCSVraw |>
  filter(Stage == "J", month(EventMonth) %in% RECRUIT_MONTHS) |>
  mutate(season = season_of_date(EventMonth),
         n_csv = HarvestedUntagged + ReturnedUntagged + MortsUntagged,
         # counted in the size split only when the mean TL clearly places them in S
         n0_csv = ifelse(!is.na(MeanTL) & MeanTL < CSV_SMALL_MAX, n_csv, 0L)) |>
  group_by(pond, season) |>
  summarise(n_csv = sum(n_csv), n0_csv = sum(n0_csv), .groups = "drop")

# A pond-season counts as "netted" (and so enters the density standardization
# with rec_total possibly 0) if it had any fall netting event in that pond
FallNetted <- Netting |>
  filter(month(EventMonth) %in% RECRUIT_MONTHS) |>
  distinct(pond, season) |>
  mutate(fall_netted = TRUE)

RecruitIndex <- expand_grid(pond = 1:N_POND, season = 1:N_SEASONS) |>
  left_join(RecruitNFWG, by = c("pond", "season")) |>
  left_join(RecruitCSV, by = c("pond", "season")) |>
  left_join(FallNetted, by = c("pond", "season")) |>
  mutate(across(c(n_tagged_juv, n0, nL, n_csv, n0_csv), ~ replace_na(as.integer(.x), 0L)),
         rec_total = n_tagged_juv + n_csv,
         n0 = n0 + n0_csv,
         netted = rec_total > 0 | coalesce(fall_netted, FALSE)) |>
  select(-fall_netted) |>
  left_join(PondTable |> select(pond, Area), by = "pond") |>
  mutate(rec_density = rec_total / Area,
         log_dens = log(rec_density + 1))

# Standardize over pond-seasons with a fall netting
dens_center <- mean(RecruitIndex$log_dens[RecruitIndex$netted])
dens_scale  <- sd(RecruitIndex$log_dens[RecruitIndex$netted])
RecruitIndex <- RecruitIndex |>
  mutate(dens_std = ifelse(netted, (log_dens - dens_center) / dens_scale, 0))

dens_std <- matrix(0, N_POND, N_SEASONS); n0_mat <- matrix(0L, N_POND, N_SEASONS); nL_mat <- n0_mat
for (r in seq_len(nrow(RecruitIndex))) {
  dens_std[RecruitIndex$pond[r], RecruitIndex$season[r]] <- RecruitIndex$dens_std[r]
  n0_mat[RecruitIndex$pond[r], RecruitIndex$season[r]]   <- RecruitIndex$n0[r]
  nL_mat[RecruitIndex$pond[r], RecruitIndex$season[r]]   <- RecruitIndex$nL[r]
}
# Size-split observations usable by the model: seasons >= 2 with fish measured
SplitObs <- RecruitIndex |> filter(season >= 2, n0 > 0) |> select(pond, season, n0, nL, dens_std)

# ---------------------------------------------------------------------------
# Mean weight by pond x season x stage (for biomass reporting)
# ---------------------------------------------------------------------------
LengthWeight <- StudyBWNFWG |>
  filter(!is.na(weight), !is.na(total_length), weight > 0, total_length > 50) |>
  group_by(species) |>
  summarise(fit = list(lm(log(weight) ~ log(total_length))), n = n(), .groups = "drop") |>
  mutate(a = sapply(fit, function(f) unname(coef(f)[1])),
         b = sapply(fit, function(f) unname(coef(f)[2])),
         sigma = sapply(fit, function(f) summary(f)$sigma)) |>
  select(species, n, a, b, sigma)

WbarObs <- StudyBWNFWG |>
  filter(event == "capture", !is.na(total_length), collection_date >= START_DATE) |>
  inner_join(PondTable |> select(location, pond, pond_species = species), by = "location") |>
  filter(species == pond_species) |>
  inner_join(LengthWeight |> select(species, a, b), by = "species") |>
  mutate(season = season_of_date(collection_date), stage = stage_of(total_length, species),
         w_g = exp(a + b * log(total_length))) |>
  filter(season >= 1, season <= N_SEASONS) |>
  group_by(pond, season, stage) |>
  summarise(Wbar_g = mean(w_g), n_w = n(), .groups = "drop")
WbarPond <- WbarObs |> group_by(pond, stage) |> summarise(Wbar_pond = weighted.mean(Wbar_g, n_w), .groups = "drop")
# Species-level fallback for pond x stage cells with no measured fish at all
# (e.g. IP4 has never handled a stage-S fish)
WbarSp <- WbarObs |>
  left_join(PondTable |> select(pond, sp), by = "pond") |>
  group_by(sp, stage) |> summarise(Wbar_sp = weighted.mean(Wbar_g, n_w), .groups = "drop")
Wbar <- array(NA_real_, dim = c(N_POND, N_SEASONS, N_STAGE))
for (j in 1:N_POND) for (y in 1:N_SEASONS) for (st in 1:N_STAGE) {
  v <- WbarObs$Wbar_g[WbarObs$pond == j & WbarObs$season == y & WbarObs$stage == st]
  if (length(v) == 0) v <- WbarPond$Wbar_pond[WbarPond$pond == j & WbarPond$stage == st]
  if (length(v) == 0) v <- WbarSp$Wbar_sp[WbarSp$sp == PondTable$sp[j] & WbarSp$stage == st]
  Wbar[j, y, st] <- if (length(v) == 0) NA_real_ else v
}

# ---------------------------------------------------------------------------
# Known die-off events (pond x calendar month indicator)
# ---------------------------------------------------------------------------
# Hand-maintained table of independently documented mass-mortality events (not
# inferred from the detection data). Applied as a fixed monthly logit-survival
# offset dE[pond] in rd3_code/rdPond3, shared across S/L/A. Under CALENDAR =
# "m7", prim_floor() collapses May-Sep dates to April, so a summer event only
# partially covered by the Oct-Apr primaries will flag the whole April
# primary; nEvt below is computed from calendar months directly (independent
# of CALENDAR) so the untagged-pool annual survival discount is exact.
DieOffsRaw <- read_csv(DIEOFF_FILE, show_col_types = FALSE) |>
  mutate(StartMonth = ym(StartMonth), EndMonth = ym(EndMonth)) |>
  inner_join(PondTable |> select(location, pond, pond_species = species), by = c("Backwater" = "location")) |>
  filter(Species == pond_species)

Evt  <- matrix(0L, N_POND, T_PRIM)
nEvt <- matrix(0L, N_POND, N_SEASONS)
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
# IP1 FY2022 adult-survival deviation term (BK request, 2026-09-10): a single
# extra logit-survival deviation cell for Pond 1 (IP1) x FY2022, layered on
# top of the existing hierarchical lphiA[j,y] and the pond-level dE[j]
# die-off event offset (see rd3_code / YCB3S_defs.R). IP1FY22[j,y] is 1 only
# at (IP1's pond index, FY2022's season index); 0 everywhere else (including
# in YCB-only builds where IP1 is not modeled) so no other pond-season is
# affected.
# ---------------------------------------------------------------------------
IP1FY22 <- matrix(0L, N_POND, N_SEASONS)
ip1_row <- which(PondTable$Backwater == "IP1")
fy2022_season <- 2022L - FIRST_FY + 1L
if (length(ip1_row) == 1 && fy2022_season >= 1 && fy2022_season <= N_SEASONS) {
  IP1FY22[ip1_row, fy2022_season] <- 1L
}

# ---------------------------------------------------------------------------
# NIMBLE constants / data
# ---------------------------------------------------------------------------
pond_off <- c(0L, cumsum(tabulate(Fish$pond, nbins = N_POND)))
stopifnot(all(diff(Fish$pond) >= 0))

RD_constants <- list(
  NPOND = N_POND, NSP = 2L, NSTAGE = N_STAGE, T = T_PRIM, NS = N_SEASONS, NMOY = N_MOY,
  NFISH = N_FISH, NEV = N_EVENTS, NSPLIT = nrow(SplitObs),
  sp = PondTable$sp, season = Primary$season, moy = Primary$moy, gap = Primary$gap,
  summer = as.integer(Primary$summer),
  tApr = Primary$t[Primary$moy == 7L],   # April is the Atot reference month in both calendars
  pond_off = pond_off,
  entry_t = Fish$entry_t, end_t = Fish$end_t, removed = as.integer(Fish$removed),
  stage0 = Fish$stage0, origin = Fish$origin,
  count_from = Fish$entry_t + as.integer(Fish$origin == 2L),
  eff = eff_std, dens = dens_std, Evt = Evt, nEvt = nEvt, IP1FY22 = IP1FY22,
  ev_pond = Netting$pond, ev_t = Netting$t_e, ev_season = Netting$season,
  ev_stage = Netting$stage, ev_n = Netting$n, ev_remU_before = Netting$remU_before,
  remU = remU_arr,
  sp_pond = SplitObs$pond, sp_season = SplitObs$season, sp_n0 = SplitObs$n0,
  Uk = as.integer(PondTable$pond_all != 1L),   # 1 = pond known empty of untagged fish at start
  Area = PondTable$Area
)

RD_data <- list(
  y = y_mat, K = K_mat, stobs = stobs_mat, post = post_mat,
  ev_m = Netting$m, sp_nL = SplitObs$nL
)

BUILD_INFO_INPUTS <- c("model/YCB3S_data.R", "model/RecruitYOYClassification.R",
                       "data/ReportingData.RData", "data/BWNettingEvents.csv",
                       "data/BWPondAreas.csv", DIEOFF_FILE)
BUILD_INFO <- list(
  built = Sys.time(),
  inputs = BUILD_INFO_INPUTS,
  input_mtime = file.info(BUILD_INFO_INPUTS)$mtime,
  MODEL_PONDS = MODEL_PONDS, StageBreaks = StageBreaks, CALENDAR = CALENDAR, DIEOFF_FILE = DIEOFF_FILE
)

save(Fish, Primary, PondTable, Netting, NettingCSV, Removals, RemovedTagged,
     SameEventRemovals, SpeciesMismatch, LengthWeight, RecruitIndex, SplitObs, Wbar, WbarObs,
     DroppedWeeks, DieOffsRaw, Evt, nEvt, IP1FY22, y_mat, K_mat, stobs_mat, post_mat, eff_wk, eff_hrs, eff_std,
     dens_std, n0_mat, nL_mat, dens_center, dens_scale,
     RD_constants, RD_data, BUILD_INFO,
     START_DATE, FIRST_FY, N_SEASONS, N_MOY, T_PRIM, N_SEC, POST_MONTHS, N_STAGE, StageBreaks,
     MODEL_PONDS, CALENDAR, DET_MONTHS, prim_floor, fy_of_date, season_of_date, stage_of,
     file = OUT_FILE)

# ---------------------------------------------------------------------------
# Reconciliation printout
# ---------------------------------------------------------------------------
cat("\nBuilt", OUT_FILE, "for ponds:", paste(PondTable$Backwater, collapse = ", "),
    "| calendar", CALENDAR, "(", N_MOY, "primaries/season, T =", T_PRIM, ")\n")
cat("\nFish by origin (1 stocked, 2 netting-tagged, 3 established) / stage at entry (1 S, 2 L, 3 A):\n")
print(Fish |> dplyr::count(pond, origin, stage0) |>
        pivot_wider(names_from = stage0, values_from = n, values_fill = 0, names_prefix = "stage"))
cat("\nRemoved tagged fish by type:\n"); print(Fish |> filter(removed) |> dplyr::count(pond, rem_type))
cat("\nSame-event tagged-and-removed fish (counted as untagged harvest):", nrow(SameEventRemovals), "\n")
cat("\nYOY vs. carryover classification of fall-netting untagged first-captures, by pond x season:\n")
print(RecruitYOYClassified |>
        group_by(pond, season) |>
        summarise(n = n(), G = first(G), n_yoy = sum(is_yoy), n_carryover = sum(!is_yoy), .groups = "drop") |>
        mutate(Backwater = PondTable$Backwater[pond], FY = FIRST_FY + season - 1) |>
        select(Backwater, FY, n, G, n_yoy, n_carryover) |> as.data.frame())
cat("\nRecruit index and size split by season (n_tagged_juv/n0/nL now YOY-classified only):\n")
print(RecruitIndex |> select(pond, season, n_tagged_juv, n_csv, rec_total, dens_std, n0, nL) |>
        mutate(pL_obs = round(nL / pmax(n0, 1), 2), dens_std = round(dens_std, 2)) |> as.data.frame())
cat("\nNetting events (stage 1 S, 2 L, 3 A):\n")
print(Netting |> select(pond, EventMonth, stage, t_e, season, m, newtag, harvU, retU_pool, n, remU) |> as.data.frame())
cat("\nDie-off events: ", nrow(DieOffsRaw), "record(s); event pond-months (Evt):", sum(Evt),
    "; event calendar-months by pond x season (nEvt):\n")
print(DieOffsRaw |> select(Backwater, Species, StartMonth, EndMonth, Cause) |> as.data.frame())
if (any(nEvt > 0)) print(which(nEvt > 0, arr.ind = TRUE) |> as.data.frame() |>
      setNames(c("pond", "season")) |> mutate(Backwater = PondTable$Backwater[pond], FY = FIRST_FY + season - 1,
                                              n_months = nEvt[cbind(pond, season)]) |>
      select(Backwater, FY, n_months) |> as.data.frame())

cat("\nDetection summary: fish-months K>0:", sum(K_mat > 0), "; y>0:", sum(y_mat > 0),
    "; stage observations:", N_STAGE_OBS, "; known-alive (code 4):", sum(stobs_mat == 4L), "\n")
cat("stobs table:\n"); print(table(stobs_mat))
cat("\nPer-pond detection summary:\n")
print(tibble(pond = 1:N_POND, Backwater = PondTable$Backwater,
             n_fish = tabulate(Fish$pond, N_POND),
             months_scanned = rowSums(eff_hrs > 0),
             scan_hrs = round(rowSums(eff_hrs)),
             fish_months_K = sapply(1:N_POND, function(j) sum(K_mat[Fish$pond == j, ] > 0)),
             fish_months_y = sapply(1:N_POND, function(j) sum(y_mat[Fish$pond == j, ] > 0)),
             stage_obs = sapply(1:N_POND, function(j) sum(stobs_mat[Fish$pond == j, ] %in% 1:3)),
             known_alive = sapply(1:N_POND, function(j) sum(stobs_mat[Fish$pond == j, ] == 4L)),
             dropped_wks = sapply(1:N_POND, function(j) sum(DroppedWeeks$pond == j))) |> as.data.frame())
if (CALENDAR == "m12") {
  sm <- Primary$summer
  cat("\nSummer (May-Sep) vs Oct-Apr, pond-months scanned:", sum(eff_hrs[, sm] > 0), "/", sum(eff_hrs[, !sm] > 0),
      "; scan hours:", round(sum(eff_hrs[, sm])), "/", round(sum(eff_hrs[, !sm])),
      "; fish-months K>0:", sum(K_mat[, sm] > 0), "/", sum(K_mat[, !sm] > 0),
      "; y>0:", sum(y_mat[, sm] > 0), "/", sum(y_mat[, !sm] > 0), "\n")
}
for (j in 1:N_POND) {
  cat("\nMean weight (g) by season x stage,", PondTable$Backwater[j], ":\n")
  print(round(Wbar[j, , ]))
}
