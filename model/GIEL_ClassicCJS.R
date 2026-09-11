# Classic (frequentist, MLE) annualized Cormack-Jolly-Seber model for Bonytail
# (GIEL) at the three IPCA ponds (IP2, IP5, IP6), fit with the `marked` package.
# A simplified comparison point to the joint hierarchical Bayesian robust-design
# model in model/BWRobustDesign_defs.R (see model/GIEL_ModelMethodology.md and
# model/GIEL_ModelResults_Summary.md for that model).
#
# Design (per user request, 2026-09-11):
#  - Occasions = annual, FY2017-FY2026 (Oct-Apr primary window each year,
#    matching the Bayesian model's primary period). Summer (May-Sep) contacts
#    are ignored entirely -- this is a much simpler model than the robust
#    design and does not use within-season secondary occasions.
#  - Tagged GIEL fish only (stocked or netting-tagged/"wild") at IP2/IP5/IP6;
#    no untagged pool, no netting join, no die-off term, no pond effect.
#  - Two "size classes" per individual, NOT a fixed group: any fish tagged at
#    < 250 mm TL is classified "young" for exactly the one annual interval
#    immediately following its own tagging, and "adult" for every interval
#    thereafter (and for its own tagging interval if tagged >= 250 mm). This
#    is implemented as an explicit id x occasion lookup (see below), not via
#    marked's built-in age-binning machinery, because that machinery's
#    auto-computed `Age` design variable was found to disagree with the
#    supplied `initial.age` for a nontrivial subset of individuals when
#    checked directly against `dp$data$initial.age` (verified via a hand
#    audit here on 2026-09-11) -- the manual approach is fully auditable.
#  - Origin (stocked vs netting-tagged/"wild") is a fixed covariate.
#  - Full 2x2x2 combination model set on Phi (time/constant x size-effect
#    yes/no x origin-effect yes/no) x 2 on p (time/constant) = 16 models,
#    ranked by AIC (AICc correction is negligible here: n ~ 2600, k <= 20).
#  - Any tagged fish with a documented removal (harvest/transfer) would be
#    dropped entirely rather than loss-on-capture-censored (per user
#    decision), but no GIEL fish in this dataset have such a flag (checked
#    below) so no removal-based exclusion was actually needed.
#  - Fish with no recorded total_length at tagging (300 of 2937, all stocked)
#    cannot be assigned a size class and are dropped.
#
# GOF: R2ucare's overall_CJS() (Program RELEASE-style tests: 3.SR transience,
# 3.Sm, 2.CT trap-dependence, 2.CL) run on the pooled (no-covariate) encounter
# histories, per the standard Lebreton et al. (1992) workflow of using these
# component tests to justify which effects belong in the global model.

source("LabFunctions.R")
packages(dplyr); packages(lubridate); packages(tidyr)
packages(marked); packages(R2ucare)

load("data/ReportingData.RData")

fy_of <- function(d) year(d) + as.integer(month(d) > 9)

FY_MIN <- 2017L; FY_MAX <- 2026L
occ_fy <- FY_MIN:FY_MAX
T_OCC  <- length(occ_fy)
occ_cols <- as.character(occ_fy)

giel_ponds <- c("IPCA (Pond 2)", "IPCA (Pond 5)", "IPCA (Pond 6)")

# ---------------------------------------------------------------------------
# Build the fish table and check for documented removals (none found for GIEL)
# ---------------------------------------------------------------------------
GielFish <- StudyBWAnalysis |>
  filter(species == "GIEL", location %in% giel_ponds, event %in% c("stocking", "capture")) |>
  transmute(PITIndex, location, total_length, first_date, event,
            entry_fy = fy_of(first_date),
            stage0 = ifelse(total_length < 250, "J", "A"),
            origin = ifelse(event == "stocking", "stocked", "wild"))

n_removed <- sum(StudyBWAnalysis$species == "GIEL" &
                  StudyBWAnalysis$location %in% giel_ponds &
                  StudyBWAnalysis$Transfer == 1)
if (n_removed > 0) {
  stop("Found ", n_removed, " GIEL fish flagged Transfer==1 -- decide on ",
       "exclusion before proceeding (script assumed none exist).")
}

n_before <- nrow(GielFish)
GielFish <- GielFish |> filter(!is.na(stage0))
cat("Dropped", n_before - nrow(GielFish),
    "fish with no recorded total_length at tagging (can't assign size class).\n")

# ---------------------------------------------------------------------------
# Annual (Oct-Apr) detection indicator per fish
# ---------------------------------------------------------------------------
GielContacts <- StudyBWContacts |>
  filter(Species == "GIEL", Backwater %in% giel_ponds, !is.na(Date)) |>
  mutate(mo = month(Date), fy = fy_of(Date)) |>
  filter(mo %in% c(10, 11, 12, 1, 2, 3, 4), fy >= FY_MIN, fy <= FY_MAX) |>
  distinct(PITIndex, fy)

DetWide <- GielContacts |>
  mutate(det = 1L) |>
  pivot_wider(names_from = fy, values_from = det, values_fill = 0L)
for (fy in occ_fy) if (!as.character(fy) %in% names(DetWide)) DetWide[[as.character(fy)]] <- 0L
DetWide <- DetWide[, c("PITIndex", occ_cols)]

giel_data <- GielFish |>
  left_join(DetWide, by = "PITIndex") |>
  mutate(across(all_of(occ_cols), ~replace_na(.x, 0L)))

# Force entry occasion = 1 (initial capture) and zero out any occasion before entry
for (i in seq_len(nrow(giel_data))) {
  ent <- giel_data$entry_fy[i]
  giel_data[i, occ_cols[occ_fy < ent]] <- 0L
  giel_data[i, as.character(ent)] <- 1L
}

giel_data$ch <- apply(giel_data[, occ_cols], 1, paste, collapse = "")
giel_data$origin <- factor(giel_data$origin, levels = c("wild", "stocked"))
giel_data$initial.age <- ifelse(giel_data$stage0 == "J", 0, 1)

cat("Final analysis dataset:", nrow(giel_data), "tagged GIEL fish,", T_OCC, "occasions (FY",
    FY_MIN, "-FY", FY_MAX, ")\n", sep = "")

# ---------------------------------------------------------------------------
# Process data in marked; build the explicit young/adult design covariate
# ---------------------------------------------------------------------------
giel_cjs <- giel_data |> select(ch, origin, initial.age, location) |> as.data.frame()

dp  <- process.data(giel_cjs, model = "CJS", begin.time = FY_MIN,
                     time.intervals = rep(1, T_OCC - 1))
ddl <- make.design.data(dp)

# marked stores the collapsed ch as comma-separated digits (verified 2026-09-11)
entry_occ <- sapply(strsplit(dp$data$ch, ","), function(x) which(x == "1")[1])
lut <- data.frame(id = dp$data$id, entry_occ = entry_occ, initial.age = dp$data$initial.age)
ddl$Phi <- merge(ddl$Phi, lut, by = "id", all.x = TRUE, sort = FALSE)
ddl$Phi$sizeclass <- factor(
  ifelse(ddl$Phi$initial.age == 0 & ddl$Phi$occ == ddl$Phi$entry_occ, "young", "adult"),
  levels = c("adult", "young")
)
ddl$Phi <- ddl$Phi[order(ddl$Phi$id, ddl$Phi$occ), ]

# ---------------------------------------------------------------------------
# GOF on the pooled (no-covariate) encounter histories
# ---------------------------------------------------------------------------
X <- as.matrix(giel_data[, occ_cols]); mode(X) <- "numeric"
freq <- rep(1, nrow(X))
gof      <- overall_CJS(X, freq)
gof_3sr  <- test3sr(X, freq)$test3sr    # transience
gof_3sm  <- test3sm(X, freq)$test3sm
gof_2ct  <- test2ct(X, freq)$test2ct    # trap-dependence, immediate
gof_2cl  <- test2cl(X, freq)$test2cl    # trap-dependence, delayed

# ---------------------------------------------------------------------------
# Full 16-model set
# ---------------------------------------------------------------------------
Phi.dot         <- list(formula = ~1)
Phi.time        <- list(formula = ~time)
Phi.size        <- list(formula = ~sizeclass)
Phi.time.size   <- list(formula = ~time + sizeclass)
Phi.origin      <- list(formula = ~origin)
Phi.time.origin <- list(formula = ~time + origin)
Phi.size.origin <- list(formula = ~sizeclass + origin)
Phi.full        <- list(formula = ~time + sizeclass + origin)
p.dot  <- list(formula = ~1)
p.time <- list(formula = ~time)

cml <- create.model.list(c("Phi", "p"))
giel_cjs_results <- crm.wrapper(cml, data = dp, ddl = ddl, external = FALSE,
                                 accumulate = FALSE, run = TRUE)
giel_cjs_table <- giel_cjs_results$model.table

# Refit the top model with the Hessian for real-scale estimates + SEs
top_formula <- as.character(giel_cjs_table$model[1])
top_model <- crm(dp, ddl, model.parameters = list(Phi = Phi.full, p = p.time),
                  hessian = TRUE, accumulate = FALSE)
giel_cjs_real <- predict(top_model, ddl = ddl, se = TRUE)

save(giel_cjs_table, giel_cjs_results, top_model, giel_cjs_real,
     gof, gof_3sr, gof_3sm, gof_2ct, gof_2cl,
     giel_data, dp, ddl,
     file = "data/GIEL_ClassicCJS.RData")

cat("\nModel selection table:\n"); print(giel_cjs_table, digits = 1)
cat("\nOverall GOF chi2/df/p:", gof$chi2, "/", gof$degree_of_freedom, "/", gof$p_value, "\n")
cat("Saved data/GIEL_ClassicCJS.RData\n")
