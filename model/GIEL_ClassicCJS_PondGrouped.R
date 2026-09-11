# Pond-grouped extension of the classic (frequentist, MLE) annualized Cormack-Jolly-Seber
# model for Bonytail (GIEL) at the three IPCA ponds (IP2, IP5, IP6), fit with the `marked`
# package. Extends model/GIEL_ClassicCJS.R (which pooled all three ponds) by adding pond
# identity as a grouping factor so Phi and p can be estimated fully independently by pond,
# per user request 2026-09-11. See model/GIEL_ClassicCJS_PondGrouped_Summary.md for results.
#
# Design (per user request):
#  - Start from the SAME data build as GIEL_ClassicCJS.R (annual Oct-Apr occasions,
#    FY2017-FY2026, 300 stocked GIEL with no recorded length dropped, no removals to
#    exclude). Re-checked here: the 300 dropped fish are ALL at IP2, but IP2 retains 1,156
#    tagged fish with a fully populated location x origin x size cross-tab (no structural
#    zero cells), so no pond is excluded on identifiability grounds.
#  - Model set: build pond ("location") interactions on top of the TOP TWO models from the
#    original (pooled) run -- Phi(~time+sizeclass+origin) with p(~time) or p(~1) -- ranging
#    from a fully independent-by-pond Phi and p (location*time + location*sizeclass +
#    location*origin on Phi; location*time on p) down to an additive pond effect only, plus
#    the two original pooled models as a baseline. 12 models total; a second exploratory
#    batch of 4 more (dropping individual location interactions) was added after the first
#    AIC pass showed 5 models within ~5 AIC of each other.
#  - Parametric bootstrap GOF (300 replicates) on the AIC-TOP pond-grouped model itself
#    (not just the null/pooled model): simulate encounter histories under the fitted
#    model's real Phi/p estimates for every actual individual's covariate path, refit the
#    identical model to each simulated dataset (warm-started from the observed fit), and
#    compute a deviance relative to a cohort-saturated reference (cohort = location x
#    origin x initial.age x entry occasion -- exactly the covariate information available
#    to differentiate individuals at release). c-hat = observed deviance / mean(bootstrap
#    deviance).

source("LabFunctions.R")
packages(dplyr); packages(marked); packages(ggplot2)

load("data/ReportingData.RData")

fy_of <- function(d) year(d) + as.integer(month(d) > 9)
FY_MIN <- 2017L; FY_MAX <- 2026L
occ_fy <- FY_MIN:FY_MAX
T_OCC  <- length(occ_fy)
occ_cols <- as.character(occ_fy)
giel_ponds <- c("IPCA (Pond 2)", "IPCA (Pond 5)", "IPCA (Pond 6)")

# ---------------------------------------------------------------------------
# Rebuild the fish table exactly as in GIEL_ClassicCJS.R (kept in sync manually)
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
if (n_removed > 0) stop("Found ", n_removed, " GIEL fish flagged Transfer==1 -- decide on exclusion.")

n_before <- nrow(GielFish)
GielFish <- GielFish |> filter(!is.na(stage0))
cat("Dropped", n_before - nrow(GielFish), "fish with no total_length at tagging.\n")

# --- Identifiability check: would any pond need exclusion? ---
missing_tl_pond <- table(
  (StudyBWAnalysis |> filter(species == "GIEL", location %in% giel_ponds,
                              event %in% c("stocking", "capture")))$location[
    is.na((StudyBWAnalysis |> filter(species == "GIEL", location %in% giel_ponds,
                                      event %in% c("stocking", "capture")))$total_length)])
cat("Dropped (missing length) fish by pond:\n"); print(missing_tl_pond)
cross_tab <- GielFish |> count(location, origin, stage0)
cat("\nlocation x origin x stage0 cross-tab (checking for empty cells):\n")
print(cross_tab)
if (any(xtabs(n ~ location + origin + stage0, data = cross_tab) == 0)) {
  warning("A location x origin x stage0 cell is empty -- consider excluding that pond.")
} else {
  cat("No empty cells: all three ponds retain a full origin x size cross-tab, so none are",
      "excluded on structural-identifiability grounds (the 300 dropped fish are all IP2",
      "stocked adults, but IP2 still has 291 measured stocked-adult and 8 stocked-young",
      "fish plus 857 netting-tagged fish).\n")
}

# ---------------------------------------------------------------------------
# Annual (Oct-Apr) detection indicator per fish (identical logic to GIEL_ClassicCJS.R)
# ---------------------------------------------------------------------------
GielContacts <- StudyBWContacts |>
  filter(Species == "GIEL", Backwater %in% giel_ponds, !is.na(Date)) |>
  mutate(mo = month(Date), fy = fy_of(Date)) |>
  filter(mo %in% c(10, 11, 12, 1, 2, 3, 4), fy >= FY_MIN, fy <= FY_MAX) |>
  distinct(PITIndex, fy)

DetWide <- GielContacts |>
  mutate(det = 1L) |>
  tidyr::pivot_wider(names_from = fy, values_from = det, values_fill = 0L)
for (fy in occ_fy) if (!as.character(fy) %in% names(DetWide)) DetWide[[as.character(fy)]] <- 0L
DetWide <- DetWide[, c("PITIndex", occ_cols)]

giel_data <- GielFish |>
  left_join(DetWide, by = "PITIndex") |>
  mutate(across(all_of(occ_cols), ~tidyr::replace_na(.x, 0L)))

for (i in seq_len(nrow(giel_data))) {
  ent <- giel_data$entry_fy[i]
  giel_data[i, occ_cols[occ_fy < ent]] <- 0L
  giel_data[i, as.character(ent)] <- 1L
}
giel_data$ch <- apply(giel_data[, occ_cols], 1, paste, collapse = "")
giel_data$origin <- factor(giel_data$origin, levels = c("wild", "stocked"))
giel_data$initial.age <- ifelse(giel_data$stage0 == "J", 0, 1)

cat("\nFinal analysis dataset:", nrow(giel_data), "tagged GIEL fish,", T_OCC, "occasions\n")
cat("By pond:\n"); print(table(giel_data$location))

# ---------------------------------------------------------------------------
# Process data with pond as a group; rebuild the sizeclass design covariate
# ---------------------------------------------------------------------------
giel_cjs2 <- giel_data |> select(ch, origin, initial.age, location) |> as.data.frame()
giel_cjs2$location <- factor(giel_cjs2$location,
                              levels = giel_ponds, labels = c("IP2", "IP5", "IP6"))

dp2  <- process.data(giel_cjs2, model = "CJS", groups = "location",
                      begin.time = FY_MIN, time.intervals = rep(1, T_OCC - 1))
ddl2 <- make.design.data(dp2)

entry_occ2 <- sapply(strsplit(dp2$data$ch, ","), function(x) which(x == "1")[1])
lut2 <- data.frame(id = dp2$data$id, entry_occ = entry_occ2, initial.age = dp2$data$initial.age)
ddl2$Phi <- merge(ddl2$Phi, lut2, by = "id", all.x = TRUE, sort = FALSE)
ddl2$Phi$sizeclass <- factor(
  ifelse(ddl2$Phi$initial.age == 0 & ddl2$Phi$occ == ddl2$Phi$entry_occ, "young", "adult"),
  levels = c("adult", "young")
)
ddl2$Phi <- ddl2$Phi[order(ddl2$Phi$id, ddl2$Phi$occ), ]

# ---------------------------------------------------------------------------
# Model set: top-2 baseline (pooled) + pond-grouped variants
# ---------------------------------------------------------------------------
Phi.full         <- list(formula = ~time + sizeclass + origin)                          # previous top-2 Phi
Phi.full.locFull <- list(formula = ~location*time + location*sizeclass + location*origin)  # fully independent
Phi.full.locTime <- list(formula = ~location*time + sizeclass + origin)                 # pond-specific time only
Phi.full.locAdd  <- list(formula = ~location + time + sizeclass + origin)               # additive pond effect
Phi.locTime.size <- list(formula = ~location*time + location*sizeclass + origin)        # drop location:origin
Phi.locTime.orig <- list(formula = ~location*time + sizeclass + location*origin)        # drop location:sizeclass

p.time         <- list(formula = ~time)
p.dot          <- list(formula = ~1)
p.time.loc     <- list(formula = ~location*time)
p.loc          <- list(formula = ~location)
p.loc.time.add <- list(formula = ~location + time)

model_specs <- list(
  M1_base_ptime        = list(Phi = Phi.full,          p = p.time),         # previous overall top model
  M2_base_pdot         = list(Phi = Phi.full,          p = p.dot),          # previous 2nd-place model
  M3_fullloc_ptl       = list(Phi = Phi.full.locFull,  p = p.time.loc),     # fully independent Phi & p
  M4_fullloc_ploc      = list(Phi = Phi.full.locFull,  p = p.loc),
  M5_loctime_ptl       = list(Phi = Phi.full.locTime,  p = p.time.loc),
  M6_locadd_ploc       = list(Phi = Phi.full.locAdd,   p = p.loc),
  M7_base_ptl          = list(Phi = Phi.full,          p = p.time.loc),
  M8_fullloc_ptime     = list(Phi = Phi.full.locFull,  p = p.time),
  M9_loctime_ptime     = list(Phi = Phi.full.locTime,  p = p.time),
  M10_fullloc_plocadd  = list(Phi = Phi.full.locFull,  p = p.loc.time.add),
  M11_dropOrig_ptime   = list(Phi = Phi.locTime.size,  p = p.time),
  M12_dropSize_ptime   = list(Phi = Phi.locTime.orig,  p = p.time)
)

fits <- list()
for (nm in names(model_specs)) {
  cat("Fitting", nm, "...\n")
  fits[[nm]] <- crm(dp2, ddl2, model.parameters = model_specs[[nm]], hessian = TRUE, accumulate = FALSE)
}

get_npar <- function(f) attr(f$results$optim.details, "npar")
model_formulas <- data.frame(
  model = names(model_specs),
  Phi = sapply(model_specs, function(s) deparse(s$Phi$formula)),
  p   = sapply(model_specs, function(s) deparse(s$p$formula))
)
aic_final <- data.frame(
  model   = names(fits),
  npar    = sapply(fits, get_npar),
  neg2lnl = sapply(fits, function(f) f$results$neg2lnl),
  AIC     = sapply(fits, function(f) f$results$AIC)
)
aic_final <- merge(aic_final, model_formulas, by = "model")
aic_final <- aic_final[order(aic_final$AIC), ]
aic_final$DeltaAIC <- aic_final$AIC - min(aic_final$AIC)
aic_final$weight   <- exp(-0.5 * aic_final$DeltaAIC)
aic_final$weight   <- aic_final$weight / sum(aic_final$weight)
rownames(aic_final) <- NULL
cat("\nFinal AIC table:\n"); print(aic_final[, c("model","Phi","p","npar","AIC","DeltaAIC","weight")], digits = 3)

top_name  <- aic_final$model[1]
top_ponded <- fits[[top_name]]
cat("\nTop model:", top_name, "\n")

preds <- predict(top_ponded, ddl = ddl2, se = FALSE)

# ---------------------------------------------------------------------------
# Parametric bootstrap GOF on the top pond-grouped model
# ---------------------------------------------------------------------------
sim_base <- giel_data |>
  transmute(location = factor(location, levels = giel_ponds, labels = c("IP2","IP5","IP6")),
            origin, initial.age, entry_occ = entry_fy - FY_MIN + 1L)
sim_base$cohort <- paste(sim_base$location, sim_base$origin, sim_base$initial.age, sim_base$entry_occ)

phi_lut <- preds$Phi |> mutate(key = paste(location, time, sizeclass, origin)) |>
  select(key, estimate) |> distinct()
phi_map <- setNames(phi_lut$estimate, phi_lut$key)
p_lut   <- preds$p |> select(time, estimate) |> distinct()
p_map   <- setNames(p_lut$estimate, as.character(p_lut$time))

simulate_ch <- function(sim_base, phi_map, p_map, T_OCC, FY_MIN) {
  n <- nrow(sim_base)
  y <- matrix(0L, n, T_OCC)
  for (i in seq_len(n)) y[i, sim_base$entry_occ[i]] <- 1L
  alive <- rep(TRUE, n)
  for (t in 1:(T_OCC - 1)) {
    elig <- which(sim_base$entry_occ <= t & alive)
    if (!length(elig)) next
    sizeclass <- ifelse(sim_base$entry_occ[elig] == t & sim_base$initial.age[elig] == 0, "young", "adult")
    key <- paste(sim_base$location[elig], FY_MIN + t - 1, sizeclass, sim_base$origin[elig])
    phi <- phi_map[key]
    if (any(is.na(phi))) stop("Missing phi lookup for: ", paste(unique(key[is.na(phi)]), collapse = ", "))
    survive <- rbinom(length(elig), 1, phi) == 1L
    alive[elig[!survive]] <- FALSE
    survivors <- elig[survive]
    if (length(survivors)) {
      pdet <- p_map[[as.character(FY_MIN + t)]]
      y[survivors, t + 1] <- rbinom(length(survivors), 1, pdet)
    }
  }
  apply(y, 1, paste, collapse = "")
}

sat_loglik <- function(ch, cohort) {
  sum(tapply(ch, cohort, function(x) { tab <- table(x); n <- sum(tab); sum(tab * log(tab / n)) }))
}

obs_ch       <- giel_data$ch
sat_ll_obs   <- sat_loglik(obs_ch, sim_base$cohort)
neg2lnl_obs  <- top_ponded$results$neg2lnl
deviance_obs <- neg2lnl_obs - (-2 * sat_ll_obs)

B <- 300
set.seed(7461)
seeds <- sample.int(1e6, B)
boot_results <- vector("list", B)
top_formula <- model_specs[[top_name]]

t0 <- Sys.time()
for (b in seq_len(B)) {
  set.seed(seeds[b])
  ch_sim <- simulate_ch(sim_base, phi_map, p_map, T_OCC, FY_MIN)
  df <- data.frame(ch = ch_sim, origin = sim_base$origin,
                    initial.age = sim_base$initial.age, location = sim_base$location)
  dp_b <- process.data(df, model = "CJS", groups = "location",
                        begin.time = FY_MIN, time.intervals = rep(1, T_OCC - 1))
  ddl_b <- make.design.data(dp_b)
  entry_occ_b <- sapply(strsplit(dp_b$data$ch, ","), function(x) which(x == "1")[1])
  lut_b <- data.frame(id = dp_b$data$id, entry_occ = entry_occ_b, initial.age = dp_b$data$initial.age)
  ddl_b$Phi <- merge(ddl_b$Phi, lut_b, by = "id", all.x = TRUE, sort = FALSE)
  ddl_b$Phi$sizeclass <- factor(
    ifelse(ddl_b$Phi$initial.age == 0 & ddl_b$Phi$occ == ddl_b$Phi$entry_occ, "young", "adult"),
    levels = c("adult", "young"))
  ddl_b$Phi <- ddl_b$Phi[order(ddl_b$Phi$id, ddl_b$Phi$occ), ]

  fit_b <- tryCatch(
    crm(dp_b, ddl_b, model.parameters = top_formula, initial = top_ponded,
        hessian = FALSE, accumulate = TRUE),
    error = function(e) NULL)
  if (is.null(fit_b)) {
    boot_results[[b]] <- list(neg2lnl = NA, deviance = NA, convergence = NA)
  } else {
    sat_ll_b <- sat_loglik(ch_sim, sim_base$cohort)
    boot_results[[b]] <- list(neg2lnl = fit_b$results$neg2lnl,
                               deviance = fit_b$results$neg2lnl - (-2 * sat_ll_b),
                               convergence = fit_b$results$convergence)
  }
  if (b %% 50 == 0) cat("bootstrap", b, "of", B, "-", round(as.numeric(Sys.time()-t0, units="mins"),1), "min\n")
}

dev_boot <- sapply(boot_results, function(x) x$deviance)
dev_boot_ok  <- dev_boot[!is.na(dev_boot)]
c_hat_mean   <- deviance_obs / mean(dev_boot_ok)
c_hat_median <- deviance_obs / median(dev_boot_ok)
boot_p_value <- mean(dev_boot_ok >= deviance_obs)

cat("\nObserved deviance:", round(deviance_obs, 2), "\n")
cat("Bootstrap deviance mean/median/sd:", round(mean(dev_boot_ok), 2), "/",
    round(median(dev_boot_ok), 2), "/", round(sd(dev_boot_ok), 2), "\n")
cat("c-hat (mean-based):", round(c_hat_mean, 3), "\n")
cat("c-hat (median-based):", round(c_hat_median, 3), "\n")
cat("Bootstrap p-value:", round(boot_p_value, 3), "\n")

p_gof <- ggplot(data.frame(deviance = dev_boot_ok), aes(x = deviance)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = deviance_obs, color = "red", linewidth = 1) +
  labs(title = "Parametric bootstrap GOF: top pond-grouped GIEL CJS model",
       subtitle = paste0("c-hat = ", round(c_hat_mean, 2), " (mean-based), bootstrap p = ", round(boot_p_value, 3)),
       x = "Deviance (relative to cohort-saturated model)", y = paste0("Bootstrap replicates (n=", B, ")"))
ggsave("model/GIEL_ClassicCJS_PondGrouped_bootGOF.png", p_gof, width = 7, height = 5, dpi = 150)

save(aic_final, fits, top_name, top_ponded, preds,
     deviance_obs, dev_boot, c_hat_mean, c_hat_median, boot_p_value,
     giel_data, dp2, ddl2, sim_base,
     file = "data/GIEL_ClassicCJS_PondGrouped.RData")

cat("\nSaved data/GIEL_ClassicCJS_PondGrouped.RData and model/GIEL_ClassicCJS_PondGrouped_bootGOF.png\n")
