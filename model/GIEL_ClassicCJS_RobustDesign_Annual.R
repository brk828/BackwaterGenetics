# Classic (frequentist, MLE) ANNUAL CJS model for Bonytail (GIEL) at IPCA
# Pond 5 (IP5) -- stocked fish only. Per user request 2026-09-11:
#  - `marked` may not be built to handle the monthly within-season "fix=1"
#    robust-design structure used in GIEL_ClassicCJS_RobustDesign.R --
#    crm() reproducibly hard-crashed (non-catchable, even under process
#    isolation via callr) on a large fraction of bootstrap-simulated
#    capture histories under that design. Root cause traced to the
#    within-season closure assumption (60 of 69 monthly transitions fixed
#    to Phi=1, survival information concentrated into just 9 annual
#    April->October gaps), likely producing degenerate/boundary outcomes
#    in the small (~18-20 fish) "young" stratum under simulation.
#  - This version collapses the 70 monthly occasions down to 10 ANNUAL
#    occasions (FY2017-FY2026), i.e. exactly the 9 between-year intervals
#    and none of the within-season ones -- there is no "fix" mechanism at
#    all now, every transition is a real, freely-estimated annual survival
#    parameter. All other model structure is unchanged from
#    GIEL_ClassicCJS_RobustDesign.R: same pond (IP5), same stocked-only
#    fish set, same "young-for-first-interval-then-adult" sizeclass logic,
#    same 3 Phi structures (const/gap/dieoff) x 2 p structures (now
#    p ~ time [annual] and p ~ 1, replacing the monthly p ~ moy since
#    there is only one detection opportunity per year to model).
#  - Detection: a fish is "detected" at an annual occasion if it had >=1
#    PIT contact anywhere in that FY's Oct-Apr window (collapsed from the
#    previous 0-4-week-detections-per-month scheme).
#  - Entry/suppression logic unchanged in spirit: a fish is released at
#    the annual occasion matching its stocking FY; any contact recorded in
#    that same release-FY window is superseded by the release marker (this
#    matters less now since there's no "later month within the release
#    season" to suppress -- occasions are annual).
#  - Same parametric-bootstrap c-hat recipe as the monthly script, but each
#    replicate's refit is isolated in its own callr subprocess (protects
#    the whole job even if crm() ever crashes again) and the simulation/
#    design-building code has no fix-column complexity to introduce a
#    crash vector in the first place.
#  - Compare the resulting c-hat, survival (esp. FY2020 die-off gap), and
#    detection against both the (incomplete) monthly attempt and the
#    existing GIEL-only Bayesian robust-design fit
#    (data/BWRobustDesign_MCMC_giel.RData).

source("LabFunctions.R")
packages(dplyr); packages(lubridate); packages(tidyr); packages(marked); packages(ggplot2); packages(callr)

load("data/ReportingData.RData")

fy_of <- function(d) year(d) + as.integer(month(d) > 9)
FY_MIN <- 2017L; FY_MAX <- 2026L
occ_fy <- FY_MIN:FY_MAX
n_fy <- length(occ_fy)
T_OCC <- n_fy  # 10 annual occasions -- 9 between-year intervals, no within-season transitions
season_months <- c(10, 11, 12, 1, 2, 3, 4)
fy_to_occ <- function(fy) fy - FY_MIN + 1L

# ---------------------------------------------------------------------------
# 1. Pond selection check (printed for the record; IP5 is used below, same
#    choice/rationale as GIEL_ClassicCJS_RobustDesign.R)
# ---------------------------------------------------------------------------
giel_ponds <- c("IPCA (Pond 2)", "IPCA (Pond 5)", "IPCA (Pond 6)")
pond_check <- StudyBWAnalysis |>
  filter(species == "GIEL", location %in% giel_ponds, event == "stocking") |>
  group_by(location) |>
  summarise(n_stocked = n(), n_with_TL = sum(!is.na(total_length)),
            n_cohorts = n_distinct(release_year), .groups = "drop")
cat("Stocked-fish TL completeness by pond (GIEL):\n"); print(pond_check)
cat("-> IP5 selected: most stocked fish with complete TL and 2 release cohorts.\n\n")

# ---------------------------------------------------------------------------
# 2. Fish table: IP5, stocked only, TL required
# ---------------------------------------------------------------------------
GielIP5 <- StudyBWAnalysis |>
  filter(species == "GIEL", location == "IPCA (Pond 5)", event == "stocking",
         !is.na(total_length)) |>
  transmute(PITIndex, total_length, first_date,
            entry_fy = fy_of(first_date),
            stage0 = ifelse(total_length < 250, "young", "adult"),
            entry_occ = fy_to_occ(entry_fy))

n_transfer <- sum(StudyBWAnalysis$species == "GIEL" &
                     StudyBWAnalysis$location == "IPCA (Pond 5)" &
                     StudyBWAnalysis$event == "stocking" & StudyBWAnalysis$Transfer == 1)
if (n_transfer > 0) stop("Found ", n_transfer, " IP5 stocked GIEL flagged Transfer==1 -- decide handling.")

cat("IP5 stocked fish:", nrow(GielIP5), "\n")
cat("By cohort/stage0:\n"); print(table(GielIP5$entry_fy, GielIP5$stage0))

# ---------------------------------------------------------------------------
# 3. Annual detection matrix and entry/ignore logic
# ---------------------------------------------------------------------------
GielIP5Contacts <- StudyBWContacts |>
  filter(Species == "GIEL", Backwater == "IPCA (Pond 5)", !is.na(Date),
         PITIndex %in% GielIP5$PITIndex) |>
  mutate(mo = month(Date), fy = fy_of(Date)) |>
  filter(mo %in% season_months, fy %in% occ_fy) |>
  distinct(PITIndex, fy) |>
  mutate(occ = fy_to_occ(fy)) |>
  distinct(PITIndex, occ)

DetWide <- GielIP5Contacts |> mutate(det = 1L) |>
  pivot_wider(names_from = occ, values_from = det, values_fill = 0L)
missing_occ <- setdiff(as.character(1:T_OCC), names(DetWide))
for (oc in missing_occ) DetWide[[oc]] <- 0L
DetWide <- DetWide[, c("PITIndex", as.character(1:T_OCC))]

giel_ip5 <- GielIP5 |> left_join(DetWide, by = "PITIndex") |>
  mutate(across(as.character(1:T_OCC), ~replace_na(.x, 0L)))

occ_cols <- as.character(1:T_OCC)
M <- as.matrix(giel_ip5[, occ_cols]); mode(M) <- "integer"
n_suppressed <- 0L
for (i in seq_len(nrow(giel_ip5))) {
  ent <- giel_ip5$entry_occ[i]
  if (ent > 1) {
    n_suppressed <- n_suppressed + sum(M[i, 1:(ent - 1)])
    M[i, 1:(ent - 1)] <- 0L
  }
  M[i, ent] <- 1L
}
giel_ip5$ch <- apply(M, 1, paste, collapse = "")
giel_ip5$stage0 <- factor(giel_ip5$stage0, levels = c("adult", "young"))
cat("\nSuppressed (ignored, pre-entry) contact-occasions:", n_suppressed, "\n")

# ---------------------------------------------------------------------------
# 4. Process data; build Phi design covariates (sizeclass, gap, dieoff) --
#    NO "fix" column: every one of the T_OCC - 1 = 9 transitions is a real,
#    freely-estimated annual survival interval.
# ---------------------------------------------------------------------------
giel_cjs_ip5 <- giel_ip5 |> select(ch, stage0) |> as.data.frame()
dp_ip5 <- process.data(giel_cjs_ip5, model = "CJS", begin.time = 1,
                        time.intervals = rep(1, T_OCC - 1))
ddl_ip5 <- make.design.data(dp_ip5)

entry_occ_lut <- sapply(strsplit(dp_ip5$data$ch, ","), function(x) which(x == "1")[1])
idmap <- data.frame(id = as.character(dp_ip5$data$id), entry_occ = entry_occ_lut)

ddl_ip5$Phi <- merge(ddl_ip5$Phi, idmap, by = "id", all.x = TRUE, sort = FALSE)
ddl_ip5$Phi <- ddl_ip5$Phi[order(ddl_ip5$Phi$id, ddl_ip5$Phi$occ), ]
occ_num <- as.numeric(as.character(ddl_ip5$Phi$occ))

ddl_ip5$Phi$sizeclass <- factor(
  ifelse(ddl_ip5$Phi$stage0 == "young" & occ_num == ddl_ip5$Phi$entry_occ, "young", "adult"),
  levels = c("adult", "young"))

gap_fy <- FY_MIN + occ_num - 1
ddl_ip5$Phi$gap <- factor(as.character(gap_fy), levels = as.character(FY_MIN:(FY_MAX - 1)))
ddl_ip5$Phi$dieoff <- factor(ifelse(gap_fy == 2020, "event", "normal"),
                               levels = c("normal", "event"))

cat("\nPhi design rows:", nrow(ddl_ip5$Phi), "(all free -- no fix mechanism)\n")
cat("sizeclass table:\n"); print(table(ddl_ip5$Phi$sizeclass))
cat("gap table:\n"); print(table(ddl_ip5$Phi$gap))
cat("dieoff table:\n"); print(table(ddl_ip5$Phi$dieoff))
cat("\np design rows:", nrow(ddl_ip5$p), "\n")

# ---------------------------------------------------------------------------
# 5. Fit the 6-model set (3 Phi x 2 p; p ~ time replaces p ~ moy)
# ---------------------------------------------------------------------------
Phi.const  <- list(formula = ~ sizeclass)
Phi.gap    <- list(formula = ~ sizeclass + gap)
Phi.dieoff <- list(formula = ~ sizeclass + dieoff)
p.time     <- list(formula = ~ time)
p.dot      <- list(formula = ~ 1)

model_specs_ip5 <- list(
  const_ptime  = list(Phi = Phi.const,  p = p.time),
  const_pdot   = list(Phi = Phi.const,  p = p.dot),
  gap_ptime    = list(Phi = Phi.gap,    p = p.time),
  gap_pdot     = list(Phi = Phi.gap,    p = p.dot),
  dieoff_ptime = list(Phi = Phi.dieoff, p = p.time),
  dieoff_pdot  = list(Phi = Phi.dieoff, p = p.dot)
)

fits_ip5 <- list()
for (nm in names(model_specs_ip5)) {
  cat("\nFitting", nm, "...\n")
  fits_ip5[[nm]] <- crm(dp_ip5, ddl_ip5, model.parameters = model_specs_ip5[[nm]],
                          hessian = TRUE, accumulate = FALSE)
  cat("  AIC:", fits_ip5[[nm]]$results$AIC, "\n")
}

get_npar <- function(f) attr(f$results$optim.details, "npar")
aic_ip5 <- data.frame(
  model = names(fits_ip5),
  npar = sapply(fits_ip5, get_npar),
  neg2lnl = sapply(fits_ip5, function(f) f$results$neg2lnl),
  AIC = sapply(fits_ip5, function(f) f$results$AIC)
)
aic_ip5 <- aic_ip5[order(aic_ip5$AIC), ]
aic_ip5$DeltaAIC <- aic_ip5$AIC - min(aic_ip5$AIC)
aic_ip5$weight <- exp(-0.5 * aic_ip5$DeltaAIC); aic_ip5$weight <- aic_ip5$weight / sum(aic_ip5$weight)
rownames(aic_ip5) <- NULL
cat("\n\nFinal AIC table:\n"); print(aic_ip5, digits = 3)

top_name_ip5 <- aic_ip5$model[1]
top_ip5 <- fits_ip5[[top_name_ip5]]
cat("\nTop model:", top_name_ip5, "\n")
preds_ip5 <- predict(top_ip5, ddl = ddl_ip5, se = TRUE)

save(pond_check, giel_ip5, dp_ip5, ddl_ip5, model_specs_ip5, fits_ip5,
     aic_ip5, top_name_ip5, top_ip5, preds_ip5, n_suppressed,
     file = "data/GIEL_ClassicCJS_RobustDesign_Annual.RData")
cat("\nSaved intermediate results to data/GIEL_ClassicCJS_RobustDesign_Annual.RData\n")

# ---------------------------------------------------------------------------
# 6. Parametric bootstrap c-hat on the top model. Each replicate's crm()
#    refit is isolated in its own callr subprocess: if it ever crashes (as
#    happened routinely under the monthly design), that surfaces as a
#    normal catchable error here instead of killing the whole job.
# ---------------------------------------------------------------------------
sim_base <- giel_ip5 |>
  transmute(entry_fy, stage0, entry_occ, cohort = paste(entry_fy, stage0))

build_get_phi <- function(preds_Phi) {
  phi_covs <- setdiff(names(preds_Phi), c("occ", "estimate", "se", "lcl", "ucl", "sizeclass"))
  function(occ, sizeclass_i) {
    row_sc <- preds_Phi$sizeclass == sizeclass_i
    if (length(phi_covs) == 0) {
      val <- preds_Phi$estimate[row_sc]
      return(val[1])
    }
    g <- as.character(FY_MIN + occ - 1)
    if ("gap" %in% phi_covs) {
      key_match <- row_sc & as.character(preds_Phi$gap) == g
      val <- preds_Phi$estimate[key_match]
      if (length(val) == 0) val <- preds_Phi$estimate[preds_Phi$sizeclass == "adult" & as.character(preds_Phi$gap) == g]
      return(val[1])
    }
    if ("dieoff" %in% phi_covs) {
      dstat <- ifelse(g == "2020", "event", "normal")
      key_match <- row_sc & as.character(preds_Phi$dieoff) == dstat
      val <- preds_Phi$estimate[key_match]
      if (length(val) == 0) val <- preds_Phi$estimate[preds_Phi$sizeclass == "adult" & as.character(preds_Phi$dieoff) == dstat]
      return(val[1])
    }
    NA
  }
}
get_phi <- build_get_phi(preds_ip5$Phi)

build_get_p <- function(preds_p) {
  p_covs <- setdiff(names(preds_p), c("occ", "estimate", "se", "lcl", "ucl"))
  if ("time" %in% p_covs) {
    p_by_time <- setNames(preds_p$estimate, as.character(preds_p$time))
    function(occ) {
      val <- p_by_time[[as.character(occ)]]
      if (is.null(val) || is.na(val)) 0 else val
    }
  } else {
    val0 <- preds_p$estimate[1]
    function(occ) if (is.na(val0)) 0 else val0
  }
}
get_p <- build_get_p(preds_ip5$p)

simulate_ch <- function(sim_base, T_OCC) {
  n <- nrow(sim_base)
  y <- matrix(0L, n, T_OCC)
  for (i in seq_len(n)) y[i, sim_base$entry_occ[i]] <- 1L
  alive <- rep(TRUE, n)
  for (t in seq_len(T_OCC - 1)) {
    elig <- which(sim_base$entry_occ <= t & alive)
    if (!length(elig)) next
    sizeclass_i <- ifelse(sim_base$stage0[elig] == "young" & sim_base$entry_occ[elig] == t,
                           "young", "adult")
    phi <- sapply(seq_along(elig), function(k) get_phi(t, sizeclass_i[k]))
    survive <- rbinom(length(elig), 1, phi) == 1L
    alive[elig[!survive]] <- FALSE
    survivors <- elig[survive]
    if (length(survivors)) {
      pdet <- get_p(t + 1)
      y[survivors, t + 1] <- rbinom(length(survivors), 1, pdet)
    }
  }
  apply(y, 1, paste, collapse = ",")
}

sat_loglik <- function(ch, cohort) {
  sum(tapply(ch, cohort, function(x) { tab <- table(x); n <- sum(tab); sum(tab * log(tab / n)) }))
}

obs_ch <- sapply(strsplit(giel_ip5$ch, ""), paste, collapse = ",")
sat_ll_obs   <- sat_loglik(obs_ch, sim_base$cohort)
neg2lnl_obs  <- top_ip5$results$neg2lnl
deviance_obs <- neg2lnl_obs - (-2 * sat_ll_obs)

# --- subprocess-isolated refit function -------------------------------------
fit_one_boot <- function(ch_sim, stage0, top_formula_ip5, top_ip5, T_OCC, FY_MIN, FY_MAX) {
  library(marked)
  df <- data.frame(ch = ch_sim, stage0 = stage0)
  dp_b <- process.data(df, model = "CJS", begin.time = 1, time.intervals = rep(1, T_OCC - 1))
  ddl_b <- make.design.data(dp_b)
  entry_occ_b <- sapply(strsplit(dp_b$data$ch, ","), function(x) which(x == "1")[1])
  idmap_b <- data.frame(id = as.character(dp_b$data$id), entry_occ = entry_occ_b)
  ddl_b$Phi <- merge(ddl_b$Phi, idmap_b, by = "id", all.x = TRUE, sort = FALSE)
  ddl_b$Phi <- ddl_b$Phi[order(ddl_b$Phi$id, ddl_b$Phi$occ), ]
  occ_num_b <- as.numeric(as.character(ddl_b$Phi$occ))
  ddl_b$Phi$sizeclass <- factor(
    ifelse(ddl_b$Phi$stage0 == "young" & occ_num_b == ddl_b$Phi$entry_occ, "young", "adult"),
    levels = c("adult", "young"))
  gap_fy_b <- FY_MIN + occ_num_b - 1
  ddl_b$Phi$gap <- factor(as.character(gap_fy_b), levels = as.character(FY_MIN:(FY_MAX - 1)))
  ddl_b$Phi$dieoff <- factor(ifelse(gap_fy_b == 2020, "event", "normal"),
                               levels = c("normal", "event"))
  fit_b <- crm(dp_b, ddl_b, model.parameters = top_formula_ip5, initial = top_ip5,
               hessian = FALSE, accumulate = TRUE)
  list(neg2lnl = fit_b$results$neg2lnl, convergence = fit_b$results$convergence)
}

top_formula_ip5 <- model_specs_ip5[[top_name_ip5]]

B <- 300
set.seed(2854)
seeds <- sample.int(1e6, B)
boot_results <- vector("list", B)
CHECKPOINT_FILE <- "data/GIEL_ClassicCJS_RobustDesign_Annual_checkpoint.RData"

t0 <- Sys.time()
for (b in seq_len(B)) {
  set.seed(seeds[b])
  ch_sim <- simulate_ch(sim_base, T_OCC)
  sat_ll_b <- sat_loglik(ch_sim, sim_base$cohort)

  res <- tryCatch(
    callr::r(fit_one_boot,
             args = list(ch_sim = ch_sim, stage0 = giel_ip5$stage0,
                          top_formula_ip5 = top_formula_ip5, top_ip5 = top_ip5,
                          T_OCC = T_OCC, FY_MIN = FY_MIN, FY_MAX = FY_MAX),
             package = TRUE, timeout = 120),
    error = function(e) { cat("  b=", b, "subprocess FAILED:", conditionMessage(e), "\n"); NULL })

  if (is.null(res)) {
    boot_results[[b]] <- list(deviance = NA, convergence = NA)
  } else {
    boot_results[[b]] <- list(deviance = res$neg2lnl - (-2 * sat_ll_b),
                               convergence = res$convergence)
  }
  save(boot_results, b, file = CHECKPOINT_FILE)
  if (b %% 25 == 0) cat("bootstrap", b, "of", B, "-",
                          round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
}

dev_boot <- sapply(boot_results, function(x) x$deviance)
dev_boot_ok  <- dev_boot[!is.na(dev_boot)]
n_crashed_or_failed <- sum(is.na(dev_boot))
c_hat_mean   <- deviance_obs / mean(dev_boot_ok)
c_hat_median <- deviance_obs / median(dev_boot_ok)
boot_p_value <- mean(dev_boot_ok >= deviance_obs)

cat("\nObserved deviance:", round(deviance_obs, 2), "\n")
cat("Bootstrap replicates: ", length(dev_boot_ok), "ok /", n_crashed_or_failed, "crashed or failed\n")
cat("Bootstrap deviance mean/median/sd:", round(mean(dev_boot_ok), 2), "/",
    round(median(dev_boot_ok), 2), "/", round(sd(dev_boot_ok), 2), "\n")
cat("c-hat (mean-based):", round(c_hat_mean, 3), "\n")
cat("c-hat (median-based):", round(c_hat_median, 3), "\n")
cat("Bootstrap p-value:", round(boot_p_value, 3), "\n")

p_gof <- ggplot(data.frame(deviance = dev_boot_ok), aes(x = deviance)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = deviance_obs, color = "red", linewidth = 1) +
  labs(title = "Parametric bootstrap GOF: IP5 annual classic GIEL CJS model",
       subtitle = paste0("c-hat = ", round(c_hat_mean, 2), " (mean-based), bootstrap p = ",
                          round(boot_p_value, 3), "; ", n_crashed_or_failed, " of ", B,
                          " replicates crashed/failed and were excluded"),
       x = "Deviance (relative to cohort-saturated model)",
       y = paste0("Bootstrap replicates (n=", length(dev_boot_ok), ")"))
ggsave("model/GIEL_ClassicCJS_RobustDesign_Annual_bootGOF.png", p_gof, width = 7, height = 5, dpi = 150)

# ---------------------------------------------------------------------------
# 7. Comparison against the GIEL-only Bayesian robust-design fit (IP5 = pond 2)
# ---------------------------------------------------------------------------
bayes_env <- new.env()
load("data/BWRobustDesign_MCMC_giel.RData", envir = bayes_env)
# rd_Fish$pond_new: IP2=1, IP5=2, IP6=3 (per AGENTS.md / GIEL_ModelMethodology.md)
samp <- do.call(rbind, bayes_env$rd_samples)
get_node <- function(pattern) samp[, grepl(pattern, colnames(samp)), drop = FALSE]

lphiA_ip5 <- get_node("^lphiA\\[2, ")
dE_ip5    <- get_node("^dE\\[2\\]")
lp_ip5    <- get_node("^lp\\[2, ")

season_fy <- FY_MIN:FY_MAX
annual_surv_bayes <- sapply(seq_len(ncol(lphiA_ip5)), function(y) {
  is_event <- (season_fy[y] == 2020)
  logit_phi <- lphiA_ip5[, y] + is_event * dE_ip5[, 1]
  phi_month <- plogis(logit_phi)
  mean(phi_month^12)
})
names(annual_surv_bayes) <- season_fy

cat("\nBayesian (GIEL-only fit) IP5 annualized adult survival by FY (posterior mean, phi^12):\n")
print(round(annual_surv_bayes, 3))

# ---------------------------------------------------------------------------
# 8. Empirical "floor" comparison for the FY2020 die-off gap (annual occasions
#    make this direct: known alive at the FY2020 occasion, re-contacted at
#    the FY2021 occasion or later)
# ---------------------------------------------------------------------------
occ_fy2020 <- fy_to_occ(2020L)
occ_fy2021 <- fy_to_occ(2021L)
M_full <- as.matrix(giel_ip5[, occ_cols]); mode(M_full) <- "integer"
alive_fy2020 <- giel_ip5$entry_occ <= occ_fy2020
seen_after_2020 <- rowSums(M_full[, occ_fy2021:T_OCC, drop = FALSE]) > 0
floor_2020 <- c(alive = sum(alive_fy2020), seen_after = sum(alive_fy2020 & seen_after_2020))
cat("\nEmpirical floor, IP5 stocked cohort known-alive FY2020 seen again FY2021+:",
    floor_2020["seen_after"], "/", floor_2020["alive"], "=",
    round(floor_2020["seen_after"] / floor_2020["alive"], 3), "\n")

save(pond_check, giel_ip5, dp_ip5, ddl_ip5, model_specs_ip5, fits_ip5,
     aic_ip5, top_name_ip5, top_ip5, preds_ip5, n_suppressed,
     deviance_obs, dev_boot, c_hat_mean, c_hat_median, boot_p_value, n_crashed_or_failed,
     annual_surv_bayes, floor_2020,
     file = "data/GIEL_ClassicCJS_RobustDesign_Annual.RData")

cat("\nSaved data/GIEL_ClassicCJS_RobustDesign_Annual.RData and",
    "model/GIEL_ClassicCJS_RobustDesign_Annual_bootGOF.png\n")
cat("\n=== SCRIPT COMPLETE ===\n")
