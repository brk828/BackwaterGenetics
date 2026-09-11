# Classic (frequentist, MLE) robust-design CJS model for Bonytail (GIEL) at
# IPCA Pond 5 (IP5) -- stocked fish only. Per user request 2026-09-11:
#  - Pick the pond with the most/best-quality stocking data (TL at release):
#    checked IP2 (599 stocked, only 299 with TL), IP5 (600 stocked, ALL with
#    TL, 2 cohorts: 2017 [280 adult+20 young] and 2022 [300 adult]), IP6 (300
#    stocked, all with TL, 1 cohort only). IP5 chosen: most fish, complete TL,
#    two independent release cohorts, and a genuine size-class contrast (20
#    young fish in 2017). Netting-tagged ("wild"/captured-tagged-released)
#    fish are excluded entirely per instructions -- stocked fish only.
#  - Occasions: 70 monthly Oct-Apr occasions (FY2017-FY2026, 7 months/season,
#    10 seasons). Within a season, consecutive monthly transitions are FIXED
#    to Phi=1 (closed session -- no mortality assumed Oct-Apr); only the
#    April(FY y) -> October(FY y+1) transition (9 such gaps) is a free
#    survival parameter -- "all mortality occurs between April and October".
#  - If a fish was stocked at any point during the Oct-Apr window, it is
#    treated as released in **April** of that FY (last occasion of that
#    season) and every contact during that same Oct-Apr window is ignored
#    (zeroed) -- the first survival opportunity is always the April->October
#    gap immediately following release ("post-release mortality is first
#    between-sessions interval").
#  - Size class: only the fish's own FIRST gap (immediately after release) is
#    classified "young" (if TL < 250 mm at tagging, the pond's existing GIEL
#    stage threshold) or "adult"; every subsequent gap is "adult" for every
#    fish, matching the classic annual CJS scripts' "young for one interval,
#    then adult" logic.
#  - Phi structure (3 candidates, sizeclass always included):
#      const  : ~ sizeclass                      (single survival for all 9 gaps)
#      gap    : ~ sizeclass + gap                 (each of 9 gaps free)
#      dieoff : ~ sizeclass + dieoff               (baseline + 1 dummy for the
#                 known IP5 die-off gap, FY2020->FY2021, per data/BWDieOffEvents.csv)
#  - p structure (2 candidates): ~ moy (7-level month-of-season factor) or ~1.
#  - 3 x 2 = 6 models, ranked by AIC.
#  - Parametric bootstrap c-hat on the top model (same recipe as
#    GIEL_ClassicCJS_PondGrouped.R): simulate every real fish's history
#    forward under the fitted model, refit, deviance vs a cohort-saturated
#    reference (cohort = entry FY x sizeclass-at-release), c-hat = observed
#    deviance / mean(bootstrap deviance).
#  - Compare survival (esp. the FY2020 die-off gap) and monthly detection
#    against the existing GIEL-only Bayesian robust-design fit
#    (data/BWRobustDesign_MCMC_giel.RData).

source("LabFunctions.R")
packages(dplyr); packages(lubridate); packages(tidyr); packages(marked); packages(ggplot2)

load("data/ReportingData.RData")

fy_of <- function(d) year(d) + as.integer(month(d) > 9)
FY_MIN <- 2017L; FY_MAX <- 2026L
occ_fy <- FY_MIN:FY_MAX
n_fy <- length(occ_fy)
T_OCC <- n_fy * 7
season_months <- c(10, 11, 12, 1, 2, 3, 4)

OccLookup <- expand.grid(local_pos = 1:7, fy_i = 1:n_fy) |>
  arrange(fy_i, local_pos) |>
  mutate(occ = row_number(), fy = occ_fy[fy_i], month = season_months[local_pos])

# ---------------------------------------------------------------------------
# 1. Pond selection check (printed for the record; IP5 is used below)
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
            stage0 = ifelse(total_length < 250, "young", "adult"))

n_transfer <- sum(StudyBWAnalysis$species == "GIEL" &
                     StudyBWAnalysis$location == "IPCA (Pond 5)" &
                     StudyBWAnalysis$event == "stocking" & StudyBWAnalysis$Transfer == 1)
if (n_transfer > 0) stop("Found ", n_transfer, " IP5 stocked GIEL flagged Transfer==1 -- decide handling.")

entry_lookup <- OccLookup |> filter(local_pos == 7) |> select(fy, entry_occ = occ)
GielIP5 <- GielIP5 |> left_join(entry_lookup, by = c("entry_fy" = "fy"))

cat("IP5 stocked fish:", nrow(GielIP5), "\n")
cat("By cohort/stage0:\n"); print(table(GielIP5$entry_fy, GielIP5$stage0))

# ---------------------------------------------------------------------------
# 3. Detection matrix (Oct-Apr months only) and entry/ignore logic
# ---------------------------------------------------------------------------
GielIP5Contacts <- StudyBWContacts |>
  filter(Species == "GIEL", Backwater == "IPCA (Pond 5)", !is.na(Date),
         PITIndex %in% GielIP5$PITIndex) |>
  mutate(mo = month(Date), fy = fy_of(Date)) |>
  filter(mo %in% season_months, fy %in% occ_fy) |>
  distinct(PITIndex, fy, mo) |>
  left_join(OccLookup |> select(fy, month, occ), by = c("fy" = "fy", "mo" = "month")) |>
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
cat("\nSuppressed (ignored, pre-April-entry) contact-occasions:", n_suppressed, "\n")

# ---------------------------------------------------------------------------
# 4. Process data; build Phi design covariates (fix, sizeclass, gap, dieoff)
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

ddl_ip5$Phi$fix <- ifelse(occ_num %% 7 == 0, NA, 1)  # within-season: Phi fixed = 1

gap_index <- ifelse(occ_num %% 7 == 0, occ_num %/% 7, NA)
gap_fy <- FY_MIN + gap_index - 1
# "gap" and "dieoff" must have no NA (fixed rows still need a valid factor
# level for the design matrix, even though their real Phi is overridden by fix)
ddl_ip5$Phi$gap <- factor(ifelse(is.na(gap_fy), "none", as.character(gap_fy)),
                            levels = c("none", as.character(FY_MIN:(FY_MAX - 1))))
ddl_ip5$Phi$dieoff <- factor(ifelse(!is.na(gap_fy) & gap_fy == 2020, "event", "normal"),
                               levels = c("normal", "event"))

occ_num_p <- as.numeric(as.character(ddl_ip5$p$occ))
local_pos_p <- ((occ_num_p - 1) %% 7) + 1
ddl_ip5$p$moy <- factor(local_pos_p, levels = 1:7,
                          labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))

cat("\nPhi design rows:", nrow(ddl_ip5$Phi), "(fixed:", sum(!is.na(ddl_ip5$Phi$fix)),
    ", free/gap:", sum(is.na(ddl_ip5$Phi$fix)), ")\n")
cat("sizeclass table:\n"); print(table(ddl_ip5$Phi$sizeclass))
cat("gap table:\n"); print(table(ddl_ip5$Phi$gap))
cat("dieoff table:\n"); print(table(ddl_ip5$Phi$dieoff))
cat("moy table:\n"); print(table(ddl_ip5$p$moy))

# ---------------------------------------------------------------------------
# 5. Fit the 6-model set
# ---------------------------------------------------------------------------
Phi.const  <- list(formula = ~ sizeclass)
Phi.gap    <- list(formula = ~ sizeclass + gap)
Phi.dieoff <- list(formula = ~ sizeclass + dieoff)
p.moy      <- list(formula = ~ moy)
p.dot      <- list(formula = ~ 1)

model_specs_ip5 <- list(
  const_pmoy  = list(Phi = Phi.const,  p = p.moy),
  const_pdot  = list(Phi = Phi.const,  p = p.dot),
  gap_pmoy    = list(Phi = Phi.gap,    p = p.moy),
  gap_pdot    = list(Phi = Phi.gap,    p = p.dot),
  dieoff_pmoy = list(Phi = Phi.dieoff, p = p.moy),
  dieoff_pdot = list(Phi = Phi.dieoff, p = p.dot)
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

save.image_partial <- function() {
  save(pond_check, giel_ip5, dp_ip5, ddl_ip5, model_specs_ip5, fits_ip5,
       aic_ip5, top_name_ip5, top_ip5, preds_ip5, n_suppressed,
       file = "data/GIEL_ClassicCJS_RobustDesign.RData")
}
save.image_partial()
cat("\nSaved intermediate results to data/GIEL_ClassicCJS_RobustDesign.RData\n")

# ---------------------------------------------------------------------------
# 6. Parametric bootstrap c-hat on the top model
# ---------------------------------------------------------------------------
# Cohort for the saturated reference = (entry_fy, stage0) -- the only
# covariate information available to differentiate individuals at release.
sim_base <- giel_ip5 |>
  transmute(entry_fy, stage0, entry_occ,
            cohort = paste(entry_fy, stage0))

# Generic Phi lookup, valid regardless of which of the 3 Phi structures
# (const/gap/dieoff) won: predict() dedupes to one row per unique *free*
# real-parameter value, plus one arbitrary representative row for the
# fixed (within-season, Phi=1, se=0) group -- se>0 excludes that
# placeholder row when matching by sizeclass alone (the "const" case).
# Same dedup issue applies to p: predict() with p ~ moy returns one row per
# unique (moy) real parameter (plus se=0 placeholder rows for the initial
# occasions before any fish are at risk), so lookup must go through moy,
# not raw occ, and must generalize to the p ~ 1 (dot) structure too.
build_get_p <- function(preds_p) {
  p_covs <- setdiff(names(preds_p), c("occ", "estimate", "se", "lcl", "ucl"))
  moy_labels <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
  if ("moy" %in% p_covs) {
    good <- preds_p$se > 0
    p_by_moy <- setNames(preds_p$estimate[good], as.character(preds_p$moy[good]))
    p_by_moy <- p_by_moy[!duplicated(names(p_by_moy))]
    function(occ) {
      local_pos <- ((occ - 1) %% 7) + 1
      val <- p_by_moy[[moy_labels[local_pos]]]
      if (is.null(val) || is.na(val)) 0 else val
    }
  } else {
    good <- preds_p$se > 0
    val0 <- preds_p$estimate[good][1]
    function(occ) if (is.na(val0)) 0 else val0
  }
}
get_p <- build_get_p(preds_ip5$p)

build_get_phi <- function(preds_Phi) {
  phi_covs <- setdiff(names(preds_Phi), c("occ", "estimate", "se", "lcl", "ucl", "sizeclass"))
  function(occ, sizeclass_i) {
    if (occ %% 7 != 0) return(1)  # within-season, fixed
    row_sc <- preds_Phi$sizeclass == sizeclass_i & preds_Phi$se > 0
    if (length(phi_covs) == 0) {
      val <- preds_Phi$estimate[row_sc]
      if (length(val) == 0) val <- preds_Phi$estimate[preds_Phi$sizeclass == "adult" & preds_Phi$se > 0]
      return(val[1])
    }
    g <- as.character(FY_MIN + occ %/% 7 - 1)
    if ("gap" %in% phi_covs) {
      key_match <- row_sc & as.character(preds_Phi$gap) == g
      val <- preds_Phi$estimate[key_match]
      if (length(val) == 0) val <- preds_Phi$estimate[preds_Phi$sizeclass == "adult" & preds_Phi$se > 0 & as.character(preds_Phi$gap) == g]
      return(val[1])
    }
    if ("dieoff" %in% phi_covs) {
      dstat <- ifelse(g == "2020", "event", "normal")
      key_match <- row_sc & as.character(preds_Phi$dieoff) == dstat
      val <- preds_Phi$estimate[key_match]
      if (length(val) == 0) val <- preds_Phi$estimate[preds_Phi$sizeclass == "adult" & preds_Phi$se > 0 & as.character(preds_Phi$dieoff) == dstat]
      return(val[1])
    }
    NA
  }
}
get_phi <- build_get_phi(preds_ip5$Phi)

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

obs_ch_marked <- dp_ip5$data$ch[match(giel_ip5$ch, apply(matrix(unlist(strsplit(dp_ip5$data$ch, ",")),
                                                                  nrow = nrow(dp_ip5$data), byrow = TRUE),
                                                          1, paste, collapse = ""))]
# Simpler/robust: just recompute comma-joined ch straight from giel_ip5$ch
obs_ch <- sapply(strsplit(giel_ip5$ch, ""), paste, collapse = ",")

sat_ll_obs   <- sat_loglik(obs_ch, sim_base$cohort)
neg2lnl_obs  <- top_ip5$results$neg2lnl
deviance_obs <- neg2lnl_obs - (-2 * sat_ll_obs)

B <- 300
set.seed(2854)
seeds <- sample.int(1e6, B)
boot_results <- vector("list", B)
top_formula_ip5 <- model_specs_ip5[[top_name_ip5]]

t0 <- Sys.time()
for (b in seq_len(B)) {
  set.seed(seeds[b])
  ch_sim <- simulate_ch(sim_base, T_OCC)
  df <- data.frame(ch = ch_sim, stage0 = giel_ip5$stage0)
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
  ddl_b$Phi$fix <- ifelse(occ_num_b %% 7 == 0, NA, 1)
  gap_index_b <- ifelse(occ_num_b %% 7 == 0, occ_num_b %/% 7, NA)
  gap_fy_b <- FY_MIN + gap_index_b - 1
  ddl_b$Phi$gap <- factor(ifelse(is.na(gap_fy_b), "none", as.character(gap_fy_b)),
                            levels = c("none", as.character(FY_MIN:(FY_MAX - 1))))
  ddl_b$Phi$dieoff <- factor(ifelse(!is.na(gap_fy_b) & gap_fy_b == 2020, "event", "normal"),
                               levels = c("normal", "event"))
  occ_num_pb <- as.numeric(as.character(ddl_b$p$occ))
  local_pos_pb <- ((occ_num_pb - 1) %% 7) + 1
  ddl_b$p$moy <- factor(local_pos_pb, levels = 1:7,
                          labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))

  fit_b <- tryCatch(
    crm(dp_b, ddl_b, model.parameters = top_formula_ip5, initial = top_ip5,
        hessian = FALSE, accumulate = TRUE),
    error = function(e) NULL)
  if (is.null(fit_b)) {
    boot_results[[b]] <- list(deviance = NA, convergence = NA)
  } else {
    sat_ll_b <- sat_loglik(ch_sim, sim_base$cohort)
    boot_results[[b]] <- list(deviance = fit_b$results$neg2lnl - (-2 * sat_ll_b),
                               convergence = fit_b$results$convergence)
  }
  if (b %% 25 == 0) cat("bootstrap", b, "of", B, "-",
                          round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
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
  labs(title = "Parametric bootstrap GOF: IP5 robust-design classic GIEL CJS model",
       subtitle = paste0("c-hat = ", round(c_hat_mean, 2), " (mean-based), bootstrap p = ",
                          round(boot_p_value, 3)),
       x = "Deviance (relative to cohort-saturated model)",
       y = paste0("Bootstrap replicates (n=", B, ")"))
ggsave("model/GIEL_ClassicCJS_RobustDesign_bootGOF.png", p_gof, width = 7, height = 5, dpi = 150)

# ---------------------------------------------------------------------------
# 7. Comparison against the GIEL-only Bayesian robust-design fit (IP5 = pond 2)
# ---------------------------------------------------------------------------
bayes_env <- new.env()
load("data/BWRobustDesign_MCMC_giel.RData", envir = bayes_env)
# rd_Fish$pond_new: IP2=1, IP5=2, IP6=3 (per AGENTS.md / GIEL_ModelMethodology.md)
samp <- do.call(rbind, bayes_env$rd_samples)
get_node <- function(pattern) samp[, grepl(pattern, colnames(samp)), drop = FALSE]

# lphiA[2, y] for y = 1..10 (IP5 season-specific adult survival, logit scale)
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
# 8. Empirical "floor" comparison for the FY2020 die-off gap
# ---------------------------------------------------------------------------
# Fish known alive (stocked cohort at IP5) at the April-2020 occasion,
# re-contacted at/after the October-2020 occasion (next season)
apr2020_occ <- OccLookup$occ[OccLookup$fy == 2020 & OccLookup$local_pos == 7]
oct2020_occ <- OccLookup$occ[OccLookup$fy == 2021 & OccLookup$local_pos == 1]
M_full <- as.matrix(giel_ip5[, occ_cols]); mode(M_full) <- "integer"
alive_apr2020 <- giel_ip5$entry_occ <= apr2020_occ
seen_after_2020 <- rowSums(M_full[, (oct2020_occ):T_OCC, drop = FALSE]) > 0
floor_2020 <- c(alive = sum(alive_apr2020), seen_after = sum(alive_apr2020 & seen_after_2020))
cat("\nEmpirical floor, IP5 stocked cohort known-alive Apr-2020 seen again Oct-2020+:",
    floor_2020["seen_after"], "/", floor_2020["alive"], "=",
    round(floor_2020["seen_after"] / floor_2020["alive"], 3), "\n")

save(pond_check, giel_ip5, dp_ip5, ddl_ip5, model_specs_ip5, fits_ip5,
     aic_ip5, top_name_ip5, top_ip5, preds_ip5, n_suppressed,
     deviance_obs, dev_boot, c_hat_mean, c_hat_median, boot_p_value,
     annual_surv_bayes, floor_2020,
     file = "data/GIEL_ClassicCJS_RobustDesign.RData")

cat("\nSaved data/GIEL_ClassicCJS_RobustDesign.RData and",
    "model/GIEL_ClassicCJS_RobustDesign_bootGOF.png\n")
cat("\n=== SCRIPT COMPLETE ===\n")
