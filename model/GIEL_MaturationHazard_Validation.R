# Validation checklist for the GIEL size-tied maturation-hazard sensitivity
# fits (Tracks F1, F2, A1) from model/BWRobustDesign_NIMBLE.R, GIEL-only
# (RUN_PONDS = c(3, 6, 7); pond_new 1 = IP2, 2 = IP5, 3 = IP6).
#
# Checks, per the plan's Step 4 / AGENTS.md "on completion" checklist:
#   1. N_tag >= known alive, for every pond-month (posterior 2.5% quantile).
#   2. lpsi_adj[GIEL] plausibility.
#   3. Hyperparameter convergence (Rhat).
#   4. Comparison across tracks and against the pre-hazard baseline fit
#      (data/BWRobustDesign_MCMC_giel.RData).
#
# Does not require ReportingData.RData or a rebuild of BWRobustDesign_data.RData
# -- everything needed (Fish table, posterior samples) is stored in each
# fit's own .RData file.

source("LabFunctions.R")
packages(dplyr)
packages(coda)
packages(tidyr)
packages(purrr)

TRACKS <- c(F1 = "data/BWRobustDesign_MCMC_giel_F1.RData",
            F2 = "data/BWRobustDesign_MCMC_giel_F2.RData",
            A1 = "data/BWRobustDesign_MCMC_giel_A1.RData")
BASELINE_FILE <- "data/BWRobustDesign_MCMC_giel.RData"

POND_LABELS <- c("IP2", "IP5", "IP6")   # pond_new 1:3
N_MOY <- 7L
FIRST_FY <- 2017L

load_fit <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e
}

# ---------------------------------------------------------------------------
# 1. Known-alive constraint check: sum(N_tag[pond,t,1:2]) >= known-alive count
# ---------------------------------------------------------------------------
known_alive_check <- function(env) {
  Fish <- env$rd_Fish
  samp <- as.matrix(env$rd_samples)
  T_PRIM <- max(Fish$end_t, na.rm = TRUE)
  # rebuild N_tag[pond,t,1:2] sum from the flat posterior matrix
  ntag_cols <- grep("^N_tag\\[", colnames(samp), value = TRUE)
  idx <- do.call(rbind, lapply(strsplit(gsub("N_tag\\[|\\]", "", ntag_cols), ", "),
                                function(x) as.integer(x)))
  colnames(idx) <- c("pond", "t", "stage")
  df_idx <- as_tibble(idx) |> mutate(col = ntag_cols)

  violations <- tibble()
  for (j in 1:3) {
    for (t in 1:T_PRIM) {
      known_alive <- Fish |>
        filter(pond_new == j, entry_t <= t, end_t >= t,
               last_alive_t >= t, count_from <= t) |>
        nrow()
      if (known_alive == 0) next
      cols <- df_idx |> filter(pond == j, t == !!t) |> pull(col)
      if (length(cols) == 0) next
      post_sum <- rowSums(samp[, cols, drop = FALSE])
      q025 <- quantile(post_sum, 0.025)
      if (q025 < known_alive) {
        violations <- bind_rows(violations,
                                 tibble(pond = j, t = t, known_alive = known_alive,
                                        post_q025 = q025, post_mean = mean(post_sum)))
      }
    }
  }
  list(n_pond_months_checked = sum(sapply(1:3, function(j)
    sum(sapply(1:T_PRIM, function(t)
      Fish |> filter(pond_new == j, entry_t <= t, end_t >= t,
                      last_alive_t >= t, count_from <= t) |> nrow() > 0)))),
    violations = violations)
}

# ---------------------------------------------------------------------------
# 2. Convergence summary (excluding known nuisance families)
# ---------------------------------------------------------------------------
NUISANCE_FAMILIES <- c("N_tag", "Tal", "U", "U0", "RecJ", "r", "Atot")

convergence_summary <- function(env) {
  rhat <- env$rhat_all[, "Point est."]
  fam <- gsub("\\[.*", "", names(rhat))
  core <- rhat[!fam %in% NUISANCE_FAMILIES]
  nuis <- rhat[fam %in% NUISANCE_FAMILIES]
  list(
    n_core = length(core), max_rhat_core = max(core, na.rm = TRUE),
    n_core_over_1.1 = sum(core > 1.1, na.rm = TRUE),
    worst_core = sort(core, decreasing = TRUE)[1:5],
    n_nuisance = length(nuis), max_rhat_nuisance = max(nuis, na.rm = TRUE),
    n_nuisance_over_1.1 = sum(nuis > 1.1, na.rm = TRUE)
  )
}

# ---------------------------------------------------------------------------
# 3. lpsi_adj[GIEL] posterior
# ---------------------------------------------------------------------------
lpsi_adj_summary <- function(env) {
  s <- env$rd_summary
  if ("lpsi_adj[1]" %in% rownames(s)) s["lpsi_adj[1]", ] else NA
}

# ---------------------------------------------------------------------------
# 4. Core hyperparameter table across fits
# ---------------------------------------------------------------------------
CORE_PARS <- c("mu_phi[1]", "sigma_phi[1]", "dJ[1]",
                "dE[1]", "dE[2]", "dE[3]",
                "b_post[1, 1]", "b_post[1, 2]", "b_post[1, 3]",
                "lpsi[1]", "lpsi_adj[1]",
                "lp[1, 1]", "lp[2, 1]", "lp[3, 1]",
                "lp[1, 2]", "lp[2, 2]", "lp[3, 2]",
                "g_eff[1]", "mu_r[1]", "sigma_r[1]", "nu[1]")

# ---------------------------------------------------------------------------
# Run for each track + baseline
# ---------------------------------------------------------------------------
envs <- lapply(TRACKS, load_fit)
env_base <- load_fit(BASELINE_FILE)

cat("=== Known-alive constraint checks ===\n")
ka_results <- list()
for (nm in names(envs)) {
  cat("\n--", nm, "--\n")
  res <- known_alive_check(envs[[nm]])
  ka_results[[nm]] <- res
  cat("Pond-months checked (known_alive > 0):", res$n_pond_months_checked, "\n")
  cat("Violations (post 2.5% quantile < known-alive):", nrow(res$violations), "\n")
  if (nrow(res$violations) > 0) {
    res$violations <- res$violations |>
      mutate(pond_label = POND_LABELS[pond],
             season = (t - 1) %/% N_MOY + 1, FY = FIRST_FY + season - 1)
    print(res$violations)
  }
}

cat("\n\n=== Convergence summary ===\n")
conv_results <- list()
for (nm in names(envs)) {
  cat("\n--", nm, "--\n")
  cs <- convergence_summary(envs[[nm]])
  conv_results[[nm]] <- cs
  cat("Core params: n =", cs$n_core, " max Rhat =", round(cs$max_rhat_core, 3),
      " n > 1.1:", cs$n_core_over_1.1, "\n")
  cat("Worst core Rhats:\n"); print(round(cs$worst_core, 3))
  cat("Nuisance params (N_tag/Tal/U/U0/RecJ/r/Atot): n =", cs$n_nuisance,
      " max Rhat =", round(cs$max_rhat_nuisance, 3),
      " n > 1.1:", cs$n_nuisance_over_1.1, "\n")
}

cat("\n\n=== lpsi_adj[GIEL] by track ===\n")
for (nm in names(envs)) {
  cat(nm, ":\n"); print(lpsi_adj_summary(envs[[nm]]))
}

cat("\n\n=== Core hyperparameter comparison (F1 / F2 / A1 / baseline) ===\n")
core_table <- map_dfr(names(envs), function(nm) {
  s <- envs[[nm]]$rd_summary
  pars <- intersect(CORE_PARS, rownames(s))
  tibble(track = nm, param = pars, mean = s[pars, "mean"],
         lo = s[pars, "2.5%"],
         hi = s[pars, "97.5%"])
})
base_pars <- intersect(setdiff(CORE_PARS, "lpsi_adj[1]"), rownames(env_base$rd_summary))
core_table_base <- tibble(track = "baseline", param = base_pars,
                            mean = env_base$rd_summary[base_pars, "mean"],
                            lo = env_base$rd_summary[base_pars, "2.5%"],
                            hi = env_base$rd_summary[base_pars, "97.5%"])
core_table_all <- bind_rows(core_table, core_table_base)
print(core_table_all |> pivot_wider(names_from = track, values_from = c(mean, lo, hi)), n = 30)

# ---------------------------------------------------------------------------
# 5. FY2022 event-adjusted realized annual survival by pond (baseline
#    comparison metric used throughout AGENTS.md)
# ---------------------------------------------------------------------------
event_windows <- list(
  IP2 = list(season = 6, months = 7),   # Jan-Jul 2022, per BWDieOffEvents.csv
  IP5 = list(season = 4, months = 12),
  IP6 = list(season = NA, months = 0)
)

annual_survival_by_pond <- function(env) {
  s <- env$rd_summary
  lphiA_rows <- grep("^lphiA\\[", rownames(s), value = TRUE)
  idx <- do.call(rbind, lapply(strsplit(gsub("lphiA\\[|\\]", "", lphiA_rows), ", "), as.integer))
  tibble(pond = idx[, 1], season = idx[, 2], row = lphiA_rows,
         lphiA_mean = s[lphiA_rows, "Mean"]) |>
    mutate(pond_label = POND_LABELS[pond],
           phi_month = plogis(lphiA_mean),
           annual_surv = phi_month^12)
}

cat("\n\n=== Baseline (no dE/nEvt in lphiA) annual adult survival by pond-season, for reference ===\n")
print(annual_survival_by_pond(env_base) |> arrange(pond, season))

save(ka_results, conv_results, core_table_all, envs, env_base,
     file = "data/GIEL_MaturationHazard_Validation.RData")
cat("\nSaved data/GIEL_MaturationHazard_Validation.RData\n")
