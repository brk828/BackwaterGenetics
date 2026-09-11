# Joint stage-structured robust-design model for backwater XYTE and GIEL
# populations, fitted in NIMBLE.
# B. Kesner / Posit Assistant, September 2026
#
# Requires data/BWRobustDesign_data.RData (from model/BWRobustDesign_data.R)
# and model/BWRobustDesign_defs.R (nimbleFunctions, model code, inits).
# Saves data/BWRobustDesign_MCMC.RData and model/BWRobustDesign_traceplots.pdf.
#
# Structure
#   Tagged fish: hidden-state process over primary periods (months Oct-Apr)
#     with states J (juvenile), A (adult), Dead; removal (harvest/transfer) is a
#     known event that ends the history with the fish alive. Each pond's fish
#     histories are marginalized by a forward-backward pass (rdPond) that also
#     returns the smoothed expected number of tagged fish alive by stage (N_tag),
#     used for the netting likelihood.
#   Untagged pool: annual per-pond juvenile / adult counts with shared survival,
#     maturation, per-adult recruitment, and known removals (tagging, harvest).
#   Netting join: previously tagged among handled fish ~ Binomial(n, T/(T+U)).
# Chains run in parallel (one R worker per chain).

source("LabFunctions.R")
packages(dplyr)
packages(nimble)
packages(coda)
packages(MCMCvis)
packages(parallel)

load("data/BWRobustDesign_data.RData")

# Stale-data guard (2026-09-10; see "Known Code Issues" in AGENTS.md): this
# script does not rebuild data/BWRobustDesign_data.RData, so a fit run after
# editing model/BWRobustDesign_data.R, data/ReportingData.RData, or
# data/BWNettingEvents.csv/BWDieOffEvents.csv without re-sourcing the data
# script would silently fit stale inputs (this happened once before, see the
# "stale-data pitfall" note in AGENTS.md). BUILD_INFO is only present in data
# builds from 2026-09-10 onward; older files skip the check with a warning.
if (exists("BUILD_INFO")) {
  cur_mtime <- file.info(BUILD_INFO$inputs)$mtime
  stale <- cur_mtime > BUILD_INFO$input_mtime
  if (any(stale, na.rm = TRUE)) {
    stop("data/BWRobustDesign_data.RData is stale relative to: ",
         paste(BUILD_INFO$inputs[which(stale)], collapse = ", "),
         ". Re-run model/BWRobustDesign_data.R before fitting.")
  }
} else {
  warning("data/BWRobustDesign_data.RData has no BUILD_INFO (pre-2026-09-10 build); ",
          "cannot verify it is current. Re-run model/BWRobustDesign_data.R if in doubt.")
}

source("model/BWRobustDesign_defs.R")

# ---------------------------------------------------------------------------
# Run settings
# ---------------------------------------------------------------------------
if (!exists("RUN_PONDS")) RUN_PONDS <- 1:7   # pond indices (1 = YCB, 2:7 = IP1..IP6); use >= 2 ponds
N_CHAINS  <- 3L
N_ITER    <- 30000L
N_BURNIN  <- 10000L
N_THIN    <- 10L
SEEDS     <- c(4127L, 8351L, 2960L)
PARALLEL  <- TRUE
RUN_TAG <- if (identical(sort(RUN_PONDS), 1:7)) "" else
  if (identical(sort(RUN_PONDS), c(3L, 6L, 7L))) "_giel" else
  paste0("_p", paste(sort(RUN_PONDS), collapse = ""))
if (!exists("OUT_FILE")) OUT_FILE <- paste0("data/BWRobustDesign_MCMC", RUN_TAG, ".RData")
if (!exists("TRACE_FILE")) TRACE_FILE <- paste0("model/BWRobustDesign_traceplots", RUN_TAG, ".pdf")

inp <- make_rd_inputs(RUN_PONDS)
cat("Ponds:", paste(PondTable$Backwater[RUN_PONDS], collapse = ", "),
    "| fish:", nrow(inp$Fish), "| unique histories:", nrow(inp$Hist),
    "| netting rows:", inp$constants$NEV, "\n")

# ---------------------------------------------------------------------------
# Run chains
# ---------------------------------------------------------------------------
t0 <- Sys.time()
if (PARALLEL && N_CHAINS > 1) {
  cl <- makeCluster(N_CHAINS)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("inp", "N_ITER", "N_BURNIN", "N_THIN"))
  clusterEvalQ(cl, {
    library(dplyr); library(nimble); library(coda)
    load("data/BWRobustDesign_data.RData")
    source("model/BWRobustDesign_defs.R")
  })
  chain_list <- parLapply(cl, SEEDS[seq_len(N_CHAINS)], function(s)
    run_rd_chain(inp, s, N_ITER, N_BURNIN, N_THIN, progress = FALSE))
  stopCluster(cl); on.exit()
} else {
  chain_list <- lapply(SEEDS[seq_len(N_CHAINS)], function(s)
    run_rd_chain(inp, s, N_ITER, N_BURNIN, N_THIN))
}
rd_samples <- as.mcmc.list(chain_list)
run_time <- difftime(Sys.time(), t0, units = "mins")
cat("Total wall time (compile + MCMC):", round(run_time, 1), "min\n")

# ---------------------------------------------------------------------------
# Diagnostics and save
# ---------------------------------------------------------------------------
core_pars <- c("mu_phi", "sigma_phi", "dJ", "dE", "b_post", "lpsi", "lp", "g_eff",
               "mu_r", "sigma_r", "nu", "U0")
rd_summary <- MCMCsummary(rd_samples, params = core_pars, round = 3)
print(rd_summary)

rhat_all <- if (N_CHAINS > 1) gelman.diag(rd_samples, multivariate = FALSE)$psrf else NULL
if (!is.null(rhat_all)) {
  cat("\nMax R-hat (all monitored):", round(max(rhat_all[, 1], na.rm = TRUE), 3), "\n")
  cat("Parameters with R-hat > 1.05:", sum(rhat_all[, 1] > 1.05, na.rm = TRUE), "\n")
}

MCMCtrace(rd_samples, params = core_pars, pdf = TRUE,
          filename = TRACE_FILE)

rd_Fish <- inp$Fish
rd_Hist <- inp$Hist
rd_Netting <- inp$Netting
rd_ponds <- RUN_PONDS
rd_settings <- list(N_CHAINS = N_CHAINS, N_ITER = N_ITER, N_BURNIN = N_BURNIN,
                    N_THIN = N_THIN, SEEDS = SEEDS, run_time = run_time)
save(rd_samples, rd_summary, rhat_all, rd_code, rd_Fish, rd_Hist, rd_Netting, rd_ponds,
     rd_settings, file = OUT_FILE)
cat("Saved", OUT_FILE, "\n")
