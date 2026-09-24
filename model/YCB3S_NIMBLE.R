# Three-stage, density-dependent robust-design model — fit driver.
# B. Kesner / Posit Assistant, September 2026
#
# Requires the data build from model/YCB3S_data.R for the chosen PONDS and
# CALENDAR, and model/YCB3S_defs.R. File names carry a pond tag ("" for YCB
# only, "_xyte" for c(1, 2, 4, 5), "_p<ponds>" otherwise) and a calendar
# suffix ("" for m7, "_m12"):
#   data/YCB3S_data<tag><sfx>.RData -> data/YCB3S_MCMC<tag><sfx>.RData,
#   model/YCB3S_traceplots<tag><sfx>.pdf
# e.g. YCB m7: data/YCB3S_MCMC.RData; all XYTE m12: data/YCB3S_MCMC_xyte_m12.RData.
# Run as a background job:
#   rstudioapi::jobRunScript("model/YCB3S_NIMBLE.R", workingDir = getwd(), importEnv = FALSE)
#
# Structure (see the plan file for rationale)
#   Tagged fish: hidden states S / L / A / Dead over Oct-Apr months, marginalized
#     per pond by rdPond3 (forward-backward), which also returns smoothed tagged
#     abundance by stage (N_tag). Removal ends a history with the fish alive.
#   Survival: logit phiA = lphiA[pond, season]; phiL adds dJ + eps_J[pond, season] (cB fixed at 0);
#     phiS adds dS on top. dens = standardized log recruit density observed at
#     the fall netting (per pond x season).
#   Growth: monthly S->L and L->A rates; fraction of recruits entering as L
#     (pLg) is a logistic function of dens with the observed netting size split
#     as its likelihood.
#   Untagged pool: annual S / L / A counts with shared rates and known removals.
#   Netting join: previously tagged among handled ~ Binomial(n, T/(T+U)) by stage.

source("LabFunctions.R")
packages(dplyr)
packages(nimble)
packages(coda)
packages(MCMCvis)
packages(parallel)

# Pond set and calendar select the data build and keep fits in separate files
if (!exists("PONDS")) PONDS <- c(1L, 2L, 4L, 5L)  # 1 = YCB; c(1, 2, 4, 5) = all XYTE (YCB, IP1, IP3, IP4)
if (!exists("CALENDAR")) CALENDAR <- "m12"         # "m7" (Oct-Apr primaries) or "m12" (all months)
# DEFS_FILE lets a diagnostic run source an alternate defs script (e.g. one
# that overrides make_3s_inits() with anchored starting values) without
# touching the default model/YCB3S_defs.R used by every other fit.
if (!exists("DEFS_FILE")) DEFS_FILE <- "model/YCB3S_defs.R"
PONDS      <- sort(as.integer(PONDS))
pond_tag   <- if (identical(PONDS, 1L)) "" else
  if (identical(PONDS, c(1L, 2L, 4L, 5L))) "_xyte" else paste0("_p", paste(PONDS, collapse = ""))
sfx        <- paste0(pond_tag, if (CALENDAR == "m12") "_m12" else "")
if (!exists("DATA_FILE"))  DATA_FILE  <- paste0("data/YCB3S_data", sfx, ".RData")
if (!exists("OUT_FILE"))   OUT_FILE   <- paste0("data/YCB3S_MCMC", sfx, ".RData")
if (!exists("TRACE_FILE")) TRACE_FILE <- paste0("model/YCB3S_traceplots", sfx, ".pdf")
load(DATA_FILE)
stopifnot(identical(BUILD_INFO$CALENDAR, CALENDAR),
          identical(sort(as.integer(BUILD_INFO$MODEL_PONDS)), PONDS))

# Refuse to fit stale inputs (see AGENTS.md, Known Code Issues)
cur_mtime <- file.info(BUILD_INFO$inputs)$mtime
stale <- which(!is.na(cur_mtime) & cur_mtime > BUILD_INFO$built)
if (length(stale) > 0) {
  stop("Inputs newer than ", DATA_FILE, ": ", paste(BUILD_INFO$inputs[stale], collapse = ", "),
       ". Re-run model/YCB3S_data.R first.")
}
source(DEFS_FILE)

# ---------------------------------------------------------------------------
# Run settings
# ---------------------------------------------------------------------------
N_CHAINS  <- 3L
N_ITER    <- 30000L
N_BURNIN  <- 10000L
N_THIN    <- 10L
SEEDS     <- c(6173L, 2948L, 8805L)
PARALLEL  <- TRUE

inp <- make_3s_inputs()
cat("Ponds:", paste(PondTable$Backwater, collapse = ", "),
    "| fish:", nrow(inp$Fish), "| unique histories:", nrow(inp$Hist),
    "| netting rows:", inp$constants$NEV, "| size-split rows:", inp$constants$NSPLIT, "\n")

# ---------------------------------------------------------------------------
# Run chains
# ---------------------------------------------------------------------------
t0 <- Sys.time()
if (PARALLEL && N_CHAINS > 1) {
  cl <- makeCluster(N_CHAINS)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("inp", "N_ITER", "N_BURNIN", "N_THIN", "DATA_FILE", "DEFS_FILE"))
  clusterEvalQ(cl, {
    library(dplyr); library(nimble); library(coda)
    load(DATA_FILE)
    source(DEFS_FILE)
  })
  chain_list <- parLapply(cl, SEEDS[seq_len(N_CHAINS)], function(s)
    run_3s_chain(inp, s, N_ITER, N_BURNIN, N_THIN, progress = FALSE))
  stopCluster(cl); on.exit()
} else {
  chain_list <- lapply(SEEDS[seq_len(N_CHAINS)], function(s)
    run_3s_chain(inp, s, N_ITER, N_BURNIN, N_THIN))
}
rd3_samples <- as.mcmc.list(chain_list)
run_time <- difftime(Sys.time(), t0, units = "mins")
cat("Total wall time (compile + MCMC):", round(run_time, 1), "min\n")

# ---------------------------------------------------------------------------
# Diagnostics and save
# ---------------------------------------------------------------------------
core_pars <- c("mu_phi", "sigma_phi", "dJ", "dS", "b_post", "sigma_J", "lpsiSL", "lpsiLA",
               "a0", "a1", "lp", "g_eff", "mu_r", "sigma_r", "nu", "U0",
               "dE", "eps_A_IP1_FY2022")
rd3_summary <- MCMCsummary(rd3_samples, params = core_pars, round = 3)
print(rd3_summary)

rhat3_all <- if (N_CHAINS > 1) gelman.diag(rd3_samples, multivariate = FALSE)$psrf else NULL
if (!is.null(rhat3_all)) {
  cat("\nMax R-hat (all monitored):", round(max(rhat3_all[, 1], na.rm = TRUE), 3), "\n")
  cat("Parameters with R-hat > 1.05:", sum(rhat3_all[, 1] > 1.05, na.rm = TRUE), "\n")
}

MCMCtrace(rd3_samples, params = core_pars, pdf = TRUE, filename = TRACE_FILE)

rd3_Fish <- inp$Fish
rd3_Hist <- inp$Hist
rd3_settings <- list(N_CHAINS = N_CHAINS, N_ITER = N_ITER, N_BURNIN = N_BURNIN,
                     N_THIN = N_THIN, SEEDS = SEEDS, run_time = run_time,
                     MODEL_PONDS = MODEL_PONDS, StageBreaks = StageBreaks, BUILD_INFO = BUILD_INFO,
                     CALENDAR = CALENDAR, PONDS = PONDS, DATA_FILE = DATA_FILE, DEFS_FILE = DEFS_FILE,
                     variant = paste(if (CALENDAR == "m12") "12-month calendar (May-Sep months are primaries);"
                                     else "Oct-Apr calendar (7 primaries/season);",
                                     "shared juvenile detection: lp[j, 2] <- lp[j, 1];",
                                     "b_postJ tied to b_post;",
                                     "season-specific juvenile deviation eps_J[j, y] on lphiL,",
                                     "centered: eps_J = eps_raw - mean(eps_raw[j, ]), eps_raw ~ N(0, sigma_J[sp]);",
                                     "no direct density effect on juvenile survival (cB fixed at 0; density acts via a1 -> pLg);",
                                     "extra IP1 x FY2022 adult-survival deviation eps_A_IP1_FY2022 ~ N(0, sd=0.5),",
                                     "added to lphiA_eff[IP1, FY2022] only (all other pond-seasons unchanged),",
                                     "dE[IP1] die-off offset still applied on top for event months;",
                                     "data built with scanner-week screening (DroppedWeeks) and real pond acreage"))
save(rd3_samples, rd3_summary, rhat3_all, rd3_code, rd3_Fish, rd3_Hist, PondTable, Netting,
     RecruitIndex, SplitObs, rd3_settings, file = OUT_FILE)
cat("Saved", OUT_FILE, "\n")
