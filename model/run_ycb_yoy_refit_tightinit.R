# Launch helper for the anchored-inits diagnostic refit (2026-09-23).
# See model/YCB3S_defs_tightinit.R for rationale.
# Run as a background job:
#   rstudioapi::jobRunScript("model/run_ycb_yoy_refit_tightinit.R", workingDir = getwd(), importEnv = FALSE)
PONDS      <- 1L
CALENDAR   <- "m12"
DEFS_FILE  <- "model/YCB3S_defs_tightinit.R"
OUT_FILE   <- "data/YCB3S_MCMC_m12_tightinit.RData"
TRACE_FILE <- "model/YCB3S_traceplots_m12_tightinit.pdf"
source("model/YCB3S_NIMBLE.R")
