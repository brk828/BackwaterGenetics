# Runner for the YCB-only m12 refit on the bulk-untagged-YOY-aware
# RecruitIndex/SplitObs data (2026-09-24; see AGENTS.md "Bulk-Harvested
# Untagged YOY vs. Individually-Tagged Recruits" section).
#
# Reuses model/YCB3S_defs_tightinit.R (the anchored-inits variant that fixed
# job BAD9373E's chain-3 mixing failure) since the anchor only touches
# species-1 survival/maturation inits, not the size-split (a0/a1) submodel
# affected by this data change.
#
# Launch as a background job:
#   rstudioapi::jobRunScript("model/run_ycb_bulkfix_refit.R", workingDir = getwd(), importEnv = FALSE)
PONDS      <- 1L
CALENDAR   <- "m12"
DEFS_FILE  <- "model/YCB3S_defs_tightinit.R"
OUT_FILE   <- "data/YCB3S_MCMC_m12_bulkfix.RData"
TRACE_FILE <- "model/YCB3S_traceplots_m12_bulkfix.pdf"
source("model/YCB3S_NIMBLE.R")
