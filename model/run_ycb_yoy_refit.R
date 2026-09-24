# Runner for the YCB-only m12 refit on YOY-classified RecruitIndex/SplitObs
# data (2026-09-23; see AGENTS.md "Young-of-year classification" section).
# Launch as a background job:
#   rstudioapi::jobRunScript("model/run_ycb_yoy_refit.R", workingDir = getwd(), importEnv = FALSE)
PONDS    <- 1L
CALENDAR <- "m12"
source("model/YCB3S_NIMBLE.R")
