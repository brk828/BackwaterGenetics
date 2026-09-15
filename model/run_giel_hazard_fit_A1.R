## Step 5 of the growth-model sensitivity plan: GIEL size-tied maturation
## hazard, Track A1 (age-based VBGF, IPCA-only).
## Assumes data/BWRobustDesign_data.RData was already (re)built with
## MATURATION_TRACK <- "A1" in model/BWRobustDesign_data.R.
RUN_PONDS <- c(3L, 6L, 7L)
OUT_FILE   <- "data/BWRobustDesign_MCMC_giel_A1.RData"
TRACE_FILE <- "model/BWRobustDesign_traceplots_giel_A1.pdf"
source("model/BWRobustDesign_NIMBLE.R")
