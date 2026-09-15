## Step 5 of the growth-model sensitivity plan: GIEL size-tied maturation
## hazard, Track F2 (Fabens increment, IPCA + Cibola).
## Assumes data/BWRobustDesign_data.RData was already (re)built with
## MATURATION_TRACK <- "F2" in model/BWRobustDesign_data.R.
RUN_PONDS <- c(3L, 6L, 7L)
OUT_FILE   <- "data/BWRobustDesign_MCMC_giel_F2.RData"
TRACE_FILE <- "model/BWRobustDesign_traceplots_giel_F2.pdf"
source("model/BWRobustDesign_NIMBLE.R")
