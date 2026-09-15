## Step 5 of the growth-model sensitivity plan: validate the GIEL size-tied
## maturation hazard mechanism on Track F1 (Fabens, IPCA-only) first.
## Assumes data/BWRobustDesign_data.RData was already (re)built with
## MATURATION_TRACK <- "F1" in model/BWRobustDesign_data.R.
RUN_PONDS <- c(3L, 6L, 7L)
OUT_FILE   <- "data/BWRobustDesign_MCMC_giel_F1.RData"
TRACE_FILE <- "model/BWRobustDesign_traceplots_giel_F1.pdf"
source("model/BWRobustDesign_NIMBLE.R")
