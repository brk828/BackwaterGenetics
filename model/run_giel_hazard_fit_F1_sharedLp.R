## GIEL sensitivity refit, Track F1 (Fabens, IPCA-only) + shared IP2
## juvenile/adult detection intercept (2026-09-15; see AGENTS.md "IP2 (Pond 2)
## adult-survival investigation" and the sharedLp mechanism added to
## model/BWRobustDesign_defs.R).
## Assumes data/BWRobustDesign_data.RData was already (re)built with
## MATURATION_TRACK <- "F1" in model/BWRobustDesign_data.R.
RUN_PONDS <- c(3L, 6L, 7L)
OUT_FILE   <- "data/BWRobustDesign_MCMC_giel_F1_sharedLp.RData"
TRACE_FILE <- "model/BWRobustDesign_traceplots_giel_F1_sharedLp.pdf"
source("model/BWRobustDesign_NIMBLE.R")
