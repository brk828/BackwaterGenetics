# Second-generation diagnostic: isolate every bootstrap replicate's crm() fit
# in its own callr subprocess. If crm() hard-crashes (segfault) inside that
# subprocess, callr::r() surfaces it to the parent as a normal, catchable R
# error instead of taking the whole job down -- this is the fix for the
# A6D7E7F2 crash (which killed the whole job at replicate b=3, ~5 sec in,
# with b=1 and b=2 also failing to produce a valid fit).
# TEST RUN: B_TEST small, to confirm the isolation approach survives past
# b=3 before committing to the full B=300 run.

source("LabFunctions.R")
packages(dplyr); packages(callr)

load("data/GIEL_ClassicCJS_RobustDesign.RData")

FY_MIN <- 2017L; FY_MAX <- 2026L
T_OCC <- 70

build_get_p <- function(preds_p) {
  p_covs <- setdiff(names(preds_p), c("occ", "estimate", "se", "lcl", "ucl"))
  moy_labels <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
  if ("moy" %in% p_covs) {
    good <- preds_p$se > 0
    p_by_moy <- setNames(preds_p$estimate[good], as.character(preds_p$moy[good]))
    p_by_moy <- p_by_moy[!duplicated(names(p_by_moy))]
    function(occ) {
      local_pos <- ((occ - 1) %% 7) + 1
      val <- p_by_moy[[moy_labels[local_pos]]]
      if (is.null(val) || is.na(val)) 0 else val
    }
  } else {
    good <- preds_p$se > 0
    val0 <- preds_p$estimate[good][1]
    function(occ) if (is.na(val0)) 0 else val0
  }
}
get_p <- build_get_p(preds_ip5$p)

build_get_phi <- function(preds_Phi) {
  phi_covs <- setdiff(names(preds_Phi), c("occ", "estimate", "se", "lcl", "ucl", "sizeclass"))
  function(occ, sizeclass_i) {
    if (occ %% 7 != 0) return(1)
    row_sc <- preds_Phi$sizeclass == sizeclass_i & preds_Phi$se > 0
    if (length(phi_covs) == 0) {
      val <- preds_Phi$estimate[row_sc]
      if (length(val) == 0) val <- preds_Phi$estimate[preds_Phi$sizeclass == "adult" & preds_Phi$se > 0]
      return(val[1])
    }
    g <- as.character(FY_MIN + occ %/% 7 - 1)
    if ("gap" %in% phi_covs) {
      key_match <- row_sc & as.character(preds_Phi$gap) == g
      val <- preds_Phi$estimate[key_match]
      if (length(val) == 0) val <- preds_Phi$estimate[preds_Phi$sizeclass == "adult" & preds_Phi$se > 0 & as.character(preds_Phi$gap) == g]
      return(val[1])
    }
    if ("dieoff" %in% phi_covs) {
      dstat <- ifelse(g == "2020", "event", "normal")
      key_match <- row_sc & as.character(preds_Phi$dieoff) == dstat
      val <- preds_Phi$estimate[key_match]
      if (length(val) == 0) val <- preds_Phi$estimate[preds_Phi$sizeclass == "adult" & preds_Phi$se > 0 & as.character(preds_Phi$dieoff) == dstat]
      return(val[1])
    }
    NA
  }
}
get_phi <- build_get_phi(preds_ip5$Phi)

sim_base <- giel_ip5 |>
  transmute(entry_fy, stage0, entry_occ, cohort = paste(entry_fy, stage0))

simulate_ch <- function(sim_base, T_OCC) {
  n <- nrow(sim_base)
  y <- matrix(0L, n, T_OCC)
  for (i in seq_len(n)) y[i, sim_base$entry_occ[i]] <- 1L
  alive <- rep(TRUE, n)
  for (t in seq_len(T_OCC - 1)) {
    elig <- which(sim_base$entry_occ <= t & alive)
    if (!length(elig)) next
    sizeclass_i <- ifelse(sim_base$stage0[elig] == "young" & sim_base$entry_occ[elig] == t,
                           "young", "adult")
    phi <- sapply(seq_along(elig), function(k) get_phi(t, sizeclass_i[k]))
    survive <- rbinom(length(elig), 1, phi) == 1L
    alive[elig[!survive]] <- FALSE
    survivors <- elig[survive]
    if (length(survivors)) {
      pdet <- get_p(t + 1)
      y[survivors, t + 1] <- rbinom(length(survivors), 1, pdet)
    }
  }
  apply(y, 1, paste, collapse = ",")
}

top_formula_ip5 <- model_specs_ip5[[top_name_ip5]]

# --- The subprocess-isolated fit function -----------------------------------
fit_one_boot <- function(ch_sim, stage0, top_formula_ip5, top_ip5, T_OCC, FY_MIN, FY_MAX) {
  library(marked)
  df <- data.frame(ch = ch_sim, stage0 = stage0)
  dp_b <- process.data(df, model = "CJS", begin.time = 1, time.intervals = rep(1, T_OCC - 1))
  ddl_b <- make.design.data(dp_b)

  entry_occ_b <- sapply(strsplit(dp_b$data$ch, ","), function(x) which(x == "1")[1])
  idmap_b <- data.frame(id = as.character(dp_b$data$id), entry_occ = entry_occ_b)
  ddl_b$Phi <- merge(ddl_b$Phi, idmap_b, by = "id", all.x = TRUE, sort = FALSE)
  ddl_b$Phi <- ddl_b$Phi[order(ddl_b$Phi$id, ddl_b$Phi$occ), ]
  occ_num_b <- as.numeric(as.character(ddl_b$Phi$occ))
  ddl_b$Phi$sizeclass <- factor(
    ifelse(ddl_b$Phi$stage0 == "young" & occ_num_b == ddl_b$Phi$entry_occ, "young", "adult"),
    levels = c("adult", "young"))
  ddl_b$Phi$fix <- ifelse(occ_num_b %% 7 == 0, NA, 1)
  gap_index_b <- ifelse(occ_num_b %% 7 == 0, occ_num_b %/% 7, NA)
  gap_fy_b <- FY_MIN + gap_index_b - 1
  ddl_b$Phi$gap <- factor(ifelse(is.na(gap_fy_b), "none", as.character(gap_fy_b)),
                            levels = c("none", as.character(FY_MIN:(FY_MAX - 1))))
  ddl_b$Phi$dieoff <- factor(ifelse(!is.na(gap_fy_b) & gap_fy_b == 2020, "event", "normal"),
                               levels = c("normal", "event"))
  occ_num_pb <- as.numeric(as.character(ddl_b$p$occ))
  local_pos_pb <- ((occ_num_pb - 1) %% 7) + 1
  ddl_b$p$moy <- factor(local_pos_pb, levels = 1:7,
                          labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))

  fit_b <- crm(dp_b, ddl_b, model.parameters = top_formula_ip5, initial = top_ip5,
               hessian = FALSE, accumulate = TRUE)
  list(neg2lnl = fit_b$results$neg2lnl, convergence = fit_b$results$convergence)
}

B_TEST <- 15
set.seed(2854)
seeds <- sample.int(1e6, B_TEST)
boot_results <- vector("list", B_TEST)
CHECKPOINT_FILE <- "data/GIEL_ClassicCJS_RobustDesign_diag2_checkpoint.RData"

t0 <- Sys.time()
for (b in seq_len(B_TEST)) {
  set.seed(seeds[b])
  ch_sim <- simulate_ch(sim_base, T_OCC)
  sat_ll_b <- sum(tapply(ch_sim, sim_base$cohort, function(x) {
    tab <- table(x); n <- sum(tab); sum(tab * log(tab / n))
  }))

  res <- tryCatch(
    callr::r(fit_one_boot,
             args = list(ch_sim = ch_sim, stage0 = giel_ip5$stage0,
                          top_formula_ip5 = top_formula_ip5, top_ip5 = top_ip5,
                          T_OCC = T_OCC, FY_MIN = FY_MIN, FY_MAX = FY_MAX),
             package = TRUE, timeout = 180),
    error = function(e) { cat("  b=", b, "subprocess FAILED:", conditionMessage(e), "\n"); NULL })

  if (is.null(res)) {
    boot_results[[b]] <- list(deviance = NA, convergence = NA, status = "crashed_or_error")
  } else {
    boot_results[[b]] <- list(deviance = res$neg2lnl - (-2 * sat_ll_b),
                               convergence = res$convergence, status = "ok")
  }
  save(boot_results, b, file = CHECKPOINT_FILE)
  cat("bootstrap", b, "of", B_TEST, "-", boot_results[[b]]$status, "-",
      round(as.numeric(Sys.time() - t0, units = "mins"), 2), "min\n")
}

cat("\n=== DIAG2 TEST COMPLETE:", B_TEST, "replicates attempted ===\n")
n_ok <- sum(sapply(boot_results, function(x) x$status == "ok"))
n_bad <- sum(sapply(boot_results, function(x) x$status == "crashed_or_error"))
cat("OK:", n_ok, " crashed/error:", n_bad, "\n")
save(boot_results, file = "data/GIEL_ClassicCJS_RobustDesign_diag2_final.RData")
