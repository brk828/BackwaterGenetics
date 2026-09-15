# Definitions for the joint stage-structured robust-design model:
# input builder, nimbleFunctions, model code, and initial values.
# Sourced by model/BWRobustDesign_NIMBLE.R (and by each parallel chain worker).
# Requires objects from data/BWRobustDesign_data.RData in the calling environment.

# ---------------------------------------------------------------------------
# Build NIMBLE inputs for a subset of ponds; identical fish histories are
# collapsed to one weighted history
# ---------------------------------------------------------------------------
make_rd_inputs <- function(ponds) {
  keep <- which(Fish$pond %in% ponds)
  FishSub <- Fish[keep, ] |>
    mutate(pond_new = match(pond, ponds),
           count_from = entry_t + as.integer(origin == 2L),
           row0 = keep)
  # TL_entry_bin / hasTL are added to the collapsing key (2026-09-14, GIEL
  # size-tied maturation hazard) so that two fish with otherwise-identical
  # detection histories but different TL bins -- and therefore different
  # size-tied maturation trajectories -- are NOT collapsed together. hasTL
  # is included too so unknown-TL fish (which all carry the same dummy
  # TL_entry_bin = 1, see model/BWRobustDesign_data.R) still collapse with
  # each other regardless of that dummy value.
  hist_key <- apply(cbind(FishSub$pond_new, FishSub$entry_t, FishSub$end_t, FishSub$stage0,
                          FishSub$origin, FishSub$removed, FishSub$count_from,
                          FishSub$hasTL, FishSub$TL_entry_bin,
                          y_mat[keep, ], K_mat[keep, ], stobs_mat[keep, ], post_mat[keep, ]),
                    1, paste, collapse = ",")
  FishSub$hist <- match(hist_key, unique(hist_key))
  Hist <- FishSub |>
    group_by(hist) |>
    summarise(pond_new = first(pond_new), entry_t = first(entry_t), end_t = first(end_t),
              stage0 = first(stage0), origin = first(origin), removed = first(removed),
              count_from = first(count_from), hasTL = first(hasTL),
              TL_entry_bin = first(TL_entry_bin),
              w = n(), row0 = first(row0), .groups = "drop") |>
    arrange(pond_new, entry_t, hist)
  NP <- length(ponds); TT <- T_PRIM
  pond_off <- c(0L, cumsum(tabulate(Hist$pond_new, nbins = NP)))
  flat <- function(m) as.vector(t(m[Hist$row0, , drop = FALSE]))

  NetSub <- Netting |> filter(pond %in% ponds) |> mutate(pond_new = match(pond, ponds))

  # Remap species to 1..n present among the chosen ponds, so a single-species
  # subset (e.g. the 3 GIEL ponds) doesn't carry a dead, prior-only species-1
  # slot the way the all-XYTE YCB3S builds do for species 2.
  sp_raw <- PondTable$sp[ponds]
  usp <- sort(unique(sp_raw))
  NSP <- length(usp)
  sp  <- match(sp_raw, usp)

  constants <- list(
    NPOND = NP, NSP = NSP, T = TT, NS = N_SEASONS, NMOY = N_MOY, NEV = nrow(NetSub),
    sp = sp, season = Primary$season, moy = Primary$moy,
    gap = Primary$gap, tApr = Primary$t[Primary$moy == N_MOY],
    pond_off = pond_off,
    K_vec = flat(K_mat), stobs_vec = flat(stobs_mat), post_vec = flat(post_mat),
    y_vec = flat(y_mat), w = Hist$w,
    entry_t = Hist$entry_t, end_t = Hist$end_t, removed = as.integer(Hist$removed),
    stage0 = Hist$stage0, origin = Hist$origin, count_from = Hist$count_from,
    hasTL = as.integer(Hist$hasTL), tlbin = as.integer(Hist$TL_entry_bin),
    eff = eff_std[ponds, , drop = FALSE],
    Evt = RD_constants$Evt[ponds, , drop = FALSE], nEvt = RD_constants$nEvt[ponds, , drop = FALSE],
    ev_pond = NetSub$pond_new, ev_t = NetSub$t_e, ev_season = NetSub$season,
    ev_stage = NetSub$stage, ev_n = NetSub$n, ev_remU_before = NetSub$remU_before,
    remU = RD_constants$remU[ponds, , , drop = FALSE],
    # GIEL size-tied maturation hazard lookup (same table for every pond
    # subset -- not pond-specific, so no [ponds, ] slicing needed)
    psi_lookup = RD_constants$psi_lookup, N_TLBIN = RD_constants$N_TLBIN,
    DELTA_T_MAX = RD_constants$DELTA_T_MAX
  )
  constants$Uk <- as.integer(PondTable$pond[ponds] != 1)   # YCB initial pool unknown
  # Shared juvenile/adult detection intercept at IP2 (2026-09-15; mirrors the
  # XYTE three-stage model's L-stage fix, see AGENTS.md "IP2 (Pond 2)
  # adult-survival investigation"): IP2 is the only GIEL pond where juvenile
  # detection is much worse than adult detection, a low-detection-as-a-
  # survival-sink failure mode that inflates apparent juvenile survival
  # (dJ[GIEL] > 0). sharedLp[j] == 1 forces lp[j,2] <- lp[j,1] at IP2; 0
  # everywhere else (a no-op for every other pond/fit, including XYTE-only
  # subsets of the joint model, which never include IP2).
  constants$sharedLp <- as.integer(PondTable$Backwater[ponds] == "IP2")
  data <- list(
    zero = rep(0, NP),
    ev_m = NetSub$m,
    okU  = array(1L, dim = c(NP, N_SEASONS, 2))
  )
  list(constants = constants, data = data, Fish = FishSub, Hist = Hist,
       Netting = NetSub, ponds = ponds)
}

# ---------------------------------------------------------------------------
# nimbleFunctions
# ---------------------------------------------------------------------------
# Forward-backward over all (weighted) histories in one pond. Returns a T x 3
# matrix: columns 1-2 = expected number of tagged fish alive in stage J / A at
# each primary period (smoothed, counted from count_from), [1, 3] = total
# log-likelihood of the detection histories.
# stobs codes: 0 none, 1 observed juvenile, 2 observed adult, 3 known alive.
rdPond <- nimbleFunction(
  run = function(y = double(1), K = double(1), stobs = double(1), post = double(1),
                 w = double(1), entry = double(1), endt = double(1), removed = double(1),
                 stage0 = double(1), origin = double(1), cfrom = double(1),
                 nfish = double(), T = double(),
                 lphiJ = double(1), lphiA = double(1), season = double(1), gap = double(1),
                 pJ = double(1), pA = double(1), bpost = double(1), lpsi = double(),
                 evt = double(1), dE = double(),
                 tlbin = double(1), hasTL = double(1), psi_lookup = double(2),
                 lpsi_adj = double(), dtmax = double()) {
    returnType(double(2))
    out  <- matrix(0, nrow = T, ncol = 3)
    phJt <- matrix(0, nrow = T, ncol = 6)   # column = post * 3 + origin
    phAt <- matrix(0, nrow = T, ncol = 6)
    pst  <- numeric(T, value = 0)
    psi_m <- ilogit(lpsi)
    for (t in 1:T) {
      g <- gap[t]
      pst[t] <- 1 - (1 - psi_m)^g
      # Known die-off event: fixed logit-survival offset dE applied over the
      # whole gap when evt[t] == 1. Assumes gap == 1 (exact) or, for the
      # 6-month Apr->Oct gap, that dE should apply to the whole gap when any
      # part of it overlaps a recorded event (see BWDieOffEvents.csv / Evt
      # construction in model/BWRobustDesign_data.R) - the same approximation
      # used for the XYTE three-stage model before its m12 calendar.
      dEt <- dE * evt[t]
      for (o in 1:3) {
        phJt[t, o]     <- ilogit(lphiJ[season[t]] + dEt)^g
        phAt[t, o]     <- ilogit(lphiA[season[t]] + dEt)^g
        phJt[t, 3 + o] <- ilogit(lphiJ[season[t]] + bpost[o] + dEt)^g
        phAt[t, 3 + o] <- ilogit(lphiA[season[t]] + bpost[o] + dEt)^g
      }
    }
    a  <- matrix(0, nrow = T, ncol = 3)
    em <- matrix(0, nrow = T, ncol = 3)
    ll <- 0
    for (i in 1:nfish) {
      off <- (i - 1) * T
      t1 <- entry[i]
      t2 <- endt[i]
      o  <- origin[i]
      a1 <- 0; a2 <- 0; a3 <- 0
      if (stage0[i] == 1) a1 <- 1 else a2 <- 1
      lli <- 0
      for (t in t1:t2) {
        yy <- y[off + t]; kk <- K[off + t]; so <- stobs[off + t]
        e1 <- 1; e2 <- 1; e3 <- 1
        if (kk > 0) {
          e1 <- dbinom(yy, kk, pJ[t])
          e2 <- dbinom(yy, kk, pA[t])
          if (yy > 0) e3 <- 0
        }
        if (so == 1) { e2 <- 0; e3 <- 0 }
        if (so == 2) { e1 <- 0; e3 <- 0 }
        if (so == 3) { e3 <- 0 }                 # known alive, stage unknown
        em[t, 1] <- e1; em[t, 2] <- e2; em[t, 3] <- e3
        a1 <- a1 * e1; a2 <- a2 * e2; a3 <- a3 * e3
        s <- a1 + a2 + a3
        if (s <= 0) s <- 1e-300
        lli <- lli + log(s)
        a1 <- a1 / s; a2 <- a2 / s; a3 <- a3 / s
        a[t, 1] <- a1; a[t, 2] <- a2; a[t, 3] <- a3
        if (t < t2) {
          idx <- post[off + t] * 3 + o
          pJs <- phJt[t, idx]; pAs <- phAt[t, idx]
          # GIEL size-tied maturation hazard (2026-09-14): for fish with a
          # known TL bin at entry (hasTL[i] == 1), replace the constant
          # per-species monthly hazard (pst[t], from lpsi[sp]) with a
          # cumulative product over the gap[t] intervening months of the
          # size- and time-since-entry-specific hazard psi_lookup[bin, dt],
          # logit-shifted by lpsi_adj[sp]. dt = months elapsed since entry
          # (dt = 1 for the month immediately after entry); beyond dtmax
          # months, maturation is treated as certain that month. Fish with
          # hasTL[i] == 0 (XYTE, or GIEL with imputed/missing TL) keep the
          # original scalar-lpsi pst[t] path unchanged.
          if (hasTL[i] == 1) {
            g <- gap[t]
            dt0 <- t - entry[i]
            surv <- 1
            for (kk2 in 1:g) {
              dtk <- dt0 + kk2
              if (dtk > dtmax) {
                psi_k <- 1
              } else {
                psi_k <- ilogit(logit(psi_lookup[tlbin[i], dtk]) + lpsi_adj)
              }
              surv <- surv * (1 - psi_k)
            }
            ps <- 1 - surv
          } else {
            ps <- pst[t]
          }
          n1 <- a1 * pJs * (1 - ps)
          n2 <- a1 * pJs * ps + a2 * pAs
          n3 <- a3 + a1 * (1 - pJs) + a2 * (1 - pAs)
          a1 <- n1; a2 <- n2; a3 <- n3
        }
      }
      if (removed[i] == 1) {
        s <- a1 + a2
        if (s <= 0) s <- 1e-300
        lli <- lli + log(s)
      }
      ll <- ll + w[i] * lli
      # backward pass for smoothed state probabilities
      b1 <- 1; b2 <- 1; b3 <- 1
      if (removed[i] == 1) b3 <- 0
      for (tt in 1:(t2 - t1 + 1)) {
        t <- t2 - tt + 1
        if (t >= cfrom[i]) {
          s1 <- a[t, 1] * b1; s2 <- a[t, 2] * b2; s3 <- a[t, 3] * b3
          s <- s1 + s2 + s3
          if (s > 0) {
            out[t, 1] <- out[t, 1] + w[i] * s1 / s
            out[t, 2] <- out[t, 2] + w[i] * s2 / s
          }
        }
        if (t > t1) {
          eb1 <- em[t, 1] * b1; eb2 <- em[t, 2] * b2; eb3 <- em[t, 3] * b3
          idx <- post[off + t - 1] * 3 + o
          pJs <- phJt[t - 1, idx]; pAs <- phAt[t - 1, idx]
          # Same size-tied override as the forward pass above, for the
          # transition starting at t - 1.
          if (hasTL[i] == 1) {
            g <- gap[t - 1]
            dt0 <- (t - 1) - entry[i]
            surv <- 1
            for (kk2 in 1:g) {
              dtk <- dt0 + kk2
              if (dtk > dtmax) {
                psi_k <- 1
              } else {
                psi_k <- ilogit(logit(psi_lookup[tlbin[i], dtk]) + lpsi_adj)
              }
              surv <- surv * (1 - psi_k)
            }
            ps <- 1 - surv
          } else {
            ps <- pst[t - 1]
          }
          nb1 <- pJs * (1 - ps) * eb1 + pJs * ps * eb2 + (1 - pJs) * eb3
          nb2 <- pAs * eb2 + (1 - pAs) * eb3
          nb3 <- eb3
          sb <- nb1 + nb2 + nb3
          if (sb > 0) { b1 <- nb1 / sb; b2 <- nb2 / sb; b3 <- nb3 / sb }
        }
      }
    }
    out[1, 3] <- ll
    return(out)
  })

# Zeros-trick distribution: log-density equals a supplied log-likelihood
dLLtrick <- nimbleFunction(
  run = function(x = double(), ll = double(), log = integer(0, default = 0)) {
    returnType(double())
    if (log) return(ll) else return(exp(ll))
  })
rLLtrick <- nimbleFunction(
  run = function(n = integer(), ll = double()) {
    returnType(double())
    return(0)
  })
registerDistributions(list(
  dLLtrick = list(BUGSdist = "dLLtrick(ll)",
                  types = c("value = double()", "ll = double()"))
), verbose = FALSE)

# ---------------------------------------------------------------------------
# Model code
# ---------------------------------------------------------------------------
rd_code <- nimbleCode({
  # Species-level parameters --------------------------------------------------
  for (s in 1:NSP) {
    mu_phi[s]    ~ dnorm(3, sd = 1.5)            # mean logit monthly adult survival
    sigma_phi[s] ~ T(dnorm(0, sd = 1), 0, )      # pond x season SD
    dJ[s]        ~ dnorm(-0.5, sd = 1)           # juvenile offset to logit survival
    for (o in 1:2) {
      b_post[s, o] ~ T(dnorm(-1, sd = 1), , 0)   # post-release offset: stocked, netting-tagged
    }
    b_post[s, 3] <- 0                            # established fish
    lpsi[s]  ~ dnorm(-2.4, sd = 1)               # logit monthly maturation J -> A (scalar fallback path)
    lpsi_adj[s] ~ dnorm(0, sd = 1)               # logit offset on the GIEL size-tied psi_lookup hazard (see AGENTS.md); prior-only for species with no known-TL juveniles (XYTE)
    g_eff[s] ~ dnorm(0, sd = 1)                  # effect of standardized log scan hours
    d_moy[s, 1] <- 0
    for (mm in 2:NMOY) { d_moy[s, mm] ~ dnorm(0, sd = 1) }
    mu_r[s]    ~ dnorm(-1.2, sd = 1.5)           # mean log recruits (juvenile stage) per adult
    sigma_r[s] ~ T(dnorm(0, sd = 1), 0, )
    nu[s]      ~ dexp(1)                         # unaccounted adult appearances per pond-season
  }

  # Pond-level structure ------------------------------------------------------
  for (j in 1:NPOND) {
    for (y in 1:NS) {
      lphiA[j, y] ~ dnorm(mu_phi[sp[j]], sd = sigma_phi[sp[j]])   # pond-season adult survival
      lphiJ[j, y] <- lphiA[j, y] + dJ[sp[j]]
    }
    for (st in 1:2) { lp_raw[j, st] ~ dnorm(1, sd = 1.5) }        # pond x stage detection
    lp[j, 1] <- lp_raw[j, 1]
    # Shared juvenile/adult detection intercept at IP2 (sharedLp[j] == 1);
    # lp_raw[j, 2] is an unused, prior-only node at IP2 (kept so every pond
    # has the same node structure) and the free juvenile/adult detection
    # difference at every other pond.
    lp[j, 2] <- lp_raw[j, 1] * sharedLp[j] + lp_raw[j, 2] * (1 - sharedLp[j])
    # Known die-off (data/BWDieOffEvents.csv): fixed, pond-specific, survival-
    # reducing-only logit offset applied to J/A alike in event months
    # (Evt[j,t] == 1). Ponds with no recorded event carry an uninformative
    # (prior-only) dE.
    dE[j] ~ T(dnorm(0, sd = 3), , 0)
    for (t in 1:T) {
      logit(pJ[j, t]) <- lp[j, 1] + g_eff[sp[j]] * eff[j, t] + d_moy[sp[j], moy[t]]
      logit(pA[j, t]) <- lp[j, 2] + g_eff[sp[j]] * eff[j, t] + d_moy[sp[j], moy[t]]
    }

    # Tagged fish: forward-backward over all histories in pond j
    RD[j, 1:T, 1:3] <- rdPond(
      y = y_vec[(pond_off[j] * T + 1):(pond_off[j + 1] * T)],
      K = K_vec[(pond_off[j] * T + 1):(pond_off[j + 1] * T)],
      stobs = stobs_vec[(pond_off[j] * T + 1):(pond_off[j + 1] * T)],
      post = post_vec[(pond_off[j] * T + 1):(pond_off[j + 1] * T)],
      w = w[(pond_off[j] + 1):pond_off[j + 1]],
      entry = entry_t[(pond_off[j] + 1):pond_off[j + 1]],
      endt = end_t[(pond_off[j] + 1):pond_off[j + 1]],
      removed = removed[(pond_off[j] + 1):pond_off[j + 1]],
      stage0 = stage0[(pond_off[j] + 1):pond_off[j + 1]],
      origin = origin[(pond_off[j] + 1):pond_off[j + 1]],
      cfrom = count_from[(pond_off[j] + 1):pond_off[j + 1]],
      nfish = pond_off[j + 1] - pond_off[j], T = T,
      lphiJ = lphiJ[j, 1:NS], lphiA = lphiA[j, 1:NS], season = season[1:T], gap = gap[1:T],
      pJ = pJ[j, 1:T], pA = pA[j, 1:T], bpost = b_post[sp[j], 1:3], lpsi = lpsi[sp[j]],
      evt = Evt[j, 1:T], dE = dE[j],
      tlbin = tlbin[(pond_off[j] + 1):pond_off[j + 1]],
      hasTL = hasTL[(pond_off[j] + 1):pond_off[j + 1]],
      psi_lookup = psi_lookup[1:N_TLBIN, 1:DELTA_T_MAX],
      lpsi_adj = lpsi_adj[sp[j]], dtmax = DELTA_T_MAX)
    zero[j] ~ dLLtrick(RD[j, 1, 3])
    for (st in 1:2) {
      for (t in 1:T) { N_tag[j, t, st] <- RD[j, t, st] }
    }

    # Untagged pool (deterministic skeleton, stochastic recruitment) ----------
    for (st in 1:2) {
      U0[j, st] ~ dunif(0, 5000)                 # only used where the initial pool is unknown
      U[j, 1, st] <- (1 - Uk[j]) * U0[j, st]     # Uk = 1: pond known empty of untagged fish
    }
    for (y in 1:(NS - 1)) {
      Atot[j, y] <- N_tag[j, tApr[y], 2] + U[j, y, 2]  # adults present at spawning (April)
      lr[j, y] ~ dnorm(mu_r[sp[j]], sd = sigma_r[sp[j]])
      r[j, y] <- exp(lr[j, y])
      RecJ[j, y] <- r[j, y] * Atot[j, y]              # recruits reaching the juvenile stage
      UJafter[j, y] <- max(0, U[j, y, 1] - remU[j, y, 1])
      UAafter[j, y] <- max(0, U[j, y, 2] - remU[j, y, 2])
      # Untagged-pool annual survival gets the same event discount as the
      # tagged model: nEvt[j,y] of the 12 calendar months in the season use
      # the event-adjusted monthly rate (dE), the rest use the baseline rate
      phiUJ[j, y] <- ilogit(lphiJ[j, y])^(12 - nEvt[j, y]) * ilogit(lphiJ[j, y] + dE[j])^nEvt[j, y]
      phiUA[j, y] <- ilogit(lphiA[j, y])^(12 - nEvt[j, y]) * ilogit(lphiA[j, y] + dE[j])^nEvt[j, y]
      psiY[j, y]  <- 1 - (1 - ilogit(lpsi[sp[j]]))^12
      U[j, y + 1, 1] <- UJafter[j, y] * phiUJ[j, y] * (1 - psiY[j, y]) + RecJ[j, y]
      U[j, y + 1, 2] <- UAafter[j, y] * phiUA[j, y] + UJafter[j, y] * phiUJ[j, y] * psiY[j, y] + nu[sp[j]]
    }
    for (y in 1:NS) {
      for (st in 1:2) {
        okU[j, y, st] ~ dconstraint(U[j, y, st] >= remU[j, y, st])
      }
    }
  }

  # Netting join: previously tagged among handled fish, by stage -------------
  for (e in 1:NEV) {
    Tal[e] <- N_tag[ev_pond[e], ev_t[e], ev_stage[e]]
    Uav[e] <- max(0, U[ev_pond[e], ev_season[e], ev_stage[e]] - ev_remU_before[e])
    pm[e]  <- (Tal[e] + 1e-9) / (Tal[e] + Uav[e] + 2e-9)
    ev_m[e] ~ dbin(pm[e], ev_n[e])
  }
})

# ---------------------------------------------------------------------------
# Initial values (must satisfy the U >= remU constraints)
# ---------------------------------------------------------------------------
make_inits <- function(inp, seed) {
  set.seed(seed)
  cst <- inp$constants
  NP <- cst$NPOND; NS <- cst$NS; NSP <- cst$NSP
  remU <- cst$remU
  U0 <- matrix(runif(NP * 2, 100, 300), NP, 2)
  for (j in 1:NP) U0[j, ] <- pmax(U0[j, ], remU[j, 1, ] + 50)
  mu_phi <- rnorm(NSP, 3, 0.3)
  # Start recruitment high so the deterministic U skeleton satisfies U >= remU;
  # chains move down from there
  mu_r <- rnorm(NSP, log(3), 0.1)
  list(
    mu_phi = mu_phi, sigma_phi = runif(NSP, 0.2, 0.6), dJ = rnorm(NSP, -0.5, 0.2),
    b_post = cbind(matrix(runif(2 * NSP, -1.5, -0.3), NSP, 2), 0),
    lpsi = rnorm(NSP, -2.4, 0.2), lpsi_adj = rnorm(NSP, 0, 0.2), g_eff = rnorm(NSP, 0.3, 0.1),
    d_moy = cbind(0, matrix(rnorm(NSP * (cst$NMOY - 1), 0, 0.2), NSP)),
    mu_r = mu_r, sigma_r = runif(NSP, 0.3, 0.8), nu = runif(NSP, 4, 6),
    lp_raw = matrix(rnorm(NP * 2, 1, 0.3), NP, 2),
    lphiA = matrix(rnorm(NP * NS, mu_phi[cst$sp], 0.2), NP, NS),
    lr = matrix(rnorm(NP * (NS - 1), mu_r[cst$sp], 0.1), NP, NS - 1),
    U0 = U0,
    dE = -runif(NP, 0.2, 1)
  )
}

# Raise recruitment in any pond-season whose deterministic untagged pool
# violates U >= remU at the initial values (rare; keeps calculate() finite)
repair_inits <- function(model, max_tries = 200) {
  for (k in seq_len(max_tries)) {
    lp <- model$calculate()
    if (is.finite(lp)) return(invisible(TRUE))
    bad <- model$getNodeNames(stochOnly = TRUE)
    bad <- bad[!is.finite(sapply(bad, function(n) model$calculate(n)))]
    idx <- regmatches(bad, regexec("okU\\[(\\d+), (\\d+), (\\d+)\\]", bad))
    idx <- do.call(rbind, lapply(idx[lengths(idx) > 0], function(v) as.integer(v[2:4])))
    if (is.null(idx)) stop("Non-finite initial log-probability not due to okU: ",
                           paste(head(bad), collapse = ", "))
    sp <- model$getConstants()$sp
    for (rr in seq_len(nrow(idx))) {
      j <- idx[rr, 1]; y <- idx[rr, 2]; st <- idx[rr, 3]
      if (y == 1) {
        model$U0[j, st] <- model$U0[j, st] + 100
      } else if (st == 1) {
        model$lr[j, y - 1] <- model$lr[j, y - 1] + 1
      } else {
        model$nu[sp[j]] <- model$nu[sp[j]] + 2       # adult pool fed by nu and maturation
        model$lr[j, y - 1] <- model$lr[j, y - 1] + 0.5
      }
    }
  }
  stop("Could not repair initial values")
}

# ---------------------------------------------------------------------------
# Build, configure, compile and run one chain (used directly or in a worker)
# ---------------------------------------------------------------------------
rd_monitors <- c("mu_phi", "sigma_phi", "dJ", "dE", "b_post", "lpsi", "lpsi_adj", "lp", "lp_raw", "g_eff",
                 "d_moy", "mu_r", "sigma_r", "nu", "U0",
                 "lphiA", "r", "U", "RecJ", "N_tag", "Tal", "Atot")

# Block samplers group every parameter that forces a forward-backward pass so
# that each pond is re-evaluated about four times per iteration
rd_configure <- function(model) {
  conf <- configureMCMC(model, monitors = rd_monitors, enableWAIC = FALSE, print = FALSE)
  NMOY <- model$getConstants()$NMOY
  NP <- model$getConstants()$NPOND; NS <- model$getConstants()$NS
  NSP <- model$getConstants()$NSP
  for (s in 1:NSP) {
    tgt_phi <- c(paste0("dJ[", s, "]"), paste0("b_post[", s, ", 1]"),
                 paste0("b_post[", s, ", 2]"), paste0("lpsi[", s, "]"), paste0("lpsi_adj[", s, "]"))
    conf$removeSamplers(tgt_phi)
    conf$addSampler(target = tgt_phi, type = "RW_block", control = list(adaptInterval = 200))
    tgt_p <- c(paste0("g_eff[", s, "]"), paste0("d_moy[", s, ", ", 2:NMOY, "]"))
    conf$removeSamplers(tgt_p)
    conf$addSampler(target = tgt_p, type = "RW_block", control = list(adaptInterval = 200))
  }
  for (j in 1:NP) {
    tgt_e <- paste0("lphiA[", j, ", ", 1:NS, "]")
    conf$removeSamplers(tgt_e)
    conf$addSampler(target = tgt_e, type = "RW_block", control = list(adaptInterval = 200))
    tgt_lp <- paste0("lp_raw[", j, ", ", 1:2, "]")
    conf$removeSamplers(tgt_lp)
    conf$addSampler(target = tgt_lp, type = "RW_block", control = list(adaptInterval = 200))
  }
  conf
}

run_rd_chain <- function(inp, seed, niter, nburnin, thin, progress = TRUE) {
  inits <- make_inits(inp, seed)
  model <- nimbleModel(rd_code, constants = inp$constants, data = inp$data,
                       inits = inits, calculate = TRUE)
  repair_inits(model)
  conf   <- rd_configure(model)
  mcmc   <- buildMCMC(conf)
  cmodel <- compileNimble(model)
  cmcmc  <- compileNimble(mcmc, project = model)
  t0 <- Sys.time()
  samples <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin,
                     nchains = 1, setSeed = seed, samplesAsCodaMCMC = TRUE,
                     progressBar = progress)
  attr(samples, "run_time") <- difftime(Sys.time(), t0, units = "mins")
  samples
}
