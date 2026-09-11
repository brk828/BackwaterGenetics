# Definitions for the three-stage, density-dependent robust-design model
# (Phase 2 of .posit/assistant/plans/2026-09-04-1335-ycb-three-stage-...md).
# Sourced by model/YCB3S_NIMBLE.R and by each parallel chain worker.
# Requires objects from data/YCB3S_data.RData in the calling environment.
#
# States: 1 = S (small juvenile), 2 = L (large juvenile), 3 = A (adult), 4 = dead.
# stobs codes: 0 none, 1/2/3 observed S/L/A, 4 known alive (stage unknown).

# ---------------------------------------------------------------------------
# Collapse identical histories to weighted rows and assemble NIMBLE inputs
# ---------------------------------------------------------------------------
make_3s_inputs <- function() {
  FishSub <- Fish |>
    mutate(count_from = entry_t + as.integer(origin == 2L), row0 = row_number())
  hist_key <- apply(cbind(FishSub$pond, FishSub$entry_t, FishSub$end_t, FishSub$stage0,
                          FishSub$origin, FishSub$removed, FishSub$count_from,
                          y_mat, K_mat, stobs_mat, post_mat), 1, paste, collapse = ",")
  FishSub$hist <- match(hist_key, unique(hist_key))
  Hist <- FishSub |>
    group_by(hist) |>
    summarise(pond = first(pond), entry_t = first(entry_t), end_t = first(end_t),
              stage0 = first(stage0), origin = first(origin), removed = first(removed),
              count_from = first(count_from), w = n(), row0 = first(row0), .groups = "drop") |>
    arrange(pond, entry_t, hist)
  NP <- RD_constants$NPOND; TT <- T_PRIM
  pond_off <- c(0L, cumsum(tabulate(Hist$pond, nbins = NP)))
  flat <- function(m) as.vector(t(m[Hist$row0, , drop = FALSE]))

  constants <- RD_constants
  constants$entry_t <- Hist$entry_t; constants$end_t <- Hist$end_t
  constants$removed <- as.integer(Hist$removed); constants$stage0 <- Hist$stage0
  constants$origin <- Hist$origin; constants$count_from <- Hist$count_from
  constants$w <- Hist$w; constants$pond_off <- pond_off
  constants$K_vec <- flat(K_mat); constants$stobs_vec <- flat(stobs_mat)
  constants$post_vec <- flat(post_mat); constants$y_vec <- flat(y_mat)
  constants$NFISH <- NULL
  constants$Area <- NULL
  # NIMBLE collapses length-1 constants to scalars and then rejects indexing
  # them (single-pond builds); pad pond-indexed vectors with an unused element
  if (NP == 1) {
    constants$sp <- c(constants$sp, 1L)
    constants$Uk <- c(constants$Uk, 1L)
  }
  data <- list(
    zero = if (NP == 1) c(0, 0) else rep(0, NP),
    ev_m = RD_data$ev_m,
    sp_nL = RD_data$sp_nL,
    okU = array(1L, dim = c(NP, N_SEASONS, N_STAGE))
  )
  list(constants = constants, data = data, Fish = FishSub, Hist = Hist)
}

# ---------------------------------------------------------------------------
# nimbleFunctions
# ---------------------------------------------------------------------------
# Forward-backward over all weighted histories in one pond with three alive
# stages. Returns a T x 4 matrix: columns 1-3 = smoothed expected tagged fish
# alive in S / L / A at each primary period (counted from cfrom), [1, 4] = total
# log-likelihood. Transition matrices over the gap g are M^g where M is the
# monthly matrix (survive, then possibly advance one stage).
rdPond3 <- nimbleFunction(
  run = function(y = double(1), K = double(1), stobs = double(1), post = double(1),
                 w = double(1), entry = double(1), endt = double(1), removed = double(1),
                 stage0 = double(1), origin = double(1), cfrom = double(1),
                 nfish = double(), T = double(),
                 lphiS = double(1), lphiL = double(1), lphiA = double(1),
                 season = double(1), gap = double(1),
                 pS = double(1), pL = double(1), pA = double(1),
                 bpost = double(1), bpostJ = double(1),
                 lpsiSL = double(), lpsiLA = double(),
                 evt = double(1), dE = double()) {
    returnType(double(2))
    out <- matrix(0, nrow = T, ncol = 4)
    # Tr[t, (idx-1)*9 + (r-1)*3 + c] = P(r -> c over the gap after t), idx = post*3 + origin
    Tr <- matrix(0, nrow = T, ncol = 54)
    M  <- matrix(0, nrow = 3, ncol = 3)
    Mg <- matrix(0, nrow = 3, ncol = 3)
    psSL <- ilogit(lpsiSL)
    psLA <- ilogit(lpsiLA)
    for (t in 1:T) {
      g <- gap[t]
      s <- season[t]
      # Known die-off event: fixed logit-survival offset dE applied to all
      # three stages in event months (evt[t] == 1). Assumes gap == 1 (m12
      # calendar); under m7, a 6-month Apr->Oct gap that partly overlaps an
      # event applies dE to the whole gap (approximation, see YCB3S_data.R).
      dEt <- dE * evt[t]
      for (o in 1:3) {
        for (pp in 1:2) {
          # post-release offsets: bpostJ for S and L, bpost for A
          bp <- 0
          bpJ <- 0
          if (pp == 2) {
            bp <- bpost[o]
            bpJ <- bpostJ[o]
          }
          phS <- ilogit(lphiS[s] + bpJ + dEt)
          phL <- ilogit(lphiL[s] + bpJ + dEt)
          phA <- ilogit(lphiA[s] + bp + dEt)
          M[1, 1] <- phS * (1 - psSL); M[1, 2] <- phS * psSL;       M[1, 3] <- 0
          M[2, 1] <- 0;                M[2, 2] <- phL * (1 - psLA); M[2, 3] <- phL * psLA
          M[3, 1] <- 0;                M[3, 2] <- 0;                M[3, 3] <- phA
          Mg <- M
          if (g > 1) {
            for (k in 2:g) Mg <- Mg %*% M
          }
          idx <- (pp - 1) * 3 + o
          for (r in 1:3) {
            for (c in 1:3) Tr[t, (idx - 1) * 9 + (r - 1) * 3 + c] <- Mg[r, c]
          }
        }
      }
    }
    a  <- matrix(0, nrow = T, ncol = 4)
    em <- matrix(0, nrow = T, ncol = 4)
    ll <- 0
    for (i in 1:nfish) {
      off <- (i - 1) * T
      t1 <- entry[i]
      t2 <- endt[i]
      o  <- origin[i]
      a1 <- 0; a2 <- 0; a3 <- 0; a4 <- 0
      if (stage0[i] == 1) a1 <- 1
      if (stage0[i] == 2) a2 <- 1
      if (stage0[i] == 3) a3 <- 1
      lli <- 0
      for (t in t1:t2) {
        yy <- y[off + t]; kk <- K[off + t]; so <- stobs[off + t]
        e1 <- 1; e2 <- 1; e3 <- 1; e4 <- 1
        if (kk > 0) {
          e1 <- dbinom(yy, kk, pS[t])
          e2 <- dbinom(yy, kk, pL[t])
          e3 <- dbinom(yy, kk, pA[t])
          if (yy > 0) e4 <- 0
        }
        if (so == 1) { e2 <- 0; e3 <- 0; e4 <- 0 }
        if (so == 2) { e1 <- 0; e3 <- 0; e4 <- 0 }
        if (so == 3) { e1 <- 0; e2 <- 0; e4 <- 0 }
        if (so == 4) { e4 <- 0 }
        em[t, 1] <- e1; em[t, 2] <- e2; em[t, 3] <- e3; em[t, 4] <- e4
        a1 <- a1 * e1; a2 <- a2 * e2; a3 <- a3 * e3; a4 <- a4 * e4
        s <- a1 + a2 + a3 + a4
        if (s <= 0) s <- 1e-300
        lli <- lli + log(s)
        a1 <- a1 / s; a2 <- a2 / s; a3 <- a3 / s; a4 <- a4 / s
        a[t, 1] <- a1; a[t, 2] <- a2; a[t, 3] <- a3; a[t, 4] <- a4
        if (t < t2) {
          b0 <- (post[off + t] * 3 + o - 1) * 9
          m11 <- Tr[t, b0 + 1]; m12 <- Tr[t, b0 + 2]; m13 <- Tr[t, b0 + 3]
          m22 <- Tr[t, b0 + 5]; m23 <- Tr[t, b0 + 6]
          m33 <- Tr[t, b0 + 9]
          n1 <- a1 * m11
          n2 <- a1 * m12 + a2 * m22
          n3 <- a1 * m13 + a2 * m23 + a3 * m33
          n4 <- a4 + a1 * (1 - m11 - m12 - m13) + a2 * (1 - m22 - m23) + a3 * (1 - m33)
          a1 <- n1; a2 <- n2; a3 <- n3; a4 <- n4
        }
      }
      if (removed[i] == 1) {
        s <- a1 + a2 + a3
        if (s <= 0) s <- 1e-300
        lli <- lli + log(s)
      }
      ll <- ll + w[i] * lli
      # backward pass for smoothed state probabilities
      b1 <- 1; b2 <- 1; b3 <- 1; b4 <- 1
      if (removed[i] == 1) b4 <- 0
      for (tt in 1:(t2 - t1 + 1)) {
        t <- t2 - tt + 1
        if (t >= cfrom[i]) {
          s1 <- a[t, 1] * b1; s2 <- a[t, 2] * b2; s3 <- a[t, 3] * b3; s4 <- a[t, 4] * b4
          s <- s1 + s2 + s3 + s4
          if (s > 0) {
            out[t, 1] <- out[t, 1] + w[i] * s1 / s
            out[t, 2] <- out[t, 2] + w[i] * s2 / s
            out[t, 3] <- out[t, 3] + w[i] * s3 / s
          }
        }
        if (t > t1) {
          eb1 <- em[t, 1] * b1; eb2 <- em[t, 2] * b2; eb3 <- em[t, 3] * b3; eb4 <- em[t, 4] * b4
          b0 <- (post[off + t - 1] * 3 + o - 1) * 9
          m11 <- Tr[t - 1, b0 + 1]; m12 <- Tr[t - 1, b0 + 2]; m13 <- Tr[t - 1, b0 + 3]
          m22 <- Tr[t - 1, b0 + 5]; m23 <- Tr[t - 1, b0 + 6]
          m33 <- Tr[t - 1, b0 + 9]
          nb1 <- m11 * eb1 + m12 * eb2 + m13 * eb3 + (1 - m11 - m12 - m13) * eb4
          nb2 <- m22 * eb2 + m23 * eb3 + (1 - m22 - m23) * eb4
          nb3 <- m33 * eb3 + (1 - m33) * eb4
          nb4 <- eb4
          sb <- nb1 + nb2 + nb3 + nb4
          if (sb > 0) { b1 <- nb1 / sb; b2 <- nb2 / sb; b3 <- nb3 / sb; b4 <- nb4 / sb }
        }
      }
    }
    out[1, 4] <- ll
    return(out)
  })

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
rd3_code <- nimbleCode({
  # Species-level parameters --------------------------------------------------
  for (s in 1:NSP) {
    mu_phi[s]    ~ dnorm(3, sd = 1.5)            # mean logit monthly adult survival
    sigma_phi[s] ~ T(dnorm(0, sd = 1), 0, )      # pond x season SD
    dJ[s]        ~ dnorm(-0.5, sd = 1)           # juvenile (L) offset to logit survival
    dS[s]        ~ T(dnorm(-1, sd = 1), , 0)     # additional small-fish (S) offset
    # Direct recruit-density effect on juvenile survival removed (2026-09-08): density
    # acts through recruit size (a1 -> pLg) and stage-specific survival; a free cB was
    # unidentifiable against eps_J with one cohort per season. Kept as a fixed 0 so
    # the lphiL expression and downstream code are unchanged.
    cB[s]        <- 0
    for (o in 1:2) {
      b_post[s, o]  ~ T(dnorm(-1, sd = 1), , 0)  # post-release offset: stocked, netting-tagged
      # Juvenile post-release offset tied to the adult one: a free b_postJ sat on an
      # unidentified ridge with dJ (almost no netting-tagged juvenile outlives the window)
      b_postJ[s, o] <- b_post[s, o]
    }
    sigma_J[s] ~ T(dnorm(0, sd = 1), 0, )        # SD of season-specific juvenile survival deviation
    b_post[s, 3]  <- 0
    b_postJ[s, 3] <- 0
    lpsiSL[s] ~ dnorm(-2.4, sd = 1)              # logit monthly S -> L
    lpsiLA[s] ~ dnorm(-2.4, sd = 1)              # logit monthly L -> A
    a0[s] ~ dnorm(-1, sd = 1.5)                  # recruit size split: intercept
    a1[s] ~ dnorm(0, sd = 1)                     # recruit size split: density slope
    g_eff[s] ~ dnorm(0, sd = 1)
    d_moy[s, 1] <- 0
    for (mm in 2:NMOY) { d_moy[s, mm] ~ dnorm(0, sd = 1) }
    mu_r[s]    ~ dnorm(-1.2, sd = 1.5)
    sigma_r[s] ~ T(dnorm(0, sd = 1), 0, )
    nu[s]      ~ dexp(1)
  }

  # Extra pond x season adult-survival deviation for IP1 x FY2022 (BK request,
  # 2026-09-10): a fixed-SD (0.5 logit) deviation added only to lphiA_eff at
  # (IP1, FY2022) via the 0/1 constant IP1FY22[j, y]; every other pond-season
  # cell has IP1FY22 == 0 so lphiA_eff == lphiA there, unchanged. This sits on
  # top of the existing hierarchical lphiA[j,y] ~ N(mu_phi, sigma_phi) draw and
  # the pond-level die-off event offset dE[j] (applied downstream in rdPond3
  # for event months) still applies as before.
  eps_A_IP1_FY2022 ~ dnorm(0, sd = 0.5)

  # Pond-level structure ------------------------------------------------------
  for (j in 1:NPOND) {
    for (y in 1:NS) {
      lphiA[j, y] ~ dnorm(mu_phi[sp[j]], sd = sigma_phi[sp[j]])
      lphiA_eff[j, y] <- lphiA[j, y] + eps_A_IP1_FY2022 * IP1FY22[j, y]
      # eps_J: season-specific juvenile deviation so a juvenile winter crash does not
      # have to be absorbed by the stage-shared lphiA. Centered (sum-to-zero) so that
      # dJ is the mean juvenile offset; the free eps_raw was confounded with dJ.
      eps_raw[j, y] ~ dnorm(0, sd = sigma_J[sp[j]])
      eps_J[j, y] <- eps_raw[j, y] - mean(eps_raw[j, 1:NS])
      lphiL[j, y] <- lphiA_eff[j, y] + dJ[sp[j]] + cB[sp[j]] * dens[j, y] + eps_J[j, y]
      lphiS[j, y] <- lphiL[j, y] + dS[sp[j]]
      logit(pLg[j, y]) <- a0[sp[j]] + a1[sp[j]] * dens[j, y]   # fraction of recruits >= T1 at netting
    }
    # Juvenile detection shared between S and L (first fit left lp[j, 2] free and it
    # collapsed to an implausibly low value with only 3 juvenile stage re-observations)
    lp[j, 1] ~ dnorm(1, sd = 1.5)
    lp[j, 2] <- lp[j, 1]
    lp[j, 3] ~ dnorm(1, sd = 1.5)
    # Known die-off (data/BWDieOffEvents.csv): fixed, pond-specific, survival-
    # reducing-only logit offset applied to S/L/A alike in event months
    # (Evt[j,t] == 1). Ponds with no recorded event carry an uninformative
    # (prior-only) dE.
    dE[j] ~ T(dnorm(0, sd = 3), , 0)
    for (t in 1:T) {
      logit(pS[j, t]) <- lp[j, 1] + g_eff[sp[j]] * eff[j, t] + d_moy[sp[j], moy[t]]
      logit(pL[j, t]) <- lp[j, 2] + g_eff[sp[j]] * eff[j, t] + d_moy[sp[j], moy[t]]
      logit(pA[j, t]) <- lp[j, 3] + g_eff[sp[j]] * eff[j, t] + d_moy[sp[j], moy[t]]
    }

    RD[j, 1:T, 1:4] <- rdPond3(
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
      lphiS = lphiS[j, 1:NS], lphiL = lphiL[j, 1:NS], lphiA = lphiA_eff[j, 1:NS],
      season = season[1:T], gap = gap[1:T],
      pS = pS[j, 1:T], pL = pL[j, 1:T], pA = pA[j, 1:T],
      bpost = b_post[sp[j], 1:3], bpostJ = b_postJ[sp[j], 1:3],
      lpsiSL = lpsiSL[sp[j]], lpsiLA = lpsiLA[sp[j]],
      evt = Evt[j, 1:T], dE = dE[j])
    zero[j] ~ dLLtrick(RD[j, 1, 4])
    for (st in 1:NSTAGE) {
      for (t in 1:T) { N_tag[j, t, st] <- RD[j, t, st] }
    }

    # Untagged pool -----------------------------------------------------------
    for (st in 1:NSTAGE) {
      U0[j, st] ~ dunif(0, 5000)
      U[j, 1, st] <- (1 - Uk[j]) * U0[j, st]
    }
    for (y in 1:(NS - 1)) {
      Atot[j, y] <- N_tag[j, tApr[y], 3] + U[j, y, 3]
      lr[j, y] ~ dnorm(mu_r[sp[j]], sd = sigma_r[sp[j]])
      r[j, y] <- exp(lr[j, y])
      RecJ[j, y] <- r[j, y] * Atot[j, y]              # recruits alive at the next fall netting
      for (st in 1:NSTAGE) { Uaft[j, y, st] <- max(0, U[j, y, st] - remU[j, y, st]) }
      # Untagged-pool annual survival gets the same event discount as the
      # tagged model: nEvt[j,y] of the 12 calendar months in the season use
      # the event-adjusted monthly rate (dE), the rest use the baseline rate
      phiUS[j, y] <- ilogit(lphiS[j, y])^(12 - nEvt[j, y]) * ilogit(lphiS[j, y] + dE[j])^nEvt[j, y]
      phiUL[j, y] <- ilogit(lphiL[j, y])^(12 - nEvt[j, y]) * ilogit(lphiL[j, y] + dE[j])^nEvt[j, y]
      phiUA[j, y] <- ilogit(lphiA_eff[j, y])^(12 - nEvt[j, y]) * ilogit(lphiA_eff[j, y] + dE[j])^nEvt[j, y]
      psiYSL[j, y] <- 1 - (1 - ilogit(lpsiSL[sp[j]]))^12
      psiYLA[j, y] <- 1 - (1 - ilogit(lpsiLA[sp[j]]))^12
      U[j, y + 1, 1] <- Uaft[j, y, 1] * phiUS[j, y] * (1 - psiYSL[j, y]) +
        RecJ[j, y] * (1 - pLg[j, y + 1])
      U[j, y + 1, 2] <- Uaft[j, y, 1] * phiUS[j, y] * psiYSL[j, y] +
        Uaft[j, y, 2] * phiUL[j, y] * (1 - psiYLA[j, y]) + RecJ[j, y] * pLg[j, y + 1]
      U[j, y + 1, 3] <- Uaft[j, y, 2] * phiUL[j, y] * psiYLA[j, y] +
        Uaft[j, y, 3] * phiUA[j, y] + nu[sp[j]]
    }
    for (y in 1:NS) {
      for (st in 1:NSTAGE) {
        okU[j, y, st] ~ dconstraint(U[j, y, st] >= remU[j, y, st])
      }
    }
  }

  # Recruit size split observed at fall nettings -----------------------------
  for (e in 1:NSPLIT) {
    sp_nL[e] ~ dbin(pLg[sp_pond[e], sp_season[e]], sp_n0[e])
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
# Initial values
# ---------------------------------------------------------------------------
make_3s_inits <- function(inp, seed) {
  set.seed(seed)
  cst <- inp$constants
  NP <- cst$NPOND; NS <- cst$NS
  U0 <- matrix(runif(NP * 3, 100, 300), NP, 3)
  for (j in 1:NP) U0[j, ] <- pmax(U0[j, ], cst$remU[j, 1, ] + 50)
  mu_phi <- rnorm(2, 3, 0.3)
  mu_r <- rnorm(2, log(3), 0.1)
  list(
    mu_phi = mu_phi, sigma_phi = runif(2, 0.2, 0.6), dJ = rnorm(2, -0.5, 0.2),
    dS = runif(2, -1.5, -0.5),
    b_post = cbind(matrix(runif(4, -1.5, -0.3), 2, 2), 0),
    sigma_J = runif(2, 0.2, 0.6),
    eps_raw = matrix(rnorm(NP * NS, 0, 0.2), NP, NS),
    lpsiSL = rnorm(2, -2.4, 0.2), lpsiLA = rnorm(2, -2.4, 0.2),
    a0 = rnorm(2, -1, 0.3), a1 = rnorm(2, -0.5, 0.2),
    g_eff = rnorm(2, 0.3, 0.1),
    d_moy = cbind(0, matrix(rnorm(2 * (cst$NMOY - 1), 0, 0.2), 2)),
    mu_r = mu_r, sigma_r = runif(2, 0.3, 0.8), nu = runif(2, 4, 6),
    lp = cbind(rnorm(NP, 1, 0.3), NA, rnorm(NP, 1, 0.3)),  # lp[, 2] is deterministic
    lphiA = matrix(rnorm(NP * NS, mu_phi[cst$sp], 0.2), NP, NS),
    lr = matrix(rnorm(NP * (NS - 1), mu_r[cst$sp], 0.1), NP, NS - 1),
    U0 = U0,
    dE = -runif(NP, 0.2, 1),
    eps_A_IP1_FY2022 = rnorm(1, 0, 0.2)
  )
}

repair_3s_inits <- function(model, max_tries = 300) {
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
      } else if (st < 3) {
        model$lr[j, y - 1] <- model$lr[j, y - 1] + 1
      } else {
        model$nu[sp[j]] <- model$nu[sp[j]] + 2
        model$lr[j, y - 1] <- model$lr[j, y - 1] + 0.5
      }
    }
  }
  stop("Could not repair initial values")
}

# ---------------------------------------------------------------------------
# MCMC configuration and chain runner
# ---------------------------------------------------------------------------
rd3_monitors <- c("mu_phi", "sigma_phi", "dJ", "dS", "dE", "b_post", "sigma_J", "eps_J", "lpsiSL", "lpsiLA",
                  "a0", "a1", "lp", "g_eff", "d_moy", "mu_r", "sigma_r", "nu", "U0",
                  "eps_A_IP1_FY2022", "lphiA", "lphiA_eff", "lphiL", "lphiS", "pLg",
                  "r", "U", "RecJ", "N_tag", "Tal", "Atot")

rd3_configure <- function(model) {
  conf <- configureMCMC(model, monitors = rd3_monitors, enableWAIC = FALSE, print = FALSE)
  cst <- model$getConstants()
  NMOY <- cst$NMOY; NP <- cst$NPOND; NS <- cst$NS
  for (s in 1:2) {
    tgt_phi <- c(paste0("dJ[", s, "]"), paste0("dS[", s, "]"),
                 paste0("b_post[", s, ", 1]"), paste0("b_post[", s, ", 2]"),
                 paste0("lpsiSL[", s, "]"), paste0("lpsiLA[", s, "]"))
    conf$removeSamplers(tgt_phi)
    conf$addSampler(target = tgt_phi, type = "RW_block", control = list(adaptInterval = 200))
    tgt_a <- c(paste0("a0[", s, "]"), paste0("a1[", s, "]"))
    conf$removeSamplers(tgt_a)
    conf$addSampler(target = tgt_a, type = "RW_block", control = list(adaptInterval = 200))
    tgt_p <- c(paste0("g_eff[", s, "]"), paste0("d_moy[", s, ", ", 2:NMOY, "]"))
    conf$removeSamplers(tgt_p)
    conf$addSampler(target = tgt_p, type = "RW_block", control = list(adaptInterval = 200))
  }
  # IP1's pond index (if modeled): identify from the IP1FY22 indicator so the
  # extra eps_A_IP1_FY2022 deviation is jointly blocked with that pond's lphiA
  ip1_j <- which(rowSums(matrix(cst$IP1FY22, nrow = NP)) > 0)
  for (j in 1:NP) {
    tgt_e <- paste0("lphiA[", j, ", ", 1:NS, "]")
    if (length(ip1_j) == 1 && j == ip1_j) tgt_e <- c(tgt_e, "eps_A_IP1_FY2022")
    conf$removeSamplers(tgt_e)
    conf$addSampler(target = tgt_e, type = "RW_block", control = list(adaptInterval = 200))
    tgt_eJ <- paste0("eps_raw[", j, ", ", 1:NS, "]")
    conf$removeSamplers(tgt_eJ)
    conf$addSampler(target = tgt_eJ, type = "RW_block", control = list(adaptInterval = 200))
    tgt_lp <- paste0("lp[", j, ", ", c(1, 3), "]")
    conf$removeSamplers(tgt_lp)
    conf$addSampler(target = tgt_lp, type = "RW_block", control = list(adaptInterval = 200))
  }
  conf
}

run_3s_chain <- function(inp, seed, niter, nburnin, thin, progress = TRUE) {
  inits <- make_3s_inits(inp, seed)
  model <- nimbleModel(rd3_code, constants = inp$constants, data = inp$data,
                       inits = inits, calculate = TRUE)
  repair_3s_inits(model)
  conf   <- rd3_configure(model)
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
