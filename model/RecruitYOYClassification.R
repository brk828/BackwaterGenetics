# ---------------------------------------------------------------------------
# Reusable young-of-year (YOY) vs. carryover classification for untagged
# fish first handled at a fall (Oct-Dec) netting event.
#
# Rationale: "recruit" counts used for density-dependent growth/survival
# analyses (model/YCBSizeDensityEDA.qmd, model/YCB3S_data.R's RecruitIndex)
# currently include ANY untagged first-capture at a fall netting, but some of
# those are older untagged fish that simply weren't captured in a prior
# year's netting, not this season's age-0 cohort. This mirrors the
# mclust::Mclust(TL, G = 1:2) approach already used for the GIEL bimodal
# reporting table in PopulationMonitoring.qmd (@tbl-giel-bimodal), generalized
# into a reusable classifier so it can drive model inputs, not just reporting
# text.
#
# classify_yoy_cohort(TL, min_n = 15, default_is_yoy = TRUE) returns, for one
# pond x fiscal-year netting cohort's vector of total lengths:
#   - is_yoy: logical, TRUE for fish assigned to the smaller-mean mixture
#     component (or all `default_is_yoy` if BIC prefers G = 1, or if n < min_n
#     -- too few fish to separate two modes reliably).
#   - G: number of components selected (1 or 2, or NA if min_n rule applied
#     without fitting).
#   - mu: sorted component mean(s) (length 1 or 2).
#   - prob_yoy: posterior probability of cohort membership in the
#     smaller-mean component (as.numeric(default_is_yoy) for every fish when
#     G == 1 or n < min_n).
#
# `default_is_yoy` (2026-09-24): when a too-small cohort can't be split, the
# old hardcoded default (always TRUE, "assume it's a small recruit sample")
# is wrong for a pond-season where the TRUE recruit cohort was already
# handled in bulk (harvested/returned untagged with only an aggregate mean
# TL, e.g. data/BWNettingEvents.csv) -- in that situation, any *individually
# tagged* leftover fish are almost certainly established/carryover fish that
# finally got caught, not recruits, because genuine YOY in a big-recruitment
# year get bulk-processed rather than individually tagged. Callers should
# pass `default_is_yoy = FALSE` (or, via classify_yoy_df's `default_col`, a
# per-group value) whenever a bulk untagged record already accounts for that
# season's recruits; otherwise the historical default (TRUE) still applies.
#
# classify_yoy_df(df, tl_col, group_vars, min_n = 15, default_col = NULL)
# applies the above to a data frame, grouped by group_vars (e.g. pond, season
# / Backwater, FY). `default_col`, if given, names a logical column in `df`
# (constant within each group) supplying that group's `default_is_yoy`;
# otherwise every group defaults to TRUE (unchanged behavior for existing
# callers). Returns the input df with is_yoy/G/mu1/mu2/prob_yoy columns
# appended. Row order is preserved.
# ---------------------------------------------------------------------------

packages(mclust)

classify_yoy_cohort <- function(TL, min_n = 15, gap_ratio_min = 5, gap_abs_min = 15,
                                 default_is_yoy = TRUE) {
  n <- length(TL)
  na_idx <- is.na(TL)
  out <- tibble::tibble(is_yoy = rep(NA, n), G = NA_integer_, mu1 = NA_real_,
                         mu2 = NA_real_, prob_yoy = NA_real_)
  if (all(na_idx)) return(out)

  tl_ok <- TL[!na_idx]
  if (length(tl_ok) < min_n) {
    out$is_yoy[!na_idx] <- default_is_yoy
    out$G[!na_idx] <- NA_integer_
    out$mu1[!na_idx] <- mean(tl_ok)
    out$prob_yoy[!na_idx] <- as.numeric(default_is_yoy)
    return(out)
  }

  m <- mclust::Mclust(tl_ok, G = 1:2, verbose = FALSE)

  # BIC alone will happily split a merely right-skewed (or otherwise
  # non-Gaussian) unimodal distribution into two Gaussians with no real
  # valley between them -- observed for two YCB cohorts (FY2019: BIC picked
  # G = 2 at means 166 vs 183 mm; FY2022: 210 vs 260 mm), both continuous,
  # gap-free declines in counts that a plain BIC-selected mixture nonetheless
  # "explained" with two narrow components. Two statistical gates were tried
  # and rejected: Ashman's D (McLachlan & Peel 2000) passed both spurious
  # cases (D = 3.0-4.7); Hartigan's dip test correctly rejected them but also
  # failed to flag a genuinely bimodal cohort with a small (8 of 58 fish),
  # sharply separated upper cluster (p = 0.95) because the dip test has low
  # power against a small, tight second cluster at modest sample sizes.
  # Gate instead on a direct, inspectable statistic: is there an actual empty
  # region of TL values between the two fitted component means, well beyond
  # the typical spacing between adjacent observed lengths in that cohort?
  # Verified across every YCB cohort with a documented upper mode: real
  # splits (FY2014, FY2024, FY2025) have gap-to-median-spacing ratios of
  # 31-60 (absolute gaps 31-62 mm); spurious BIC splits (FY2019, FY2022) have
  # ratios of 1-3 (gaps 1-3 mm). Default thresholds (ratio >= 5, gap >= 15 mm)
  # sit well inside that separation.
  gap_between_modes <- function(TL, mu) {
    s <- sort(unique(TL))
    d <- diff(s)
    if (length(d) == 0) return(c(maxgap = 0, ratio = 0))
    mids <- s[-length(s)] + d / 2
    in_range <- mids > min(mu) & mids < max(mu)
    if (!any(in_range)) return(c(maxgap = 0, ratio = 0))
    maxgap <- max(d[in_range])
    c(maxgap = maxgap, ratio = maxgap / median(d))
  }
  gap_ok <- !is.null(m) && m$G == 2 && {
    g <- gap_between_modes(tl_ok, m$parameters$mean)
    g["ratio"] >= gap_ratio_min && g["maxgap"] >= gap_abs_min
  }
  spurious_split <- !is.null(m) && m$G == 2 && !gap_ok

  if (is.null(m) || m$G == 1 || spurious_split) {
    out$is_yoy[!na_idx] <- TRUE
    out$G[!na_idx] <- if (is.null(m)) NA_integer_ else 1L
    out$mu1[!na_idx] <- if (is.null(m)) mean(tl_ok) else mean(tl_ok)
    out$prob_yoy[!na_idx] <- 1
    return(out)
  }

  mu_order <- order(m$parameters$mean)          # 1 = smaller mean (YOY)
  yoy_comp <- mu_order[1]
  cls <- m$classification                        # MAP assignment, 1 or 2
  z <- m$z[, yoy_comp]                            # posterior prob of YOY component

  out$is_yoy[!na_idx] <- cls == yoy_comp
  out$G[!na_idx] <- 2L
  out$mu1[!na_idx] <- m$parameters$mean[mu_order[1]]
  out$mu2[!na_idx] <- m$parameters$mean[mu_order[2]]
  out$prob_yoy[!na_idx] <- z
  out
}

classify_yoy_df <- function(df, tl_col, group_vars, min_n = 15, default_col = NULL) {
  tl_col <- rlang::ensym(tl_col)
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::group_modify(~ {
      default_yoy <- if (is.null(default_col)) TRUE else .x[[default_col]][1]
      dplyr::bind_cols(.x, classify_yoy_cohort(.x[[rlang::as_string(tl_col)]], min_n = min_n,
                                                default_is_yoy = default_yoy))
    }) |>
    dplyr::ungroup()
}
