## GIEL_GrowthTrackComparison.R
##
## Three-track bonytail (GIEL) growth-model comparison (2026-09-14), Step 2 of
## the growth-model sensitivity plan feeding the joint robust-design model's
## size-tied maturation hazard (see AGENTS.md "GIEL Size-Tied Maturation
## Hazard" and the plan file .posit/assistant/plans/2026-09-14-1444-giel-
## growth-model-tracks-feeding-robust-design-maturation-hazard.md).
##
## Tracks compared (Track A2, age-based + Cibola, is explicitly out of scope):
##   F1 -- Fabens increment model, IPCA-only         (GIEL_GrowthSeasonThreshold.R)
##   F2 -- Fabens increment model, IPCA + Cibola      (GIEL_GrowthSeasonThreshold_withCibola.R)
##   A1 -- Age-based VBGF, IPCA-only                  (GIEL_VBGF_AgeLength.R)
##
## Self-contained: loads all three tracks' saved .RData fresh (does not rely
## on session state), so it can be re-run standalone once the three source
## scripts have been run at least once.

library(tidyverse)

## ---- 1. Load each track's saved fit + derived probability table ----------

f1 <- new.env(); load("data/GIEL_GrowthSeasonThreshold.RData", envir = f1)
f2 <- new.env(); load("data/GIEL_GrowthSeasonThreshold_withCibola.RData", envir = f2)
a1 <- new.env(); load("data/GIEL_VBGF_AgeLength.RData", envir = a1)

## ---- 2. Combined P(TL >= 250mm after one growth season | TL1) table ------

tab_f1 <- f1$growth_season_threshold_probs |>
  transmute(TL1, F1_expected_increment_mm = Expected_increment_mm,
            F1_P_cross_250 = P_cross_250_one_growth_season)
tab_f2 <- f2$growth_season_threshold_probs_combined |>
  transmute(TL1, F2_expected_increment_mm = Expected_increment_mm,
            F2_P_cross_250 = P_cross_250_one_growth_season)
tab_a1 <- a1$growth_season_threshold_probs_A1 |>
  transmute(TL1, A1_expected_increment_mm = Expected_increment_mm,
            A1_P_cross_250 = P_cross_250_one_growth_season)

GrowthTrackComparisonTable <- tab_f1 |>
  inner_join(tab_f2, by = "TL1") |>
  inner_join(tab_a1, by = "TL1")

cat("Three-track P(TL >= 250mm after one growth season | TL1) comparison:\n")
print(GrowthTrackComparisonTable)

## Simple sensitivity summary: range and max spread across tracks at each TL1
GrowthTrackSpread <- GrowthTrackComparisonTable |>
  rowwise() |>
  mutate(
    min_P = min(F1_P_cross_250, F2_P_cross_250, A1_P_cross_250),
    max_P = max(F1_P_cross_250, F2_P_cross_250, A1_P_cross_250),
    spread_P = max_P - min_P
  ) |>
  ungroup() |>
  select(TL1, min_P, max_P, spread_P)

cat("\nSpread across tracks (max - min) at each TL1:\n")
print(GrowthTrackSpread)

## ---- 3. Posterior VBGF curves for the three tracks (common TL0 = 100mm) --

pred_curve <- function(Linf_draws, K_draws, t_grid, TL0) {
  TL_mat <- sapply(t_grid, function(t) Linf_draws - (Linf_draws - TL0) * exp(-K_draws * t))
  data.frame(t = t_grid, median = apply(TL_mat, 2, median),
             lwr = apply(TL_mat, 2, quantile, 0.025), upr = apply(TL_mat, 2, quantile, 0.975))
}

t_grid <- seq(0, 12, by = 0.25)
TL0_plot <- 100

## F1 (backwater-only Fabens)
samps_f1 <- do.call(rbind, f1$gs_winner$samples)
Linf_f1 <- samps_f1[, "Linf"]
set.seed(1301)
K_f1 <- exp(samps_f1[, "logK[1]"] + rnorm(nrow(samps_f1), 0, samps_f1[, "sigma_v"]))
curve_f1 <- pred_curve(Linf_f1, K_f1, t_grid, TL0_plot) |> mutate(Track = "F1: Fabens, IPCA-only")

## F2 (backwater + Cibola Fabens)
samps_f2 <- do.call(rbind, f2$comb_winner$samples)
Linf_f2 <- samps_f2[, "Linf"]
set.seed(1302)
sigma_v_f2 <- if (f2$winner_use_group_comb) samps_f2[, "sigma_v"] else rep(0, nrow(samps_f2))
K_f2 <- exp(samps_f2[, "logK[1]"] + rnorm(nrow(samps_f2), 0, sigma_v_f2))
curve_f2 <- pred_curve(Linf_f2, K_f2, t_grid, TL0_plot) |> mutate(Track = "F2: Fabens, IPCA + Cibola")

## A1 (age-based VBGF, IPCA-only)
samps_a1 <- do.call(rbind, a1$al_final$samples)
Linf_a1 <- samps_a1[, "Linf"]
set.seed(1303)
K_a1 <- exp(samps_a1[, "logK[1]"] + rnorm(nrow(samps_a1), 0, samps_a1[, "sigma_v"]))
curve_a1 <- pred_curve(Linf_a1, K_a1, t_grid, TL0_plot) |> mutate(Track = "A1: Age-based VBGF, IPCA-only")

curves_all <- bind_rows(curve_f1, curve_f2, curve_a1)

linf_lines <- data.frame(
  Track = c("F1: Fabens, IPCA-only", "F2: Fabens, IPCA + Cibola", "A1: Age-based VBGF, IPCA-only"),
  Linf = c(median(Linf_f1), median(Linf_f2), median(Linf_a1))
)

track_colors <- c("F1: Fabens, IPCA-only" = "#2c7fb8",
                   "F2: Fabens, IPCA + Cibola" = "#33a02c",
                   "A1: Age-based VBGF, IPCA-only" = "#e6550d")

p_tracks <- ggplot(curves_all, aes(t, median, color = Track, fill = Track)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(data = linf_lines, aes(yintercept = Linf, color = Track),
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  scale_color_manual(values = track_colors) +
  scale_fill_manual(values = track_colors) +
  labs(title = "Bonytail VBGF: three growth-model tracks",
       subtitle = paste0("Common start TL0 = ", TL0_plot, " mm; dashed = posterior median Linf"),
       x = "Growing-season years elapsed", y = "Total length (mm)",
       caption = paste0(
         "F1/F2 fit on paired growth increments (Fabens); A1 fit on absolute age-length pairs.\n",
         "All three tracks' Linf sits at or below each dataset's own observed max TL -- see\n",
         "GIEL_GrowthTrackComparison_Summary.md for the shared identifiability caveat.")) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

ggsave("model/GIEL_GrowthTrackComparison.png", p_tracks, width = 7.5, height = 5.5, dpi = 150)

## ---- 4. P(TL >= 250mm) overlay plot ---------------------------------------

prob_long <- GrowthTrackComparisonTable |>
  select(TL1, F1_P_cross_250, F2_P_cross_250, A1_P_cross_250) |>
  pivot_longer(-TL1, names_to = "Track", values_to = "P_cross_250") |>
  mutate(Track = recode(Track,
    F1_P_cross_250 = "F1: Fabens, IPCA-only",
    F2_P_cross_250 = "F2: Fabens, IPCA + Cibola",
    A1_P_cross_250 = "A1: Age-based VBGF, IPCA-only"
  ))

p_prob <- ggplot(prob_long, aes(TL1, P_cross_250, color = Track)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = track_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "P(TL >= 250mm after one Apr-Sep growth season | TL1)",
       subtitle = "Three growth-model tracks compared on the same TL1 grid",
       x = "TL at start of growth season (mm)", y = "P(cross 250mm)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

ggsave("model/GIEL_GrowthTrackComparison_ProbCurves.png", p_prob, width = 7, height = 5, dpi = 150)

## ---- 5. Save results -------------------------------------------------------

save(
  GrowthTrackComparisonTable, GrowthTrackSpread,
  curve_f1, curve_f2, curve_a1, curves_all, linf_lines,
  file = "data/GIEL_GrowthTrackComparison.RData"
)

cat("\nSaved: data/GIEL_GrowthTrackComparison.RData, model/GIEL_GrowthTrackComparison.png,",
    "model/GIEL_GrowthTrackComparison_ProbCurves.png\n")
