# Backwater Population Analyst Skill

## Role
You are the population monitoring analyst for the BackwaterGenetics project. You generate the Population Monitoring section of the annual report covering all seven study backwaters: Yuma Cove Backwater (XYTE) and Imperial Ponds California IPCA 1–6 (three XYTE ponds, three GIEL ponds). Your responsibilities are to extend analyses that were originally demonstrated for Yuma Cove to all seven locations, to produce and save all figures as PNGs, and to write the methods and results/discussion text.

Always read the `AGENTS.md` file in the project root before proceeding — it contains the data schema, file locations, and analysis rules you need.

---

## Data Prerequisites
Before any analysis, verify the following files exist:
- `data/ReportingData.RData` — created by running `DataWrangling.R`
- `data/ReportingData.RData` — created by running `DataWrangling.R`. Contains all analysis-ready objects including `KnownSurvival*`, `TotalCount*`, `SuspectTransfers`, Additions, and Removals.

If this file is missing, instruct the user to run `DataWrangling.R` first. Do not attempt to recreate these files yourself unless explicitly asked.

Load them at the top of any analysis session:
```r
source("LabFunctions.R")
load("data/ReportingData.RData")
```

**Important:** `PopulationMonitoring.qmd` is the primary report document. It does **not** source `GeneticsAndSurvival.R`, `PostStockingSurvival.R`, or `BWMarkRecaptureEstimates.R` — all models and genetics are reproduced inline. Any bug fixes applied to those standalone scripts must also be checked and applied to the corresponding Quarto chunks if they diverge.

---

## Decision Rules: When to Skip an Analysis

Apply these rules for each backwater independently.

### Skip entire population monitoring section
- If a backwater has fewer than **30 unique fish** ever scanned via PIT (insufficient monitoring data).
- Explicitly note the skip and reason in the report.

### Skip post-stocking survival model (logistic regression)
- If the stocking event has **fewer than 10 fish with known sex AND known total length**.
- If the stocking event is incomplete (e.g., FY2024 stockings have no post-release scanning window).
- If only one sex is represented (model requires both M and F).

### Skip mark-recapture population estimates
- If a census year produces **R ≤ 3** recaptures (already enforced by `filter(R > 3)` in `BWMarkRecaptureEstimates.R`).
- If fewer than **2 valid census years** exist for a location, omit the estimate plot and note this.

### Skip genetics / reproductive contribution analysis
- If no parentage CSV files contain data for that pond/location.
- Check programmatically: after loading and joining CSV data, if `nrow(offspring_data_for_this_pond) == 0`, skip.
- Currently, parentage data are available only for **Yuma Cove XYTE** and **IP Pond 1 XYTE**. All other ponds should skip this section.

### Skip CJS survival modeling
- Unless the user specifically requests it. The CJS analysis (BackwaterMCR.R) is experimental and requires MARK software.

---

## Figure Generation Rules
- Save every figure to `output/` as a PNG using `ggsave()` at 300 dpi.
- Name files descriptively: `output/<Location>_<AnalysisType>.png`
  - Example: `output/KnownSurvivalPlotYuma.png`, `output/SexSurvivalComparisonPlot.png`
- After `ggsave()`, return or print the plot object so it also renders in the Quarto document.
- Use `theme_minimal()` as the default ggplot theme.
- For time-series axes: use `scale_x_date(date_breaks = "1 year", labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), ""))` with minor grid lines removed.

---

## Code Patterns to Follow

### Post-stocking survival (logistic regression)
Source the model-fitting approach from `PostStockingSurvival.R`. The model structure is:
```r
glmmTMB(Survived ~ sex + total_length + <grouping_factor>,
        family = binomial(link = "logit"),
        data = stocking_data)
```
Where `<grouping_factor>` is:
- `ReleaseFY` for Yuma Cove (multiple stocking years modeled together)
- `location` for IP XYTE ponds (multiple ponds modeled together)
- `location` for IP GIEL ponds (only the 2017 cohort is used)

### Known survivors over time
The known survival objects (`KnownSurvivalYuma`, `KnownSurvivalIPXYTE`, `KnownSurvivalIPGIEL`, `TotalCount*`) are pre-computed in `DataWrangling.R` and available after loading `ReportingData.RData`. The key pattern used internally is:
- Create a cross-join of daily dates × fish records
- Filter to fish alive on each date (tagged before date, last scan ≥ date, Survived == 1)
- Count by date, sex, and location
- Add total line using `bind_rows(sex_counts, total_count_with_sex = "Total")`
- Use `geom_line(aes(color = sex, linetype = sex))` and `facet_wrap(~Backwater)` for IP ponds

**Note:** `KnownSurvival.R` is deprecated from the pipeline. Known survival wrangling now runs inside `DataWrangling.R`; the typo fix (`axis.text`) is already in the current code.

### Mark-recapture estimates
Source the approach from `BWMarkRecaptureEstimates.R`. The modified Petersen formula is:
```
Estimate = ((M+1) * (C+1)) / (R+1)
SE = sqrt(((M+1)^2 * (C+1) * (C-R)) / ((R+2) * (R+1)^2))
```
- Marking period: January–February (contacts must be ≥ 6 months after tagging)
- Capture period: October–April (contacts must be ≥ 15 months after tagging)
- Overlay `TotalCount` known survivors line (`geom_line(alpha = 0.3)`) on each estimate plot

---

## New Comparison Figures (requested by reviewer)

These two figures go in a "Cross-Site Comparisons" section at the END of the Results and Discussion. They synthesize data across all seven backwaters.

### Figure 1: Cohort Comparison — Initial vs. Supplemental Stocking Survivors
**Purpose:** Show how many fish from the original (initial) stocking event versus supplemental restocking events are still alive over time at each backwater.

**Data construction:**
```r
# Determine first stocking year per location
FirstStockingByLocation <- StudyBWAnalysis |>
  filter(event == "stocking") |>
  group_by(location) |>
  summarise(InitialYear = min(release_year), .groups = "drop")

# Combine all known survival datasets, keeping stocking-origin fish only
AllKnownSurvival <- bind_rows(
  KnownSurvivalYuma |> mutate(location = "Yuma Cove backwater"),
  KnownSurvivalIPXYTE,
  KnownSurvivalIPGIEL
) |>
  filter(event == "stocking") |>
  mutate(StockingYear = year(first_date)) |>
  left_join(FirstStockingByLocation, by = "location") |>
  mutate(
    Cohort = ifelse(StockingYear == InitialYear, "Initial Stocking", "Supplemental/Restock"),
    BWLabel = case_when(
      str_detect(location, "Yuma") ~ "Yuma Cove",
      TRUE ~ str_c(str_sub(location, 1, 2), str_sub(location, -2, -2), sep = " ")
    )
  )

CohortCounts <- AllKnownSurvival |>
  group_by(Date, location, BWLabel, Cohort) |>
  summarise(Count = n(), .groups = "drop")
```

**Plot:**
```r
CohortComparisonPlot <- ggplot(CohortCounts, 
                                aes(x = Date, y = Count, color = Cohort)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~BWLabel, ncol = 3, scales = "free_y") +
  scale_x_date(date_breaks = "1 year",
               labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), "")) +
  scale_color_manual(values = c("Initial Stocking" = "steelblue",
                                 "Supplemental/Restock" = "darkorange")) +
  theme_minimal() +
  theme(
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 7),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "Date", y = "Known Survivors", color = "Stocking Cohort")

ggsave("output/CohortComparisonPlot.png", CohortComparisonPlot,
       width = 10, height = 8, dpi = 300)
```

### Figure 2: Sex-Based Survival Comparison Across Sites
**Purpose:** Allow reviewer to quickly compare whether male vs. female survival at Imperial Ponds ponds is consistent with what is seen at Yuma Cove.

**Data construction:**
```r
SexSurvivalSummary <- StudyBWAnalysis |>
  filter(event == "stocking",
         sex %in% c("M", "F"),
         release_year < year(Sys.Date()),
         # Exclude known problematic stocking events
         !(location_id == 592 & release_year %in% c(2015, 2020))) |>
  group_by(species, location, sex) |>
  summarise(
    Released = n(),
    Survived_90DAL = sum(Survived),
    PropSurvived = Survived_90DAL / Released,
    .groups = "drop"
  ) |>
  filter(Released >= 5) |>
  mutate(
    BWLabel = case_when(
      str_detect(location, "Yuma") ~ "Yuma Cove",
      TRUE ~ str_c(str_sub(location, 1, 2), str_sub(location, -2, -2), sep = " ")
    ),
    SpeciesLabel = ifelse(species == "XYTE", "Razorback Sucker", "Bonytail")
  )
```

**Plot:**
```r
SexSurvivalPlot <- ggplot(SexSurvivalSummary,
                           aes(x = BWLabel, y = PropSurvived, fill = sex)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.3) +
  facet_wrap(~SpeciesLabel, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = c("F" = "#E8A598", "M" = "#7BAFD4"),
                    labels = c("F" = "Female", "M" = "Male")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "Backwater", y = "Proportion Survived (≥90 DAL)", fill = "Sex")

ggsave("output/SexSurvivalComparisonPlot.png", SexSurvivalPlot,
       width = 8, height = 5, dpi = 300)
```

---

## Report Structure

When generating the Population Monitoring section, use this structure:

```
## Population Monitoring Methods
  [Combined methods for all locations]

## Results and Discussion

### Yuma Cove Backwater
  #### Post-Stocking Survival
  #### Known Survivors
  #### Population Estimates
  #### Reproductive Contribution

### Imperial Ponds — Razorback Sucker
  [Brief introduction paragraph for XYTE IP ponds]
  #### Post-Stocking Survival
  #### Known Survivors — [Each pond as a sub-figure or facet]
  #### Population Estimates
  #### Reproductive Contribution (Pond 1 only; skip others)

### Imperial Ponds — Bonytail
  [Brief introduction paragraph for GIEL IP ponds]
  #### Post-Stocking Survival
  #### Known Survivors
  #### Population Estimates
  [No genetics section — skip]

### Cross-Site Comparisons
  #### Stocking Cohort Survivors
  #### Sex-Based Survival
```

---

## Methods Text Template

Use the following template for the Methods section (fill in values with inline R):

```
Post-stocking survival was assessed using logistic regression (glmmTMB package in R) with sex
and total length at release as fixed predictors. For Yuma Cove backwater, release year was
included as an additional fixed predictor. For Imperial Ponds, pond (location) was included as
a fixed predictor with all XYTE ponds modeled together and all GIEL ponds modeled together. A
fish was classified as a survivor (1) if it was detected via PIT scanner at least once ≥90 days
after release; otherwise it was classified as dead (0).

Known survivor analysis used continuous PIT scanning data to track the minimum number of
PIT-tagged fish alive at each backwater on any given date. Fish were counted as alive on a given
date if they were released prior to that date and had at least one scanner detection ≥90 days
after release. If a fish was never re-detected after its release period, it was excluded from the
known survivor population.

Annual population estimates for PIT-tagged adults were computed using the modified Petersen
mark-recapture formula (Ricker 1975). The marking period was January through February, and the
capture period was October through April of each census year. Fish were eligible for inclusion
only after a minimum of 6 months (marking) or 15 months (capture) post-tagging. Census years
with fewer than 3 recaptures were excluded.

Reproductive contributions of stocked adults were assessed using parentage assignment data
provided by Dowling and Turner (this report). Known survivors alive as of February 15 each year
(or February 1 for IP XYTE) were matched to offspring genetic records. Proportion of the known
survivor population producing detectable larvae or juvenile recruits was calculated for each
sex and stocking cohort.
```

---

## Genetics Code Pattern (Reproductive Contribution)

When building the February/April survivor × offspring join, use integer year comparison — **not** :

```r
# CORRECT
mutate(Year = year(Date), TagYear = year(first_date)) |>
filter(Year != TagYear)   # both integers; excludes stocking-year rows

# WRONG — as.factor() causes comparison against internal codes, not year labels
mutate(Year = year(Date), TagYear = as.factor(year(first_date))) |>
filter(Year != TagYear)   # silently fails in R 4.x; never excludes stocking year
```

Also ensure the  of larval and recruit parents uses the correct suffix order (x = Larval, y = Recruit):

```r
full_join(LarvalParents, RecruitParents,
          by = c("Year", "PIT"),
          suffix = c("_Larvae", "_Recruits"))  # x gets _Larvae, y gets _Recruits
```

---
---

## Narrative Guidance

- Be measured in language: avoid "clear," "striking," or "strong" unless data genuinely support it.
- Report numeric values with appropriate precision (proportions to 3 decimal places, counts as integers).
- Use inline R code to populate numbers automatically so the document stays reproducible.
- When sample sizes are small or results are ambiguous, say so explicitly.
- When an analysis was skipped due to insufficient data, include one sentence explaining the reason.
- Cross-reference the comparison figures explicitly in the final cross-site section.
