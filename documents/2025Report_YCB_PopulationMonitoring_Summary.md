# Summary: 2025Report_YCB_PopulationMonitoring.docx

## Description
Population monitoring narrative section for **Yuma Cove Backwater (YCB) Razorback Sucker (RASU)**, covering data through 2025. This appears to be an extracted chapter from the larger 2025 project report, focused on demographic monitoring rather than genetics. Prepared by Marsh & Associates (B.R. Kesner).

---

## Content Overview

This short document (~11 paragraphs plus figure captions) covers three interrelated topics:

### Post-Stocking Survival and Detection Probability
- Detection probability (post-stocking survival) for 2013 and 2014 stockings was significantly dependent on **total length (TL), sex, and release year**.
- YCB was stocked with 100M / 100F in 2013; only 24 males survived the 90-day post-stocking period vs. 80 females — a strongly skewed sex ratio that persisted across all subsequent stocking years.
- Of 119 males stocked 2014–2020, only 2 survived the post-stocking period.
- 2020 stocking (100 fish: 21M, 79 unknown sex): no males survived; 28 unknown-sex fish survived. Overall apparent survival was higher than 2014/2015 but lower than 2013.

### Known Survivors and Recruitment
- Sudden increases in the known survivor count correspond to netting/trapping events that harvested surplus RASU and PIT-tagged naturally recruited fish.
- Most recruit cohorts experienced significant early post-release mortality (not shown in known survival plots as they didn't survive 90 days).
- Exception: ~100 RASU recruits survived the 90-day period in 2018 but disappeared in early 2019, possibly indicating **overcrowding**. Management shifted post-2018 toward harvesting (transferring) fish to reduce biomass while continuing supplemental stockings.

### Mark-Recapture Estimates
- Initial 2013 PIT-tagged population estimated at ~100 fish — consistent with known survivor analysis.
- Population exceeded 200 PIT-tagged individuals in four estimates from 2016–2019 (after recruit PIT tagging began).
- Mark-recapture estimates closely aligned with known survivor counts, validating the known survivor approach.
- The majority of the current YCB population consists of **naturally recruited fish**, not stocked adults.

### Reproductive Contribution (genetics-based)
- A majority of surviving male and female RASU stocked in 2013 contributed to offspring in the first five years post-stocking.
- Proportion contributing from the 2013 cohort dropped below 50% for females in 2018 and males in 2021, likely reflecting increasing reproductive contribution from untagged recruits.
- RASU stocked in 2014 and 2015 had near-zero contribution in their stocking year (likely sexually immature at stocking size). Both single surviving males from 2014 and 2015 cohorts contributed offspring at least once in later years.

---

## Relevance to BackwaterGenetics Project
- Narrative text and figure interpretations for the YCB section of `PopulationMonitoring.qmd` and the 2025 annual report.
- Confirms the key finding that male post-stocking survival is dramatically lower than female survival at YCB — a consistent pattern across all stocking years.
- Background for the known survivor and mark-recapture analyses in `DataWrangling.R` and `BWMarkRecaptureEstimates.R`.
- The 2018 recruit die-off and subsequent harvest-based management shift are important context for interpreting population trends.
