# BackwaterGenetics Project — Agent Memory

## Project Overview
Fisheries population monitoring for PIT-tagged Razorback Sucker (XYTE) and Bonytail (GIEL) in two sets of backwater habitats along the lower Colorado River: Yuma Cove Backwater and Imperial Ponds California (IPCA 1–6). Analyses cover post-stocking survival, minimum known survivor counts, annual mark-recapture population estimates, and reproductive contribution from parentage genetics (Dowling lab). Managed by Brian Kesner, Marsh & Associates.

---

## Species Codes
| Code | Common Name |
|------|-------------|
| XYTE | Razorback Sucker (RASU) |
| GIEL | Bonytail (BONY) |

---

## Study Locations
| Location Name | LID | Species |
|---|---|---|
| Yuma Cove backwater | 592 | XYTE |
| IPCA (Pond 1) | 1043 | XYTE |
| IPCA (Pond 2) | 1044 | GIEL |
| IPCA (Pond 3) | 1045 | XYTE |
| IPCA (Pond 4) | 1046 | XYTE |
| IPCA (Pond 5) | 1047 | GIEL |
| IPCA (Pond 6) | 1048 | GIEL |

> **Note:** Pond assignments verified from source data (July 2026); Ponds 1, 3, 4 = XYTE; Ponds 2, 5, 6 = GIEL.

---

## Key Analysis Parameters (set in DataWrangling.R)
- `SurvivalDAL = 90` — minimum days at large post-release to be classified a "survivor"
- `SizeClass2 = 350` mm TL
- `SizeClass3 = 500` mm TL
- Fiscal year (FY): October 1 – September 30; `ScanFY = year + 1` when `month > 9`
- Yuma Cove data filtered from `2013-01-01`; IPCA ponds from `2016-01-01`

---

## Analysis Pipeline (run in order)
| Step | Script | Output |
|---|---|---|
| 1 | `DataWrangling.R` | `data/ReportingData.RData` (includes all known survival objects) |
| — | `KnownSurvival.R` | *(Deprecated — wrangling moved to DataWrangling.R; script retained for standalone PNG generation only)* |
| 3 | `PostStockingSurvival.R` | 3 PNG survival curve figures |
| 4 | `BWMarkRecaptureEstimates.R` | 3 PNG estimate figures |
| 5 | `BackwaterGrowth.R` | Growth summaries (no saved PNGs currently) |
| 6 | `GeneticsAndSurvival.R` | 2 PNG offspring contribution figures |
| 7 | `BackwaterMCR.R` | CJS model results (Pond 1 only; requires MARK software) |
| 8 | `PopulationMonitoring.qmd` | Full 7-location report (HTML or DOCX); regenerates all PNGs and adds 2 cross-site comparison figures. Requires steps 1-2 first. Does **not** source scripts 3-6; models and genetics run inline. |
| — | `ScanningEffectiveness.qmd` | Standalone PIT scanner effectiveness report (HTML). Requires step 1 only. Covers effort reliability filtering, monthly detection probability, seasonal patterns, per-sub-effort detection rates, and a diagnostic section identifying tagged fish with zero passive antenna contacts. |

---

## Key Data Files
| File | Contents |
|---|---|
| `data/BWScanningIndex.RData` | PIT scanner contacts (`BWContacts`, `BWEffort`); downloaded from ncreased.net |
| `data/NFWGAnalysis.RData` | NFWG database captures, stockings, transfers, mortality; downloaded from ncreased.net |
| `data/ReportingData.RData` | All processed analysis-ready objects (created by DataWrangling.R) |
| `data/ReportingData.RData` | All processed analysis-ready objects including known survival datasets (created by DataWrangling.R) |
| `data/Yuma_larvae_mothers_repro_counts_PIT.csv` | Yuma larval parentage – females |
| `data/Yuma_larvae_fathers_repro_counts_PIT.csv` | Yuma larval parentage – males |
| `data/Yuma_recruits_mothers_repro_counts_PIT.csv` | Yuma recruit parentage – females |
| `data/Yuma_recruits_fathers_repro_counts_PIT.csv` | Yuma recruit parentage – males |
| `data/IP_larvae_mothers_repro_counts_PIT.csv` | IP larval parentage – females |
| `data/IP_larvae_fathers_repro_counts_PIT.csv` | IP larval parentage – males |
| `data/IP_recruits_mothers_repro_counts_PIT.csv` | IP recruit parentage – females |
| `data/IP_recruits_fathers_repro_counts_PIT.csv` | IP recruit parentage – males |

---

## Key R Objects (loaded via `data/ReportingData.RData`)
| Object | Description |
|---|---|
| `StudyBWAnalysis` | One row per PIT-tagged fish (first record in each backwater). Contains: `location`, `location_id`, `species`, `sex`, `total_length`, `first_date`, `event` (stocking/capture), `release_year`, `release_month`, `MaxScanDate`, `MaxDAL`, `Survived` (0/1), `SurvivedFY24`, `SurvivedFY25`, `Transfer` (0/1) |
| `StudyBWNFWG` | All NFWG records for study backwaters. Contains `FirstRecord` flag, `event`, `collection_date`, `PITIndex`, `disposition`, `species`, `sex`, `total_length` |
| `StudyBWContacts` | PIT scanner contacts summarized per PIT-Date-Hour. Contains `Backwater`, `LID`, `Species`, `PITIndex`, `Date`, `ScanHr`, `ScanMonth`, `ScanFY`, `available_date` |
| `StudyBWTransfersAnalysis` | Transfer records with associated capture data |
| `StudyBWCaptureUniques` | Unique capture records per fish per month per backwater |
| `StudyBWGrowth` | Growth data for fish with ≥1 year between tagging and recapture |
| `SurvivalDAL`, `SizeClass2`, `SizeClass3` | Analysis parameters |

## Key R Objects (loaded via `data/ReportingData.RData` — Known Survival additions)
| Object | Description |
|---|---|
| `KnownSurvivalYuma` | Daily known survivor records for Yuma Cove XYTE. Columns: `Date`, `PITIndex`, `first_date`, `sex`, `total_length`, `event`, `MaxScanDate`, `MaxDAL` |
| `KnownSurvivalIPXYTE` | Daily known survivor records for IP XYTE ponds. Adds `location` column |
| `KnownSurvivalIPGIEL` | Daily known survivor records for IP GIEL ponds. Adds `location` column |
| `TotalCountYuma` | Daily total known survivor count (Yuma, all sex combined) |
| `TotalCountIPXYTE` | Daily total known survivor count by IP XYTE pond |
| `TotalCountIPGIEL` | Daily total known survivor count by IP GIEL pond |
| `SuspectTransfers` | Fish with NFWG transfer records flagged as suspect because PIT scanner contacts at the same backwater continued across >1 calendar month after the listed transfer date. Columns: `PITIndex`, `Backwater`, `TransferDate`, `TransferLocation`, `TransferYear`, `PostTransferMonths`, `LastPostTransferContact`. These fish are retained in the known survival analysis despite `Transfer == 1`. Also written to `output/TransfersWithContinuedScanning.csv` by the QMD. |

---

## Output PNG Figures (in `output/`)
| File | Script | Description |
|---|---|---|
| `KnownSurvivalPlotYuma.png` | PopulationMonitoring.qmd | Known survivors over time – Yuma, by sex |
| `KnownSurvivalPlotIPXYTE.png` | PopulationMonitoring.qmd | Known survivors over time – IP XYTE, faceted by pond |
| `KnownSurvivalPlotIPGIEL.png` | PopulationMonitoring.qmd | Known survivors over time – IP GIEL, faceted by pond |
| `YumaEstimatePlot.png` | BWMarkRecaptureEstimates.R | Annual mark-recapture estimates – Yuma |
| `IPXYTEEstimatePlot.png` | BWMarkRecaptureEstimates.R | Annual mark-recapture estimates – IP XYTE, faceted |
| `IPGIELEstimatePlot.png` | BWMarkRecaptureEstimates.R | Annual mark-recapture estimates – IP GIEL, faceted |
| `YumaCoveSurvival90DaysPostStocking.png` | PostStockingSurvival.R / PopulationMonitoring.qmd | Post-stocking survival curves – Yuma; exported at 6×4 in |
| `IPXYTESurvival90DaysPostStocking.png` | PostStockingSurvival.R / PopulationMonitoring.qmd | Post-stocking survival curves – IP XYTE; exported at 6×4 in |
| `IPGIEL2017Survival90DaysPostStocking.png` | PostStockingSurvival.R / PopulationMonitoring.qmd | Post-stocking survival curves – IP GIEL (2017); exported at 6×4 in |
| `YumaOffspringPlot.png` | GeneticsAndSurvival.R | Reproductive contribution by sex/year – Yuma XYTE |
| `IPCA1OffSpringPlot.png` | GeneticsAndSurvival.R | Reproductive contribution – IP Pond 1 XYTE |
| `CohortComparisonXYTEPlot.png` | PopulationMonitoring.qmd | Cross-site initial vs supplemental stocking survivors — Razorback Sucker |
| `CohortComparisonGIELPlot.png` | PopulationMonitoring.qmd | Cross-site initial vs supplemental stocking survivors — Bonytail |
| `SexSurvivalXYTEPlot.png` | PopulationMonitoring.qmd | Sex-based survival comparison — Razorback Sucker |
| `SexSurvivalGIELPlot.png` | PopulationMonitoring.qmd | Sex-based survival comparison — Bonytail |
| `TransfersWithContinuedScanning.csv` | PopulationMonitoring.qmd | Records for fish with suspect transfer entries retained in the known survival analysis |

---

## Analysis Rules and Caveats
- **Post-stocking survival model:** Requires ≥5 stocked fish with known sex and total length. Excludes Yuma 2020 (no sex data) and Yuma 2015 (no male contacts). IP GIEL model only uses 2017 stocking cohort. Sex "J" (juvenile) is recoded to "U" before analysis.
- **Mark-recapture estimates:** Computed only when recaptures R > 3. Marking period = Jan–Feb, capture period = Oct–Apr. Fish must be tagged ≥6 months before marking period.
- **Known survivors:** A fish is counted as "alive" on a date if it was tagged before that date AND was scanned at least once ≥90 DAL.
- **Genetics/offspring:** Only Yuma Cove XYTE and IP Pond 1 XYTE currently have parentage data from the Dowling lab. Skip this section for all other ponds.
- **CJS model:** Only implemented for IP Pond 1 in BackwaterMCR.R. Requires MARK software (RMark package).
- **Pond-to-pond transfers:** One known transfer (003C06F28D) has been documented; `FirstRecord` is tracked per backwater.
- **Transferred fish exclusion:** `DataWrangling.R` filters out fish with `Transfer == 1` from the known survivor groups via `SuspectTransfers`. Exception: fish in `SuspectTransfers` are retained because their post-transfer scanning record (>1 month of contacts after the listed transfer date) indicates they remained in or returned to the backwater. The cohort comparison figure in the QMD also excludes confirmed transfers via `filter(event == "stocking")` combined with the upstream exclusion in the `KnownSurvival*` objects.

---

## Utility Functions (`LabFunctions.R`)
- `packages(pkg)` — auto-installs and loads a package
- `no_na(x)` / `no_na_df(x)` — replaces NA with 0 (scalar / vector)
- `download_backwater("data")` / `download_nfwg("data")` — refresh data from ncreased.net
- `euclid(p1, q1, p2, q2)` — Euclidean distance between UTM coordinates
- `split_hourly(df, id, start, end)` — expands start/end records to hourly rows

---

## Known Code Issues
- `BackwaterMCR.R` references `StudyBWNFWGAnalysis` which does not exist in ReportingData.RData; likely should be `StudyBWAnalysis`. Also references column `Location` which should be `Backwater`. This script is incomplete/experimental.

## Post-Edit File Format Checks

After editing any `.qmd` or `.R` file in this project, verify:

1. **No blank lines within code chunks** — Blank lines between every line of code is a known artifact in this project. After making edits, check that no blank line between every code line was inadvertently introduced. Quick diagnostic in R:
   ```r
   lines <- readLines("PopulationMonitoring.qmd")  # or other file
   in_chunk <- FALSE
   blank_in_chunk <- 0
   for (line in lines) {
     if (grepl("^```\\{r", line)) in_chunk <- TRUE
     else if (in_chunk && grepl("^```\\s*$", line)) in_chunk <- FALSE
     else if (in_chunk && trimws(line) == "") blank_in_chunk <- blank_in_chunk + 1
   }
   cat("Blank lines inside code chunks:", blank_in_chunk, "\n")
   ```
   If the count is high (more than a handful of intentional separators), strip them with the loop in `Known Code Issues` below or re-apply the cleanup used on 2026-07-24.

2. **YAML front matter has no blank lines** — Blank lines inside nested YAML blocks (`format:`, `execute:`) break parsing and cause the front matter to print as raw text. The YAML in `.qmd` files must be compact with no blank lines between keys.

3. **Global chunk options present** — `PopulationMonitoring.qmd` and `ScanningEffectiveness.qmd` rely on `execute: echo: false / warning: false / message: false` in the YAML. Confirm these are intact after edits.

---

## Fixed Code Issues (resolved)
- `DataWrangling.R` (recent, formerly `KnownSurvival.R`): Added `SuspectTransfers` detection — fish flagged as transferred but with >1 month of post-transfer PIT contacts at the same backwater are retained in the known survival analysis. Group filters use `Transfer == 0 | PITIndex %in% SuspectTransfers$PITIndex`.
- `PopulationMonitoring.qmd` cohort comparison (recent): Cohort comparison now correctly excludes transferred fish because `KnownSurvival*` objects upstream exclude them; the existing `filter(event == "stocking")` additionally ensures captured-and-tagged fish never appear in either cohort line.
- `PopulationMonitoring.qmd` (recent): Added `suspect-transfers` chunk that writes `output/TransfersWithContinuedScanning.csv` and emits a bold note paragraph listing flagged tags when `SuspectTransfers` is non-empty.
- `ScanningEffectiveness.qmd` (recent): Added "Tagged Fish Never Contacted by Passive Antennas" section with a summary table of never-scanned fish by location and a long-absence table listing any fish physically re-encountered ≥365 days after tagging with zero intervening passive contacts.
- `GeneticsAndSurvival.R` / `PopulationMonitoring.qmd` (recent): Cohort comparison `filter(event == "stocking")` confirmed present and correct; no captured-and-tagged fish flow into the Supplemental/Restock cohort lines.
- `KnownSurvival.R` line 144: `ax is.text` typo fixed to `axis.text`. *(Script is deprecated from pipeline; plotting now done in PopulationMonitoring.qmd.)*
- `DataWrangling.R` line 313: `writeData(..., ContactSummaryDuplicateLocations)` would error when no pond-hoppers exist (variable was `rm()`'d). Wrapped in `if (exists(...))` guard.
- `GeneticsAndSurvival.R` lines 133/155/175: `TagYear = as.factor(year(first_date))` caused numeric-vs-factor comparison in `filter(Year != TagYear)` to fail silently — stocking-year records were never excluded, leading to double-counting in `PropOffspring` for stocking years. Fixed by removing `as.factor()`.
- `GeneticsAndSurvival.R` lines 65/99: `suffix = c("_Recruits", "_Larvae")` was backwards on `full_join(LarvalParents, RecruitParents)`. Results were numerically correct due to symmetric logic, but column names were misleading. Fixed to `c("_Larvae", "_Recruits")`.
- `PopulationMonitoring.qmd` (Yuma offspring chunk line ~420, IP XYTE offspring chunk line ~719): same `TagYear = as.factor()` bug copied from original script. Fixed to `TagYear = year(first_date)` (integer comparison).
- `DataWrangling.R` / `KnownSurvival.R` (recent): All known survival data wrangling (`SuspectTransfers`, `KnownSurvival*`, `TotalCount*`, Additions, Removals) moved from `KnownSurvival.R` into `DataWrangling.R` and saved in `ReportingData.RData`. `KnownSurvivalCounts.RData` is no longer generated or loaded by any script. `KnownSurvival.R` is retained but deprecated from the pipeline.
