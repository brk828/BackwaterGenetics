# BackwaterGenetics Project — Agent Memory

> **AGENTS.md File Maintenance:** This file must use **LF-only line endings** (Unix-style). Never write CRLF () line endings to this file. When editing via tools that write raw bytes (e.g., Python scripts, bash), always open/write in text mode or explicitly convert  → . The Error 0x80070006: The handle is invalid. tool may silently introduce CRLF if the file already contains them — if you must use the edit tool and then verify, run: `python3 -c "print(open('AGENTS.md','rb').read().count(b'
'))"` and strip any introduced CRLF immediately.


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
| 8 | `PopulationMonitoring.qmd` | Full 7-location report (HTML or DOCX); regenerates all PNGs and adds 2 cross-site comparison figures. Requires steps 1-2 first. Does **not** source scripts 3-6; models and genetics run inline. Includes an **Introduction** section (species conservation status, OCH concept, site histories, BONY pond sizes) and a **Data Sources and Monitoring Infrastructure** subsection in Methods. Narrative paragraphs for each location draw on `2025Report_IP_PopulationMonitoring.md`, `Marsh_et_al_2024_Summary.md`, and related document summaries. |
| — | `ScanningEffectiveness.qmd` | Standalone PIT scanner effectiveness report (HTML). Requires step 1 only. Covers effort reliability filtering, monthly detection probability, seasonal patterns, per-sub-effort detection rates, and a diagnostic section identifying tagged fish with zero passive antenna contacts. |

---

## Reference Document Lookup Guide

> **Standing instruction:** Whenever a document in the `documents/` folder is converted to a markdown summary at the user's request, automatically: (1) add a row to the **Key Reference Documents** table below with the PDF/DOCX filename, summary filename, and a brief description, and (2) add a row to this **Reference Document Lookup Guide** table with one or more trigger topics and the summary filename. Both updates must be made together before ending the response.

When a request touches one of the topics below, **read the corresponding markdown file before responding**. These summaries contain the key facts, formulas, and citations needed to give accurate answers without having to re-parse the original PDFs or DOCX files.

| Topic / Request Type | File(s) to read |
|---|---|
| Chapman estimator formula, mark-recapture methods, population estimates, wildlife abundance estimation | `documents/PierceEtAl_ChapterAbundance_Summary.md` |
| Bootstrap confidence intervals for Chapman, `ciChapman()` behavior, CI coverage at small R | `documents/Sangnawakij_et_al_2026_Summary.md` |
| Both of the above together (population estimate methods section of report) | Both `PierceEtAl_ChapterAbundance_Summary.md` and `Sangnawakij_et_al_2026_Summary.md` |
| OCH concept, IPCA site history, BONY pond sizes, Marsh et al. 2024 citation | `documents/Marsh_et_al_2024_Summary.md` |
| YCB population narrative, sex-biased survival history, 2013 stocking cohort | `documents/2025Report_YCB_PopulationMonitoring_Summary.md` |
| IPCA population narrative, recruit counts by pond, IP RASU or BONY narrative text | `documents/2025Report_IP_PopulationMonitoring.md` |
| IPCA site design rationale, 2005 design workshop, original pond failures, six-pond layout | `documents/IP_DesignWorkshop2005_Summary.md` |
| 2014–2015 renovation, gap before 2017 stocking, well-water conversion | `documents/IPRenovationPlan2014_Summary.md` |
| FY2018 annual report, 2017 stocking results, first-year survival at IPCA | `documents/IPCA_AnnualReport2018_Summary.md` |
| YCB stocking history, harvest events, fish returned or removed by year | `documents/YumaBW_SummaryOfEvents_Summary.md` |
| GT-seq panel, CKMR methods, genetic database (lcrgenetics.net), MSB tissue archive | `documents/2025FinalReport_Summary.md` |
| 2019/2020 genetics report, microsatellite methods, AJ/DAN ephemeral ponds, early IPCA data | `documents/BackwaterFinalReport2020_Summary.md` |
| Pierce et al. full citation or page range | `documents/PierceEtAl_ChapterAbundance_Summary.md` |
| Sangnawakij et al. full citation or DOI | `documents/Sangnawakij_et_al_2026_Summary.md` |

---

## Key Reference Documents

| File | Description |
|---|---|
| `documents/2025FinalReport.docx` | Full 2025 Dowling lab genetics and demographic final report for BONY and RASU (mainstem and OCH). GT-seq panel, CKMR methods, genetic database, tissue archive. See `2025FinalReport_Summary.md`. |
| `documents/2025FinalReport_Summary.md` | Summary of `2025FinalReport.docx` — GT-seq panel results for BONY (238 loci) and RASU, CKMR abundance estimates, database (http://lcrgenetics.net), and MSB tissue archive through 2025. |
| `documents/2025Report_YCB_PopulationMonitoring.docx` | Population monitoring narrative section for Yuma Cove Backwater RASU through 2025. Covers post-stocking survival (sex-biased), known survivors, mark-recapture estimates, and reproductive contribution. See `2025Report_YCB_PopulationMonitoring_Summary.md`. |
| `documents/2025Report_YCB_PopulationMonitoring_Summary.md` | Summary of `2025Report_YCB_PopulationMonitoring.docx`. Key facts: 24M vs 80F survived 2013 stocking; 2018 recruit cohort die-off preceded management shift to harvesting; majority of current YCB population is naturally recruited. |
| `documents/2025Report_IP_PopulationMonitoring.md` | Draft population monitoring narrative for Imperial Ponds RASU (Ponds 1, 3, 4) and BONY (Ponds 2, 5, 6). Includes recruit counts by pond and year. |
| `documents/Marsh_et_al_2024_Summary.md` | Summary of Marsh et al. (2024) — the primary reference for Imperial Ponds site history, design rationale, OCH concept, and BONY pond sizes. See file for full citation. |
| `documents/backwater final report draft 16 Jan 2020.docx` | Dowling lab 2019 final report (draft Jan 2020) on RASU genetics and demographics in OCH including ephemeral ponds (AJ, DAN), YCB, and early IPCA data. Used 14-locus microsatellites. Predecessor to `2025FinalReport.docx`. See `BackwaterFinalReport2020_Summary.md`. |
| `documents/BackwaterFinalReport2020_Summary.md` | Summary of `backwater final report draft 16 Jan 2020.docx`. Covers AJ and DAN ephemeral ponds (not in current project), YCB, and early IPCA. Established sex-biased survival pattern; identified need for SNP parentage. |
| `documents/IP Design Workshop 2005.pdf` | Bureau of Reclamation design workshop report for reconstructing the DU2 Ponds at INWR into what became the IPCA. Documents why original ponds failed and six-pond layout design goals. See `IP_DesignWorkshop2005_Summary.md`. |
| `documents/IP_DesignWorkshop2005_Summary.md` | Summary of `IP Design Workshop 2005.pdf`. Foundational document for IPCA site design and rationale for river-disconnected, well-water-supplied ponds. |
| `documents/Pierce et al. - Chpt 11 Estimating Animal Abundance.pdf` | Chapter 11 from *The Wildlife Techniques Manual* (7th ed., 2012), Pierce, Lopez, and Silvy. Covers the full range of abundance estimation methods including the Chapman (1951) mark-recapture estimator. Primary reference for Chapman estimator context and formula in wildlife practice. See `PierceEtAl_ChapterAbundance_Summary.md`. |
| `documents/PierceEtAl_ChapterAbundance_Summary.md` | Summary of Pierce et al. Chapter 11. Includes the Chapman estimator formula N̂ = (M+1)(C+1)/(R+1) − 1, its variance, assumptions, and notation mapping to the backwater project. |
| `documents/Sanwnawakij et al 2026 - Two-way capture-recapture methods.pdf` | Sangnawakij et al. (2026), *Statistical Methods & Applications* 35:1–22. Primary reference supporting bootstrap confidence intervals for the Chapman estimator; compares imputed bootstrap (preferred), double bootstrap, and simple bootstrap methods. See `Sangnawakij_et_al_2026_Summary.md`. |
| `documents/Sangnawakij_et_al_2026_Summary.md` | Summary of Sangnawakij et al. (2026). Covers Chapman estimator formula and bias, three bootstrap CI methods, simulation findings (imputed bootstrap best), and relevance to `recapr::ciChapman()`. Includes note on the special case R = C where Chapman formula simplifies to M. |
| `documents/LCR MSCP 2021 - Imperial Ponds Conservation Area Annual Report 2018.pdf` | FY2018 IPCA annual report (Swatzell et al. 2020). First full monitoring year after 2017 initial stocking. Contains stocking table, FY18 population estimates (60–66% RASU survival at 21 months), and mass BONY spawning observation. See `IPCA_AnnualReport2018_Summary.md`. |
| `documents/IPCA_AnnualReport2018_Summary.md` | Summary of the FY2018 IPCA annual report. Key reference for 2017 stocking conditions and first-year survival. |
| `documents/LCR MSCP Imperial Ponds Renovation Plan 2014.pdf` | LCR MSCP plan (Finnegan 2014) for rotenone renovation of all six IPCA ponds in Dec 2014–Jan 2015 prior to 2017 stocking. See `IPRenovationPlan2014_Summary.md`. |
| `documents/IPRenovationPlan2014_Summary.md` | Summary of the 2014 renovation plan. Explains the IPCA gap before 2017 stocking and the conversion to 100% well-water supply. |
| `documents/Yuma BW summary of events.xlsx` | Reference spreadsheet of major management events at YCB 2013–2025: stockings (2013, 2014, 2015, 2020, 2025) and netting/harvest events with fish counts. See `YumaBW_SummaryOfEvents_Summary.md`. |
| `documents/YumaBW_SummaryOfEvents_Summary.md` | Summary of `Yuma BW summary of events.xlsx`. Stocking and netting events at YCB with counts of tagged/untagged fish returned and harvested. |

> **Marsh et al. (2024)** is the canonical reference for the OCH program context, sex-biased survival patterns, and Imperial Ponds history. Cite as: Marsh, P.C., T.E. Dowling, T.F. Turner, M.J. Osborne, and B.R. Kesner. 2024. Monographs of the Western North American Naturalist, Vol. 15, Article 1.

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
| `PopulationEstimates_<date>.xlsx` | PopulationMonitoring.qmd | All-site mark-recapture estimates (M, C, R, N̂, 95% bootstrap CI) exported each time the QMD is rendered. Filename includes the render date (e.g. `PopulationEstimates_07Aug2026.xlsx`); re-rendering on the same day overwrites that day's file. Sheet: "Population Estimates". |

---

## Analysis Rules and Caveats
- **Post-stocking survival model:** Requires ≥5 stocked fish with known sex and total length. Excludes Yuma 2020 (no sex data) and Yuma 2015 (no male contacts). IP GIEL model only uses 2017 stocking cohort. Sex "J" (juvenile) is recoded to "U" before analysis.
- **Mark-recapture estimates:** Computed only when recaptures R > 3. Marking period = Jan–Feb, capture period = Oct–Apr. Fish must be tagged ≥6 months before marking period. **Chapman (1951) estimator:** N̂ = (M+1)(C+1)/(R+1) − 1, the nearly unbiased modification of the Lincoln-Petersen estimator (see `PierceEtAl_ChapterAbundance_Summary.md`). **95% CI method:** Parametric bootstrap via `ciChapman(n1 = M, n2 = C, m2 = R, method = "boot", bootreps = 10000)` from the `recapr` R package (Tyers 2021). Recaptures are resampled from Binomial(n = C, p = R/C); bootstrap CIs outperform normal-approximation methods especially at low R. See `Sangnawakij_et_al_2026_Summary.md` for the theoretical basis and comparison of bootstrap methods. **Special case:** When R = C (all captured fish were marked), the Chapman formula simplifies algebraically to M, identical to the simple ratio MC/R = M — the −1 is present and correct; it is a mathematical identity. Applied consistently in `BWMarkRecaptureEstimates.R` and `PopulationMonitoring.qmd`.
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

## Line Ending Policy

All `.R` and `.qmd` files in this project **must use LF-only line endings** (Unix-style, `\n`). The files originated with Windows CRLF (`\r\n`) endings and were fully normalized to LF on 2026-08-07. The RStudio edit tool can silently reintroduce CRLF on Windows — if a string-replacement edit fails with "String not found" on any `.R` or `.qmd` file, CRLF re-contamination is the most likely cause.

**After any editing session**, run this one-liner from the project root to re-normalize all `.R` and `.qmd` files:

```python
# Run in a bash or terminal pane from C:\GIT\BackwaterGenetics
python3 -c "
import glob
for ext in ('*.R', '*.qmd'):
    for path in glob.glob(ext):
        data = open(path, 'rb').read()
        if b'\r\n' in data:
            open(path, 'wb').write(data.replace(b'\r\n', b'\n'))
            print('Normalized:', path)
"
```

To verify a single file is clean:
```python
python3 -c "data=open('FILE.R','rb').read(); print('CRLF count:', data.count(b'\r\n'))"
```

---

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
- `PopulationMonitoring.qmd` (ip-pond1-offspring chunk, line ~863): `filter(Pond == "IPCA ( Pond 1)", ...)` had an extra space before "Pond", returning 0 rows and causing `facet_wrap` to error with "Faceting variables must have at least one value." Fixed to `"IPCA (Pond 1)"`.
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
