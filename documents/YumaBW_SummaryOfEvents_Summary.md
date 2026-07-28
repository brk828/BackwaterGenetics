# Summary: Yuma BW summary of events.xlsx

## Description
A reference spreadsheet (single sheet: "Yuma Cove BW") recording major management events at Yuma Cove Backwater from 2013 through 2025. Prepared by Marsh & Associates. Contains 24 rows of stocking and netting events.

---

## Columns
| Column | Description |
|--------|-------------|
| Date | Date of event |
| Location | Always "Yuma BW" |
| Event | "Stocking" or "Netting" |
| Number stocked | Adults stocked (stocking events only) |
| Previously tagged Returned to Pond | Fish with prior PIT tags returned to pond after netting |
| New tag returned to pond | Fish newly PIT-tagged and returned during the event |
| Returned to Pond untagged | Untagged fish returned to pond |
| Harvested Untagged | Untagged fish removed/harvested |
| Harvested Previously Tagged | Previously tagged fish removed/harvested |
| Harvested New Tags | Newly tagged fish removed/harvested |
| Avg. TL of YoY | Average total length of young-of-year fish captured |
| Morts | Mortalities recorded at event |
| File | Reference to associated data file |
| Notes | Event-specific notes |

---

## Event Summary

| Date | Event | Notes |
|------|-------|-------|
| 2013-02-11 | Stocking | 200 adults (initial stocking) |
| 2013-11-13/14 | Netting | 79 previous tags returned; first post-stocking sampling |
| 2014-01-02 | Netting | |
| 2014-01-29 | Stocking | 98 adults (supplemental) |
| 2015-01-29 | Stocking | 100 adults (supplemental) |
| 2015-05-14 | Netting | |
| 2015-05-28 | Netting | |
| 2016-05-03 | Netting | |
| 2016-05-18 | Netting | |
| 2016-10-13 | Netting | |
| 2017-05-09 | Netting | |
| 2017-10-12 | Netting | 36 previous tags returned |
| 2018-10-12 | Netting | |
| 2019-11-05/06 | Netting (×3) | Multiple events; population declining |
| 2020-02-12 | Stocking | 100 adults (supplemental) |
| 2020-11-03 | Netting | 56 previous tags returned |
| 2021-10-29 | Netting | 40 previous tags returned |
| 2022-11-02 | Netting | 17 previous tags returned |
| 2023-11-14 | Netting | 13 previous tags returned |
| 2024-10-22 | Netting | 21 previous tags returned |
| 2025-01-17 | Stocking | 100 adults (supplemental) |

**Stocking events:** 2013 (200), 2014 (98), 2015 (100), 2020 (100), 2025 (100)

---

## Relevance to BackwaterGenetics Project
- Cross-reference for stocking dates and counts used in `DataWrangling.R` and `StudyBWAnalysis`.
- Netting events correspond to capture events that appear in the NFWG database (`data/NFWGAnalysis.RData`) and the `StudyBWNFWG` object.
- The YoY average TL data can contextualize recruit size-at-age estimates in `BackwaterGrowth.R`.
- Harvest and return counts provide management context for interpreting known survivor trajectories.
