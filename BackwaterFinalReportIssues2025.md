---
output:
  pdf_document: default
  html_document: default
---
# Backwater Issues for 2025 Final Report

Initial summary statistics for stocking, transfer, and capture data were sent to Jeff Lantow on September 10, 2025. These analysis are included in a Word Document entitled "Backwater Demographics September 2025" in the FinaReport2025 folder of this repository.

One issue identified early was a record of 66 fish that were harvested out of Yuma Cove backwater in 2013. These 66 razorback sucker were captured, tagged, and released into Big Bend Conservation Area, but in the NFWG database they only have a “transfer release” record. Normally with transfers from a backwater that is used as a long-term backwater and not just grow-out there will be a capture and transfer record. This is because capture events in these backwaters are used for monitoring as well as transfers and the two events must be considered separately with all transfer captures marked as "Capture - Retained." Still, these fish can be identified with their single transfer-released record and the DataWrangling script was adapted to accomodate this issue.

Overall backwater issue for backwaters, don't rely on min_date from PITIndex for any analysis as this is the first record overall and not the first record in the backwater. Backwater specific min-date should be in backwater index dataframe.

Individual BW Tag issues (do not impact IP or Yuma Backwater)

424E023F68 - First captured in 2001 but 2001 record listed as recap "Y" and same trip "Y", only other record is from 2005 where it is listed as recap "N". Both records seem to indicate it escaped and was recaptured within the pond.

424F193A7D - recap "Y" and same trip "Y" for a 2002 record, but a 2003 record has recap "N".

1C2D06B06A - Raceway study fish were first record in 2009 has recap "Y" and second record in 2010 has recap "N". However, there are two 2010 records for this fish one of which has no TL or WT, the other with a TL and weight that is less than the original 2009 record. Definitely not the same fish as the 2009 record as the fish was still less than the original 2009 record in TL and weight when captured again in November 2010. Release record into Yuma Cove backwater in March 2011 has no TL or weight recorded.
