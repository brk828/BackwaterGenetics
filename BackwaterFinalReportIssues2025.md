# Backwater Issues for 2025 Final Report

Initial summary statistics for stocking, transfer, and capture data were sent to Jeff Lantow on September 10, 2025. These analysis are included in a Word Document entitled "Backwater Demographics September 2025" in the FinaReport2025 folder of this repository.

One issue identified early was a record of 66 fish that were harvested out of Yuma Cove backwater in 2013. These 66 razorback sucker were captured, tagged, and released into Big Bend Conservation Area, but in the NFWG database they only have a “transfer release” record. Normally with transfers from a backwater that is used as a long-term backwater and not just grow-out there will be a capture and transfer record. This is because capture events in these backwaters are used for monitoring as well as transfers and the two events must be considered separately with all transfer captures marked as "Capture - Retained." Still, these fish can be identified with their single transfer-released record and the DataWrangling script was adapted to accomodate this issue.

Overall backwater issue for backwaters, don't rely on min_date from PITIndex for any analysis as this is the first record overall and not the first record in the backwater.
