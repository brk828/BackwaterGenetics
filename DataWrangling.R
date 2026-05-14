# Data Wrangling for genetic backwater analysis IP 1-6 and Yuma Cove BW
# B. Kesner September 2025

# Set parameters
SurvivalDAL <- 90 # minimum days a fish is at large after tagging to be considered "survived"
SizeClass2 <- 350 # TL (mm) cutoff for Size Class 2
SizeClass3 <- 500 # TL (mm) cutoff for Size Class 3

# Load useful lab functions
source("LabFunctions.R")

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(readxl) # import Excel spreadsheets
packages(data.table) # faster at indexing than grouping in dplyr
packages(openxlsx) # package openxlsx is required to create to Excel files

# Download and load workspaces
download_backwater("data")
load("data/BWScanningIndex.RData")

download_nfwg("data")
load("data/NFWGAnalysis.RData")

rm(split_hourly, download_nfwg, download_genetics, download_PITindex, download_basin, download_backwater, euclid)
  
# Restrict PITindex dataframe to study backwaters
StudyBWContacts <- BWContacts %>%
  filter(LID > 1042 & LID < 1049| LID == 592) %>%
  mutate(ScanHr = hour(DateTime),
         ScanMonth = month(DateTime)) %>% 
  select(EID, PIT, PITIndex, Date, ScanHr, ScanMonth, Backwater, LID, Species = species, available_date) 

# Use data.tables to reduce contacts to summary per PIT-Date-Hour
StudyBWContactsdt <- as.data.table(StudyBWContacts)

StudyBWContacts <- StudyBWContactsdt[
  , .(count = .N), 
  by = .(EID, Backwater, LID, Species, PIT, PITIndex, Date, ScanHr, ScanMonth, available_date)] %>%
  as.data.frame() %>%
  mutate(ScanFY = as.integer(ifelse(ScanMonth > 9, year(Date)+1, year(Date))),
         ScanMonthName = format(Date, "%b"))

rm(StudyBWContactsdt)

# All study backwater records for relevant time period
StudyBWNFWG <- NFWGAnalysis %>%
  filter(collection_date > as.Date("2013-01-01") & 
           location_id == 592|
           collection_date > as.Date("2016-01-01") & 
           location_id > 1042 & 
           location_id < 1049) %>%
  mutate(CollectionYear = year(collection_date),
         CollectionMonth = month(collection_date), 
         SizeClass = ifelse(is.na(total_length), 0,
                            ifelse(total_length >= SizeClass3, 3,
                                   ifelse(total_length >= SizeClass2, 2, 1)))) %>%
  arrange(PITIndex, collection_date) %>%
  group_by(PITIndex, location) %>%
  mutate(FirstRecord = ifelse(row_number() > 1, "no", "yes")) %>%
  ungroup() %>%
  mutate(SizeClassText = case_when(
    SizeClass == 1 ~ "Size Class One",
    SizeClass == 2 ~ "Size Class Two",
    SizeClass == 3 ~ "Size Class Three",
    SizeClass == 0 ~ "No TL Provided"))

# Check for any fish contacted after study start date that were released prior to start date
StudyBWHoldovers <- NFWGAnalysis %>%
  filter(collection_date < as.Date("2013-01-01") & 
           location_id == 592|
           collection_date < as.Date("2016-01-01") & 
           location_id > 1042 & 
           location_id < 1049) %>%
  inner_join(StudyBWContacts %>% select(PITIndex, LID, Date, ScanHr, ScanMonth, count) %>%
               filter(Date > as.Date("2013-01-01") & 
                        LID == 592|
                        Date > as.Date("2016-01-01") & 
                        LID > 1042 & 
                        LID < 1049), by = "PITIndex")
  
# Study backwater captures used for growth and survival 
StudyBWCaptures <- StudyBWNFWG %>%
    filter(event == "capture") %>%
  mutate(Transfer = ifelse(disposition == "retained", 1, 0),
         Mortality = ifelse(disposition == "recovered", 1, 0))

# Fish transferred from study backwaters to mainstem will have 
# transfer record as well. Some transfers will not have a capture
# record if PIT tagged during transfer. Also include IP pond to pond
# transfers
StudyBWTransfers <- NFWGAnalysis %>%
  filter(event == "transfer") %>%
  filter(collection_date > as.Date("2013-01-01") & 
           transfer_from == "Yuma Cove backwater"|
           collection_date > as.Date("2016-01-01") & 
           grepl("^IPCA", transfer_from)|
           collection_date > as.Date("2016-01-01") & 
           grepl("^IPCA", location)) %>%
  mutate(Backwater = transfer_from, TransferLocation = location,
         TransferDate = as.Date(collection_date),
         TransferYear = year(collection_date),
         TransferMonth = month(collection_date),
         SizeClass = ifelse(is.na(total_length), 0,
                            ifelse(total_length >= SizeClass3, 3,
                                   ifelse(total_length >= SizeClass2, 2, 1))),
         SizeClassText = case_when(
                  SizeClass == 1 ~ "Size Class One",
                  SizeClass == 2 ~ "Size Class Two",
                  SizeClass == 3 ~ "Size Class Three",
                  SizeClass == 0 ~ "No TL Provided"))

# These fish have a capture retained record
# which will include backwater to backwater transfers
StudyBWRetained <- StudyBWCaptures %>%
  filter(Transfer == 1)

# Add retained data when available for transfer records for analysis of transferred fish
# This data will not include fish harvested and not PIT tagged. Those must be acquired
# from summary Excel tables.
StudyBWTransfersAnalysis <- StudyBWTransfers %>%
  select(location, PITIndex, Backwater, TransferDate, TransferLocation, TransferTL = total_length,
         TransferWeight = weight, SizeClass, SizeClassText, species) %>%
  left_join(StudyBWRetained %>%
               group_by(PITIndex) %>%
               slice_max(order_by = collection_date, n = 1, with_ties = FALSE) %>%
               ungroup() %>%
               select(PITIndex, CaptureDisposition = disposition, CaptureLocation = location, CaptureDate = collection_date, 
                      CaptureTL = total_length, Captureweight = weight, TaggedAtCapture = pit1_new), 
             by = "PITIndex") %>%
  mutate(Transfer = 1, 
         Mortality = 0,
         CaptureDate = as.Date(CaptureDate),
         TransferYear = year(TransferDate), 
         TransferMonth = month(TransferDate),
         TransferDays = as.integer(TransferDate - CaptureDate),
         TLDiff = as.integer(CaptureTL - TransferTL),
         TaggedAtCapture = ifelse(is.na(TaggedAtCapture), "yes", TaggedAtCapture)) %>%
  filter(TransferDays < 4|is.na(TransferDays)) # in case a fish is transferred more than once


if(nrow(StudyBWTransfersAnalysis) != n_distinct(StudyBWTransfersAnalysis$PITIndex)) {
  warning("Duplicates in Transfer analysis data frame.")
}

# Transfers without a capture will need to be manually added to BWCaptures data
StudyBWTransfersNoCapture <- StudyBWTransfersAnalysis %>% 
  filter(is.na(CaptureLocation)) %>%
  mutate(Contacts = 1,
         FirstRecord = 1) %>%
  select(Backwater, CollectionYear = TransferYear, CollectionMonth = TransferMonth,
         PITIndex, Contacts, FirstRecord, SizeClass, total_length = TransferTL, 
         Transfer, Mortality, SizeClassText, species)

# Ensure only one record per BW capture records before binding to Transfer records without capture
# record for a complete tagged capture dataframe
StudyBWCaptureUniques <- StudyBWCaptures %>%
  arrange(location, CollectionYear, CollectionMonth, PITIndex, collection_date) %>%
  group_by(Backwater, CollectionYear, CollectionMonth, PITIndex, species) %>%
  summarise(Contacts = n(), 
            FirstRecord = max(ifelse(FirstRecord == "yes", 1, 0)),
            SizeClass = max(SizeClass),
            total_length = max(total_length),
            Mortality = max(Mortality)) %>%
  ungroup() %>%
  left_join(StudyBWTransfersAnalysis %>%
              select(TransferYear, TransferMonth, PITIndex, Transfer), 
            by = c("CollectionYear" = "TransferYear",
                   "CollectionMonth" = "TransferMonth",
                   "PITIndex" = "PITIndex")) %>%
  mutate(Transfer = if_else(is.na(Transfer), 0, 1),
         SizeClassText = case_when(
           SizeClass == 1 ~ "Size Class One",
           SizeClass == 2 ~ "Size Class Two",
           SizeClass == 3 ~ "Size Class Three",
           SizeClass == 0 ~ "No TL Provided")) %>%
  rbind(StudyBWTransfersNoCapture)

# cleanup
rm(StudyBWTransfersNoCapture, StudyBWTransfers, StudyBWRetained)

# One pond to pond transfer 003C06F28D
# No long issue as First record is per backwater
StudyBWNFWGTaggingN <- StudyBWNFWG %>%
  filter(pit1_new == "yes", FirstRecord == "no")

if(nrow(StudyBWNFWGTaggingN) > 0){
  warning("NFWG data for study backwaters includes a new pit1 tag that is not the
          first record. Please see StudyBWNFGTaggingN for more information")
}
# Any first record that doesn't have a pit1_new "yes" should have a tagging_date
# prior to first record collection date. Except for pond to pond transfer if 
# if tagging and capture transfer happened on same day.
FirstRecordNotTagging <- StudyBWNFWG %>%
  filter(FirstRecord == "yes", collection_date > as.Date("2013-01-01"),
         pit1_new == "no", tagging_date >= collection_date)

StudyBWEffort <- BWEffort %>%
  filter(LID > 1042 & LID < 1049| LID == 592) %>%
  mutate(ScanMonth = month(MidDate),
         ScanMonthName = format(MidDate, "%b"),
         DeployedHrs = as.numeric(difftime(Retrieve, Deploy, units = "hours"))) %>%
  select(EID, Location, Deploy, Retrieve, Issue, UnitType, MidDate, 
         Comments, ScanTimeHrs, DeployedHrs, ScanMonth, ScanMonthName) %>%
  left_join(BWContacts %>% 
              select(EID, PIT, Date) %>% 
              group_by(EID) %>%
              summarise(Contacts = n(), MinScan = min(Date), MaxScan = max(Date)) %>%
              ungroup(), by = "EID") %>%
  mutate(Contacts = no_na_df(Contacts)) %>%
  filter(Contacts > 0, !is.na(MinScan)) %>%
  mutate(EffectiveScanHrs = as.numeric(difftime(MaxScan, MinScan, units = "hours")),
         ScanFY = as.integer(ifelse(ScanMonth > 9, year(MidDate)+1, year(MidDate))))

EffectiveTimeLong <- StudyBWEffort %>%
  filter(EffectiveScanHrs - DeployedHrs > 24)

# Identify contact PIT records without a NFWG entry (no PITIndex)
# Ignore short-term contacts and recent releases that have not been added to database
ContactsNoNFWG <- StudyBWContacts %>%
  filter(is.na(PITIndex), Date < as.Date("2024-10-01")) %>%
  group_by(PIT) %>%
  summarise(Backwater = min(Backwater), MinScanDate = min(Date), MaxScanDate = max(Date), contacts = sum(count)) %>%
  ungroup() %>%
  mutate(DAL = as.numeric(difftime(MaxScanDate, MinScanDate, units = "days"))) %>%
  filter(DAL > SurvivalDAL)


rm(BWContacts) # cleanup

# Summary of contact information for tags with a NFWG record
ContactsSummary <- StudyBWContacts %>%
  filter(!is.na(PITIndex)) %>%
  group_by(PITIndex, Backwater) %>%
  summarise(MinScanDate = min(Date), MaxScanDate = max(Date), contacts = sum(count)) %>%
  ungroup() %>%
  mutate(ScanDAL = as.numeric(difftime(MaxScanDate, MinScanDate, units = "days"))) 

# Fish hopping ponds or improperly allocated scan data records
ContactSummaryDuplicateLocations <- ContactsSummary %>%
  group_by(PITIndex) %>%
  summarise(count = n(), FirstLocation = min(Backwater), SecondLocation = max(Backwater)) %>%
  ungroup() %>%
  filter(count > 1)

if(nrow(ContactSummaryDuplicateLocations) > 0) {
  ContactSummaryDuplicateLocations <- ContactSummaryDuplicateLocations %>%
    left_join(StudyBWContacts %>%
                select(PITIndex, FirstEID = EID, Backwater, FirstDate = Date), 
              by = c("PITIndex", "FirstLocation" = "Backwater")) %>%
    group_by(PITIndex, FirstLocation, SecondLocation) %>%
    summarise(FirstContacts = n(), FirstDate = min(FirstDate), FirstEID = min(FirstEID)) %>%
    ungroup() %>%
    left_join(StudyBWContacts %>% select(PITIndex, SecondEID = EID, Backwater, SecondDate = Date), 
              by = c("PITIndex", "SecondLocation" = "Backwater")) %>%
    group_by(PITIndex, FirstLocation, SecondLocation, FirstEID, FirstContacts, FirstDate) %>%
    summarise(SecondContacts = n(), SecondDate = min(SecondDate), SecondEID = min(SecondEID)) %>%
    ungroup()
} else{rm(ContactSummaryDuplicateLocations)}


# Dataframe of first record fish in backwater adding contact data
# Cannot rely on pit1_new alone as some fish were tagged in hatchery
# Add summary contact data and calculate tagging DAL
StudyBWAnalysis <- StudyBWNFWG %>% 
  filter(FirstRecord == "yes") %>%
  select(first_date = collection_date, location, disposition, event, fin_clip, primary_method,
         species, PITIndex, sex, total_length, location_id) %>%
  left_join(ContactsSummary, by = c("PITIndex", "location" = "Backwater")) %>%
  mutate(release_month = month(first_date),
         release_year = year(first_date),
         MaxDAL = ifelse(!is.na(MaxScanDate), 
                as.numeric(difftime(MaxScanDate, first_date, units = "days")),
                0), 
         contacts = no_na_df(contacts)) %>%
  select(-ScanDAL) %>%
  mutate(Survived = ifelse(MaxDAL > SurvivalDAL, 1, 0),
         SurvivedFY24 = ifelse(!is.na(MaxScanDate) & 
                                 MaxScanDate > as.Date("2024-09-30"), 1, 0),
         SurvivedFY25 = ifelse(!is.na(MaxScanDate) & 
                                 MaxScanDate > as.Date("2025-09-30"), 1, 0),
         MaxScanDate = if_else(is.na(MaxScanDate), first_date, MaxScanDate)) %>%
  left_join(StudyBWTransfersAnalysis %>% select(PITIndex, Transfer, TransferDate), by = "PITIndex") %>%
  mutate(Transfer = ifelse(is.na(Transfer), 0, 1))

StudyBWGrowth <- StudyBWAnalysis %>%
  inner_join(StudyBWCaptures %>%
               filter(FirstRecord == "no", !is.na(total_length)) %>%
               select(PITIndex, RecaptureDate = collection_date, RecaptureTL = total_length),
             by = "PITIndex") %>%
  select(location, PITIndex, FirstDate = first_date, disposition, species, PITIndex, sex, FirstTL = total_length,
         RecaptureDate, RecaptureTL) %>%
  mutate(DAL = as.integer(RecaptureDate - FirstDate),
         YAL = DAL/365,
         DeltaTL = RecaptureTL - FirstTL,
         Growth_mm_year = DeltaTL/YAL) %>%
  filter(YAL >=1)

# Look for Tags with an NFWG record but no tagging record
ContactNoTagging <- ContactsSummary %>%
  anti_join(StudyBWAnalysis, by = "PITIndex") %>%
  filter(ScanDAL > SurvivalDAL)

# Create workbook for contacts with NO PITIndex
wb <- createWorkbook() # creates object to hold workbook sheets
addWorksheet(wb, "ContactsNoNFWG") # add worksheet
writeData(wb, "ContactsNoNFWG", ContactsNoNFWG) # write dataframe

addWorksheet(wb, "BadScanUploads") # add worksheet
writeData(wb, "BadScanUploads", EffectiveTimeLong) # write dataframe

addWorksheet(wb, "PondHoppers") # add worksheet
writeData(wb, "PondHoppers", ContactSummaryDuplicateLocations) # write dataframe

addWorksheet(wb, "ScannedWithoutTagging") # add worksheet
writeData(wb, "ScannedWithoutTagging", ContactNoTagging) # write dataframe

saveWorkbook(wb, paste0("output/BackwaterIssues",
                        format(Sys.time(), "%Y%m%d"), ".xlsx"), overwrite = TRUE)

# Summarize all tagged fish
# Need to remove transferred fish from survival analysis
BackwaterSurvivalSummary <- StudyBWAnalysis %>% # remove transferred fish from survival
  group_by(species, location, release_year, release_month, event) %>%
  summarise(count = n(), transferred = sum(Transfer), meanTL = as.integer(mean(total_length)), 
            minTL = min(total_length), maxTL = max(total_length),
            survivedDAL = sum(Survived), 
            survivedFY24 = sum(SurvivedFY24),
            survivedFY25 = sum(SurvivedFY25)) %>%
  ungroup() %>%
  mutate(PropSurvivedDAL = round(as.numeric(survivedDAL/count), 3),
         PorpSurvivedFY24 = round(as.numeric(survivedFY24/(count-transferred)), 3),
         PorpSurvivedFY25 = round(as.numeric(survivedFY25/(count-transferred)), 3)) %>%
  arrange(location, release_year, release_month) %>%
  filter(release_year < year(Sys.Date()))

# Summary for stockings only
StockingBackwaterSummary <- BackwaterSurvivalSummary %>%
  filter(event == "stocking") 

# Summary for captures only
CaptureBackwaterSummary <- BackwaterSurvivalSummary %>%
  filter(event == "capture") 

TaggedBWTransferSummary <- StudyBWTransfersAnalysis %>%
  group_by(species, Backwater, TransferYear, TransferMonth, TransferLocation, TaggedAtCapture) %>%
  summarise(UniqueFish = n_distinct(PITIndex),
            Mortalities = sum(Mortality),
            minTL = min(TransferTL),
            maxTL = max(TransferTL)) %>%
  ungroup()

# Summary to compare with Reclamation records
StudyBWCaptureSummary <- StudyBWCaptureUniques %>%
  filter(Transfer == 0) %>%
  mutate(TaggedAtCapture = ifelse(FirstRecord == 1, "yes", "no")) %>%
  group_by(species, Backwater, CaptureYear = CollectionYear, CaptureMonth = CollectionMonth, TaggedAtCapture) %>%
  summarise(UniqueFish = n_distinct(PITIndex),
            Mortalities = sum(Mortality),
            minTL = min(total_length),
            maxTL = max(total_length)) %>%
  ungroup() 

ScanningBackwaterSummary <- StudyBWContacts %>%
  arrange(Species, Backwater, ScanFY) %>%
  group_by(Species, Backwater, ScanFY) %>%
  summarise(Contacts = n(), Uniques = n_distinct(PITIndex)) %>%
  ungroup() 

# Create workbook for contacts with NO PITIndex
wb <- createWorkbook() # creates object to hold workbook sheets
addWorksheet(wb, "StockingBackwaterSummary") # add worksheet
writeData(wb, "StockingBackwaterSummary", StockingBackwaterSummary) # write dataframe

addWorksheet(wb, "CaptureBackwaterSummary") # add worksheet
writeData(wb, "CaptureBackwaterSummary", CaptureBackwaterSummary) # write dataframe

addWorksheet(wb, "TaggedBWTransferSummary") # add worksheet
writeData(wb, "TaggedBWTransferSummary", TaggedBWTransferSummary) # write dataframe

addWorksheet(wb, "StudyBWCaptureSummary") # add worksheet
writeData(wb, "StudyBWCaptureSummary", StudyBWCaptureSummary) # write dataframe

addWorksheet(wb, "ScanningBackwaterSummary") # add worksheet
writeData(wb, "ScanningBackwaterSummary", ScanningBackwaterSummary) # write dataframe

saveWorkbook(wb, paste0("output/BackwaterSummary",
                        format(Sys.time(), "%Y%m%d"), ".xlsx"), overwrite = TRUE)

save(StudyBWNFWG, StudyBWAnalysis, StudyBWEffort, StudyBWContacts, SurvivalDAL, 
     SizeClass2, SizeClass3, no_na, no_na_df, packages, StudyBWCaptureUniques, 
     StudyBWTransfersAnalysis, StudyBWGrowth, file = "data/ReportingData.RData")

