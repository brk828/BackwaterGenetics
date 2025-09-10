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

# Load BWScanning data workspace or downlod and load if more than 7 days old
if(file.exists("data/BWScanningIndex.RData")){
  data_info <- file.info("data/BWScanningIndex.RData")
  data_date <- as.Date(data_info$mtime)
  if(data_date>Sys.Date() - 7){
    load("data/BWScanningIndex.RData")
  } else {
    download_backwater("data")
    load("data/BWScanningIndex.RData")
  }
} else {
  download_backwater("data")
  load("data/BWScanningIndex.RData")
}

# Load NFWG data workspace or downlod and load if more than 7 days old
if(file.exists("data/NFWGAnalysis.RData")){
  data_info <- file.info("data/NFWGAnalysis.RData")
  data_date <- as.Date(data_info$mtime)
  if(data_date>Sys.Date() - 7){
    load("data/NFWGAnalysis.RData")
  } else {
    download_nfwg("data")
    load("data/NFWGAnalysis.RData")
  }
} else {
  download_nfwg("data")
  load("data/NFWGAnalysis.RData")
}

rm(split_hourly, download_nfwg, download_basin, download_backwater)
  
# Restrict PITindex dataframe to study backwaters
StudyBWContacts <- BWContacts %>%
  filter(LID > 1042 & LID < 1049| LID == 592) %>%
  mutate(ScanHr = hour(DateTime),
         ScanMonth = month(DateTime)) %>% 
  select(EID, PIT, PITIndex, Date, ScanHr, ScanMonth, Location, Species = species, tagging_date) 

# Use data.tables to reduce contacts to summary per PIT-Date-Hour
StudyBWContactsdt <- as.data.table(StudyBWContacts)

StudyBWContacts <- StudyBWContactsdt[
  , .(count = .N), 
  by = .(EID, Location, Species, PIT, PITIndex, Date, ScanHr, ScanMonth, tagging_date)] %>%
  as.data.frame() %>%
  mutate(ScanFY = as.integer(ifelse(ScanMonth > 9, year(Date)+1, year(Date))),
         ScanMonthName = format(Date, "%b"))

rm(StudyBWContactsdt)

StudyBWNFWG <- NFWGAnalysis %>%
  filter(location_id > 1042 & location_id < 1049| location_id == 592) 

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
  summarise(Backwater = min(Location), MinScanDate = min(Date), MaxScanDate = max(Date), contacts = sum(count)) %>%
  ungroup() %>%
  mutate(DAL = as.numeric(difftime(MaxScanDate, MinScanDate, units = "days"))) %>%
  filter(DAL > SurvivalDAL)


rm(BWContacts, NFWGAnalysis) # cleanup

# Summary of contact information for tags with a NFWG record
ContactsSummary <- StudyBWContacts %>%
  filter(!is.na(PITIndex)) %>%
  group_by(PITIndex, Location) %>%
  summarise(MinScanDate = min(Date), MaxScanDate = max(Date), contacts = sum(count)) %>%
  ungroup() %>%
  mutate(ScanDAL = as.numeric(difftime(MaxScanDate, MinScanDate, units = "days"))) 

# Fish hopping ponds or improperly allocated scan data records
ContactSummaryDuplicateLocations <- ContactsSummary %>%
  group_by(PITIndex) %>%
  summarise(count = n(), FirstLocation = min(Location), SecondLocation = max(Location)) %>%
  ungroup() %>%
  filter(count > 1) %>%
  left_join(StudyBWContacts %>%
              select(PITIndex, FirstEID = EID, Location, FirstDate = Date), 
            by = c("PITIndex", "FirstLocation" = "Location")) %>%
  group_by(PITIndex, FirstLocation, SecondLocation) %>%
  summarise(FirstContacts = n(), FirstDate = min(FirstDate), FirstEID = min(FirstEID)) %>%
  ungroup() %>%
  left_join(StudyBWContacts %>% select(PITIndex, SecondEID = EID, Location, SecondDate = Date), 
            by = c("PITIndex", "SecondLocation" = "Location")) %>%
  group_by(PITIndex, FirstLocation, SecondLocation, FirstEID, FirstContacts, FirstDate) %>%
  summarise(SecondContacts = n(), SecondDate = min(SecondDate), SecondEID = min(SecondEID)) %>%
  ungroup()

# Dataframe of first tagging (capture tagged or stocked)
# Cannot rely on pit1_new alone as some fish were tagged in hatchery
# Add summary contact data and calculate tagging DAL
StudyBWNFWGTagging <- StudyBWNFWG %>% 
  filter(pit1_new == "yes"| event == "stocking") %>%
  select(collection_date, location, disposition, event, fin_clip, primary_method,
         species, PITIndex, sex, total_length, tagging_date, location_id) %>%
  left_join(ContactsSummary, by = c("PITIndex", "location" = "Location")) %>%
  mutate(MaxDAL = ifelse(!is.na(MaxScanDate), 
                as.numeric(difftime(MaxScanDate, collection_date, units = "days")),
                0), 
         contacts = no_na_df(contacts)) %>%
  select(-ScanDAL) %>%
  mutate(Survived = ifelse(MaxDAL > SurvivalDAL, 1, 0),
         SurvivedFY24 = ifelse(!is.na(MaxScanDate) & 
                                 MaxScanDate > as.Date("2024-09-30"), 1, 0))

# Check for duplicates, there should be none
TaggingDuplicates <- StudyBWNFWGTagging %>%
  group_by(PITIndex) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(count >1)

# Look for Tags with an NFWG record but no tagging record
ContactNoTagging <- ContactsSummary %>%
  anti_join(StudyBWNFWGTagging, by = "PITIndex") %>%
  filter(ScanDAL > SurvivalDAL)

# Fish holdovers scanned during study period but released prior
StudyBWNFWGTaggingHoldover <- StudyBWNFWGTagging %>%
  filter(collection_date < as.Date("2013-01-01"),
         MaxScanDate > as.Date("2013-01-01"))

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

addWorksheet(wb, "Holdovers") # add worksheet
writeData(wb, "Holdovers", StudyBWNFWGTaggingHoldover) # write dataframe

saveWorkbook(wb, paste0("output/BWIssues",
                        format(Sys.time(), "%Y%m%d"), ".xlsx"), overwrite = TRUE)


# filter out older tagging records
StudyBWNFWGTagging <- StudyBWNFWGTagging %>%
  filter(collection_date > as.Date("2013-01-01") & 
           location_id == 592|
           collection_date > as.Date("2016-01-01") & 
           location_id > 1042 & 
           location_id < 1049)

# Summarize all tagged fish
BackwaterSummary <- StudyBWNFWGTagging %>%
  group_by(species, location, collection_date, event, disposition, sex) %>%
  summarise(count = n(), meanTL = as.integer(mean(total_length)), 
            minTL = min(total_length), maxTL = max(total_length),
            survivedDAL = sum(Survived), survivedFY24 = sum(SurvivedFY24)) %>%
  ungroup() %>%
  mutate(PropSurvivedDAL = round(as.numeric(survivedDAL/count), 3))

# Summary for stockings only
StockingBackwaterSummary <- BackwaterSummary %>%
  filter(event == "stocking", collection_date < as.Date("2024-10-01")) %>%
  rename(release_date = collection_date)

ScanningBackwaterSummary <- StudyBWContacts %>%
  arrange(Location, ScanFY) %>%
  group_by(Location, ScanFY) %>%
  summarise(Contacts = n(), Uniques = n_distinct(PITIndex)) %>%
  ungroup() 

# Create workbook for contacts with NO PITIndex
wb <- createWorkbook() # creates object to hold workbook sheets
addWorksheet(wb, "StockingBackwaterSummary") # add worksheet
writeData(wb, "StockingBackwaterSummary", StockingBackwaterSummary) # write dataframe

addWorksheet(wb, "ScanningBackwaterSummary") # add worksheet
writeData(wb, "ScanningBackwaterSummary", ScanningBackwaterSummary) # write dataframe

saveWorkbook(wb, paste0("output/StockingBackwaterSummary",
                        format(Sys.time(), "%Y%m%d"), ".xlsx"), overwrite = TRUE)

save(StudyBWNFWG, StudyBWNFWGTagging,StudyBWEffort, StudyBWContacts, SurvivalDAL, 
     SizeClass2, SizeClass3, no_na, no_na_df, packages,
     file = "data/ReportingData.RData")

