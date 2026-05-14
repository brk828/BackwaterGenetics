# Script to look at Genetic data from Dowling et al. and known survivors 
# to calculate proportion of known survivors that contributed to larvae and adults
# Only data is for Yuma Cove backwater and IP Pond 1 (XYTE)
# B. Kesner May 13, 2026

load("data/KnownSurvivalCounts.RData")
load("data/ReportingData.RData")

YumaStockings <- StudyBWAnalysis %>% filter(event == "stocking",
                                            location_id == 592,
                                            (release_year == 2013|
                                               release_year == 2015))

IPStockings <- StudyBWAnalysis %>% filter(event == "stocking",
                                          location_id == 1043,
                                          release_year == 2016)

rm(StudyBWGrowth, StudyBWEffort, StudyBWCaptureUniques, StudyBWTransfersAnalysis)

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(ggplot2) # Plotting
packages(tidyr) # replace_na for dataframes function plus others
packages(openxlsx) # package openxlsx is required to create to Excel files


# Bring in larval and recruit counts by PIT tag from Dowling data for Yuma and IP1 stored as csv files in
# data directory
YumaLarvalMoms <- read.csv("C:/GIT/BackwaterGenetics/data/Yuma_larvae_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "F")

YumaLarvalDads <- read.csv("C:/GIT/BackwaterGenetics/data/Yuma_larvae_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "M")

YumaRecruitsMoms <- read.csv("C:/GIT/BackwaterGenetics/data/Yuma_recruits_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "F")

YumaRecruitsDads <- read.csv("C:/GIT/BackwaterGenetics/data/Yuma_recruits_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "M")

YumaRecruitParents <- rbind(YumaRecruitsMoms, YumaRecruitsDads)

YumaLarvalParents <- rbind(YumaLarvalMoms, YumaLarvalDads) 

rm(YumaLarvalMoms, YumaLarvalDads, YumaRecruitsMoms, YumaRecruitsDads)  

IPLarvalMoms <- read.csv("C:/GIT/BackwaterGenetics/data/IP_larvae_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "F")

IPLarvalDads <- read.csv("C:/GIT/BackwaterGenetics/data/IP_larvae_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "M")

IPRecruitsMoms <- read.csv("C:/GIT/BackwaterGenetics/data/IP_recruits_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "F")

IPRecruitsDads <- read.csv("C:/GIT/BackwaterGenetics/data/IP_recruits_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "M")

IPRecruitParents <- rbind(IPRecruitsMoms, IPRecruitsDads)

IPLarvalParents <- rbind(IPLarvalMoms, IPLarvalDads) 

rm(IPLarvalMoms, IPLarvalDads, IPRecruitsMoms, IPRecruitsDads)  


FebYumaSurvivors <- KnownSurvivalYuma %>%
  filter(month(Date) == 2, day(Date) == 15) %>%
  mutate(Species = "XYTE", Year = year(Date), TagYear = as.factor(year(first_date)),
         location = "Yuma Cove backwater") %>% 
  arrange(Year, PITIndex) %>%
  select(Species, location, Year, TagYear, PITIndex, sex, total_length, event, MaxDAL, MaxScanDate)

FebYumaSurvivorsStocked <- FebYumaSurvivors %>%
  filter(event == "stocking")

FebYumaCounts <- FebYumaSurvivors %>%
  group_by(Species, location, Year) %>%
  summarise(Count = n()) %>%
  ungroup()

FebIPXYTESurvivors <- KnownSurvivalIPXYTE %>%
  filter(month(Date) == 2, day(Date) == 15) %>%
  mutate(Species = "XYTE", Year = year(Date), TagYear = as.factor(year(first_date))) %>% 
  arrange(location, Year, PITIndex) %>%
  select(Species, location, Year, TagYear, PITIndex, sex, total_length, event, MaxDAL, MaxScanDate)

FebIPXYTECounts <- FebIPXYTESurvivors %>%
  group_by(Species, location, Year) %>%
  summarise(Count = n()) %>%
  ungroup()

AprIPGIELSurvivors <- KnownSurvivalIPGIEL %>%
  filter(month(Date) == 4, day(Date) == 15) %>%
  mutate(Species = "GIEL", Year = year(Date), TagYear = as.factor(year(first_date))) %>% 
  arrange(location, Year, PITIndex) %>%
  select(Species, location, Year, TagYear, PITIndex, sex, total_length, event, MaxDAL, MaxScanDate)

AprIPGIELCounts <- AprIPGIELSurvivors %>%
  group_by(Species, location, Year) %>%
  summarise(Count = n()) %>%
  ungroup()

SpawnCounts <- rbind(FebYumaCounts, FebIPXYTECounts, AprIPGIELCounts)
# Create workbook for contacts with NO PITIndex
wb <- createWorkbook() # creates object to hold workbook sheets
addWorksheet(wb, "YumaFebSurvivors") # add worksheet
writeData(wb, "YumaFebSurvivors", FebYumaSurvivors) # write dataframe

addWorksheet(wb, "IPXYTEFebSurvivors") # add worksheet
writeData(wb, "IPXYTEFebSurvivors", FebIPXYTESurvivors) # write dataframe

addWorksheet(wb, "IPGIELAprSurvivors") # add worksheet
writeData(wb, "IPGIELAprSurvivors", AprIPGIELSurvivors) # write dataframe

addWorksheet(wb, "SpawnCounts") # add worksheet
writeData(wb, "SpawnCounts", SpawnCounts) # write dataframe

saveWorkbook(wb, paste0("output/BWSpawnSurvivors",
                        format(Sys.time(), "%Y%m%d"), ".xlsx"), overwrite = TRUE)

IPXYTESurvivorsHistogram <- ggplot(FebIPXYTESurvivors %>%
                                     filter(Year >= 2021, !is.na(total_length)), 
                                   aes(x = total_length, fill = TagYear)) +
  geom_histogram(binwidth = 20, position = "stack", color = "black") +
  facet_wrap(~ location + Year, ncol = 5) +
  labs(
    x = "Total Length (mm)",
    y = "Count") +
  theme_minimal()

IPXYTESurvivorsHistogram

IPGIELSurvivorsHistogram <- ggplot(AprIPGIELSurvivors %>%
                                     filter(Year >= 2021, !is.na(total_length)), 
                                   aes(x = total_length, fill = TagYear)) +
  geom_histogram(binwidth = 20, position = "stack", color = "black") +
  facet_wrap(~ location + Year, ncol = 4) +
  labs(
    x = "Total Length (mm)",
    y = "Count") +
  theme_minimal()

IPGIELSurvivorsHistogram

AprIPGIELSurvivorsSummary <- AprIPGIELSurvivors %>%
  mutate(TLCM = as.integer(total_length*.10)) %>%
  group_by(location, Year, TagYear, TLCM) %>%
  summarise(Count = n()) %>%
  ungroup()