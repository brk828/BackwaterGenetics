# Script to look at Genetic data from Dowling et al. and known survivors 
# to calculate proportion of known survivors that contributed to larvae and adults
# Only data is for Yuma Cove backwater and IP Pond 1 (XYTE)
# B. Kesner May 13, 2026

load("data/KnownSurvivalCounts.RData")
load("data/ReportingData.RData")

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(readxl) # import Excel spreadsheets
packages(data.table) # faster at indexing than grouping in dplyr
packages(openxlsx) # package openxlsx is required to create to Excel files
packages(tidyr) # needed for some pivot functions

YumaStockings <- StudyBWAnalysis %>% filter(event == "stocking",
                                            location_id == 592,
                                            release_year < 2025) %>%
  select(Species = species, PITIndex, StockingYear = release_year, StockingMonth = release_month,
         StockingSex = sex, StockingTL = total_length, Survived, MaxDAL, )

IP1Stockings <- StudyBWAnalysis %>% filter(event == "stocking",
                                        location_id == 1043,
                                        release_year == 2016) %>%
  select(Species = species, PITIndex, StockingYear = release_year, StockingMonth = release_month,
         StockingSex = sex, StockingTL = total_length, Survived, MaxDAL, )

rm(StudyBWGrowth, StudyBWEffort, StudyBWCaptureUniques, StudyBWTransfersAnalysis)

# Bring in larval and recruit counts by PIT tag from Dowling data for Yuma and IP1 stored as csv files in
# data directory
YumaLarvalMoms <- read.csv("C:/GIT/BackwaterGenetics/data/Yuma_larvae_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "F", PIT = toupper(PIT))

YumaLarvalDads <- read.csv("C:/GIT/BackwaterGenetics/data/Yuma_larvae_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "M", PIT = toupper(PIT))

YumaRecruitsMoms <- read.csv("C:/GIT/BackwaterGenetics/data/Yuma_recruits_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "F", PIT = toupper(PIT))

YumaRecruitsDads <- read.csv("C:/GIT/BackwaterGenetics/data/Yuma_recruits_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "M", PIT = toupper(PIT))

YumaRecruitParents <- rbind(YumaRecruitsMoms, YumaRecruitsDads)

YumaLarvalParents <- rbind(YumaLarvalMoms, YumaLarvalDads) 

YumaOffspring <- full_join(YumaLarvalParents, YumaRecruitParents, by = c("Year", "PIT"), suffix = c("_Recruits", "_Larvae")) %>%
  mutate(
    Recruits = no_na_df(Recruits),
    Larvae   = no_na_df(Larvae)
  ) %>%
  mutate( GSex = case_when(
    is.na(GSex_Recruits) ~ GSex_Larvae,
    is.na(GSex_Larvae) ~ GSex_Recruits,
    GSex_Recruits == GSex_Larvae ~ GSex_Recruits,
    TRUE ~ "U")) %>% 
  select(-GSex_Recruits, -GSex_Larvae)

rm(YumaLarvalMoms, YumaLarvalDads, YumaRecruitsMoms, YumaRecruitsDads, YumaLarvalParents)  

IPLarvalMoms <- read.csv("C:/GIT/BackwaterGenetics/data/IP_larvae_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "F", PIT = toupper(PIT))

IPLarvalDads <- read.csv("C:/GIT/BackwaterGenetics/data/IP_larvae_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "M", PIT = toupper(PIT))

IPRecruitsMoms <- read.csv("C:/GIT/BackwaterGenetics/data/IP_recruits_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "F", PIT = toupper(PIT))

IPRecruitsDads <- read.csv("C:/GIT/BackwaterGenetics/data/IP_recruits_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "M", PIT = toupper(PIT))

IPRecruitParents <- rbind(IPRecruitsMoms, IPRecruitsDads)

IPLarvalParents <- rbind(IPLarvalMoms, IPLarvalDads) 

IPOffspring <- full_join(IPLarvalParents, IPRecruitParents, by = c("Year", "PIT"), suffix = c("_Recruits", "_Larvae")) %>%
  mutate(
    Recruits = no_na_df(Recruits),
    Larvae   = no_na_df(Larvae)
  ) %>%
  mutate( GSex = case_when(
    is.na(GSex_Recruits) ~ GSex_Larvae,
    is.na(GSex_Larvae) ~ GSex_Recruits,
    GSex_Recruits == GSex_Larvae ~ GSex_Recruits,
    TRUE ~ "U")) %>% 
  select(-GSex_Recruits, -GSex_Larvae)

rm(IPLarvalMoms, IPLarvalDads, IPRecruitsMoms, IPRecruitsDads, IPLarvalParents)  

YumaOffspringNoStocking <- YumaOffspring %>% 
  anti_join(YumaStockings, by = c("PIT" = "PITIndex"))

IPOffspringNoStocking <- IPOffspring %>% 
  anti_join(IP1Stockings %>% filter(), by = c("PIT" = "PITIndex"))


if(sum(nrow(YumaOffspringNoStocking) + nrow(IPOffspringNoStocking)) > 0){
  stop("You have genetic data that doesn't match stocking records, please review the 
       NoStocking dataframes for details.")
}else{rm(YumaOffspringNoStocking, IPOffspringNoStocking)}

FebYumaSurvivors <- KnownSurvivalYuma %>%
  filter(month(Date) == 2, day(Date) == 15) %>%
  mutate(Species = "XYTE", Year = year(Date), TagYear = as.factor(year(first_date)),
         location = "Yuma Cove backwater") %>% 
  arrange(Year, PITIndex) %>%
  select(PITIndex, Year, TagYear) %>% 
  filter(Year != TagYear) %>%
  select(-TagYear) %>%
  left_join(YumaStockings, by = "PITIndex") %>%
  filter(!is.na(StockingYear)) %>%
  rbind(YumaStockings %>% mutate(Year = StockingYear)) %>%
  left_join(YumaOffspring %>%
              select(-GSex)
            , by = c("PITIndex" = "PIT", "Year")) %>%
  mutate(Recruits = no_na_df(Recruits),
         Larvae = no_na_df(Larvae)) %>%
  mutate(Offspring = Recruits + Larvae)


FebYumaCounts <- FebYumaSurvivors %>%
  group_by(Species, StockingSex, StockingYear, Year) %>%
  summarise(Count = n(),
            LarvaePos = sum(Larvae > 0),
            RecruitsPos = sum(Recruits > 0),
            Offspring = sum(Offspring > 0)) %>%
  ungroup() %>%
  mutate(PropOffspring = round(Offspring/Count, 3))

FebIP1Survivors <- KnownSurvivalIPXYTE %>%
  filter(month(Date) == 2, day(Date) == 15, location == "IPCA (Pond 1)") %>%
  mutate(Species = "XYTE", Year = year(Date), TagYear = as.factor(year(first_date))) %>% 
  arrange(Year, PITIndex) %>%
  select(PITIndex, Year, TagYear) %>% 
  filter(Year != TagYear) %>%
  select(-TagYear) %>%
  left_join(IP1Stockings, by = "PITIndex") %>%
  filter(!is.na(StockingYear)) %>%
  rbind(IP1Stockings %>% mutate(Year = StockingYear)) %>%
  left_join(IPOffspring %>%
              select(-GSex)
            , by = c("PITIndex" = "PIT", "Year")) %>%
  mutate(Recruits = no_na_df(Recruits),
         Larvae = no_na_df(Larvae)) %>%
  mutate(Offspring = Recruits + Larvae)

FebIP1Counts <- FebIP1Survivors %>%
  group_by(Species, StockingSex, StockingYear, Year) %>%
  summarise(Count = n(),
            LarvaePos = sum(Larvae > 0),
            RecruitsPos = sum(Recruits > 0),
            Offspring = sum(Offspring > 0)) %>%
  ungroup() %>%
  mutate(PropOffspring = round(Offspring/Count, 3))

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