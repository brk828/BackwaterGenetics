# Script to look at Genetic data from Dowling et al. and known survivors 
# to calculate proportion of known survivors that contributed to larvae and adults
# Only data is for Yuma Cove backwater and IP Pond 1 (XYTE)
# B. Kesner May 13, 2026

load("data/ReportingData.RData")

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(readxl) # import Excel spreadsheets
packages(data.table) # faster at indexing than grouping in dplyr
packages(openxlsx) # package openxlsx is required to create to Excel files
packages(tidyr) # needed for some pivot functions
packages(stringr) # string queries
packages(ggplot2) # plotting

# Pull stocking only records from Yuma (XYTE), IP for XYTE, and IP for GIEL. Avoiding any
# stockings after 2024
YumaStockings <- StudyBWAnalysis %>% filter(event == "stocking",
                                            location_id == 592,
                                            release_year < 2025) %>%
  select(Species = species, PITIndex, StockingYear = release_year, StockingMonth = release_month,
         StockingSex = sex, StockingTL = total_length, Survived, MaxDAL)

IPXYTEStockings <- StudyBWAnalysis %>% filter(event == "stocking",
                                              species == "XYTE",
                                              str_starts(location, "IPCA"),
                                              release_year < 2025) %>%
  select(Pond = location, Species = species, PITIndex, StockingYear = release_year, 
         StockingMonth = release_month,StockingSex = sex, StockingTL = total_length, Survived, MaxDAL)

IPGIELStockings <- StudyBWAnalysis %>% filter(event == "stocking",
                                              str_starts(location, "IPCA"),
                                              species == "GIEL",
                                              release_year < 2025) %>%
  select(Pond = location, Species = species, PITIndex, StockingYear = release_year, 
         StockingMonth = release_month,StockingSex = sex, StockingTL = total_length, Survived, MaxDAL)

rm(StudyBWGrowth, StudyBWEffort, StudyBWCaptureUniques, StudyBWTransfersAnalysis)

# Bring in larval and recruit counts by PIT tag from Dowling data for Yuma and IP1 stored as csv files in
# data directory Combine larvae and recruits per backwater and create a single dataframe for all sexes
# and all reproductive output
YumaLarvalMoms <- read.csv("data/Yuma_larvae_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "F", PIT = toupper(PIT))

YumaLarvalDads <- read.csv("data/Yuma_larvae_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "M", PIT = toupper(PIT))

YumaRecruitsMoms <- read.csv("data/Yuma_recruits_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "F", PIT = toupper(PIT))

YumaRecruitsDads <- read.csv("data/Yuma_recruits_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "M", PIT = toupper(PIT))

YumaRecruitParents <- rbind(YumaRecruitsMoms, YumaRecruitsDads)

YumaLarvalParents <- rbind(YumaLarvalMoms, YumaLarvalDads) 

YumaOffspring <- full_join(YumaLarvalParents, YumaRecruitParents, by = c("Year", "PIT"), suffix = c("_Larvae", "_Recruits")) %>%
  mutate(
    Recruits = no_na_df(Recruits),
    Larvae   = no_na_df(Larvae)
  ) %>%
  mutate( GSex = case_when(
    is.na(GSex_Larvae) ~ GSex_Recruits,
    is.na(GSex_Recruits) ~ GSex_Larvae,
    GSex_Larvae == GSex_Recruits ~ GSex_Larvae,
    TRUE ~ "U")) %>% 
  select(-GSex_Larvae, -GSex_Recruits)

rm(YumaLarvalMoms, YumaLarvalDads, YumaRecruitsMoms, YumaRecruitsDads, YumaLarvalParents)  

IPLarvalMoms <- read.csv("data/IP_larvae_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "F", PIT = toupper(PIT))

IPLarvalDads <- read.csv("data/IP_larvae_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Larvae = offspring_count) %>% select(Year, PIT, Larvae) %>%
  mutate(GSex = "M", PIT = toupper(PIT))

IPRecruitsMoms <- read.csv("data/IP_recruits_mothers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "F", PIT = toupper(PIT))

IPRecruitsDads <- read.csv("data/IP_recruits_fathers_repro_counts_PIT.csv") %>%
  rename(PIT = PIT.TAG, Recruits = offspring_count) %>% select(Year, PIT, Recruits) %>%
  mutate(GSex = "M", PIT = toupper(PIT))

IPRecruitParents <- rbind(IPRecruitsMoms, IPRecruitsDads)

IPLarvalParents <- rbind(IPLarvalMoms, IPLarvalDads) 

IPOffspring <- full_join(IPLarvalParents, IPRecruitParents, by = c("Year", "PIT"), suffix = c("_Larvae", "_Recruits")) %>%
  mutate(
    Recruits = no_na_df(Recruits),
    Larvae   = no_na_df(Larvae)
  ) %>%
  mutate( GSex = case_when(
    is.na(GSex_Larvae) ~ GSex_Recruits,
    is.na(GSex_Recruits) ~ GSex_Larvae,
    GSex_Larvae == GSex_Recruits ~ GSex_Larvae,
    TRUE ~ "U")) %>% 
  select(-GSex_Larvae, -GSex_Recruits)

rm(IPLarvalMoms, IPLarvalDads, IPRecruitsMoms, IPRecruitsDads, IPLarvalParents)  

# Check for any genetic results that aren't in the stocking dataframes
YumaOffspringNoStocking <- YumaOffspring %>% 
  anti_join(YumaStockings, by = c("PIT" = "PITIndex"))

IPOffspringNoStocking <- IPOffspring %>% 
  anti_join(IPXYTEStockings %>% filter(), by = c("PIT" = "PITIndex"))

# if there are any genetic results without stocking data, stop process
if(sum(nrow(YumaOffspringNoStocking) + nrow(IPOffspringNoStocking)) > 0){
  stop("You have genetic data that doesn't match stocking records, please review the 
       NoStocking dataframes for details.")
}else{rm(YumaOffspringNoStocking, IPOffspringNoStocking)}

# Start with Known survivors dataframe, look for fish that were known to be alive 
# as of February 15 each year. These fish are considered "available" for reproduction
# except for stocking year, which we consider all stocked fish to be available. 
# Create a dataframe with stocking information and known survivors for each year for
# all stocked fish, then add offspring data to this dataframe. 
FebYumaSurvivors <- KnownSurvivalYuma %>%
  filter(month(Date) == 2, day(Date) == 15) %>%
  mutate(Species = "XYTE", Year = year(Date), TagYear = year(first_date),
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

# Same process for all XYTE ponds at IP. Only Pond 1 currently has offspring data
# but future proofed to include all ponds and all stocking years.
FebIPXYTESurvivors <- KnownSurvivalIPXYTE %>%
  filter(month(Date) == 2, day(Date) == 1) %>%
  mutate(Species = "XYTE", Year = year(Date), 
         Pond = location,
         TagYear = year(first_date)) %>% 
  arrange(Pond, PITIndex, Year) %>%
  select(Pond, Species, PITIndex, Year, TagYear) %>% 
  filter(Year != TagYear) %>%
  select(-TagYear) %>%
  left_join(IPXYTEStockings, by = c("PITIndex","Pond", "Species")) %>%
  filter(!is.na(StockingYear)) %>%
  rbind(IPXYTEStockings %>% mutate(Year = StockingYear)) %>%
  left_join(IPOffspring %>%
              select(-GSex)
            , by = c("PITIndex" = "PIT", "Year")) %>%
  mutate(Recruits = no_na_df(Recruits),
         Larvae = no_na_df(Larvae)) %>%
  mutate(Offspring = Recruits + Larvae)

# Same process for all GIEL at IP. No offspring data currently, but can be added later
AprIPGIELSurvivors <- KnownSurvivalIPGIEL %>%
  filter(month(Date) == 4, day(Date) == 15) %>%
  mutate(Species = "GIEL", Year = year(Date), 
         Pond = location,
         TagYear = year(first_date)) %>% 
  arrange(Pond, PITIndex, Year) %>%
  select(Pond, Species, PITIndex, Year, TagYear) %>% 
  filter(Year != TagYear) %>%
  select(-TagYear) %>%
  left_join(IPGIELStockings, by = c("PITIndex","Pond", "Species")) %>%
  filter(!is.na(StockingYear)) %>%
  rbind(IPGIELStockings %>% mutate(Year = StockingYear)) %>%
  left_join(IPOffspring %>%
              select(-GSex)
            , by = c("PITIndex" = "PIT", "Year")) %>%
  mutate(Recruits = no_na_df(Recruits),
         Larvae = no_na_df(Larvae)) %>%
  mutate(Offspring = Recruits + Larvae)

# Reduce data to counts for survivors, fish that produced larvae, recruits, or either
# Calculate proportion of known population that produced offspring
FebYumaCounts <- FebYumaSurvivors %>%
  group_by(Species, StockingSex, StockingYear, Year) %>%
  summarise(Count = n(),
            LarvaePos = sum(Larvae > 0),
            RecruitsPos = sum(Recruits > 0),
            Offspring = sum(Offspring > 0)) %>%
  ungroup() %>%
  mutate(PropOffspring = round(Offspring/Count, 3))

FebIPXYTECounts <- FebIPXYTESurvivors %>%
  group_by(Pond, StockingSex, StockingYear, Year) %>%
  summarise(Count = n(),
            Larvae = sum(Larvae > 0),
            Recruits = sum(Recruits > 0),
            Offspring = sum(Offspring > 0)) %>%
  ungroup() %>%
  mutate(PropOffspring = round(Offspring/Count, 3))


AprIPGIELCounts <- AprIPGIELSurvivors %>%
  group_by(Pond, StockingSex, StockingYear, Year) %>%
  summarise(Count = n(),
            Larvae = sum(Larvae > 0),
            Recruits = sum(Recruits > 0),
            Offspring = sum(Offspring > 0)) %>%
  ungroup() %>%
  mutate(PropOffspring = round(Offspring/Count, 3))

# Create workbook for contacts with NO PITIndex
wb <- createWorkbook() # creates object to hold workbook sheets
addWorksheet(wb, "FebYumaCounts") # add worksheet
writeData(wb, "FebYumaCounts", FebYumaCounts) # write dataframe

addWorksheet(wb, "FebIPXYTECounts") # add worksheet
writeData(wb, "FebIPXYTECounts", FebIPXYTECounts) # write dataframe

addWorksheet(wb, "AprIPGIELCounts") # add worksheet
writeData(wb, "AprIPGIELCounts", AprIPGIELCounts) # write dataframe

saveWorkbook(wb, paste0("output/BWSpawnSurvivors",
                        format(Sys.time(), "%Y%m%d"), ".xlsx"), overwrite = TRUE)

IPXYTESurvivorsHistogram <- ggplot(FebIPXYTESurvivors %>%
                                     filter(Year >= 2021, !is.na(StockingTL)), 
                                   aes(x = StockingTL, fill = as.factor(StockingYear))) +
  geom_histogram(binwidth = 20, position = "stack", color = "black") +
  facet_wrap(~ Pond + Year, ncol = 5) +
  labs(
    x = "Stocking Total Length (mm)",
    y = "Count",
    fill = "Stocking Year") +
  theme_minimal()

IPXYTESurvivorsHistogram

IPGIELSurvivorsHistogram <- ggplot(AprIPGIELSurvivors %>%
                                     filter(Year < 2021, !is.na(StockingTL)), 
                                   aes(x = StockingTL, fill = as.factor(StockingYear))) +
  geom_histogram(binwidth = 20, position = "stack", color = "black") +
  facet_wrap(~ Pond + Year, ncol = 4) +
  labs(
    x = "Stocking Total Length (mm)",
    y = "Count",
    fill = "Stocking Year") +
  theme_minimal()

IPGIELSurvivorsHistogram

AprIPGIELSurvivorsSummary <- AprIPGIELSurvivors %>%
  mutate(TLCM = as.integer(StockingTL*.10)) %>%
  group_by(Pond, Year, StockingYear, TLCM) %>%
  summarise(Count = n()) %>%
  ungroup()

pos <- position_dodge(width = 0.4)

YumaOffSpringPlot <- ggplot(FebYumaCounts %>% filter(StockingYear != 2020),
                            aes(x = as.factor(Year),
                                y = PropOffspring,
                                color = as.factor(StockingYear),
                                group = as.factor(StockingYear))) +
  geom_point(position = pos) +
  geom_line(position = pos) +
  geom_text(aes(label = Count),
            vjust = -0.6,
            size = 2.5,
            color = "black",
            show.legend = FALSE) +
  facet_wrap(~ StockingSex, ncol = 1) +
  scale_color_brewer(palette = "Dark2", name = "Stocking Year") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Year",
    y = "Proportion") +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

png(paste0("output/YumaOffspringPlot.png"), width = 6, height = 4, units = 'in', res = 300)   
YumaOffSpringPlot
dev.off()

YumaOffSpringPlot

IPCA1OffSpringPlot <- ggplot(FebIPXYTECounts %>% filter(StockingYear == 2016, Year < 2023, Pond == "IPCA (Pond 1)"),
                            aes(x = as.factor(Year),
                                y = PropOffspring,
                                color = as.factor(StockingYear),
                                group = as.factor(StockingYear))) +
  geom_point(position = pos) +
  geom_line(position = pos) +
  geom_text(aes(label = Count),
            vjust = -0.6,
            size = 2.5,
            color = "black",
            show.legend = FALSE) +
  facet_wrap(~ StockingSex, ncol = 1) +
  scale_color_brewer(palette = "Dark2", name = "Stocking Year") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Year",
    y = "Proportion") +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

IPCA1OffSpringPlot

png(paste0("output/IPCA1OffSpringPlot.png"), width = 6, height = 4, units = 'in', res = 300)   
IPCA1OffSpringPlot
dev.off()
