# Minimum known survivors for IP 1-6 and Yuma Cove backwater
# B. Kesner 11 Sep 2025

load("data/ReportingData.RData")

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(ggplot2) # Plotting

# If MaxScanDay is null, make it the release or tagging date of the fish
# change juvenile to unknown sex to simplify sexes
StudyBWNFWGAnalysis <- StudyBWNFWGAnalysis %>%
  mutate(sex = if_else(sex == "J", "U", sex))

# separate out the 3 study groups for analysis
StudyBWXYTEIP <- StudyBWNFWGAnalysis %>% 
  filter(location_id != 592, species == "XYTE", MaxDAL > SurvivalDAL)

StudyBWGIELIP <- StudyBWNFWGAnalysis %>% 
  filter(location_id != 592, species == "GIEL", MaxDAL > SurvivalDAL)

StudyBWYuma <- StudyBWNFWGAnalysis %>% 
  filter(location_id == 592, MaxDAL > SurvivalDAL)

# Determine start date for each study group and create a dataframe 
# with each day as a value from starting date to 2 x SurvivalDAL
YumaMinReleaseDate <- as.Date(min(StudyBWYuma$first_date))
IPXYTEMinReleaseDate <- as.Date(min(StudyBWXYTEIP$first_date))
IPGIELMinReleaseDate <- as.Date(min(StudyBWGIELIP$first_date))

SurvDaysYuma <- data.frame(Date = seq(YumaMinReleaseDate, Sys.Date() - SurvivalDAL*2,
                                      by =  "day"))     
SurvDaysXYTEIP <- data.frame(Date = seq(IPXYTEMinReleaseDate, Sys.Date() - SurvivalDAL*2, 
                                      by =  "day"))   
SurvDaysGIELIP <- data.frame(Date = seq(IPGIELMinReleaseDate, Sys.Date() - SurvivalDAL*2, 
                                        by =  "day"))   
# Create a crossdf of the days and the BW data, this duplicates all the rows of the 
# study fish data for all dates in the SurvDays dataframe
CrossDFYuma <- SurvDaysYuma %>%
  mutate(key = 1) %>%
  full_join(StudyBWYuma %>% mutate(key = 1), 
            by = "key", relationship = "many-to-many") %>%
  select(-key)

CrossDFIPXYTE <- SurvDaysXYTEIP %>%
  mutate(key = 1) %>%
  full_join(StudyBWXYTEIP %>% mutate(key = 1), 
            by = "key", relationship = "many-to-many") %>%
  select(-key)

CrossDFIPGIEL <- SurvDaysGIELIP %>%
  mutate(key = 1) %>%
  full_join(StudyBWXYTEIP %>% mutate(key = 1), 
            by = "key", relationship = "many-to-many") %>%
  select(-key)

# filter out for each day the fish that were tagged or released prior to the Date
# their MaxScanDate is after the Date, and they were know to survive past the SurvivalDAL
KnownSurvivalYuma <- CrossDFYuma %>% 
  filter(first_date <= Date &
           MaxScanDate >= Date, 
         Survived == 1) %>%
  select(Date, PITIndex, first_date, sex, total_length, event, MaxScanDate, MaxDAL) 

KnownSurvivalIPXYTE <- CrossDFIPXYTE %>% 
  filter(first_date <= Date &
           MaxScanDate >= Date, 
         Survived == 1) %>%
  select(location, Date, PITIndex, first_date, sex, total_length, event, MaxScanDate, MaxDAL) 

KnownSurvivalIPGIEL <- CrossDFIPGIEL %>% 
  filter(first_date <= Date &
           MaxScanDate >= Date, 
         Survived == 1) %>%
  select(location, Date, PITIndex, first_date, sex, total_length, event, MaxScanDate, MaxDAL) 

KnownSurvivalPlotYuma <- ggplot(KnownSurvivalYuma, aes(x = Date)) +
  geom_line(aes(color = sex, linetype = sex), stat = "count", width = 1) + 
  scale_x_date(date_breaks = "1 month", 
               labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), "")) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "black")) +
  labs(x = "Date", y = "Count")

KnownSurvivalPlotYuma

KnownSurvivalPlotIPXYTE <- ggplot(KnownSurvivalIPXYTE, aes(x = Date)) +
  geom_line(aes(color = sex, linetype = sex), stat = "count", width = 1) + 
  scale_x_date(date_breaks = "1 month", 
               labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), "")) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "black")) +
  labs(x = "Date", y = "Count") +
  facet_wrap(~location, nrow=3) 

KnownSurvivalPlotIPXYTE

KnownSurvivalPlotIPGIEL <- ggplot(KnownSurvivalIPGIEL, aes(x = Date)) +
  geom_line(aes(color = sex, linetype = sex), stat = "count", width = 1) + 
  scale_x_date(date_breaks = "1 month", 
               labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), "")) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "black")) +
  labs(x = "Date", y = "Count") +
  facet_wrap(~location, nrow=3) 

KnownSurvivalPlotIPGIEL
