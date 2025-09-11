# Minimum known survivors for IP 1-6 and Yuma Cove backwater
# B. Kesner 11 Sep 2025

load("data/ReportingData.RData")

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(ggplot2) # Plotting

StudyBWNFWGAnalysis <- StudyBWNFWGAnalysis %>%
  mutate(MaxScanDate = if_else(is.na(MaxScanDate), first_date, MaxScanDate))

StudyBWXYTEIP <- StudyBWNFWGAnalysis %>% 
  filter(location_id != 592, species == "XYTE", MaxDAL > SurvivalDAL)

StudyBWGIELIP <- StudyBWNFWGAnalysis %>% 
  filter(location_id != 592, species == "GIEL", MaxDAL > SurvivalDAL)

StudyBWYuma <- StudyBWNFWGAnalysis %>% 
  filter(location_id == 592, MaxDAL > SurvivalDAL)

YumaMinReleaseDate <- as.Date(min(StudyBWYuma$first_date))
IPXYTEMinReleaseDate <- as.Date(min(StudyBWXYTEIP$first_date))
IPGIELMinReleaseDate <- as.Date(min(StudyBWGIELIP$first_date))

SurvDaysYuma <- data.frame(Date = seq(YumaMinReleaseDate, Sys.Date() - SurvivalDAL*2,
                                      by =  "day"))     

SurvDaysIPXYTE <- data.frame(Date = seq(IPXYTEMinReleaseDate, Sys.Date() - SurvivalDAL*2, 
                                      by =  "day"))   

SurvDaysIPGIEL <- data.frame(Date = seq(IPGIELMinReleaseDate, Sys.Date() - SurvivalDAL*2, 
                                        by =  "day"))   

CrossDFYuma <- SurvDaysYuma %>%
  mutate(key = 1) %>%
  full_join(StudyBWYuma %>% mutate(key = 1), 
            by = "key", relationship = "many-to-many") %>%
  select(-key)

KnownSurvivalYuma <- CrossDFYuma %>% 
  filter(first_date <= Date &
           MaxScanDate >= Date) %>%
  select(Date, PITIndex, first_date, sex, total_length, event, MaxScanDate, MaxDAL) 

KnownSurvivalPlotYuma <- ggplot(KnownSurvivalYuma, aes(x = Date)) +
  geom_line(stat = "count", aes(linetype = sex), width = 1) + 
  scale_x_date(date_breaks = "1 month", 
               labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), "")) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "black")) +
  labs(x = "Date", y = "Count", color = "sex")

KnownSurvivalPlotYuma
