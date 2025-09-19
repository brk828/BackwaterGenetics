# Annual mark-recapture tagged population estimates for IP 1-6 and Yuma Cove Backwater
# B. Kesner September 2025

load("data/ReportingData.RData")
load("data/KnownSurvivalCounts.RData")

TotalCountIPXYTE <- TotalCountIPXYTE %>%
  rename(Location = location)

TotalCountIPGIEL <- TotalCountIPGIEL %>%
  rename(Location = location)

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(ggplot2) # Plotting
packages(stringr) # string manipulations and function

# Razorback sucker stocked in Winter, bonytail in spring. 

# Marks are restricted to January or February of the census year (month < 3)
# tagging must have occurred at least six months prior to the first scannning month
BWMarks <- StudyBWContacts %>% 
  filter(!is.na(tagging_date)) %>%
  filter(Date >= tagging_date %m+% months(6)) %>%
  select(Location, Species, Date, PITIndex, tagging_date) %>%
  filter(month(Date) < 3) %>%
  mutate(CensusYear = year(Date)) %>%
  group_by(CensusYear, Species, Location, PITIndex, tagging_date) %>%
  summarise(Contacts = n(), FirstScan = min(Date), LastScan = max(Date)) %>%
  ungroup()

# Captures are restricted to October through April (month > 9 or month < 5)
# The census year is the year in which the October scanning begins but a year less
# than scans for January through April, Fiscal Year is a year ahead of this schedule 
# Filtered for 15 months so as same tagging_date cutoff as for marks (10 - 1 = 9)
BWCaptures <- StudyBWContacts %>% 
  filter(!is.na(tagging_date)) %>%
  filter(Date >= tagging_date %m+% months(15)) %>%
  select(Location, Species, Date, PITIndex, ScanFY, tagging_date)  %>%
  filter(month(Date) > 9| month(Date) < 5) %>%
  mutate(CensusYear = ScanFY - 1) %>%
  group_by(CensusYear, Species, Location, PITIndex, tagging_date) %>%
  summarise(Contacts = n(), FirstScan = min(Date), LastScan = max(Date)) %>%
  ungroup()

BWRecaptures <- BWMarks %>%
  select(Location, Species, CensusYear, PITIndex) %>%
  inner_join(BWCaptures %>%
               select(Location, Species, CensusYear, PITIndex), 
             by = c("Location" = "Location",
                    "Species" = "Species",
                    "CensusYear" = "CensusYear", 
                    "PITIndex" = "PITIndex"))

BWMark <-BWMarks %>%
  group_by(Location, Species, CensusYear) %>%
  summarise(M = n_distinct(PITIndex)) %>%
  ungroup()

BWCapture <- BWCaptures %>%
  group_by(Location, Species, CensusYear) %>%
  summarise(C = n_distinct(PITIndex)) %>%
  ungroup()

BWRecapture <- BWRecaptures %>%
  group_by(Location, Species, CensusYear) %>%
  summarise(R = n_distinct(PITIndex)) %>%
  ungroup()

BWEstimates <- BWMark %>%
  inner_join(BWCapture, by = c("Species" = "Species",
                              "Location" = "Location",
                              "CensusYear" = "CensusYear")) %>%
  inner_join(BWRecapture, by =  c("Species" = "Species",
                                  "Location" = "Location",
                                  "CensusYear" = "CensusYear")) %>%
  filter(R>3) %>%
  arrange(Location, Species, CensusYear) %>%
  mutate(LowerBoundR = qpois(0.025, R),
         UpperBoundR = qpois(0.975, R),
         Estimate = as.integer(((M + 1) * (C + 1))/(R + 1))) %>%
  mutate(LowerQP95CI = as.integer(((M + 1) * (C + 1))/(UpperBoundR + 1)),
         UpperQP95CI = as.integer(((M + 1) * (C + 1))/(LowerBoundR + 1)),
         SE = as.integer((sqrt(((M+1)^2*(C+1)*(C-R)/((R+2)*(R+1)^2))))),
         LowerN95CI = as.integer(Estimate-(1.96*SE)),
         UpperN95CI = as.integer(Estimate+(1.96*SE))) %>%
  na.omit() %>%
  mutate(CensusDate = as.Date(paste0(CensusYear, "-01-01")))

YumaEstimatePlot <- ggplot(BWEstimates %>%
                             filter(Location == "Yuma Cove backwater"), 
                           aes(x = CensusDate, y = Estimate)) +
  geom_point(shape = 21, fill = "steelblue", color = "black", size = 3, stroke = 0.5) +
  geom_errorbar(aes(ymin = LowerN95CI, ymax = UpperN95CI), width = 10, color = "steelblue") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "Estimate") +
  theme_minimal() +
  theme(axis.text = element_text(size = 11),
        axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank()) +
  geom_line(data = TotalCountYuma,
            aes(x = Date, y = Count),
            linewidth = 1.2, alpha = 0.3)


YumaEstimatePlot

IPXYTEEstimatePlot <- ggplot(BWEstimates %>%
                             filter(str_starts(Location, "IPCA"), Species == "XYTE"), 
                           aes(x = CensusDate, y = Estimate)) +
  geom_point(shape = 21, fill = "steelblue", color = "black", size = 3, stroke = 0.5) +
  geom_errorbar(aes(ymin = LowerN95CI, ymax = UpperN95CI), width = 10, color = "steelblue") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "Estimate") +
  theme_minimal() +
  theme(axis.text = element_text(size = 11),
        axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Location, nrow=3) +
  geom_line(data = TotalCountIPXYTE,
            aes(x = Date, y = Count),
            linewidth = 1.2, alpha = 0.3)

IPXYTEEstimatePlot

IPGIELEstimatePlot <- ggplot(BWEstimates %>%
                               filter(str_starts(Location, "IPCA"), Species == "GIEL"), 
                             aes(x = CensusDate, y = Estimate)) +
  geom_point(shape = 21, fill = "steelblue", color = "black", size = 3, stroke = 0.5) +
  geom_errorbar(aes(ymin = LowerN95CI, ymax = UpperN95CI), width = 10, color = "steelblue") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "Estimate") +
  theme_minimal() +
  theme(axis.text = element_text(size = 11),
        axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Location, nrow=3) +
  geom_line(data = TotalCountIPGIEL,
            aes(x = Date, y = Count),
            linewidth = 1.2, alpha = 0.3)

IPGIELEstimatePlot

png(paste0("output/YumaEstimatePlot.png"), width = 6, height = 4, units = 'in', res = 300)   
YumaEstimatePlot
dev.off()

png(paste0("output/IPXYTEEstimatePlot.png"), width = 6, height = 4, units = 'in', res = 300)   
IPXYTEEstimatePlot
dev.off()

png(paste0("output/IPGIELEstimatePlot.png"), width = 6, height = 4, units = 'in', res = 300)   
IPGIELEstimatePlot
dev.off()
