# Annual mark-recapture tagged population estimates for IP 1-6 and Yuma Cove Backwater
# B. Kesner September 2025

load("data/ReportingData.RData")
packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(ggplot2) # Plotting
packages(stringr) # string manipulations and function

# Razorback sucker stocked in Winter, bonytail in spring. 

# Marks are restricted to January or February of the census year (month < 3)
# The census year is equal to the year of scanning
BWMarks <- StudyBWContacts %>% 
  filter(year(tagging_date) < year(Date)) %>%
  select(Location, Species, Date, PITIndex, tagging_date) %>%
  filter(month(Date) < 3) %>%
  mutate(CensusYear = year(Date)) %>%
  group_by(CensusYear, Species, Location, PITIndex, tagging_date) %>%
  summarise(Contacts = n(), FirstScan = min(Date), LastScan = max(Date)) %>%
  ungroup()

# Captures are restricted to October through April (month > 9 or month < 5)
# The census year is the year in which the October scanning begins but a year less
# than scans for January through April, Fiscal Year is a year ahead of this schedule 
BWCaptures <- StudyBWContacts %>% 
  filter(year(tagging_date) < ScanFY - 1) %>% 
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
  na.omit()
