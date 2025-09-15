
# CJS survival analysis for IP 1-6 and Yuma Cove Backwater
# B. Kesner September 2025

load("data/ReportingData.RData")
packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(ggplot2) # Plotting
packages(tidyr) # pivot table
packages(RMark) # MCR analysis
packages(stringr) # string manipulations and function

# Razorback sucker stocked in Winter, bonytail in spring. First survival period is from
# stocking to fall-winter of next FY. bonytail and razorback sucker stocked in FY2017 - IP.
# First session FY2018 months Oct through January (ScanMonthFY 1-4). Survival determined annually 
# based on Scanning through those months and simple CMS model for stocked fish

StockingBW <- StudyBWNFWGAnalysis %>%
  filter(event == "stocking", first_date < as.Date("2024-10-01")) %>%
  rename(ReleaseDate = first_date) %>%
  mutate(ReleaseFY = as.factor(ifelse(month(ReleaseDate) > 9, year(ReleaseDate)+1, year(ReleaseDate))),
         ReleaseMonth = as.integer(month(ReleaseDate)),
         ReleaseMonthName = format(ReleaseDate, "%b"),
         ReleaseMonthFY = ifelse(ReleaseMonth>9, ReleaseMonth-9, ReleaseMonth+3),
         sex = trimws(sex),
         location = as.factor(location),
         TLClass = as.integer(total_length*.1)*10) 

ContactsFYMCR <- StudyBWContacts %>%
  mutate(ScanMonthFY = ifelse(ScanMonth>9, ScanMonth-9, ScanMonth+3)) %>%
  filter(ScanMonthFY >= 1, ScanMonthFY <= 4) %>%
  arrange(ScanFY, Location, PITIndex) %>%
  group_by(Location, LID, ScanFY, PITIndex) %>%
  summarise(Contacts = n()) %>%
  ungroup()

ContactsFYMCRSummary <- ContactsFYMCR %>%
  group_by(Location, ScanFY) %>%
  summarise(Contacts = sum(Contacts), Uniques = n_distinct(PITIndex)) %>%
  ungroup()

FirstBWStockingIP <- StockingBW %>%
  filter(ReleaseFY == 2017) %>%
  select(ReleaseFY, ReleaseMonthFY, PITIndex, Location = location, location_id, species, sex, total_length, TLClass) %>%
  mutate(PostReleaseMonths = 13 - ReleaseMonthFY) %>%
  filter(!is.na(sex), !is.na(total_length))

FirstBWStockingYuma <- StockingBW %>%
  filter(ReleaseFY == 2013) %>%
  select(ReleaseFY, ReleaseMonthFY, PITIndex, Location = location, location_id, species, sex, total_length, TLClass) %>%
  mutate(PostReleaseMonths = 13 - ReleaseMonthFY) %>%
  filter(!is.na(sex), !is.na(total_length))

ContactsFYMCRPivotIP <- ContactsFYMCR %>%
  filter(LID != 592) %>%
  mutate(Contacts = 1) %>%
  pivot_wider(
    names_from = ScanFY,
    values_from = Contacts,
    values_fill = list(Contacts = 0))


ContactsFYMCRPivotYuma <- ContactsFYMCR %>%
  filter(LID == 592) %>%
  mutate(Contacts = 1) %>%
  pivot_wider(
    names_from = ScanFY,
    values_from = Contacts,
    values_fill = list(Contacts = 0))


MCRDataIP <- FirstBWStockingIP %>%
  mutate(Release = 1) %>%
  left_join(ContactsFYMCRPivotIP, by = c("PITIndex", "Location")) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

MCRDataYuma <- FirstBWStockingYuma %>%
  mutate(Release = 1) %>%
  left_join(ContactsFYMCRPivotYuma, by = c("PITIndex", "Location")) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))


rm(ContactsFYMCRPivotYuma, ContactsFYMCRPivotIP, FirstBWStockingIP, FirstBWStockingYuma, ContactsFYMCR)


IPHistory <- unite(MCRDataIP %>% select(PITIndex, Release, starts_with("20")), ch, -1, sep="", remove = TRUE) %>%
  inner_join(MCRDataIP %>% select(PITIndex, Location, TL = total_length, sex), by = "PITIndex") %>%
  mutate(ch = as.character(ch))

YumaHistory <- unite(MCRDataYuma %>% select(PITIndex, Release, starts_with("20")), ch, -1, sep="", remove = TRUE) %>%
  inner_join(MCRDataYuma %>% select(PITIndex, Location, TL = total_length, sex), by = "PITIndex") %>%
  mutate(ch = as.character(ch))

# Wait to factor sex until individual pond selection as some ponds have U for sex
IPHistory1 <- IPHistory %>% filter(Location == "IPCA (Pond 1)") %>%
  mutate(sex = as.factor(sex))

dpPond1 = RMark::process.data(IPHistory1, model="CJS", groups = "sex")
ddlPond1 = RMark::make.design.data(dpPond1)

# Building a max 2 age structure
ddlPond1$Phi$age2 <- ifelse(ddlPond1$Phi$Age >2,3,ddlPond1$Phi$age)
ddlPond1$Phi$age2 <- as.factor(ddlPond1$Phi$age2)
ddlPond1$Phi$first_year <- ifelse(ddlPond1$Phi$Age == 1, 1, 0)

# Create two time sample structure as earlier ages have smaller scan windows
# ddl$p$time2 <- ifelse(ddl$p$Age > 2, 3, 1)
# ddl$p$time2 <- as.factor(ddl$p$time2)

Phi.age2TL = list(formula=~age2 + first_year:TL) 
Phi.timeTL = list(formula=~time + sex + first_year:TL) 
p.dot = list(formula=~1)
p.time = list(formula=~time)

Pond1TLAgeModel <- mark(data=dpPond1,ddl=ddlPond1,model.parameters=list(Phi=Phi.age2TL,p=p.time))
Pond1TLTimeModel <- mark(data=dpPond1,ddl=ddlPond1,model.parameters=list(Phi=Phi.timeTL,p=p.dot))

Pond1TLAgeResults <- Pond1TLAgeModel$results$real %>%
  tibble::rownames_to_column(var = "parameter")  # Store row names in a column

Pond1TLTimeResults <- Pond1TLTimeModel$results$real %>%
  tibble::rownames_to_column(var = "parameter")  # Store row names in a column

# Extract survival (S) estimates by age
Pond1TimeSurvival <- Pond1TLTimeResults %>%
  filter(str_starts(parameter, "Phi")) %>%
  mutate(age = row_number(),  # Assign sequential age numbers
         parameter = "Survival")  # Rename for clarity

# Extract recapture (p) estimates by time
Pond1TimeRecapture <- Pond1TLTimeResults %>%
  filter(str_starts(parameter, "p")) %>%
  mutate(time = row_number(),  # Assign sequential time numbers
         parameter = "Recapture")  # Rename for clarity



NewPond1 <- expand.grid(
  sex = unique(IPHistory1$sex),
  TL = as.integer(seq(min(IPHistory1$TL, na.rm = TRUE),
                      max(IPHistory1$TL, na.rm = TRUE),
                      length.out = 50)),
  time = factor(1),  # First-year only
  first_year = 1     # Indicator for first-year survival
)

TLPredictions <- covariate.predictions(Pond1TLTimeModel, data = NewPond1, parameter = "Phi")
# NewPhiPond1 <- predict_real(Pond1TLTimeModel, ddl$Phi[1,,drop=FALSE], data = NewPond1, parameter = "Phi")

TLPredictions$estimates
