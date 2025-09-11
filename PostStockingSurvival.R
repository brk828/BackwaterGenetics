# Post-stocking survival analysis for IP 1-6 and Yuma Cove Backwater
# B. Kesner September 2025

# Mark-recapture based on Months since release with FY of release, Sex, Size, and Backwater as factors
# Logistic regression based on SurvivalDAL as well

load("data/ReportingData.RData")
packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(glmmTMB) # General linear mixed model analysis built on TMB automatic differentiation engine
packages(ggplot2) # Plotting
packages(tidyr) # pivot table
packages(RMark) # MCR analysis
packages(stringr) # string manipulations and function

StockingBW <- StudyBWNFWGTagging %>%
  filter(event == "stocking", collection_date < as.Date("2024-10-01"), sex != "U") %>%
  rename(ReleaseDate = collection_date) %>%
  mutate(ReleaseFY = as.factor(ifelse(month(ReleaseDate) > 9, year(ReleaseDate)+1, year(ReleaseDate))),
         ReleaseMonth = as.integer(month(ReleaseDate)),
         ReleaseMonthName = format(ReleaseDate, "%b"),
         ReleaseMonthFY = ifelse(ReleaseMonth>9, ReleaseMonth-9, ReleaseMonth+3),
         sex = as.factor(trimws(sex)),
         location = as.factor(location),
         TLClass = as.integer(total_length*.1)*10) 

StockingBWSurvivalSummary <- StockingBW %>%
  group_by(location, location_id, ReleaseFY, species, sex, TLClass) %>%
  summarise(Released = n(), Survivors = sum(Survived), ContactedProp = Survivors/Released) %>%
  ungroup() %>%
  filter(Released >= 5)

# Yuma stocking in 2020 did not have proper sexing done
YumaStocking <- StockingBW %>%
  filter(location_id == 592, ReleaseFY != 2020) %>%
  select(PITIndex, ReleaseFY, sex, total_length, Survived, SurvivedFY24)

YumaReleaseFY <- unique(YumaStocking$ReleaseFY)
#zero inflated model for Yuma Cove backwater.
DALModelYuma <- glmmTMB(Survived ~ sex * total_length * ReleaseFY,
                    family = binomial(link = 'logit'), 
                    data = YumaStocking)

#creating a predictions data.frame from model for graphing
PredictorsYuma <- expand.grid(total_length = as.integer(seq(350, 500, by = 10)), 
                          ReleaseFY =YumaReleaseFY,
                          sex = as.factor(c("M", "F")))

PredictionsYuma <- predict(DALModelYuma, newdata = PredictorsYuma, type = "response", se.fit = TRUE)
lowerCI <- PredictionsYuma$fit - (1.96 * PredictionsYuma$se.fit)
upperCI <- PredictionsYuma$fit + (1.96 * PredictionsYuma$se.fit)
fit <- PredictionsYuma$fit
PredictedYuma <- cbind(PredictorsYuma, fit, lowerCI, upperCI)

# cleanup
rm(PredictorsYuma, PredictionsYuma, lowerCI, upperCI, fit)

# Yuma Figure
YumaDALSurvival <- ggplot(PredictedYuma %>% filter(ReleaseFY != 2015), 
                     aes(x = total_length, y = fit, color = sex)) + 
  geom_line(linewidth = 1) + 
  labs(x = 'Total Length (mm)', y = 'Detection Probability', color = 'sex') + 
  theme(plot.margin = margin(.75,.75,.75,.75, unit = 'cm'), 
        axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 5)) +
  geom_ribbon(aes(x = total_length, ymin = lowerCI, ymax = upperCI, fill = sex), alpha = 0.3) + 
  facet_wrap(~ReleaseFY, nrow=3) +
  geom_point(data = StockingBWSurvivalSummary %>%
               filter(location_id == 592, ReleaseFY != 2020, 
                      ReleaseFY != 2015, TLClass >=350), 
             aes(x = TLClass, y = ContactedProp)) 


YumaDALSurvival

# IP XYTE stockings 
IPXYTEStocking <- StockingBW %>%
  filter(location_id != 592, species == "XYTE") %>%
  select(location, PITIndex, ReleaseFY, sex, total_length, Survived, SurvivedFY24)

IPXYTELocations <- unique(IPXYTEStocking$location)

#zero inflated model for IP XYTE stockings
DALModelIPXYTE <- glmmTMB(Survived ~ sex * total_length * location,
                        family = binomial(link = 'logit'), 
                        data = IPXYTEStocking)

PredictorsIPXYTE <- expand.grid(total_length = as.integer(seq(350, 500, by = 10)), 
                              location = IPXYTELocations,
                              sex = as.factor(c("M", "F")))

PredictionsIPXYTE <- predict(DALModelIPXYTE, newdata = PredictorsIPXYTE, type = "response", se.fit = TRUE)
lowerCI <- PredictionsIPXYTE$fit - (1.96 * PredictionsIPXYTE$se.fit)
upperCI <- PredictionsIPXYTE$fit + (1.96 * PredictionsIPXYTE$se.fit)
fit <- PredictionsIPXYTE$fit
PredictedIPXYTE <- cbind(PredictorsIPXYTE, fit, lowerCI, upperCI)

# cleanup
rm(PredictionsIPXYTE, PredictorsIPXYTE, lowerCI, upperCI, fit)

# IPXTYE figure
IPXYTEDALSurvival <- ggplot(PredictedIPXYTE, 
                          aes(x = total_length, y = fit, color = sex)) + 
  geom_line(linewidth = 1) + 
  labs(x = 'Total Length (mm)', y = 'Detection Probability', color = 'sex') + 
  theme(plot.margin = margin(.75,.75,.75,.75, unit = 'cm'), 
        axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 5)) +
  geom_ribbon(aes(x = total_length, ymin = lowerCI, ymax = upperCI, fill = sex), alpha = 0.3) + 
  facet_wrap(~location, nrow=3) +
  geom_point(data = StockingBWSurvivalSummary %>%
               filter(location_id != 592, species == "XYTE", TLClass >=350), 
             aes(x = TLClass, y = ContactedProp))

IPXYTEDALSurvival

# IP GIEL 2017 stockings 
IPGIEL2017Stocking <- StockingBW %>%
  filter(location_id != 592, species == "GIEL", ReleaseFY == 2017) %>%
  select(location, PITIndex, ReleaseFY, sex, total_length, Survived, SurvivedFY24)

IPGIEL2017Locations <- unique(IPGIEL2017Stocking$location)

#zero inflated model for IP GIEL 2017
DALModelIPGIEL2017 <- glmmTMB(Survived ~ sex * total_length * location,
                          family = binomial(link = 'logit'), 
                          data = IPGIEL2017Stocking)

PredictorsIPGIEL2017 <- expand.grid(total_length = as.integer(seq(230, 290, by = 5)), 
                                location = IPGIEL2017Locations,
                                sex = as.factor(c("M", "F")))

PredictionsIPGIEL2017 <- predict(DALModelIPGIEL2017, newdata = PredictorsIPGIEL2017, type = "response", se.fit = TRUE)
lowerCI <- PredictionsIPGIEL2017$fit - (1.96 * PredictionsIPGIEL2017$se.fit)
upperCI <- PredictionsIPGIEL2017$fit + (1.96 * PredictionsIPGIEL2017$se.fit)
fit <- PredictionsIPGIEL2017$fit
PredictedIPGIEL2017 <- cbind(PredictorsIPGIEL2017, fit, lowerCI, upperCI)

# cleanup
rm(PredictionsIPGIEL2017, PredictorsIPGIEL2017, lowerCI, upperCI, fit)

# IPGIEL 2017 figure
IPGIEL2017DALSurvival <- ggplot(PredictedIPGIEL2017, 
                            aes(x = total_length, y = fit, color = sex)) + 
  geom_line(linewidth = 1) + 
  labs(x = 'Total Length (mm)', y = 'Detection Probability', color = 'sex') + 
  theme(plot.margin = margin(.75,.75,.75,.75, unit = 'cm'), 
        axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 5)) +
  geom_ribbon(aes(x = total_length, ymin = lowerCI, ymax = upperCI, fill = sex), alpha = 0.3) + 
  facet_wrap(~location, nrow=3) +
  geom_point(data = StockingBWSurvivalSummary %>%
               filter(location_id != 592, species == "GIEL", ReleaseFY == 2017), 
             aes(x = TLClass, y = ContactedProp))

IPGIEL2017DALSurvival


png(paste0("output/YumaCoveSurvival", SurvivalDAL, "DaysPostStocking.png"), width = 6, height = 4, units = 'in', res = 300)   
YumaDALSurvival
dev.off()

png(paste0("output/IPXYTESurvival", SurvivalDAL, "DaysPostStocking.png"), width = 6, height = 4, units = 'in', res = 300)   
IPXYTEDALSurvival
dev.off()

png(paste0("output/IPGIEL2017Survival", SurvivalDAL, "DaysPostStocking.png"), width = 6, height = 4, units = 'in', res = 300)   
IPGIEL2017DALSurvival
dev.off()


# Razorback sucker stocked in Winter, bonytail in spring. First survival period is from
# stocking to fall-winter of next FY. Bonytail and razorbak sucker stocked in FY2017.
# First session FY2018 months Oct through January (ScanMonthFY 1-4). Survival determined annually 
# based on Scanning through those months and simple CMS model for stocked fish

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
  mutate(ch = as.character(ch), sex = as.factor(sex), Location = as.factor(Location))


dp=RMark::process.data(IPHistory,model="CJS")
ddl=RMark::make.design.data(dp)

# Building a max 2 age structure
ddl$Phi$age2 <- ifelse(ddl$Phi$Age >2,3,ddl$Phi$age)
ddl$Phi$age2 <- as.factor(ddl$Phi$age2)

# Create two time sample structure as earlier ages have smaller scan windows
ddl$p$time2 <- ifelse(ddl$p$Age > 2, 3, 1)
ddl$p$time2 <- as.factor(ddl$p$time2)

Phi.age2 = list(formula=~age2) 
p.dot = list(formula=~1)
p.time2 = list(formula=~time2)

AgeModel <- mark(data=dp,ddl=ddl,model.parameters=list(Phi=Phi.age2,p=p.time2))

AgeResults <- AgeModel$results$real %>%
  tibble::rownames_to_column(var = "parameter")  # Store row names in a column

# Extract survival (S) estimates by age
AgeSurvival <- AgeResults %>%
  filter(str_starts(parameter, "Phi")) %>%
  mutate(age = row_number(),  # Assign sequential age numbers
         parameter = "Survival")  # Rename for clarity

# Extract recapture (p) estimates by time
TimeRecapture <- AgeResults %>%
  filter(str_starts(parameter, "p")) %>%
  mutate(time = row_number(),  # Assign sequential time numbers
         parameter = "Recapture")  # Rename for clarity

