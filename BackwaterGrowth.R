load("data/ReportingData.RData")

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(ggplot2) # Plotting
packages(tidyr) # replace_na for dataframes function plus others
packages(openxlsx) # package openxlsx is required to create to Excel files

StudyBWRecaptures <- StudyBWNFWGAnalysis %>%
  select(first_date, location, disposition, event, species, PITIndex, sex, total_length,
         location_id, MaxScanDate, MaxDAL, Survived) %>%
  inner_join(StudyBWNFWG %>% 
               filter(FirstRecord == "no", !is.na(total_length)) %>%
               select(collection_date, PITIndex, location_id, RecapTL = total_length),
             by = c("PITIndex", "location_id")) %>%
  mutate(RecapDAL = as.numeric(difftime(collection_date, first_date, units = "days")),
         RecapGrowthmm = RecapTL - total_length) %>%
  filter(RecapDAL > 365)

StudyBWGrowthSummary <- StudyBWGrowth %>%
  group_by(location, species) %>%
  summarise(Count = n(), 
            Growthmin = as.integer(min(Growth_mm_year)), 
            Growthmean = as.integer(mean(Growth_mm_year)), 
            meanYAL = round(mean(YAL), 3)) %>%
  ungroup()

yoyHistogramDataIP <- StudyBWAnalysis %>%
  filter(event == "capture", !is.na(total_length), location_id != 592) %>%
  mutate(CaptureYear = year(first_date))

yoyHistogram <- ggplot(yoyHistogramDataIP, aes(x = total_length, fill = factor(CaptureYear))) +
  geom_histogram(binwidth = 10, position = "stack", color = "black") +
  facet_wrap(~ location) +
  scale_fill_viridis_d(name = "Year") +
  labs(
    x = "Total Length (mm)",
    y = "Count") +
  coord_cartesian(ylim = c(0, 130)) +
  theme_minimal()

yoyHistogram

yoyHistogramS <- ggplot(yoyHistogramDataIP %>%
                          filter(Survived == 1), aes(x = total_length, fill = factor(CaptureYear))) +
  geom_histogram(binwidth = 10, position = "stack", color = "black") +
  facet_wrap(~ location) +
  scale_fill_viridis_d(name = "Year") +
  labs(
    x = "Total Length (mm)",
    y = "Count") +
  coord_cartesian(ylim = c(0, 130)) +
  theme_minimal()

yoyHistogramS

yoyHistogramS2024 <- ggplot(yoyHistogramDataIP %>%
                          filter(SurvivedFY24 == 1), aes(x = total_length, fill = factor(CaptureYear))) +
  geom_histogram(binwidth = 10, position = "stack", color = "black") +
  facet_wrap(~ location) +
  scale_fill_viridis_d(name = "Year") +
  labs(
    x = "Total Length (mm)",
    y = "Count") +
  coord_cartesian(ylim = c(0, 130)) +
  theme_minimal()
# Histogram for all captured fish for each location and 
# color coded by harvest year.
yoyHistogramS2024
