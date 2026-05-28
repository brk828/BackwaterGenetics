# Minimum known survivors for IP 1-6 and Yuma Cove backwater
# B. Kesner 11 Sep 2025

load("data/ReportingData.RData")

packages(dplyr)     # data manipulation
packages(lubridate) # date and time manipulation
packages(ggplot2) # Plotting
packages(tidyr) # replace_na for dataframes function plus others
packages(openxlsx) # package openxlsx is required to create to Excel files

# If MaxScanDay is null, make it the release or tagging date of the fish
# change juvenile to unknown sex to simplify sexes
StudyBWAnalysis <- StudyBWAnalysis %>%
  mutate(sex = if_else(sex == "J", "U", sex))

# separate out the 3 study groups for analysis
StudyBWXYTEIP <- StudyBWAnalysis %>% 
  filter(location_id != 592, species == "XYTE", MaxDAL > SurvivalDAL)

StudyBWGIELIP <- StudyBWAnalysis %>% 
  filter(location_id != 592, species == "GIEL", MaxDAL > SurvivalDAL)

StudyBWYuma <- StudyBWAnalysis %>% 
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
  full_join(StudyBWGIELIP %>% mutate(key = 1), 
            by = "key", relationship = "many-to-many") %>%
  select(-key)

# filter out for each day the fish that were tagged or released prior to the Date
# their MaxScanDate is after the Date, and they were know to survive past the SurvivalDAL
KnownSurvivalYuma <- CrossDFYuma %>% 
  filter(first_date <= Date &
           MaxScanDate >= Date, 
         Survived == 1) %>%
  select(Date, PITIndex, first_date, sex, total_length, event, MaxScanDate, MaxDAL) 

SurvivorCountsYuma <- expand_grid(
  Date = SurvDaysYuma$Date,
  sex = unique(KnownSurvivalYuma$sex)) %>% 
  left_join(KnownSurvivalYuma %>%
      count(Date, sex, name = "Count"),
    by = c("Date", "sex")) %>%
  mutate(Count = replace_na(Count, 0))

KnownSurvivalIPXYTE <- CrossDFIPXYTE %>% 
  filter(first_date <= Date &
           MaxScanDate >= Date, 
         Survived == 1) %>%
  select(location, Date, PITIndex, first_date, sex, total_length, event, MaxScanDate, MaxDAL) 

SurvivorCountsIPXYTE <- expand_grid(
  Date = SurvDaysXYTEIP$Date,
  sex = unique(KnownSurvivalIPXYTE$sex),
  location = unique(KnownSurvivalIPXYTE$location)) %>% 
  left_join(KnownSurvivalIPXYTE %>%
      count(Date, sex, location, name = "Count"),
    by = c("Date", "sex", "location")) %>%
  mutate(Count = replace_na(Count, 0))

KnownSurvivalIPGIEL <- CrossDFIPGIEL %>% 
  filter(first_date <= Date &
           MaxScanDate >= Date, 
         Survived == 1) %>%
  select(location, Date, PITIndex, first_date, sex, total_length, event, MaxScanDate, MaxDAL) 

SurvivorCountsIPGIEL <- expand_grid(
  Date = SurvDaysGIELIP$Date,
  sex = unique(KnownSurvivalIPGIEL$sex),
  location = unique(KnownSurvivalIPGIEL$location)) %>%
  left_join(KnownSurvivalIPGIEL %>%
      count(Date, sex, location, name = "Count"),
    by = c("Date", "sex", "location")) %>%
  mutate(Count = replace_na(Count, 0))

TotalCountYuma <- SurvivorCountsYuma %>%
  group_by(Date) %>%
  summarise(Count = sum(Count, na.rm = TRUE)) %>%
  mutate(sex = "Total") %>%
  ungroup()

TotalCountIPXYTE <- SurvivorCountsIPXYTE %>%
  group_by(Date, location) %>%
  summarise(Count = sum(Count, na.rm = TRUE)) %>%
  mutate(sex = "Total") %>%
  ungroup()

TotalCountIPGIEL <- SurvivorCountsIPGIEL %>%
  group_by(Date, location) %>%
  summarise(Count = sum(Count, na.rm = TRUE)) %>%
  mutate(sex = "Total") %>%
  ungroup()

SurvivorCountsYuma$sex <- factor(SurvivorCountsYuma$sex, levels = c("Total", "F", "M", "U"))
TotalCountYuma$sex <- factor(TotalCountYuma$sex, levels = c("Total", "F", "M", "U"))

SurvivorCountsIPXYTE$sex <- factor(SurvivorCountsIPXYTE$sex, levels = c("Total", "F", "M", "U"))
TotalCountIPXYTE$sex <- factor(TotalCountIPXYTE$sex, levels = c("Total", "F", "M", "U"))

SurvivorCountsIPGIEL$sex <- factor(SurvivorCountsIPGIEL$sex, levels = c("Total", "F", "M", "U"))
TotalCountIPGIEL$sex <- factor(TotalCountIPGIEL$sex, levels = c("Total", "F", "M", "U"))

KnownSurvivalPlotYuma <- ggplot(SurvivorCountsYuma, aes(x = Date, y = Count)) +
  geom_line(aes(color = sex, linetype = sex)) +
  scale_x_date(date_breaks = "1 month",
               labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), "")) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "black")) +
  labs(x = "Date", y = "Count") +
  geom_line(data = TotalCountYuma,
            aes(x = Date, y = Count, color = sex, linetype = sex),
            linewidth = 1.2, alpha = 0.4)

KnownSurvivalPlotIPXYTE <- ggplot(SurvivorCountsIPXYTE, aes(x = Date, y = Count)) +
  geom_line(aes(color = sex, linetype = sex)) +
  scale_x_date(date_breaks = "1 month",
               labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), "")) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "black")) +
  labs(x = "Date", y = "Count") +
  facet_wrap(~location, nrow=3) +
  geom_line(data = TotalCountIPXYTE,
            aes(x = Date, y = Count, color = sex, linetype = sex),
            linewidth = 1.2, alpha = 0.4)

KnownSurvivalPlotIPGIEL <- ggplot(SurvivorCountsIPGIEL, aes(x = Date, y = Count)) +
  geom_line(aes(color = sex, linetype = sex)) +
  scale_x_date(date_breaks = "1 month",
               labels = function(x) ifelse(month(x) == 1, format(x, "%Y"), "")) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "black")) +
  labs(x = "Date", y = "Count") +
  facet_wrap(~location, nrow=3) +
  geom_line(data = TotalCountIPGIEL,
            aes(x = Date, y = Count, color = sex, linetype = sex),
            linewidth = 1.2, alpha = 0.4)

KnownSurvivalPlotIPGIEL

# --- Overlay annotations: new fish additions (+) and harvested fish removals (-) ---
# Additions: captured, tagged, and released back (newly entered tagged population)
# Removals: previously tagged fish captured and retained (FirstRecord == "no" excludes fish
#           tagged at the time of harvest, which were never part of the tagged population)

# Additions: captured-and-released fish + all stockings after the initial one per location.
# The initial stocking is excluded by filtering to dates after the minimum stocking date.
YumaAdditions <- bind_rows(
  StudyBWNFWG |>
    filter(location_id == 592, event == "capture", disposition == "released") |>
    count(Date = collection_date, name = "n"),
  StudyBWNFWG |>
    filter(location_id == 592, event == "stocking") |>
    filter(collection_date > min(collection_date)) |>
    count(Date = collection_date, name = "n")
) |>
  group_by(Date) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(label = paste0("+", n))

IPXYTEAdditions <- bind_rows(
  StudyBWNFWG |>
    filter(location_id != 592, species == "XYTE", event == "capture", disposition == "released") |>
    count(Date = collection_date, location, name = "n"),
  StudyBWNFWG |>
    filter(location_id != 592, species == "XYTE", event == "stocking") |>
    group_by(location) |>
    filter(collection_date > min(collection_date)) |>
    ungroup() |>
    count(Date = collection_date, location, name = "n")
) |>
  group_by(Date, location) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(label = paste0("+", n))

IPGIELAdditions <- bind_rows(
  StudyBWNFWG |>
    filter(location_id != 592, species == "GIEL", event == "capture", disposition == "released") |>
    count(Date = collection_date, location, name = "n"),
  StudyBWNFWG |>
    filter(location_id != 592, species == "GIEL", event == "stocking") |>
    group_by(location) |>
    filter(collection_date > min(collection_date)) |>
    ungroup() |>
    count(Date = collection_date, location, name = "n")
) |>
  group_by(Date, location) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(label = paste0("+", n))

# Removals: previously tagged fish captured and retained (not tagged at time of harvest)
YumaRemovals <- StudyBWNFWG |>
  filter(location_id == 592, event == "capture", disposition == "retained", FirstRecord == "no") |>
  count(Date = collection_date, name = "n") |>
  mutate(label = paste0("-", n))

IPXYTERemovals <- StudyBWNFWG |>
  filter(location_id != 592, species == "XYTE", event == "capture", disposition == "retained", FirstRecord == "no") |>
  count(Date = collection_date, location, name = "n") |>
  mutate(label = paste0("-", n))

IPGIELRemovals <- StudyBWNFWG |>
  filter(location_id != 592, species == "GIEL", event == "capture", disposition == "retained", FirstRecord == "no") |>
  count(Date = collection_date, location, name = "n") |>
  mutate(label = paste0("-", n))

KnownSurvivalPlotYuma <- KnownSurvivalPlotYuma +
  geom_text(data = YumaAdditions, aes(x = Date, y = Inf, label = label),
            vjust = 1.5, size = 2.5, color = "black", inherit.aes = FALSE) +
  geom_text(data = YumaRemovals, aes(x = Date, y = Inf, label = label),
            vjust = 3.5, size = 2.5, color = "red", inherit.aes = FALSE)

KnownSurvivalPlotIPXYTE <- KnownSurvivalPlotIPXYTE +
  geom_text(data = IPXYTEAdditions, aes(x = Date, y = Inf, label = label),
            vjust = 1.5, size = 2.5, color = "black", inherit.aes = FALSE) +
  geom_text(data = IPXYTERemovals, aes(x = Date, y = Inf, label = label),
            vjust = 3.5, size = 2.5, color = "red", inherit.aes = FALSE)

KnownSurvivalPlotIPGIEL <- KnownSurvivalPlotIPGIEL +
  geom_text(data = IPGIELAdditions, aes(x = Date, y = Inf, label = label),
            vjust = 1.5, size = 2.5, color = "black", inherit.aes = FALSE) +
  geom_text(data = IPGIELRemovals, aes(x = Date, y = Inf, label = label),
            vjust = 3.5, size = 2.5, color = "red", inherit.aes = FALSE)

KnownSurvivalPlotYuma
KnownSurvivalPlotIPXYTE
KnownSurvivalPlotIPGIEL

rm(SurvivorCountsIPXYTE, SurvivorCountsIPGIEL, SurvivorCountsYuma, CrossDFYuma, CrossDFIPXYTE, 
   CrossDFIPGIEL, SurvDaysXYTEIP, SurvDaysGIELIP, SurvDaysYuma)

png(paste0("output/KnownSurvivalPlotYuma.png"), width = 8, height = 5, units = 'in', res = 300)   
KnownSurvivalPlotYuma
dev.off()

png(paste0("output/KnownSurvivalPlotIPXYTE.png"), width = 8, height = 5, units = 'in', res = 300)   
KnownSurvivalPlotIPXYTE
dev.off()

png(paste0("output/KnownSurvivalPlotIPGIEL.png"), width = 8, height = 5, units = 'in', res = 300)   
KnownSurvivalPlotIPGIEL
dev.off()

save(KnownSurvivalYuma, KnownSurvivalIPXYTE, KnownSurvivalIPGIEL,
     TotalCountIPGIEL, TotalCountIPXYTE, TotalCountYuma, packages, 
     file = "data/KnownSurvivalCounts.RData")



