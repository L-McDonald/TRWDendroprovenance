library(tidyverse)
library(sf)
library(dplR)

source("Code/00_Functions.R")


# LOAD DATA---------------------------------------------------------------------
crn <- readRDS("Outputs/08_CalculatedChronologies_Step1.Rdata")
corr <- read.csv("Outputs/09_RawCorrelations_Step1.csv")


# IDENTIFY VERY HIGH CORRELATING SITES------------------------------------------
corr %>%
  filter(cor >= 0.70 | dist_m <= 10000) %>%
  arrange(desc(cor))


# COMBINE SITES-----------------------------------------------------------------
clusters <- tibble(
  cluster = seq(1:3),
  members = list(
    c("FRAN003", "FRAN049"),
    c("GERM170", "GERM168_GERM169"),
    c("DEU314", "DEU306_DEU312")
  )
)

merged_rows <- clusters %>%
  mutate(
    members_df = map(members, ~ crn %>% filter(site_code %in% .x)),
    merged = map(members_df, merge_cluster)
  ) %>%
  pull(merged) %>%
  bind_rows() %>%
  st_as_sf(coords = c("x_3035", "y_3035"), crs = 3035, remove = FALSE)

merged_codes <- unlist(clusters$members)

crn_merged_sites <- crn %>%
  filter(!site_code %in% merged_codes) %>%
  mutate(merged_from = NA_character_) %>%
  bind_rows(merged_rows) %>%
  arrange(source, site_code)


# REMOVE POLA006----------------------------------------------------------------
## POLA006 is a regional chronology including many historical sites. POLA005 is
## most of the forest chronologies from POLA006 so retain.

crn_filtered <- crn_merged_sites %>% 
  filter(site_code != "POLA006")


# CHECK FOR PROBLEM SERIES------------------------------------------------------
## removes problem series and flags chronologies which have problems with more
## than half of their series
crn_checked <- remove_problem_series(crn_filtered) %>% 
    filter(problem != TRUE) %>%
    select(-problem, -series_removed, -merged_from)


# SAVE CLEANED CHRONOLOGIES-----------------------------------------------------
saveRDS(crn_checked, "Outputs/10_CleanedChronologies_Step2.Rdata")
