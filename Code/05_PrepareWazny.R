library(tidyverse)
library(dplR)
library(terra)
library(sf)

source("Code/00_Functions.R")


# LOAD DATA---------------------------------------------------------------------
meta_wazny <- read.csv("Data/metadata_Wazny.csv")
range <- read_sf("Outputs/EnvironmentalVariables/Region_4326.gpkg")


# READ RWL FILES----------------------------------------------------------------
rwl_files <- list.files("Data/rwl_Wazny/", full.names = TRUE, pattern = "\\.rwl")

wazny <- map_dfr(rwl_files, rwl_tibble) %>%
  left_join(meta_wazny, by = "site_code") %>%
  select(
    site_code,
    site_name,
    start,
    end,
    length,
    n_series,
    species,
    lat,
    lon,
    elev,
    rwl
  ) %>% 
  mutate(source = "Wazny", .before = everything())


# SAVE CHRONOLOGIES-------------------------------------------------------------
saveRDS(wazny, "Outputs/03_InitialChronologies_Wazny.Rdata")
