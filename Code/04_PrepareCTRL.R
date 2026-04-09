library(tidyverse)
library(dplR)
library(terra)
library(sf)

source("Code/00_Functions.R")


# LOAD DATA---------------------------------------------------------------------

# load worldclim elevation and reproject to EPSG 3035
# elev <- terra::rast("Data/wc2.1_30s_elev/wc2.1_30s_elev.tif") %>%
#   terra::crop(c(xmin = -25, xmax = 90, ymin = 10, ymax = 90)) %>%
#   terra::project("EPSG:3035")

# writeRaster(elev, "Data/wc2.1_30s_elev/wc2.1_30s_elev_reprojected3035.tif")

elev <- terra::rast("Data/wc2.1_30s_elev/wc2.1_30s_elev.tif")
meta_ctrl <- read.csv("Data/metadata_CTRL.csv")
range <- read_sf("Data/SpeciesDistributionRanges/QURO_QUPE_QUCE_QUFR.gpkg")


# READ RWL FILES----------------------------------------------------------------
rwl_files <- list.files("Data/rwl_CTRL/", full.names = TRUE, pattern = "\\.rwl")

ctrl <- map_dfr(rwl_files, rwl_tibble) %>%
  left_join(meta_ctrl, by = "site_code") %>%
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
    elev_final,
    elev_initial,
    rwl
  )


# EXTRACT ELEVATION IF UNKNOWN--------------------------------------------------
ctrl_elev <- ctrl %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  mutate(
    lat = st_coordinates(.)[, "Y"],
    lon = st_coordinates(.)[, "X"],
    .after = species
  ) %>%
  mutate(
    elev = if_else(
      is.na(elev_final),
      terra::extract(.env$elev, vect(.), method = "bilinear")[, 2],
      elev_final
    ),
    .after = lon
  ) %>% 
  select(-elev_initial, -elev_final) %>% 
  st_drop_geometry() %>% 
  mutate(source = "CTRL", .before = everything())


# SAVE CHRONOLOGIES
saveRDS(ctrl_elev, "Data/InitialChronologies_CTRL.Rdata")
