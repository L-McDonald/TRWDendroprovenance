library(tidyverse)
library(sf)
library(terra)


# LOAD DATA---------------------------------------------------------------------
crn <- readRDS("Outputs/12_CalculatedChronologies_Step2.Rdata")
cor_all <- read.csv("Outputs/13_RawCorrelations_Step2.csv")

elev <- terra::rast("Outputs/EnvironmentalVariables/Elevation_3035.tif")

cru_files <- list.files(
  "Outputs/EnvironmentalVariables/",
  full.names = TRUE,
  pattern = "CRU.*_3035\\.tif$"
)
cru <- terra::rast(cru_files)

soil <- terra::rast("Outputs/EnvironmentalVariables/Soil_3035.tif")


# UPSCALE SOIL RASTERS----------------------------------------------------------
## native resolution of the soil raster is approx 625 m but site coordinates are
## unlikely to be this precise. Upscale by a factor of 6. This gives a resolution
## of approx 3750 which is approx 1 tenth of the prediction grid resolution.
soil_up <- aggregate(soil, fact = 10, fun = "mean", na.rm = TRUE)


# GET ENVIRONMENTAL VALUES AT EACH SITE-----------------------------------------
coords <- crn %>%
  select(x = x_3035, y = y_3035)

elev_df <- crn %>%
  select(site_code, x = x_3035, y = y_3035, elev)

cru_df <- terra::extract(cru, coords, ID = FALSE)

soil_df <- map_dfc(as.list(soil_up), \(l) {
  terra::extract(l, coords, ID = FALSE)
})


# BIND INTO TRAINING TABLE AND SAVE---------------------------------------------
env <- bind_cols(elev_df, cru_df, soil_df)

write.csv(env, "Outputs/17_SiteEnvironmentalVariables.csv", row.names = FALSE)
