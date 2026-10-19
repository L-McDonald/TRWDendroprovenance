library(tidyverse)
library(sf)
library(terra)


# LOAD DATA---------------------------------------------------------------------
grid_df <- read.csv("Outputs/15_PredictionGrid_Coordinates.csv")
elev <- terra::rast("Outputs/EnvironmentalVariables/Elevation_3035.tif")

cru_files <- list.files(
  "Outputs/EnvironmentalVariables/",
  full.names = TRUE,
  pattern = "CRU.*_3035\\.tif$"
)
cru <- terra::rast(cru_files)

soil <- terra::rast("Outputs/EnvironmentalVariables/Soil_3035.tif")


# EXTRACT VALUES----------------------------------------------------------------
## search radius over which to look when a point falls on a missing value. Set to
## the resolution of each raster, i.e. can look one grid square over.
search_rad_elev <- res(elev)[[1]] * 2
search_rad_cru <- res(cru)[[1]] * 2
search_rad_soil <- res(soil)[[1]] * 2

grid_elev <- terra::extract(
  elev,
  grid_df,
  ID = FALSE,
  search_radius = search_rad_elev
) %>%
  select(-all_of(c("distance", "cell")))

grid_cru <- map_dfc(as.list(cru), \(l) {
  terra::extract(l, grid_df, ID = FALSE, search_radius = search_rad_cru) %>%
    select(-all_of(c("distance", "cell")))
})

grid_soil <- map_dfc(as.list(soil), \(l) {
  terra::extract(l, grid_df, ID = FALSE, search_radius = search_rad_soil) %>%
    select(-all_of(c("distance", "cell")))
})


# COMBINE INTO PREDICTION TABLE AND SAVE----------------------------------------
pred_table <- bind_cols(grid_df, grid_elev, grid_cru, grid_soil)

write.csv(
  pred_table,
  "Outputs/16_PredictionGrid_EnvironmentalVariables.csv",
  row.names = FALSE
)
