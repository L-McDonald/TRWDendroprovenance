library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)


# LOAD STUDY REGION AND COAST---------------------------------------------------
region <- st_read("Outputs/EnvironmentalVariables/Region_3035.gpkg")
coast <- st_read("Outputs/EnvironmentalVariables/Coastline_3035.gpkg")


# GENERATE GRID OVER SPECIES DISTRIBUTION RANGE---------------------------------
## 37.5 Km resolution is approximately equal to 0.5 degrees
region_v <- vect(region)
grid_r <- rast(ext(region_v), resolution = 37500, crs = "EPSG:3035") %>%
  rasterize(region_v, y = ., field = 1)

grid_df <- as.data.frame(grid_r, xy = TRUE, na.rm = TRUE) %>%
  select(x, y)


# SAVE GRIDS--------------------------------------------------------------------
terra::writeRaster(grid_r, "Outputs/PredictionGrid/PredictionGrid_Raster.tif")
write.csv(
  grid_df,
  "Outputs/PredictionGrid/PredictionGrid_DataFrame_Coordinates.csv"
)