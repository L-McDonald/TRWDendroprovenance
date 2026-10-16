library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)


# LOAD STUDY REGION AND COAST---------------------------------------------------
region <- st_read(
  "Data/SpeciesDistributionRanges/QURO_QUPE_QUCE_QUFR_QUPU_3035.gpkg"
)

coast <- ne_countries(
  scale = 10,
  continent = c("Africa", "Europe", "Asia")
) %>%
  st_make_valid() %>%
  st_crop(c(xmin = -25, xmax = 90, ymin = 10, ymax = 90)) %>%
  st_geometry() %>%
  st_transform(crs = 3035)


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