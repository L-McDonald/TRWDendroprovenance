library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)


# LOAD REGION-------------------------------------------------------------------
region <- read_sf(
  "Data/SpeciesDistributionRanges/QURO_QUPE_QUCE_QUFR_QUPU.gpkg"
)


# COASTLINE---------------------------------------------------------------------
## load from natural earth, crop to approximate region and convert to EPSG:3035
coast <- ne_countries(
  scale = 10,
  continent = c("Africa", "Europe", "Asia")
) %>%
  st_make_valid() %>%
  st_crop(c(xmin = -15, xmax = 65, ymin = 20, ymax = 70)) %>%
  st_geometry()

coast_3035 <- st_transform(coast, crs = 3035)

st_write(coast, "Outputs/EnvironmentalVariables/Coastline_4326.gpkg")
st_write(coast_3035, "Outputs/EnvironmentalVariables/Coastline_3035.gpkg")


# ELEVATION---------------------------------------------------------------------
## load worldclim, crop to region with buffer and convert to EPSG:3035
elev <- terra::rast(
  "Data/EnvironmentalVariables/wc2.1_30s_elev/wc2.1_30s_elev.tif"
) %>%
  terra::crop(st_buffer(region, dist = 5))

elev_3035 <- terra::project(elev, "EPSG:3035")

writeRaster(elev, "Outputs/EnvironmentalVariables/Elevation_4326.tif")
writeRaster(elev_3035, "Outputs/EnvironmentalVariables/Elevation_3035.tif")


# CRU---------------------------------------------------------------------------
## load, crop to buffered region, select relevant layers, calculate monthly mean,
## convert to EPSG:3035

##Precipitation
pre <- terra::rast(
  "Data/EnvironmentalVariables/CRU_Unzipped/cru_ts4.10.1901.2025.pre.dat.nc"
) %>%
  terra::crop(st_buffer(region, dist = 5)) %>% 
  select(starts_with("pre_")) %>% 
  tapp(index = "months", fun = mean, na.rm = TRUE) %>%
  setNames(paste0("pre_", month.abb))

pre_3035 <- terra::project(pre, "EPSG:3035")

writeRaster(pre, "Outputs/EnvironmentalVariables/CRU_Precipitation_4326.tif")
writeRaster(pre_3035, "Outputs/EnvironmentalVariables/CRU_Precipitation_3035.tif")

## Minimum Temperature
tmn <- terra::rast(
  "Data/EnvironmentalVariables/CRU_Unzipped/cru_ts4.10.1901.2025.tmn.dat.nc"
) %>%
  terra::crop(st_buffer(region, dist = 5)) %>% 
  select(starts_with("tmn_")) %>% 
  tapp(index = "months", fun = mean, na.rm = TRUE) %>%
  setNames(paste0("tmn_", month.abb))

tmn_3035 <- terra::project(tmn, "EPSG:3035")

writeRaster(tmn, "Outputs/EnvironmentalVariables/CRU_MinimumTemperature_4326.tif")
writeRaster(tmn_3035, "Outputs/EnvironmentalVariables/CRU_MinimumTemperature_3035.tif")

## Average Temperature
tmp <- terra::rast(
  "Data/EnvironmentalVariables/CRU_Unzipped/cru_ts4.10.1901.2025.tmp.dat.nc"
) %>%
  terra::crop(st_buffer(region, dist = 5)) %>% 
  select(starts_with("tmp_")) %>% 
  tapp(index = "months", fun = mean, na.rm = TRUE) %>%
  setNames(paste0("tmp_", month.abb))

tmp_3035 <- terra::project(tmp, "EPSG:3035")

writeRaster(tmp, "Outputs/EnvironmentalVariables/CRU_AverageTemperature_4326.tif")
writeRaster(tmp_3035, "Outputs/EnvironmentalVariables/CRU_AverageTemperature_3035.tif")

## Maximum Temperature
tmx <- terra::rast(
  "Data/EnvironmentalVariables/CRU_Unzipped/cru_ts4.10.1901.2025.tmx.dat.nc"
) %>%
  terra::crop(st_buffer(region, dist = 5)) %>% 
  select(starts_with("tmx_")) %>% 
  tapp(index = "months", fun = mean, na.rm = TRUE) %>%
  setNames(paste0("tmx_", month.abb))

tmx_3035 <- terra::project(tmx, "EPSG:3035")

writeRaster(tmx, "Outputs/EnvironmentalVariables/CRU_MaximumTemperature_4326.tif")
writeRaster(tmx_3035, "Outputs/EnvironmentalVariables/CRU_MaximumTemperature_3035.tif")


# SOIL GRIDS--------------------------------------------------------------------
files <- list.files("Data/EnvironmentalVariables/SoilGrids/", full.names = TRUE)

soil <- terra::rast(files) %>% 
  terra::crop(st_buffer(region, dist = 5)) %>% 
  setNames(str_extract(names(.), "^[^_]+"))

soil_3035 <- terra::project(soil, "EPSG:3035")

writeRaster(soil, "Outputs/EnvironmentalVariables/Soil_4326.tif")
writeRaster(soil_3035, "Outputs/EnvironmentalVariables/Soil_3035.tif")