library(tidyverse)
library(sf)
library(terra)


# LOAD STUDY REGION AND COAST---------------------------------------------------
region <- st_read("Outputs/EnvironmentalVariables/Region_3035.gpkg")
coast <- st_read("Outputs/EnvironmentalVariables/Coastline_3035.gpkg")

sites <- readRDS("Outputs/12_CalculatedChronologies_Step2.Rdata") %>%
  st_as_sf(coords = c("x_3035", "y_3035"), remove = FALSE, crs = 3035)


# GENERATE GRID OVER SPECIES DISTRIBUTION RANGE---------------------------------
## 37.5 Km resolution is approximately equal to 0.5 degrees
region_v <- vect(region)
grid_r <- rast(ext(region_v), resolution = 37500, crs = "EPSG:3035") %>%
  rasterize(region_v, y = ., field = 1)


# MAKE SURE THAT EVERY SITE IS WITHIN THE GRID----------------------------------
## because the grid is quite coarse some sites fall outside of it even though 
## they are within the region. This adds cells to the raster to make sure that 
## all sites are within the grid.
sites_v <- vect(st_transform(sites, crs(grid_r)))
grid_r <- cover(grid_r, rasterize(sites_v, grid_r, field = 1))

## RUS328's cell is isolated; connect it to the nearest cell with a one-cell
## buffered bridge
## RUS328 is 33.45 km from the region and the cell added under it ends up as an
## island. This connects it to the nearest cell with a one-cell bridge.
site_rus328 <- sites_v[sites_v$site_code == "RUS328", ]
main_r <- grid_r
main_r[cellFromXY(grid_r, crds(site_rus328))] <- NA
dist_r <- distance(main_r, site_rus328)
dist_r[is.na(main_r)] <- NA
nearest_xy <- xyFromCell(dist_r, which.min(values(dist_r)))
bridge <- vect(
  rbind(crds(site_rus328), nearest_xy),
  type = "lines",
  crs = crs(grid_r)
)
grid_r <- cover(
  grid_r,
  rasterize(buffer(bridge, 37500), grid_r, field = 1, touches = TRUE)
)

## the one-cell bridge is artificial s smooth the edge around RUS328: add cells 
# where a majority of the 7x7 neighbourhood is filled. Restricted to +/- 8 cells
# of the site and only ever adds cells, never removes them
smooth_r <- grid_r
smooth_r[is.na(smooth_r)] <- 0
smooth_r <- focal(smooth_r, w = 7, fun = "mean", na.rm = TRUE) > 0.5
site_xy <- crds(site_rus328)
smooth_win <- ext(
  site_xy[1] - 8 * 37500, site_xy[1] + 8 * 37500,
  site_xy[2] - 8 * 37500, site_xy[2] + 8 * 37500
)
smooth_r <- extend(crop(smooth_r, smooth_win), grid_r)
grid_r <- cover(grid_r, ifel(smooth_r == 1, 1, NA))


# GRID COORDINATES--------------------------------------------------------------
grid_df <- as.data.frame(grid_r, xy = TRUE, na.rm = TRUE) %>%
  select(x, y)


# SAVE GRIDS--------------------------------------------------------------------
terra::writeRaster(grid_r, "Outputs/14_PredictionGrid_Raster.tif", overwrite = TRUE)
write.csv(
  grid_df,
  "Outputs/15_PredictionGrid_Coordinates.csv",
  row.names = FALSE
)
