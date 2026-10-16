library(tidyverse)
library(sf)


source("Code/00_Functions.R")


# LOAD CHRONOLOGIES-------------------------------------------------------------
itrdb <- readRDS("Outputs/04_HistoricalSeriesRemoved_ITRDB.Rdata")
ctrl <- readRDS("Outputs/02_InitialChronologies_CTRL.Rdata")
wazny <- readRDS("Outputs/03_InitialChronologies_Wazny.Rdata")


# MERGE INTO ONE DATAFRAME------------------------------------------------------
chron <- itrdb %>% 
  bind_rows(ctrl) %>% 
  bind_rows(wazny)


# ADD EPSG 3035 COORDINATES-----------------------------------------------------
coords_3035 <- chron %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(3035) %>%
  st_coordinates()

chron <- chron %>%
  mutate(x_3035 = coords_3035[, "X"], y_3035 = coords_3035[, "Y"], .after = lon)


# SAVE--------------------------------------------------------------------------
saveRDS(chron, "Outputs/05_InitialChronologies_All.Rdata")
