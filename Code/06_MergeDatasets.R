library(tidyverse)


source("Code/00_Functions.R")


# LOAD CHRONOLOGIES-------------------------------------------------------------
itrdb <- readRDS("Data/04_HistoricalSeriesRemoved_ITRDB.Rdata")
ctrl <- readRDS("Data/02_InitialChronologies_CTRL.Rdata")
wazny <- readRDS("Data/03_InitialChronologies_Wazny.Rdata")


# MERGE INTO ONE DATAFRAME AND SAVE---------------------------------------------
chron <- itrdb %>% 
  bind_rows(ctrl) %>% 
  bind_rows(wazny)

saveRDS(chron, "Data/05_InitialChronologies_All.Rdata")
