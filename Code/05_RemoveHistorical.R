library(tidyverse)
library(dplR)

source("Code/00_Functions.R")


# LOAD DATA---------------------------------------------------------------------
itrdb <- readRDS("Data/01_InitialChronologies_ITRDB.Rdata")


# PLOT INDIVIDUAL SERIES FROM EACH CHRONOLOGY-----------------------------------
for (i in 1:nrow(itrdb)){
  rwl <- itrdb %>% 
    slice(i) %>% 
    pluck("rwl", 1)

  title <- itrdb$site_code[i]

  png(filename = paste0("Figures/ITRDB_CheckHistorical/", title, ".png"))
  plot(rwl, main = title)
  dev.off()
}


# REMOVE HISTORICAL CHRONOLOGIES------------------------------------------------
## Method is to remove series where the most recent year is before 1900,
## then trim any leading/trailing years left with no data across all
## remaining series. See remove_historical() in Code/00_Functions.R.

itrdb_clean <- itrdb %>%
  mutate(
    rwl = map(rwl, remove_historical),
    start = map_int(rwl, ~ min(as.numeric(rownames(.x)))),
    end = map_int(rwl, ~ max(as.numeric(rownames(.x)))),
    length = end - start + 1,
    n_series = map_int(rwl, ncol)
  )

saveRDS(itrdb_clean, "Data/04_CleanedChronologies_ITRDB.Rdata")


# REPLOT INDIVIDUAL SERIES FROM EACH CHRONOLOGY---------------------------------
for (i in 1:nrow(itrdb_clean)){
  rwl <- itrdb_clean %>% 
    slice(i) %>% 
    pluck("rwl", 1)

  title <- itrdb_clean$site_code[i]

  png(filename = paste0("Figures/ITRDB_CheckHistorical_ForestOnly/", title, ".png"))
  plot(rwl, main = title)
  dev.off()
}
