library(tidyverse)
library(dplR)

source("Code/00_Functions.R")


# LOAD DATA---------------------------------------------------------------------
rwl <- readRDS("Data/06_CleanedChronologies.Rdata")


# DETREND CHRONOLOGIES----------------------------------------------------------
rwi <- rwl %>%
  mutate(
    rwi = map(rwl, ~ detrend(.x, method = "Spline", make.plot = FALSE))
  )

saveRDS(rwi, "Data/07_DetrendedChronologies.Rdata")


# CALCULATE CHRONOLOGIES--------------------------------------------------------
rwi <- readRDS("Data/07_DetrendedChronologies.Rdata")

crn <- rwi %>% 
  mutate(
    chronologies = map(
      rwi,
      ~ chron.ars.mod(.x, verbose = FALSE, prewhitenMethod = "ar.yw")
    ),
    crn = map(chronologies, "chronologies"),
    crn = map(crn, ~ {
      .x[] <- lapply(.x, function(col) replace(col, is.nan(col), NA))
      .x
    }),
    forcedAR = map_lgl(chronologies, "forcedAR")
  ) %>% 
  select(-chronologies)

saveRDS(crn, "Data/08_CalculatedChronologies.Rdata")
