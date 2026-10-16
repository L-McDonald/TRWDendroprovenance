library(tidyverse)
library(dplR)

source("Code/00_Functions.R")


# LOAD DATA---------------------------------------------------------------------
rwl <- readRDS("Outputs/06_CleanedChronologies_Step1.Rdata")


# DETREND CHRONOLOGIES----------------------------------------------------------
rwi <- rwl %>%
  mutate(
    rwi = map(rwl, ~ detrend(.x, method = "Spline", make.plot = FALSE))
  )

saveRDS(rwi, "Outputs/07_DetrendedChronologies_Step1.Rdata")


# CALCULATE CHRONOLOGIES--------------------------------------------------------
rwi <- readRDS("Outputs/07_DetrendedChronologies_Step1.Rdata")

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

saveRDS(crn, "Outputs/08_CalculatedChronologies_Step1.Rdata")
