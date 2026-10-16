library(tidyverse)
library(dplR)

source("Code/00_Functions.R")


# LOAD DATA---------------------------------------------------------------------
rwl <- readRDS("Outputs/10_CleanedChronologies_Step2.Rdata")


# DETREND CHRONOLOGIES----------------------------------------------------------
rwi <- rwl %>%
  mutate(
    rwi = map(rwl, ~ detrend(.x, method = "Spline", make.plot = FALSE))
  )

saveRDS(rwi, "Outputs/11_DetrendedChronologies_Step2.Rdata")


# CALCULATE CHRONOLOGIES--------------------------------------------------------
rwi <- readRDS("Outputs/11_DetrendedChronologies_Step2.Rdata")

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

saveRDS(crn, "Outputs/12_CalculatedChronologies_Step2.Rdata")
