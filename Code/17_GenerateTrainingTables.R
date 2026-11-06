library(tidyverse)


# LOAD DATA---------------------------------------------------------------------
crn <- readRDS("Outputs/12_CalculatedChronologies_Step2.Rdata")
site_env <- read.csv("Outputs/17_SiteEnvironmentalVariables.csv")
cor_all <- read.csv("Outputs/13_RawCorrelations_Step2.csv")


# TRAINING TABLES---------------------------------------------------------------
train <- crn %>%
  select(site_code, x_3035, y_3035, elev) %>%
  crossing(type = c("std", "res", "ars")) %>%
  mutate(
    train = map2(
      site_code,
      type,
      ~ cor_all %>%
        filter(site_1 == .x | site_2 == .x) %>%
        filter(type == .y) %>%
        mutate(
          map_site = .x,
          test_site = if_else(site_1 == .x, site_2, site_1),
          dist_km = dist_m / 1000
        ) %>%
        select(map_site, test_site, cor, dist_km) %>%
        left_join(site_env, by = join_by(test_site == site_code)) %>% 
        relocate(dist_km, .after = y)
    )
  )


# SAVE TRAINING TABLES----------------------------------------------------------
saveRDS(train, "Outputs/18_TrainingTables.Rdata")

train %>%
  select(site_code, type, train) %>%
  pwalk(\(site_code, type, train) {
    write.csv(
      train,
      paste0("Outputs/TrainingTables/", site_code, "_", type, ".csv"),
      row.names = FALSE
    )
  })
