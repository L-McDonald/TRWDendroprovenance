library(tidyverse)
library(sf)
library(ranger)
library(automap)
library(furrr)


# LOAD DATA---------------------------------------------------------------------
sites <- readRDS("Outputs/18_TrainingTables.Rdata") %>%
  select(site_code, type, x_3035, y_3035)

pred_grid <- read.csv("Outputs/16_PredictionGrid_EnvironmentalVariables.csv")


# FUNCTION FOR LOOCV------------------------------------------------------------
run_loocv <- function(site, type, centre_x, centre_y, pred_grid) {
  path_in <- paste0("Outputs/TrainingTables/", site, "_", type, ".csv")
  path_out <- paste0("Outputs/LOOCV_Models/", site, "_", type, ".Rdata")

  # skip finished sites so the run can be resumed
  if (file.exists(path_out)) {
    return(invisible(path_out))
  }

  # generate site specific training table
  train_df <- read.csv(path_in)
  train_sf <- train_df %>%
    st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)

  # generate site specific prediction tables
  grid_df <- pred_grid %>%
    mutate(dist_km = sqrt((centre_x - x)^2 + (centre_y - y)^2) / 1000) %>%
    relocate(dist_km, .after = y)
  grid_sf <- grid_df %>%
    st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)

  # perform loocv for each site and each model
  loocv <- train_df %>%
    select(map_site, test_site, cor, x, y, dist_km, elev) %>%
    mutate(type = type, .after = test_site) %>% 
    filter(map_site != test_site) %>%
    crossing(model = c("RF", "UK1", "UK2")) %>%
    mutate(
      model_obj = map2(model, test_site, \(m, test) {
        if (m == "RF") {
          train0 <- train_df %>%
            filter(test_site != test) %>%
            select(-map_site, -test_site, -x, -y)

          ranger(cor ~ ., data = train0, quantreg = TRUE, num.threads = 1)
        } else if (m == "UK1") {
          train0 <- train_sf %>%
            filter(test_site != test) %>%
            select(cor, x, y, dist_km)

          autoKrige(
            cor ~ x + y + dist_km,
            input_data = train0,
            new_data = grid_sf
          )
        } else if (m == "UK2") {
          train0 <- train_sf %>% filter(test_site != test) %>% select(cor, x, y)

          autoKrige(
            cor ~ x + y + I(x^2) + I(y^2) + I(x * y),
            input_data = train0,
            new_data = grid_sf
          )
        }
      })
    )

  saveRDS(loocv, path_out)
  invisible(path_out)
}


# PERFMORM LOOCV----------------------------------------------------------------
plan(multisession, workers = 9)

future_pmap(
  list(sites$site_code, sites$type, sites$x_3035, sites$y_3035),
  \(s, t, cx, cy) run_loocv(s, t, cx, cy, pred_grid),
  .options = furrr_options(
    seed = TRUE,
    packages = c("tidyverse", "sf", "ranger", "automap")
  ),
  .progress = TRUE
)

plan(sequential)