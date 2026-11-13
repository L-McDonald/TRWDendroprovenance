library(tidyverse)
library(terra)
library(sf)
library(ranger)
library(furrr)


# LOAD DATA---------------------------------------------------------------------
sites <- readRDS("Outputs/18_TrainingTables.Rdata") %>%
  select(site_code, type, x_3035, y_3035)

pred_grid <- read.csv("Outputs/16_PredictionGrid_EnvironmentalVariables.csv")


# FUNCTION TO EXTRACT PREDICTIONS AND SAVE RASTERS------------------------------
extract_predictions <- function(site, type, centre_x, centre_y, pred_grid) {
  path_in <- paste0("Outputs/LOOCV_Models/", site, "_", type, ".Rdata")
  path_out <- paste0("Outputs/LOOCV_Predictions/", site, "_", type, ".Rdata")

  # generate site specific prediction tables
  grid_df <- pred_grid %>%
    mutate(dist_km = sqrt((centre_x - x)^2 + (centre_y - y)^2) / 1000) %>%
    relocate(dist_km, .after = y)
  grid_sf <- grid_df %>%
    st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)

  # load models
  models <- readRDS(path_in) %>%
    mutate(
      predictions = map2(model, model_obj, \(m, obj) {
        if (grepl("UK", m)) {
          obj %>%
            pluck("krige_output") %>%
            select(pred = var1.pred, se = var1.stdev) %>%
            st_drop_geometry() %>%
            mutate(x = grid_df$x, y = grid_df$y, .before = pred)
        } else if (m == "RF") {
          pred <- predict(obj, grid_df)$predictions

          se <- predict(
            obj,
            grid_df,
            type = "quantiles",
            quantiles = c(0.159, 0.841)
          )$predictions

          grid_df %>%
            select(x, y) %>%
            mutate(pred = pred) %>%
            bind_cols(as.data.frame(se) %>% setNames(c("lwr", "upr"))) %>%
            mutate(se = (upr - lwr) / 2) %>%
            select(-lwr, -upr)
        }
      })
    ) %>%
    select(-model_obj) %>%
    mutate(
      rast = map(
        predictions,
        ~ terra::rast(.x, type = "xyz", crs = "EPSG:3035")
      )
    ) %T>%
    pwalk(\(rast, map_site, test_site, type, model, ...) {
      dir.create(
        file.path("Outputs/LOOCV_Rasters", map_site),
        recursive = TRUE,
        showWarnings = FALSE
      )

      terra::writeRaster(
        rast,
        paste0(
          "Outputs/LOOCV_Rasters/",
          map_site,
          "/",
          map_site,
          "_",
          test_site,
          "excl",
          "_",
          type,
          "_",
          model,
          ".tif"
        )
      )
    }) %>%
    select(-rast)

  saveRDS(models, path_out)
}

# PERFMORM PREDICTIONS----------------------------------------------------------
plan(multisession, workers = 9)

future_pmap(
  list(sites$site_code, sites$type, sites$x_3035, sites$y_3035),
  \(s, t, cx, cy) extract_predictions(s, t, cx, cy, pred_grid),
  .options = furrr_options(
    seed = TRUE,
    packages = c("tidyverse", "sf", "ranger", "terra")
  ),
  .progress = TRUE
)

plan(sequential)