library(tidyverse)
library(terra)
library(assignR)
library(furrr)


# LOAD DATA---------------------------------------------------------------------
sites <- readRDS("Outputs/18_TrainingTables.Rdata") %>%
  select(site_code, type, x_3035, y_3035)

make_assignments <- function(site, type) {
  path_in <- paste0("Outputs/LOOCV_Predictions/", site, "_", type, ".Rdata")
  path_out <- paste0("Outputs/LOOCV_Assignments/", site, "_", type, ".Rdata")

  assign <- readRDS(path_in) %>%
    select(-predictions) %>%
    mutate(
      rast_path = paste0(
        "Outputs/LOOCV_Rasters/",
        map_site,
        "_",
        test_site,
        "excl_",
        type,
        "_",
        model,
        ".tif"
      ),
      location = map2(x, y, ~ terra::vect(cbind(.x, .y), crs = "EPSG:3035")),
      assignment_df = vector("list", nrow(assign)),
      accuracy_68 = NA,
      accuracy_95 = NA,
      area_68_count = NA,
      area_95_count = NA,
      precision_68 = NA,
      precision_95 = NA
    )

  for (i in 1:nrow(assign)) {
    rast_in <- assign$rast_path[i]
    rast <- terra::rast(rast_in)

    location <- assign$location[[i]]

    cor <- assign %>%
      slice(i) %>%
      select(test_site, cor) %>%
      as.data.frame()

    prob_rast <- pdRaster(rast, cor, genplot = FALSE)

    quant_68_rast <- qtlRaster(
      prob_rast,
      0.68,
      thresholdType = "prob",
      genplot = FALSE
    )
    quant_95_rast <- qtlRaster(
      prob_rast,
      0.95,
      thresholdType = "prob",
      genplot = FALSE
    )

    assign[i, "accuracy_68"] <- terra::extract(quant_68_rast, location)[[2]]
    assign[i, "accuracy_95"] <- terra::extract(quant_95_rast, location)[[2]]

    assign[i, "area_68_count"] <- sum(values(quant_68_rast), na.rm = TRUE)
    assign[i, "area_95_count"] <- sum(values(quant_95_rast), na.rm = TRUE)

    assign[i, "precision_68"] <- (1 -
      mean(values(quant_68_rast), na.rm = TRUE)) *
      100
    assign[i, "precision_95"] <- (1 -
      mean(values(quant_95_rast), na.rm = TRUE)) *
      100

    df <- terra::as.data.frame(prob_rast, xy = TRUE) %>%
      bind_cols(terra::as.data.frame(quant_68_rast)) %>%
      bind_cols(terra::as.data.frame(quant_95_rast)) %>%
      setNames(c("x", "y", "prob", "quant_68", "quant_95"))

    assign$assignment_df[[i]] <- df
  }
  saveRDS(assign, path_out)
}


# PERFMORM ASSIGNMENTS----------------------------------------------------------
plan(multisession, workers = 9)

future_pmap(
  list(sites$site_code, sites$type),
  \(s, t) make_assignments(s, t),
  .options = furrr_options(
    seed = TRUE,
    packages = c("tidyverse", "sf", "terra")
  ),
  .progress = TRUE
)

plan(sequential)