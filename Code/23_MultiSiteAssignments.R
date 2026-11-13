library(tidyverse)
library(terra)
library(assignR)
library(furrr)


# LOAD DATA---------------------------------------------------------------------
sites <- readRDS("Outputs/18_TrainingTables.Rdata") %>%
  select(site_code, type, x_3035, y_3035)

pairwise_cor <- read.csv("Outputs/13_RawCorrelations_Step2.csv")


# FUNCTION FOR MULTI-CHRONOLOGY ASSIGNMENTS-------------------------------------
assign_multi_site <- function(test_site, test_x, test_y) {
  output <- vector("list", 9)

  chron_types <- c("ars", "res", "std")
  model_types <- c("RF", "UK1", "UK2")

  count <- 0

  location <- vect(cbind(test_x, test_y), crs = "EPSG:3035")

  for (i in 1:length(chron_types)) {
    for (j in 1:length(model_types)) {
      count <- count + 1

      assign <- pairwise_cor %>%
        filter(site_1 == test_site | site_2 == test_site) %>%
        filter(site_1 != site_2) %>%
        filter(type == chron_types[i]) %>%
        mutate(
          test_site = test_site,
          map_site = if_else(site_1 == test_site, site_2, site_1)
        ) %>%
        select(test_site, map_site, type, dist_m, n, cor) %>%
        arrange(desc(cor)) %>%
        mutate(
          n_sites = row_number(),
          model = model_types[j],
          rast_path = paste0(
            "Outputs/LOOCV_Rasters/",
            map_site,
            "/",
            map_site,
            "_",
            test_site,
            "excl_",
            chron_types[i],
            "_",
            model,
            ".tif"
          ),
          assignment_df = list(NULL),
          accuracy_68 = NA,
          accuracy_95 = NA,
          area_68_count = NA,
          area_95_count = NA,
          precision_68 = NA,
          precision_95 = NA
        )

      for (k in 1:nrow(assign)) {
        rast_list <- lapply(assign$rast_path[1:k], terra::rast) %>%
          setNames(assign$map_site[1:k])

        rast <- isoStack(rast_list)

        test_cor <- assign %>%
          slice(1:k) %>%
          select(map_site, cor) %>%
          pivot_wider(names_from = map_site, values_from = cor) %>%
          mutate(ID = test_site, .before = everything()) %>%
          as.data.frame()

        prob_rast <- pdRaster(rast, test_cor, genplot = FALSE)

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

        assign[k, "accuracy_68"] <- terra::extract(quant_68_rast, location)[[2]]
        assign[k, "accuracy_95"] <- terra::extract(quant_95_rast, location)[[2]]

        assign[k, "area_68_count"] <- sum(values(quant_68_rast), na.rm = TRUE)
        assign[k, "area_95_count"] <- sum(values(quant_95_rast), na.rm = TRUE)

        assign[k, "precision_68"] <- (1 -
          mean(values(quant_68_rast), na.rm = TRUE)) *
          100
        assign[k, "precision_95"] <- (1 -
          mean(values(quant_95_rast), na.rm = TRUE)) *
          100

        df <- terra::as.data.frame(prob_rast, xy = TRUE) %>%
          bind_cols(
            terra::as.data.frame(quant_68_rast),
            terra::as.data.frame(quant_95_rast),
            .name_repair = "minimal"
          ) %>%
          setNames(c("x", "y", "prob", "quant_68", "quant_95"))

        assign$assignment_df[[k]] <- df
      }

      output[[count]] <- assign
    }
  }

  output_long <- bind_rows(output)

  saveRDS(
    output_long,
    paste0("Outputs/LOOCV_MultiSiteAssignments/", test_site, ".Rdata")
  )
}


# PERFMORM ASSIGNMENTS----------------------------------------------------------
plan(multisession, workers = 9)

future_pmap(
  list(sites$site_code, sites$x_3035, sites$y_3035),
  \(s, x, y) assign_multi_site(s, x, y),
  .options = furrr_options(
    seed = TRUE,
    packages = c("tidyverse", "sf", "terra")
  ),
  .progress = TRUE
)

plan(sequential)