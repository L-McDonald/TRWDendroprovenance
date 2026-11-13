library(tidyverse)
library(terra)
library(assignR)
library(furrr)


# LOAD DATA---------------------------------------------------------------------
sites <- readRDS("Outputs/18_TrainingTables.Rdata") %>%
  select(site_code, type, x_3035, y_3035)

pairwise_cor <- read.csv("Outputs/13_RawCorrelations_Step2.csv")


# FUNCTION FOR MULTI-CHRONOLOGY ASSIGNMENTS-------------------------------------
assign_multi_chron <- function(site) {
  path_in <- paste0("Outputs/LOOCV_Predictions/", site, "_ars.Rdata")
  path_out <- paste0("Outputs/LOOCV_MultiChronAssignments/", site, ".Rdata")

  assign <- readRDS(path_in) %>%
    select(map_site, test_site, x, y, dist_km, elev, model) %>%
    mutate(
      location = map2(x, y, ~ terra::vect(cbind(.x, .y), crs = "EPSG:3035")),
      results = list(NULL)
    )

  for (i in 1:nrow(assign)) {
    # set parameters
    map_site <- assign$map_site[i]
    test_site <- assign$test_site[i]
    model <- assign$model[i]

    # load rasters and name layers
    rast_paths <- sprintf(
      "Outputs/LOOCV_Rasters/%s/%s_%sexcl_%s_%s.tif",
      map_site,
      map_site,
      test_site,
      c("ars", "res", "std"),
      model
    )
    rast <- terra::rast(rast_paths)
    names(rast) <- paste(
      c("pred", "se"),
      rep(c("ars", "res", "std"), each = 2),
      sep = "_"
    )

    # select relevant pairwise correlations
    site_cor <- pairwise_cor %>%
      filter(
        (site_1 == map_site & site_2 == test_site) |
          (site_1 == test_site & site_2 == map_site)
      )

    # generate output tibble
    output <- tibble(
      type = c("res_std", "ars_res", "ars_std", "ars_res_std")
    ) %>%
      mutate(
        type_vect = strsplit(type, "_"),
        location = assign$location[i],
        correlation = map(
          type_vect,
          ~ site_cor %>%
            filter(type %in% .x) %>%
            select(type, cor) %>%
            pivot_wider(names_from = "type", values_from = "cor") %>%
            select(all_of(.x)) %>%
            mutate(ID = test_site, .before = everything()) %>%
            as.data.frame()
        ),
        raster = map(
          type_vect,
          ~ rast[[grepl(paste(.x, collapse = "|"), names(rast))]] %>%
            split(., sub(".*_", "", names(.))) %>%
            isoStack()
        ),
        prob = map2(
          raster,
          correlation,
          ~ suppressWarnings(pdRaster(.x, .y, genplot = FALSE))
        ),
        quant_68 = map(
          prob,
          ~ qtlRaster(.x, 0.68, thresholdType = "prob", genplot = FALSE)
        ),
        quant_95 = map(
          prob,
          ~ qtlRaster(.x, 0.95, thresholdType = "prob", genplot = FALSE)
        ),
        assignment_df = pmap(list(prob, quant_68, quant_95), \(p, q68, q95) {
          as.data.frame(p, xy = TRUE) %>%
            bind_cols(
              as.data.frame(q68),
              as.data.frame(q95),
              .name_repair = "minimal"
            )
        }),
        accuracy_68 = map2_lgl(
          quant_68,
          location,
          ~ terra::extract(.x, .y)[[2]]
        ),
        accuracy_95 = map2_lgl(
          quant_95,
          location,
          ~ terra::extract(.x, .y)[[2]]
        ),
        area_68_count = map_int(quant_68, ~ sum(values(.x), na.rm = TRUE)),
        area_95_count = map_int(quant_95, ~ sum(values(.x), na.rm = TRUE)),
        precision_68 = map_dbl(
          quant_68,
          ~ (1 - mean(values(.x), na.rm = TRUE)) * 100
        ),
        precision_95 = map_dbl(
          quant_95,
          ~ (1 - mean(values(.x), na.rm = TRUE)) * 100
        )
      ) %>%
      select(
        type,
        correlation,
        assignment_df,
        starts_with("accuracy"),
        starts_with("area"),
        starts_with("precision")
      )

    assign$results[[i]] <- output
  }

  assign_out <- assign %>% unnest("results")
  
  saveRDS()
}


# PERFMORM ASSIGNMENTS----------------------------------------------------------
plan(multisession, workers = 9)

future_pmap(
  list(sites$site_code),
  \(s) assign_multi_chron(s),
  .options = furrr_options(
    seed = TRUE,
    packages = c("tidyverse", "sf", "terra")
  ),
  .progress = TRUE
)

plan(sequential)