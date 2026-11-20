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

  map_site <- site

  # check for and skip finished sites so the run can be resumed
  if (file.exists(path_out)) {
    return(invisible(path_out))
  }

  assign <- readRDS(path_in) %>%
    select(map_site, test_site, x, y, dist_km, elev, model) %>%
    mutate(
      location = map2(x, y, ~ terra::vect(cbind(.x, .y), crs = "EPSG:3035")),
      results = list(NULL)
    )

  for (i in 1:nrow(assign)) {
    # set parameters
    test_site <- assign$test_site[i]
    model <- assign$model[i]

    # check rasters exist, load and name layers
    rast_paths <- sprintf(
      "Outputs/LOOCV_Rasters/%s/%s_%sexcl_%s_%s.tif",
      map_site,
      map_site,
      test_site,
      c("ars", "res", "std"),
      model
    )

    files_exist <- file.exists(rast_paths)

    if (sum(files_exist) < 2) {
      next
    }

    rast_paths <- rast_paths[files_exist]

    rast <- terra::rast(rast_paths)
    chrons <- str_extract(sources(rast), "ars|res|std")
    names(rast) <- paste(
      c("pred", "se"),
      rep(chrons, each = 2),
      sep = "_"
    )

    # select relevant pairwise correlations
    site_cor <- pairwise_cor %>%
      filter(type %in% chrons) %>%
      filter(
        (site_1 == map_site & site_2 == test_site) |
          (site_1 == test_site & site_2 == map_site)
      )

    # generate output tibble
    output <- tibble(
      type = unlist(lapply(2:length(chrons), function(k) {
        combn(chrons, k, paste, collapse = "_")
      }))
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
            split(., sub(".*_", "", rep(.x, each = 2))) %>%
            isoStack()
        ),
        prob = map2(
          raster,
          correlation,
          ~ suppressWarnings(pdRaster(.x, .y, genplot = FALSE))
        ),
        prob_fail = map_lgl(prob, ~ all(is.na(values(.x)))),
        quant_68 = map2(
          prob,
          prob_fail,
          \(prob, run) {
            if (!run) {
              qtlRaster(prob, 0.68, thresholdType = "prob", genplot = FALSE)
            } else {
              NULL
            }
          }
        ),
        quant_95 = map2(
          prob,
          prob_fail,
          \(prob, run) {
            if (!run) {
              qtlRaster(prob, 0.95, thresholdType = "prob", genplot = FALSE)
            } else {
              NULL
            }
          }
        ),
        assignment_df = pmap(
          list(prob, quant_68, quant_95, prob_fail),
          \(p, q68, q95, run) {
            if (!run) {
              as.data.frame(p, xy = TRUE) %>%
                bind_cols(
                  as.data.frame(q68),
                  as.data.frame(q95),
                  .name_repair = "minimal"
                )
            } else {
              NULL
            }
          }
        ),
        accuracy_68 = pmap_lgl(
          list(quant_68, location, prob_fail),
          \(quant, loc, run) {
            if (!run) {
              terra::extract(quant, loc)[[2]]
            } else {
              NA
            }
          }
        ),
        accuracy_95 = pmap_lgl(
          list(quant_95, location, prob_fail),
          \(quant, loc, run) {
            if (!run) {
              terra::extract(quant, loc)[[2]]
            } else {
              NA
            }
          }
        ),
        area_68_count = map2_int(quant_68, prob_fail, \(quant, run) {
          if (!run) {
            sum(values(quant), na.rm = TRUE)
          } else {
            NA_integer_
          }
        }),
        area_95_count = map2_int(quant_95, prob_fail, \(quant, run) {
          if (!run) {
            sum(values(quant), na.rm = TRUE)
          } else {
            NA_integer_
          }
        }),
        precision_68 = map2_dbl(
          quant_68,
          prob_fail,
          \(quant, run) {
            if (!run) {
              (1 - mean(values(quant), na.rm = TRUE)) * 100
            } else {
              NA
            }
          }
        ),
        precision_95 = map2_dbl(
          quant_95,
          prob_fail,
          \(quant, run) {
            if (!run) {
              (1 - mean(values(quant), na.rm = TRUE)) * 100
            } else {
              NA
            }
          }
        )
      ) %>%
      select(
        type,
        correlation,
        prob_fail,
        assignment_df,
        starts_with("accuracy"),
        starts_with("area"),
        starts_with("precision")
      )

    assign$results[[i]] <- output
  }

  assign_out <- assign %>% unnest("results")

  saveRDS(assign_out, path_out)
  invisible(path_out)
}


# PERFMORM ASSIGNMENTS----------------------------------------------------------
plan(multisession, workers = 10)

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
