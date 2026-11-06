library(tidyverse)
library(automap)
library(ranger)
library(RandomForestsGLS)
library(sf)
library(terra)


# LOAD DATA---------------------------------------------------------------------
train_init <- readRDS("Outputs/18_TrainingTables.Rdata")
pred_df <- read.csv("Outputs/16_PredictionGrid_EnvironmentalVariables.csv")

pred_sf <- pred_df %>%
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)

coast <- read_sf("Outputs/EnvironmentalVariables/Coastline_3035.gpkg")
region <- read_sf("Outputs/EnvironmentalVariables/Region_3035.gpkg")


# SET UP TRAINING AND PREDICTION DATA-------------------------------------------
train <- train_init %>%
  crossing(model = c("UK_1", "UK_2", "UK_3", "UK_4", "UK_5", "UK_6", "RF")) %>%
  relocate(model, .before = train) %>%
  mutate(
    train = map2(train, model, \(df, m) {
      if (grepl("UK", m)) {
        df %>%
          select(cor, x, y, dist_km, elev) %>%
          st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)
      } else {
        df %>% select(-map_site, -test_site) %>% relocate(dist_km, .after = y)
      }
    }),
    pred_grid = map(
      model,
      ~ if (grepl("UK", .x)) {
        pred_sf
      } else {
        pred_df
      }
    ),
    pred_grid = pmap(
      list(x_3035, y_3035, pred_grid),
      \(x0, y0, grid) {
        grid %>%
          mutate(dist_km = sqrt((x0 - x)^2 + (y0 - y)^2) / 1000) %>%
          relocate(dist_km, .after = y)
      }
    )
  )

saveRDS(train, "tmp/train.RDS")


# RUN MODELS--------------------------------------------------------------------
models <- train %>%
  mutate(
    model_obj = pmap(list(train, model, pred_grid), \(df, m, grid) {
      if (m == "UK_1") {
        autoKrige(
          cor ~ x + y + I(x^2) + I(y^2) + I(x * y),
          input_data = df,
          new_data = grid
        )
      } else if (m == "UK_2") {
        autoKrige(
          cor ~ x + y + I(x^2) + I(y^2) + I(x * y) + elev,
          input_data = df,
          new_data = grid
        )
      } else if (m == "UK_3") {
        autoKrige(
          cor ~ x + y + I(x^2) + I(y^2) + I(x * y) + dist_km,
          input_data = df,
          new_data = grid
        )
      } else if (m == "UK_4") {
        autoKrige(
          cor ~ x + y + I(x^2) + I(y^2) + I(x * y) + elev + dist_km,
          input_data = df,
          new_data = grid
        )
      } else if (m == "UK_5") {
        autoKrige(
          cor ~ x + y + dist_km,
          input_data = df,
          new_data = grid
        )
      } else if (m == "UK_6") {
        autoKrige(
          cor ~ x + y + dist_km + elev,
          input_data = df,
          new_data = grid
        )
      } else if (m == "RF") {
        ranger(cor ~ ., data = df, quantreg = TRUE, num.trees = 500)
      }
    })
  )

saveRDS(models, "tmp/models.RDS")

# PREDICTIONS AND RESIDUALS-----------------------------------------------------
pred <- models %>%
  mutate(
    predictions = pmap(list(model_obj, model, pred_grid), \(obj, m, grid) {
      if (grepl("UK", m)) {
        obj[["krige_output"]] %>%
          mutate(x = st_coordinates(.)[, "X"], y = st_coordinates(.)[, "Y"]) %>%
          select(x, y, pred = var1.pred, var = var1.var, stdev = var1.stdev) %>%
          st_drop_geometry()
      } else if (m == "RF") {
        pred_df %>%
          select(x, y) %>%
          bind_cols(
            predict(obj, grid) %>%
              pluck("predictions") %>%
              as.data.frame() %>%
              setNames("pred")
          ) %>%
          bind_cols(
            predict(
              obj,
              grid,
              type = "quantiles",
              quantiles = c(0.159, 0.841)
            ) %>%
              pluck("predictions") %>%
              as.data.frame() %>%
              setNames(c("quant_lwr", "quant_upr"))
          )
      }
    }),
    residuals = pmap(
      list(model_obj, train, model, predictions),
      \(obj, train, m, pred) {
        if (grepl("UK", m)) {
          raster <- terra::rast(pred[, 1:3], type = "xyz", crs = "EPSG:3035")

          train %>%
            select(x, y, cor) %>%
            st_drop_geometry() %>%
            bind_cols(
              train %>%
                select(x, y) %>%
                st_as_sf(coords = c("x", "y"), crs = 3035) %>%
                terra::extract(
                  raster,
                  .,
                  ID = FALSE,
                  search_radius = 100000
                ) %>%
                select(pred_ib = pred) %>%
                st_drop_geometry()
            )
        } else if (m == "RF") {
          raster <- terra::rast(pred[, 1:3], type = "xyz", crs = "EPSG:3035")

          train %>%
            select(x, y, cor) %>%
            bind_cols(
              obj %>%
                pluck("predictions") %>%
                as.data.frame() %>%
                setNames("pred_oob")
            ) %>%
            bind_cols(
              train %>%
                select(x, y) %>%
                st_as_sf(coords = c("x", "y"), crs = 3035) %>%
                terra::extract(
                  raster,
                  .,
                  ID = FALSE,
                  search_radius = 100000
                ) %>%
                select(pred_ib = pred) %>%
                st_drop_geometry()
            )
        }
      }
    )
  )

saveRDS(pred, "tmp/pred.RDS")

# CALCULATE R2 AND RMSE THEN SAVE RESULTS---------------------------------------
results <- pred %>%
  mutate(
    rmse = map_dbl(residuals, ~ sqrt(sum((.x$pred_ib - .x$cor)^2) / nrow(.x))),
    rsq = map_dbl(
      residuals,
      ~ 1 - sum((.x$cor - .x$pred_ib)^2) / sum((.x$cor - mean(.x$cor))^2)
    ),
    rmse_oob = map2_dbl(
      model,
      residuals,
      ~ if (.x == "RF") {
        sqrt(sum((.y$pred_oob - .y$cor)^2) / nrow(.y))
      } else {
        NA_real_
      }
    ),
    rsq_oob = map2_dbl(
      model,
      residuals,
      ~ if (.x == "RF") {
        1 - sum((.y$cor - .y$pred_oob)^2) / sum((.y$cor - mean(.y$cor))^2)
      } else {
        NA_real_
      }
    )
  )

calc <- results %>%
  group_by(model, type) %>%
  summarise(
    mean_rmse = mean(rmse),
    mean_rsq = mean(rsq),
    mean_rmse_oob = mean(rmse_oob),
    mean_rsq_oob = mean(rsq_oob)
  )

saveRDS(results, "Outputs/19_ModelComparisons.Rdata")
write.csv(calc, "Outputs/20_ModelComparisonCalculations.csv", row.names = FALSE)

# PLOT MAPS---------------------------------------------------------------------
exp_fct <- 100000

plot_site <- function(site, results, exp_fct = 100000) {
  site_res <- results %>% filter(site_code == site)

  test <- site_res %>%
    select(site_code, x_3035, y_3035, type, model, predictions) %>%
    unnest(predictions)

  labels <- site_res %>%
    select(site_code, type, model, rmse, rsq, rmse_oob, rsq_oob) %>%
    mutate(across(where(is.numeric), \(x) round(x, 3)))

  point <- site_res %>%
    select(site_code, x_3035, y_3035) %>%
    distinct()

  plt <- ggplot() +
    geom_sf(data = coast, fill = NA) +
    geom_tile(data = test, aes(x = x, y = y, fill = pred)) +
    annotate(
      "point",
      x = point$x_3035,
      y = point$y_3035,
      shape = 4,
      size = 1.5,
      stroke = 1,
      colour = "firebrick"
    ) +
    geom_label(
      data = labels,
      aes(
        x = I(0.025),
        y = I(0.975),
        label = paste0("R2 = ", rsq, "\nRMSE = ", rmse)
      ),
      hjust = 0,
      vjust = 1,
      family = "serif",
      size = 6,
      size.unit = "pt",
      lineheight = 1,
      fill = "white",
      border.colour = "grey40"
    ) +
    geom_label(
      data = labels %>% filter(grepl("RF", model)),
      aes(
        x = I(0.975),
        y = I(0.025),
        label = paste0("OOB\nR2 = ", rsq_oob, "\nRMSE = ", rmse_oob)
      ),
      hjust = 1,
      vjust = 0,
      family = "serif",
      size = 6,
      size.unit = "pt",
      lineheight = 1,
      fill = "white",
      border.colour = "grey40"
    ) +
    facet_grid(type ~ model, switch = "y") +
    scico::scale_fill_scico(
      palette = "batlow",
      name = "Correlation (r)",
      guide = guide_colourbar(frame.colour = "black", frame.linewidth = 0.25)
    ) +
    coord_sf(
      xlim = c(st_bbox(region)[1] - exp_fct, st_bbox(region)[3] + exp_fct),
      ylim = c(st_bbox(region)[2] - exp_fct, st_bbox(region)[4] + exp_fct),
      expand = FALSE
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "serif"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(family = "serif", face = "bold", size = 12),
      legend.key.height = unit(1, "null"),
      legend.ticks = element_line(colour = "black", linewidth = 0.25),
      legend.title.position = "right",
      legend.title = element_text(angle = 270, hjust = 0.5),
      legend.text = element_text(size = 10)
    )

  ggsave(
    paste0("Figures/ModelComparisons/", site, ".png"),
    plot = plt,
    width = 56,
    height = 15,
    units = "cm",
    dpi = 300
  )
}

purrr::walk(unique(results$site_code), \(s) plot_site(s, results, exp_fct))


# PLOT STATISTICS---------------------------------------------------------------
results %>%
  mutate(plot_x = paste(model, type, sep = "_")) %>%
  ggplot() +
  geom_jitter(aes(x = plot_x, y = rmse, colour = type)) +
  scico::scale_colour_scico_d(palette = "bamako") +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.justification.inside = c(0, 1),
    legend.position.inside = c(0, 1),
    legend.background = element_blank(),
    legend.box = element_blank()
  )

results %>%
  mutate(plot_x = paste(model, type, sep = "_")) %>%
  ggplot() +
  geom_jitter(aes(x = plot_x, y = rsq, colour = type)) +
  scico::scale_colour_scico_d(palette = "bamako") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.justification.inside = c(0, 1),
    legend.position.inside = c(0, 1),
    legend.background = element_blank(),
    legend.box = element_blank()
  )

results %>%
  select(type, model, rmse_oob, rsq_oob) %>%
  filter(model == "RF") %>%
  pivot_longer(cols = c("rmse_oob", "rsq_oob")) %>%
  ggplot() +
  geom_jitter(aes(x = type, y = value, colour = type)) +
  facet_wrap(~name, scales = "free_y") +
  scico::scale_colour_scico_d(palette = "bamako") +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.justification.inside = c(0, 1),
    legend.position.inside = c(0, 1),
    legend.background = element_blank(),
    legend.box = element_blank()
  )


results %>%
  select(type, model, rmse, rsq) %>%
  filter(str_detect(model, "UK")) %>%
  ggplot() +
  geom_jitter(aes(x = model, y = rmse, colour = model)) +
  facet_wrap(~type) +
  scico::scale_colour_scico_d(palette = "bamako") +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.justification.inside = c(0, 1),
    legend.position.inside = c(0, 1),
    legend.background = element_blank(),
    legend.box = element_blank()
  )

results %>%
  select(type, model, rmse, rsq) %>%
  filter(str_detect(model, "UK")) %>%
  ggplot() +
  geom_jitter(aes(x = model, y = rsq, colour = model)) +
  facet_wrap(~type) +
  scico::scale_colour_scico_d(palette = "bamako") +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.justification.inside = c(0, 0),
    legend.position.inside = c(0, 0),
    legend.background = element_blank(),
    legend.box = element_blank()
  )

calc %>%
  select(model, type, mean_rsq) %>%
  pivot_wider(names_from = type, values_from = mean_rsq)


# PROPORTION IMPROVEMENT--------------------------------------------------------
## for the UK models, examine how often the addition of dist_km and/or elevation
## leads to improvement over the base model

results %>%
  select(site_code, type, model, rsq) %>%
  filter(str_detect(model, "UK")) %>%
  pivot_wider(names_from = model, values_from = rsq) %>%
  mutate(
    comp_1_2 = UK_2 > UK_1,
    comp_1_3 = UK_3 > UK_1,
    comp_1_4 = UK_4 > UK_1,
    comp_1_5 = UK_5 > UK_1,
    comp_5_6 = UK_6 > UK_5
  ) %>% 
  group_by(type) %>% 
  summarise(across(starts_with("comp"), ~mean(.x, na.rm = TRUE)))


results %>%
  select(site_code, type, model, rmse) %>%
  filter(str_detect(model, "UK")) %>%
  pivot_wider(names_from = model, values_from = rmse) %>%
  mutate(
    comp_1_2 = UK_2 < UK_1,
    comp_1_3 = UK_3 < UK_1,
    comp_1_4 = UK_4 < UK_1,
    comp_1_5 = UK_5 < UK_1,
    comp_5_6 = UK_6 < UK_5
  ) %>% 
  group_by(type) %>% 
  summarise(across(starts_with("comp"), ~mean(.x, na.rm = TRUE)))
