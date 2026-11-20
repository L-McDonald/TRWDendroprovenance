library(tidyverse)
library(scico)
library(ggnewscale)
library(sf)

coast <- read_sf("Outputs/EnvironmentalVariables/Coastline_3035.gpkg")
region <- read_sf("Outputs/EnvironmentalVariables/Region_3035.gpkg")

bbox <- st_bbox(region)


pred <- readRDS("Outputs/LOOCV_Assignments/AVRA_res.Rdata") %>%
  filter(test_site == "VIZ") %>%
  rename(test_x = x, test_y = y) %>%
  unnest("assignment_df") %>%
  pivot_longer(
    cols = c("prob", "quant_68", "quant_95"),
    names_to = "raster",
    values_to = "value"
  ) %>%
  mutate(rast_type = ifelse(raster == "prob", "prob", "quant"))

test_x <- unique(pred$test_x)
test_y <- unique(pred$test_y)

ggplot() +
  geom_sf(data = coast, fill = "white") +
  geom_tile(
    data = filter(pred, raster == "prob" & model == "RF"),
    aes(x = x, y = y, fill = value)
  ) +
  scale_fill_scico(palette = "davos", direction = -1) +
  # scale_fill_gradient(low = "white", high = "steelblue") +
  new_scale_fill() +
  geom_tile(
    data = filter(pred, raster == "prob" & model == "UK1"),
    aes(x = x, y = y, fill = value)
  ) +
  scale_fill_scico(palette = "batlowW", direction = -1) +
  new_scale_fill() +
  geom_tile(
    data = filter(pred, raster == "prob" & model == "UK2"),
    aes(x = x, y = y, fill = value)
  ) +
  scale_fill_scico(palette = "oslo", direction = -1) +
  new_scale_fill() +
  geom_tile(
    data = filter(pred, raster == "quant_95"),
    aes(x = x, y = y, fill = as.logical(value))
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = NA,
      "TRUE" = "#E7EABF"
    ),
    na.value = NA
  ) +
  new_scale_fill() +
  geom_tile(
    data = filter(pred, raster == "quant_68"),
    aes(x = x, y = y, fill = as.logical(value))
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = NA,
      "TRUE" = "#77958D"
    ),
    na.value = NA
  ) +
  geom_sf(data = coast, fill = NA) +
  annotate(
    "point",
    x = test_x,
    y = test_y,
    shape = 4,
    colour = "firebrick",
    stroke = 1.5,
    size = 2
  ) +
  facet_grid(model ~ rast_type) +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  theme_bw()
