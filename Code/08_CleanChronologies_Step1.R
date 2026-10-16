library(rnaturalearth)
library(tidyverse)
library(sf)
library(dplR)


source("Code/00_Functions.R")


# LOAD DATA---------------------------------------------------------------------
rwl <- readRDS("Outputs/05_InitialChronologies_All.Rdata")

coast <- ne_countries(scale = 10, continent = c("Africa", "Europe", "Asia")) %>%
    st_make_valid() %>%
    st_crop(c(xmin = -25, xmax = 90, ymin = 10, ymax = 90)) %>%
    st_geometry() %>%
    st_transform(3035)


# CHECK FOR AND REMOVE CHRONOLOGIES LOCATED AT SEA------------------------------
rwl_sea <- rwl %>%
  st_as_sf(coords = c("x_3035", "y_3035"), crs = 3035, remove = FALSE) %>%
  filter(apply(st_intersects(., coast, sparse = FALSE), 1, any))


# COMBINE CHRONOLOGIES WITHIN 1 KM OF EACH OTHER (ANY SPECIES)------------------

## Identify sites within 1 km of each other
dist_mat <- st_distance(rwl_sea)
diag(dist_mat) <- NA

close_idx <- which(as.matrix(dist_mat) < units::set_units(1000, "m"), arr.ind = TRUE)
close_idx <- close_idx[close_idx[, 1] < close_idx[, 2], , drop = FALSE]

close_pairs <- tibble(
  i = close_idx[, 1],
  j = close_idx[, 2],
  site_a = rwl_sea$site_code[i],
  site_b = rwl_sea$site_code[j],
  species_a = rwl_sea$species[i],
  species_b = rwl_sea$species[j],
  same_species = species_a == species_b,
  dist_m = as.numeric(dist_mat[close_idx])
)

## Use tidygraph to identify clusters of chronologies that are close together,
## regardless of species.
clusters <- close_pairs %>% 
  select(site_a, site_b) %>%
  tidygraph::as_tbl_graph(directed = FALSE) %>%
  tidygraph::activate(nodes) %>%
  tidygraph::mutate(cluster = tidygraph::group_components()) %>%
  as_tibble() %>%
  rename(site_code = name) %>%
  group_by(cluster) %>%
  summarise(members = list(site_code)) %>% 
  ungroup()

## merge each cluster into a single chronology
merged_rows <- clusters %>%
  mutate(
    members_df = map(members, ~ rwl_sea %>% filter(site_code %in% .x)),
    merged = map(members_df, merge_cluster)
  ) %>%
  pull(merged) %>%
  bind_rows() %>%
  st_as_sf(coords = c("x_3035", "y_3035"), crs = 3035, remove = FALSE)

merged_codes <- unlist(clusters$members)

rwl_merged_sites <- rwl_sea %>%
  filter(!site_code %in% merged_codes) %>%
  mutate(merged_from = NA_character_) %>%
  bind_rows(merged_rows) %>%
  arrange(source, site_code)


# CHECK FOR PROBLEM SERIES------------------------------------------------------
## removes problem series and flags chronologies which have problems with more
## than half of their series
rwl_checked <- remove_problem_series(rwl_merged_sites) %>% 
    filter(problem != TRUE) %>%
    select(-problem, -series_removed, -merged_from)


# FINAL CLEANING----------------------------------------------------------------
## remove chronologies shorter than 40 years (cf Fowler and Bridge 2017 on 
## insensitivity to series length). Also remove specific series which cannot be
## detrended using spline
rwl_clean <- rwl_checked %>% 
  filter(length >= 40) %>% 
  mutate(rwl = map(rwl, ~ remove_series(.x, c(
    "BAF-4A", "KDR0110A", "KDR0110B", "KDR0228A",
    "ZALIOJ05", "PAGEGI10", "Q20A", "Q20B"
  )))) %>%
  mutate(
    n_series = map_int(rwl, ncol),
    length = map_int(rwl, nrow),
    start = map_int(rwl, ~ min(as.numeric(rownames(.x)))),
    end = map_int(rwl, ~ max(as.numeric(rownames(.x))))
  ) 



# SAVE CLEANED CHRONOLOGIES-----------------------------------------------------
saveRDS(rwl_clean, "Outputs/06_CleanedChronologies_Step1.Rdata")
