library(tidyverse)
library(sf)

source("Code/00_Functions.R")

# LOAD DATA---------------------------------------------------------------------
crn <- readRDS("Data/08_CalculatedChronologies.Rdata")


# EXPORT SITE LOCATIONS-----------------------------------------------------------
# site_code/lat/lon (+ identifying metadata), used e.g. by Code/10_DownloadElevation.py
# to sample FABDEM elevation via Google Earth Engine.
crn %>%
  select(site_code, site_name, source, species, lat, lon, elev) %>%
  distinct() %>%
  write.csv("Data/site_locations.csv", row.names = FALSE)


# BUILD WIDE-FORMAT CHRONOLOGY MATRICES-----------------------------------------
# extracts one chronology type (std/res/ars) from the crn list-column into a
# year x site_code matrix, needed for pairwise correlation.
extract_chron_wide <- function(crn, type) {
  crn %>%
    select(site_code, crn) %>%
    mutate(
      series = map(crn, ~ tibble(year = as.integer(rownames(.x)), value = .x[[type]]))
    ) %>%
    select(site_code, series) %>%
    unnest(series) %>%
    pivot_wider(names_from = site_code, values_from = value) %>%
    arrange(year)
}

std_wide <- extract_chron_wide(crn, "std")
res_wide <- extract_chron_wide(crn, "res")
ars_wide <- extract_chron_wide(crn, "ars")


# PAIRWISE CORRELATIONS----------------------------------------------------------
# computes pairwise Pearson correlations (and number of overlapping years) between
# all chronologies of a given type, using pairwise complete observations. Returned
# in long/tidy format with one row per unique site pair.
pairwise_cor_long <- function(wide, value_name) {
  mat <- as.matrix(select(wide, -year))

  cor_mat <- cor(mat, use = "pairwise.complete.obs")
  n_mat <- crossprod(!is.na(mat)) # number of overlapping years per pair

  cor_long <- as_tibble(cor_mat, rownames = "site_1") %>%
    pivot_longer(-site_1, names_to = "site_2", values_to = value_name)

  n_long <- as_tibble(n_mat, rownames = "site_1") %>%
    pivot_longer(-site_1, names_to = "site_2", values_to = paste0("n_", value_name))

  left_join(cor_long, n_long, by = c("site_1", "site_2")) %>%
    filter(site_1 < site_2) # keep each unique pair once
}

cor_std <- pairwise_cor_long(std_wide, "cor_std")
cor_res <- pairwise_cor_long(res_wide, "cor_res")
cor_ars <- pairwise_cor_long(ars_wide, "cor_ars")


# PAIRWISE DISTANCES--------------------------------------------------------------
# computes pairwise Euclidean distances (m) between chronology locations using
# EPSG:3035 coordinates.
sites_sf <- crn %>%
  select(site_code, x_3035, y_3035) %>%
  st_as_sf(coords = c("x_3035", "y_3035"), crs = 3035, remove = FALSE)

dist_mat <- st_distance(sites_sf)
dist_mat <- matrix(as.numeric(dist_mat), nrow = nrow(dist_mat))
rownames(dist_mat) <- sites_sf$site_code
colnames(dist_mat) <- sites_sf$site_code

dist_long <- as_tibble(dist_mat, rownames = "site_1") %>%
  pivot_longer(-site_1, names_to = "site_2", values_to = "dist_m") %>%
  filter(site_1 < site_2)


# COMBINE---------------------------------------------------------------------
# long format with one row per unique chronology pair per type (std/res/ars),
# giving site1, site2, dist, type, n (overlapping years), and cor.
pairwise_all <- dist_long %>%
  left_join(cor_std, by = c("site_1", "site_2")) %>%
  left_join(cor_res, by = c("site_1", "site_2")) %>%
  left_join(cor_ars, by = c("site_1", "site_2")) %>%
  rename(
    std_cor = cor_std, std_n = n_cor_std,
    res_cor = cor_res, res_n = n_cor_res,
    ars_cor = cor_ars, ars_n = n_cor_ars
  ) %>%
  pivot_longer(
    cols = c(std_cor, std_n, res_cor, res_n, ars_cor, ars_n),
    names_to = c("type", ".value"),
    names_pattern = "(std|res|ars)_(cor|n)"
  ) %>% 
  select(site_1, site_2, type, dist_m, n, cor)

# NOTE: 37 of 19,306 pairs have zero overlapping years (cor = NA), and 64 pairs
# overlap in fewer than 10 years, so those correlations are unreliable.

write.csv(pairwise_all, "Data/09_RawCorrelations.csv", row.names = FALSE)


# NEST CORRELATIONS--------------------------------------------------------------
# nests each site's pairwise correlations with every other site (across std,
# res, and ars) as a data frame in a new list-column of crn.
site_correlations <- pairwise_all %>%
  transmute(site_code = site_1, other_site = site_2, dist_m, type, n, cor) %>%
  bind_rows(
    pairwise_all %>%
      transmute(site_code = site_2, other_site = site_1, dist_m, type, n, cor)
  ) %>%
  arrange(site_code, other_site, type) %>%
  nest(correlations = c(other_site, dist_m, type, n, cor))

crn <- crn %>%
  left_join(site_correlations, by = "site_code") %>%
  relocate(correlations, .after = crn)

#saveRDS(crn, "Data/09_ChronologiesWithCorrelations.Rdata")


crn_long <- crn %>%
  mutate(type = list(c("std", "res", "ars"))) %>%
  unnest(type) %>%
  mutate(
    correlations = map2(correlations, type, ~ filter(.x, type == .y) %>% select(-type))
  ) %>%
  relocate(type, .after = site_code)

dim(crn)

