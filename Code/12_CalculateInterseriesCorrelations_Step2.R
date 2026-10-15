library(tidyverse)
library(sf)

source("Code/00_Functions.R")

# LOAD DATA---------------------------------------------------------------------
crn <- readRDS("Data/12_CalculatedChronologies_Step2.Rdata")


# EXTRACT CHRONOLOGIES OF EACH TYPE IN WIDE FORMAT------------------------------
std_wide <- extract_chron_wide(crn, "std")
res_wide <- extract_chron_wide(crn, "res")
ars_wide <- extract_chron_wide(crn, "ars")


# CALCULATE PAIRWISE CORRELATIONS FOR EACH CHRONOLOGY TYPE----------------------
cor_std <- pairwise_cor_long(std_wide, "cor_std")
cor_res <- pairwise_cor_long(res_wide, "cor_res")
cor_ars <- pairwise_cor_long(ars_wide, "cor_ars")


# PAIRWISE DISTANCES------------------------------------------------------------
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
    std_cor = cor_std,
    std_n = n_cor_std,
    res_cor = cor_res,
    res_n = n_cor_res,
    ars_cor = cor_ars,
    ars_n = n_cor_ars
  ) %>%
  pivot_longer(
    cols = c(std_cor, std_n, res_cor, res_n, ars_cor, ars_n),
    names_to = c("type", ".value"),
    names_pattern = "(std|res|ars)_(cor|n)"
  ) %>%
  select(site_1, site_2, type, dist_m, n, cor) %>%
  filter(n >= 40)


# CALCULATE SELF-CORRELATIONS---------------------------------------------------
## self correlation is the maximum of the rbar and the maximum correlation that
## the site has with any other site. It is the same for std, res, and ars.
self_cor <- crn %>%
  mutate(
    rbar = map_dbl(
      rwl,
      ~ .x %>% interseries.cor(method = "pearson") %>% pull(res.cor) %>% mean()
    ),
    max_cor = map_dbl(
      site_code,
      ~ pairwise_all %>%
        filter(site_1 == .x | site_2 == .x) %>%
        pull(cor) %>%
        max()
    ),
    self_cor = map2_dbl(rbar, max_cor, max),
  ) %>%
  select(site_code, self_cor, crn) %>%
  crossing(type = c("std", "res", "ars")) %>%
  mutate(
    site_1 = site_code,
    site_2 = site_code,
    dist_m = 0,
    n = map2_int(crn, type, ~ .x %>% pull(.y) %>% na.omit() %>% length())
  ) %>% 
  select(site_1, site_2, type, dist_m, n, cor = self_cor)


# ADD SELF-CORRELATIONS TO PAIRWISE CORRELATIONS AND SAVE-----------------------
cor_all <- pairwise_all %>% 
  bind_rows(self_cor)


# SAVE CORRELATIONS-------------------------------------------------------------
write.csv(pairwise_all, "Data/13_RawCorrelations_Step2.csv", row.names = FALSE)
