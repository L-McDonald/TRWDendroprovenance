library(tidyverse)
library(sf)

source("Code/00_Functions.R")

# LOAD DATA---------------------------------------------------------------------
crn <- readRDS("Data/08_CalculatedChronologies_Step1.Rdata")


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
    std_cor = cor_std, std_n = n_cor_std,
    res_cor = cor_res, res_n = n_cor_res,
    ars_cor = cor_ars, ars_n = n_cor_ars
  ) %>%
  pivot_longer(
    cols = c(std_cor, std_n, res_cor, res_n, ars_cor, ars_n),
    names_to = c("type", ".value"),
    names_pattern = "(std|res|ars)_(cor|n)"
  ) %>% 
  select(site_1, site_2, type, dist_m, n, cor) %>% 
  filter(n >= 40)

write.csv(pairwise_all, "Data/09_RawCorrelations_Step1.csv", row.names = FALSE)