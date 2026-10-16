# ORIGINAL SCRIPT FROM: https://github.com/chguiterman/read_ITRDB_rawData/tree/v1.1.0

library(tidyverse)
library(dplR)
library(httr)
library(jsonlite)
library(lubridate)
library(furrr)
library(sf)
library(rnaturalearth)


# LOAD COAST and SPECIES DISTRIBUTION RANGE-------------------------------------
coast <- ne_countries(scale = 10, continent = c("Africa", "Europe", "Asia")) %>%
  st_make_valid() %>%
  st_crop(c(xmin = -25, xmax = 90, ymin = 10, ymax = 90)) %>%
  st_geometry()

range <- read_sf("Data/SpeciesDistributionRanges/QURO_QUPE_QUCE_QUFR_QUPU.gpkg")


# ITRDB SEARCH PARAMETERS-------------------------------------------------------
base_url <- "https://www.ncei.noaa.gov/paleo-search/study/search.json?"

search_params <- list(
  dataPublisher = "NOAA",
  dataTypeId = 18,
  minLat = 33,
  maxLat = 62,
  minLon = -10,
  maxLon = 60,
  species = "QUCE|QUFR|QUPE|QURO|QUPU",
  metadataOnly = "true"
)

# Ping WDS-Paleo API------------------------------------------------------------
resp <- GET(base_url, query = search_params)

parsed <- fromJSON(
  content(resp, "text", encoding = "UTF-8"),
  simplifyDataFrame = TRUE
)

p_dat <- structure(list(
  url = url,
  content = parsed,
  response = resp
))

in_list <- unnest(p_dat$content$study, site)

# Filter for ITRDB sites
itrdb_list <- in_list %>%
  filter(str_detect(studyName, "ITRDB"))


# 1. Create ITRDB site metadata table -----------------------------------------

# Helper function
get_ITRDB_Meta <- function(itrdb_list) {
  out_df <- itrdb_list %>%
    transmute(
      NOAAStudyId,
      siteName,
      studyCode,
      investigators,
      first_year = earliestYearCE,
      last_year = mostRecentYearCE
    )
  out_noaa_params <- itrdb_list %>%
    transmute(
      NOAAStudyId,
      doi,
      contr_year = year(contributionDate),
      NOAASiteId,
      url = onlineResourceLink
    )
  # Extract Lat-Long coordinates
  out_coords <- out_df %>%
    select(NOAAStudyId) %>%
    mutate(
      coords = pluck(itrdb_list, "geo", "geometry", "coordinates"),
      latitude = map_dbl(coords, ~ as.numeric(.x[1])),
      longitude = map_dbl(coords, ~ as.numeric(.x[2])),
      elevation = as.numeric(pluck(
        itrdb_list,
        "geo",
        "properties",
        "minElevationMeters"
      ))
    ) %>%
    select(-coords)
  # Extract species codes
  out_spp <- out_df %>%
    select(NOAAStudyId) %>%
    mutate(
      species_list = map(itrdb_list[["paleoData"]], "species"),
      species_df = map(species_list, ~ .x[[1]]$speciesCode),
      species = map_chr(species_df, str_c, collapse = ", ")
    ) %>%
    select(NOAAStudyId, species)
  out_age <- out_df %>%
    select(NOAAStudyId) %>%
    mutate(date_type = map_chr(itrdb_list[["paleoData"]], "timeUnit"))
  # Extract published references
  pub_list <- unnest(itrdb_list, publication)
  if (nrow(pub_list) > 0) {
    out_pub <- pub_list %>%
      select(NOAAStudyId, citation) %>%
      group_by(NOAAStudyId) %>%
      summarize(reference = str_c(citation, collapse = "; "))
  } else {
    out_pub <- data.frame(NOAAStudyId = NA, citation = NA)
  }
  # Combine
  itrdb_site_meta <- out_df %>%
    left_join(out_age, by = "NOAAStudyId") %>%
    left_join(out_coords, by = "NOAAStudyId") %>%
    left_join(out_spp, by = "NOAAStudyId") %>%
    left_join(out_pub, by = "NOAAStudyId") %>%
    left_join(out_noaa_params, by = "NOAAStudyId") %>%
    filter(!is.na(studyCode))

  return(itrdb_site_meta)
}

itrdb_site_meta <- get_ITRDB_Meta(itrdb_list)


# 2. ITRDB list of raw measurement files -----------------------------------------------

itrdb_rawmeas_files <- itrdb_list %>%
  select(NOAAStudyId, first_year = earliestYearCE) %>%
  mutate(
    data = map(itrdb_list[["paleoData"]], "dataFile"),
    files = map(data, ~ .x[[1]][c("linkText", "urlDescription", "fileUrl")]),
    vars = map(data, ~ .x[[1]]$variables)
  ) %>%
  select(-data) %>%
  unnest(cols = c(files, vars)) %>%
  filter(
    urlDescription == "Raw Measurements",
    !str_detect(linkText, "-noaa")
  ) %>%
  unnest(vars) %>%
  filter(cvWhat != "age variable>age") %>%
  select(NOAAStudyId, first_year, linkText, fileUrl, cvWhat, cvUnit) %>%
  filter(str_detect(cvWhat, "total ring width")) # whole ring width only


# 3. ITRDB RWL files ------------------------------------------------------

# Helper function
read_itrdb_rwl <- function(first_year, fileUrl) {
  require(dplR)
  if (any(c(length(first_year) > 1, length(fileUrl) > 1))) {
    abort("the function can only handle one case at a time")
  }
  # read file
  if (first_year < -999) {
    read.rwl(fileUrl, long = TRUE, format = "tucson")
  } else {
    read.rwl(fileUrl, format = "tucson")
  }
}


# clean up
rm(in_list, itrdb_list, resp, p_dat, parsed)

# Read rwls from the HTTPS server (this may take some time)
all_rwl <- itrdb_rawmeas_files %>%
  transmute(
    NOAAStudyId,
    RWL = future_map2(
      first_year,
      fileUrl,
      ~ try(read_itrdb_rwl(.x, .y)),
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )
  )

# Sometimes sites are loaded on the ITRDB before the data are released,
# triggering access errors in the script above. This check_files will identify
# those so they can be removed form further use, as below
check_files <- all_rwl %>%
  mutate(
    check = map_chr(
      RWL,
      ~ {
        # if (class(.x) %in% "try-error") {
        #   out <- .x[[1]]
        # }
        if ("rwl" %in% class(.x)) {
          out <- "Good"
        } else {
          out = "Error"
        }
        out
      }
    )
  ) %>%
  filter(check == "Error")


# 4. Output data -------------------------------------------------------------

# generate output tibble and filter to end after 1900 
out_dat <- all_rwl %>%
  filter(!NOAAStudyId %in% check_files$NOAAStudyId) %>%
  left_join(itrdb_site_meta, by = "NOAAStudyId") %>%
  select(
    -NOAAStudyId,
    investigators,
    -date_type,
    -reference,
    -doi,
    -contr_year,
    -NOAASiteId,
    -url
  ) %>%
  mutate(
    length = as.integer((last_year - first_year) + 1),
    n_series = map_int(RWL, ncol)
  ) %>%
  select(
    site_code = studyCode,
    site_name = siteName,
    start = first_year,
    end = last_year,
    length,
    n_series,
    species,
    lat = latitude,
    lon = longitude,
    elev = elevation,
    rwl = RWL
  ) %>%
  filter(end >= 1900)


# walk2(
#   out_dat$RWL,
#   out_dat$linkText,
#   ~ write.tucson(.x, str_glue("Data/RWL_ITRDB/{.y}"))
# )

# Save metadata files
# write_csv(itrdb_rawmeas_files, "Data/ITRDB_raw_measurement_files.csv")
# write_csv(itrdb_site_meta, "Data/ITRDB_site_metadata.csv")

# PLOT CHRONOLOGIES-------------------------------------------------------------
plot_limits <- st_bbox(range) %>%
  st_as_sfc() %>%
  st_transform(3035) %>%
  st_bbox()

plot_crn <- out_dat %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

ggplot() +
  geom_sf(data = coast, fill = "grey95") +
  geom_sf(data = range, fill = "grey80", linetype = "longdash") +
  geom_sf(data = coast, fill = NA) +
  geom_sf(data = plot_crn, aes(colour = species)) +
  coord_sf(
    crs = 3035,
    xlim = c(plot_limits["xmin"], plot_limits["xmax"]),
    ylim = c(plot_limits["ymin"], plot_limits["ymax"]),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "azure1"),
    legend.position = "inside",
    legend.justification.inside = c(1, 1),
    legend.position.inside = c(1, 1),
    legend.background = element_rect(fill = NA),
    legend.box.background = element_rect(colour = NA, fill = NA),
    legend.key = element_rect(fill = NA)
  )


# FILTER AND SAVE SELECTED CHRONOLOGIES-----------------------------------------

# remove chronologies which are also in CTRL dataset

save_dat <- out_dat %>% 
  filter(site_code != "GREE006") %>% 
  filter(site_code != "TURK028") %>% 
  mutate(source = "ITRDB", .before = everything())

saveRDS(save_dat, "Outputs/01_InitialChronologies_ITRDB.Rdata")
