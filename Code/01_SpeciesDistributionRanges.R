library(tidyverse)
library(sf)
library(rnaturalearth)


# UNZIP EUFORGEN SPECIES DISTRIBUTION RANGES------------------------------------
unzip(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset_20260126.zip",
    exdir = "Data/SpeciesDistributionRanges"
)


# LOAD COAST--------------------------------------------------------------------
coast <- ne_countries(scale = 10, continent = c("Africa", "Europe", "Asia")) %>%
    st_make_valid() %>%
    st_crop(c(xmin = -15, xmax = 90, ymin = 20, ymax = 90)) %>%
    st_geometry()


# LOAD SPECIES DISTRIBUTION RANGES----------------------------------------------

# Quercus petraea
pe <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus petraea/shapefiles/",
    layer = "quercus_petraea_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    mutate(species = "Quercus petraea")

# Quercus robur
ro <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus robur/shapefiles/",
    layer = "quercus_robur_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    mutate(species = "Quercus robur")

# Quercus cerris
ce <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus cerris/shapefiles/",
    layer = "quercus_cerris_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    mutate(species = "Quercus cerris")

# Quercus frainetto
fr <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus frainetto/shapefiles/",
    layer = "quercus_frainetto_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    mutate(species = "Quercus frainetto")

# merge all SDRs
all_species <- bind_rows(pe, ro, ce, fr)

merged <- st_union(c(
    st_geometry(pe),
    st_geometry(ro),
    st_geometry(ce),
    st_geometry(fr)
))

st_write(merged, "Data/SpeciesDistributionRanges/QURO_QUPE_QUCE_QUFR.gpkg")


# BOUNDING BOX------------------------------------------------------------------
st_bbox(merged)

ggplot() +
    geom_sf(data = coast) +
    geom_sf(data = x, fill = "steelblue")
