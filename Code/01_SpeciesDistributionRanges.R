library(tidyverse)
library(sf)
library(rnaturalearth)


# UNZIP EUFORGEN SPECIES DISTRIBUTION RANGES------------------------------------
# unzip(
#     "Data/SpeciesDistributionRanges/chorological_maps_dataset_20260126.zip",
#     exdir = "Data/SpeciesDistributionRanges"
# )


# LOAD COAST--------------------------------------------------------------------
coast <- ne_countries(scale = 10, continent = c("Africa", "Europe", "Asia")) %>%
    st_make_valid() %>%
    st_crop(c(xmin = -25, xmax = 90, ymin = 10, ymax = 90)) %>%
    st_geometry()


# LOAD SPECIES DISTRIBUTION RANGES----------------------------------------------

# Quercus petraea
pe <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus petraea/shapefiles/",
    layer = "quercus_petraea_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    st_union() %>%
    st_sf() %>%
    mutate(species = "Quercus petraea")

# Quercus robur
ro <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus robur/shapefiles/",
    layer = "quercus_robur_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    st_union() %>%
    st_sf() %>%
    mutate(species = "Quercus robur")

# Quercus cerris
ce <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus cerris/shapefiles/",
    layer = "quercus_cerris_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    st_union() %>%
    st_sf() %>%
    mutate(species = "Quercus cerris")

# Quercus frainetto
fr <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus frainetto/shapefiles/",
    layer = "quercus_frainetto_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    st_union() %>%
    st_sf() %>%
    mutate(species = "Quercus frainetto")

# Quercus pubescens
pu <- read_sf(
    "Data/SpeciesDistributionRanges/chorological_maps_dataset/Quercus pubescens/shapefiles/",
    layer = "quercus_pubescens_plg"
) %>%
    st_intersection(coast) %>%
    st_make_valid() %>%
    st_union() %>%
    st_sf() %>%
    mutate(species = "Quercus pubescens")

# merge all SDRs
all_species <- bind_rows(pe, ro, ce, fr, pu)

merged <- st_union(c(
    st_geometry(pe),
    st_geometry(ro),
    st_geometry(ce),
    st_geometry(fr),
    st_geometry(pu)
)) %>% 
    st_make_valid()

merged_3035 <- st_transform(merged, crs = 3035)

st_write(merged, "Data/SpeciesDistributionRanges/QURO_QUPE_QUCE_QUFR_QUPU.gpkg")
st_write(merged_3035, "Data/SpeciesDistributionRanges/QURO_QUPE_QUCE_QUFR_QUPU_3035.gpkg")

st_write(merged, "Outputs/EnvironmentalVariables/Region_4326.gpkg")
st_write(merged_3035, "Outputs/EnvironmentalVariables/Region_3035.gpkg")


# BOUNDING BOX------------------------------------------------------------------
st_bbox(merged)


# PLOT SPECIES DISTRIBUTION RANGES----------------------------------------------
plot_limits <- st_bbox(merged) %>%
    st_as_sfc() %>%
    st_transform(3035) %>%
    st_bbox()

ggplot() +
    geom_sf(data = coast, fill = "grey95", colour = NA) +
    geom_sf(data = merged, linetype = "longdash", fill = "grey80", stroke = 2) +
    geom_sf(data = all_species, fill = "steelblue") +
    geom_sf(data = coast, fill = NA) +
    geom_text(
        data = all_species,
        aes(x = I(0.95), y = I(0.95), label = species),
        hjust = 1,
        vjust = 1,
        family = "serif",
        fontface = "bold.italic",
        size = 7,
        size.unit = "pt"
    ) +
    facet_wrap(~species, ncol = 2) +
    theme_bw() +
    theme(
        text = element_text(family = "serif"),
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 6)
    ) +
    coord_sf(
        crs = 3035,
        xlim = c(plot_limits["xmin"], plot_limits["xmax"]),
        ylim = c(plot_limits["ymin"], plot_limits["ymax"]),
        expand = FALSE
    )

ggsave(
    "Figures/SpeciesDistributionRanges.png",
    dpi = 300,
    width = 190,
    height = 170,
    units = "mm"
)
