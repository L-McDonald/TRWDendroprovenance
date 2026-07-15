# RWL_TIBBLE--------------------------------------------------------------------
# reads rwl files and stores them as a column in a tibble. Also adds some basic
# metadata calculated from the rwl file.
rwl_tibble <- function(path) {
  # Read the rwl file
  rwl_data <- read.rwl(path)

  #extract metadata
  metadata <- tibble(
    site_code = tools::file_path_sans_ext(basename(path)),
    n_series = ncol(rwl_data),
    length = nrow(rwl_data),
    start = min(as.numeric(row.names(rwl_data))),
    end = max(as.numeric(row.names(rwl_data)))
  )

  #combine metadata and rwl data into a list-column tibble
  metadata %>%
    mutate(rwl = list(rwl_data))
}


# REMOVE_HISTORICAL---------------------------------------------------------------
# removes series whose most recent ring is before min_last_year, then trims any
# leading/trailing years left with no data across all remaining series.
remove_historical <- function(rwl, min_last_year = 1900) {
  # identify series where most recent year is before min_last_year and remove
  keep <- summary(rwl)$last >= min_last_year
  rwl <- rwl[, keep, drop = FALSE]

  # trim rows with no data across all remaining series
  has_data <- rowSums(!is.na(rwl)) > 0
  rwl[min(which(has_data)):max(which(has_data)), , drop = FALSE]
}


# REMOVE_SERIES---------------------------------------------------------------
# drops named series (columns) from an rwl data frame, then trims any
# leading/trailing years left with no data across the remaining series.
remove_series <- function(rwl, series) {
  keep <- !colnames(rwl) %in% series
  rwl <- rwl[, keep, drop = FALSE]

  has_data <- rowSums(!is.na(rwl)) > 0
  rwl <- rwl[min(which(has_data)):max(which(has_data)), , drop = FALSE]

  class(rwl) <- c("rwl", "data.frame")
  rwl
}


# MERGE_RWL_LIST--------------------------------------------------------------
# combines a list of rwl data frames into one, aligning by year via dplR's
# combine.rwl(), with no additional filtering or duplicate checking.
merge_rwl_list <- function(rwl_list) {
  combined <- combine.rwl(rwl_list)
  combined <- combined[order(as.numeric(rownames(combined))), , drop = FALSE]
  class(combined) <- c("rwl", "data.frame")
  combined
}


# MERGE_CLUSTER-----------------------------------------------------------------
# combines a cluster of sites into a single metadata row + merged rwl data frame;
# site_code/site_name/source are concatenated across members, and lat/lon/elev
# are averaged (members are already within 1 km of each other); x_3035/y_3035
# are recomputed from the averaged lat/lon rather than also being averaged.
merge_cluster <- function(members) {
  rwl_combined <- merge_rwl_list(members$rwl)

  mean_lat <- mean(members$lat)
  mean_lon <- mean(members$lon)
  mean_coords_3035 <- tibble(lon = mean_lon, lat = mean_lat) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(3035) %>%
    st_coordinates()

  tibble(
    source = paste(unique(members$source), collapse = "_"),
    site_code = paste(members$site_code, collapse = "_"),
    site_name = paste(unique(members$site_name), collapse = "_"),
    start = min(as.numeric(rownames(rwl_combined))),
    end = max(as.numeric(rownames(rwl_combined))),
    length = end - start + 1,
    n_series = ncol(rwl_combined),
    species = unique(members$species),
    lat = mean_lat,
    lon = mean_lon,
    x_3035 = mean_coords_3035[, "X"],
    y_3035 = mean_coords_3035[, "Y"],
    elev = mean(members$elev),
    rwl = list(rwl_combined),
    merged_from = paste(members$site_code, collapse = ", ")
  )
}


# REMOVE_PROBLEM_SERIES---------------------------------------------------------
# checks for and removes problem series within a chronology and flags chronologies
# which have overall problems according to the method of van Sluijs et al.
remove_problem_series <- function(oaks) {
  cleaned_rows <- list() # Store cleaned data row-by-row
  site_codes <- oaks$site_code

  for (i in seq_along(site_codes)) {
    site <- site_codes[i]
    rwl <- filter(oaks, site_code == site) %>% pluck("rwl", 1)

    # Ensure rwl is not empty before proceeding
    if (is.null(rwl) || ncol(rwl) == 0) {
      message(paste0(
        site,
        " has no valid data and will be marked as problematic"
      ))
      problematic <- TRUE
      n.removed <- NA
    } else {
      n.all <- ncol(rwl)
      n.half <- floor(n.all / 2)

      c <- interseries.cor(rwl)
      max.p <- max(c$p.val)

      while (max.p > 0.01 & ncol(rwl) >= n.half) {
        if (ncol(rwl) < n.half) {
          message(paste0(
            site,
            " has internal issues but will still be included in the output"
          ))
          break
        }

        worst.series <- row.names(c)[which.max(c$p.val)]
        rwl <- rwl %>% select(-all_of(worst.series))
        c <- interseries.cor(rwl)
        max.p <- max(c$p.val)
      }

      n.removed <- n.all - ncol(rwl)
      problematic <- ifelse(ncol(rwl) < n.half, TRUE, FALSE)
    }

    # Create row for both good and problem sites
    new.row <- update_metadata(rwl) %>%
      mutate(
        site_code = site,
        source = filter(oaks, site_code == site) %>% pull("source"),
        site_name = filter(oaks, site_code == site) %>% pull("site_name"),
        species = filter(oaks, site_code == site) %>% pull("species"),
        lat = filter(oaks, site_code == site) %>% pull("lat"),
        lon = filter(oaks, site_code == site) %>% pull("lon"),
        x_3035 = filter(oaks, site_code == site) %>% pull("x_3035"),
        y_3035 = filter(oaks, site_code == site) %>% pull("y_3035"),
        elev = filter(oaks, site_code == site) %>% pull("elev"),
        merged_from = filter(oaks, site_code == site) %>% pull("merged_from"),
        series_removed = n.removed,
        problem = problematic # Now always assigned
      ) %>%
      relocate(
        source,
        site_code,
        site_name,
        start,
        end,
        length,
        n_series = n.series,
        species,
        lat,
        lon,
        x_3035,
        y_3035,
        elev,
        merged_from,
        problem,
        series_removed,
        rwl
      )

    cleaned_rows[[length(cleaned_rows) + 1]] <- new.row
  }

  # Combine cleaned rows into a final cleaned data frame
  oaks.cln <- bind_rows(cleaned_rows)

  return(oaks.cln)
}


# UPDATE_METADATA--------------------------------------------------------------
# removes empty rows from a rwl dataframe after series have been removed and 
# updates the basic metadata.
update_metadata <- function(rwl) {
  # trim rows with no data across all remaining series
  has_data <- rowSums(!is.na(rwl)) > 0
  rwl[min(which(has_data)):max(which(has_data)), , drop = FALSE]

  #extract metadata
  metadata <- tibble(
    n.series = ncol(rwl),
    length = nrow(rwl),
    start = min(as.numeric(row.names(rwl))),
    end = max(as.numeric(row.names(rwl)))
  )

  #combine metadata and rwl data into a list-column tibble
  metadata %>%
    mutate(rwl = list(rwl))
}
