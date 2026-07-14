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
