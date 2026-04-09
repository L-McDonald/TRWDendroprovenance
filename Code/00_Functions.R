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
