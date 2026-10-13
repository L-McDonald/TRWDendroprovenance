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
# species are combined into a single " & "-separated label if members differ.
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
    # split existing combined labels (e.g. "QUPE & QUCE") so species aren't repeated
    species = members$species %>%
      str_split(" & ") %>%
      unlist() %>%
      unique() %>%
      sort() %>%
      paste(collapse = " & "),
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
  rwl <- rwl[min(which(has_data)):max(which(has_data)), , drop = FALSE]

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


# CHRON.ARS.MOD-----------------------------------------------------------------
# modifies chron.ars from dplR to force an AR model with order 1 when the AIC
# selected model is order 0. Returns a list with the calculated chronologies and
# TRUE/FALSE indicating whether the model order was forced or not.
chron.ars.mod <- function (x, biweight = TRUE, maxLag = 10, firstAICmin = TRUE, 
    verbose = TRUE, prewhitenMethod = c("ar.yw", "arima.CSS-ML")) 
{
    pooledAR <- function(x, maxLag = maxLag, firstAICmin = TRUE) {
        nYrs <- dim(x)[1]
        nSeries <- dim(x)[2]
        x <- scale(x, center = TRUE, scale = FALSE)
        productSumMat <- matrix(0, maxLag + 1, 1)
        forcedAR <- FALSE
        for (k in 0:maxLag) {
            for (j in 1:nSeries) {
                for (i in 1:nSeries) {
                  series <- cbind(x[, j], x[, i])
                  mask <- which(rowSums(!is.na(series)) == 2)
                  if (length(mask) == 0) {
                    break
                  }
                  series1 <- x[mask, j]
                  series2 <- x[mask, i]
                  if (length(series2) <= k) {
                    break
                  }
                  series2 <- series2[1:(length(series2) - k)]
                  series2 <- c(rep(NA, k), series2)
                  mask2 <- which(rowSums(!is.na(cbind(series1, 
                    series2))) == 2)
                  if (length(mask2) == 0) {
                    break
                  }
                  productSumMat[k + 1, 1] = productSumMat[k + 
                    1, 1] + t(series1[mask2]) %*% series2[mask2]
                }
            }
        }
        productSumMat <- t(productSumMat)
        r0 <- productSumMat[1, 1]
        xAcf <- productSumMat/r0
        A <- acf2AR(xAcf) * -1
        vp_shortcut <- c(r0, c(r0 - -1 * A %*% matrix(t(productSumMat)[-1, 
            ])))
        aic2 <- nYrs * log(vp_shortcut) + 2 * 1:(maxLag + 1)
        names(aic2) <- paste0("ar(", 0:maxLag, ")")
        xAcf <- xAcf[1, ]
        names(xAcf) <- paste0("ar(", 0:maxLag, ")")
        getFirstMin <- function(y) {
            if (y[1] < y[2]) 
                out <- 1
            else out <- min(which(diff(y) > 0))
            return(out)
        }
        if (firstAICmin) {
            orderOut <- getFirstMin(aic2) - 1
        }
        else {
            orderOut <- which.min(aic2) - 1
        }
        if (orderOut == 0){
            orderOut <- 1
            forcedAR <- TRUE
        }
        res <- list(ACF = xAcf, ARcoefs = A, aic = aic2, orderOut = orderOut, forcedAR = forcedAR)
        return(res)
    }
    postAR <- function(x, phi) {
        mask <- is.na(x)
        x0 <- x
        x <- x[!mask]
        xMean <- mean(x)
        x <- x - xMean
        nx <- length(x)
        nPhi <- length(phi)
        xRev <- rev(x)
        xRevPad <- c(rep(0, nPhi), xRev)
        for (i in 1:nx) {
            for (j in 1:nPhi) {
                xRevPad[i + nPhi] <- xRevPad[i + nPhi] + phi[j] * 
                  xRevPad[i + nPhi - j]
            }
        }
        xRevInit <- xRevPad[-c(1:nPhi)]
        for (i in 1:nPhi) {
            cntr <- nPhi - i + 1
            xRevInit[cntr] <- 0
            for (j in 1:nPhi) {
                xRevInit[cntr] <- xRevInit[cntr] + phi[j] * xRevInit[cntr + 
                  j]
            }
        }
        xAR <- c(xRevInit[1:nPhi], x)
        for (i in 1:length(xRevInit)) {
            for (j in 1:nPhi) {
                xAR[i + nPhi] <- xAR[i + nPhi] + phi[j] * xAR[i + 
                  nPhi - j]
            }
        }
        xAR <- xAR[-c(1:nPhi)]
        x0[!mask] <- xAR + xMean
        return(x0)
    }
    prewhitenAR <- function(x, p) {
        mask <- is.na(x)
        x2 <- x[!mask]
        out <- ar(x2, aic = FALSE, order.max = p, demean = TRUE)
        x[!mask] <- out$resid + mean(x2)
        return(x)
    }
    prewhitenARIMA <- function(x, p) {
        mask <- is.na(x)
        x2 <- x[!mask]
        out <- arima(x2, order = c(p, 0, 0), method = "CSS-ML", 
            optim.control = list(maxit = 1000))
        x[!mask] <- out$residuals + mean(x2)
        return(x)
    }
    known.prewhitenMethods <- c("ar.yw", "arima.CSS-ML")
    prewhitenMethod2 <- match.arg(arg = prewhitenMethod, choices = known.prewhitenMethods, 
        several.ok = FALSE)
    samps <- rowSums(!is.na(x))
    if (!biweight) {
        stdCrn <- rowMeans(x, na.rm = TRUE)
    }
    else {
        stdCrn <- apply(x, 1, tbrm, C = 9)
    }
    outAR <- pooledAR(x, firstAICmin = firstAICmin, maxLag = maxLag)
    p <- outAR$orderOut
    if (verbose) {
        cat("Pooled AR Summary\n")
        cat("ACF\n")
        print(outAR$ACF)
        cat("AR Coefs\n")
        print(outAR$ARcoefs)
        cat("AIC\n")
        print(outAR$aic, digits = 4)
        cat("Selected Order\n")
        print(outAR$orderOut)
    }
    if (prewhitenMethod2 == "ar.yw") {
        RWIclean <- apply(x, 2, prewhitenAR, p = p)
    }
    if (prewhitenMethod2 == "arima.CSS-ML") {
        RWIclean <- apply(x, 2, prewhitenARIMA, p = p)
    }
    if (!biweight) {
        if (prewhitenMethod2 == "ar.yw") {
            resCrn <- prewhitenAR(rowMeans(RWIclean, na.rm = TRUE), 
                p = p)
        }
        if (prewhitenMethod2 == "arima.CSS-ML") {
            resCrn <- prewhitenARIMA(rowMeans(RWIclean, na.rm = TRUE), 
                p = p)
        }
    }
    else {
        if (prewhitenMethod2 == "ar.yw") {
            resCrn <- prewhitenAR(apply(RWIclean, 1, tbrm, C = 9), 
                p = p)
        }
        if (prewhitenMethod2 == "arima.CSS-ML") {
            resCrn <- prewhitenARIMA(apply(RWIclean, 1, tbrm, 
                C = 9), p = p)
        }
    }
    phi <- outAR$ARcoefs[p, ][1:p]
    RWIar <- apply(RWIclean, 2, postAR, -phi)
    if (!biweight) {
        arsCrn <- rowMeans(RWIar, na.rm = TRUE)
    }
    else {
        arsCrn <- apply(RWIar, 1, tbrm, C = 9)
    }
    out <- data.frame(std = stdCrn, res = resCrn, ars = arsCrn, 
        samp.depth = samps)
    row.names(out) <- row.names(x)
    class(out) <- c("crn", "data.frame")
    out_final <- list(chronologies = out, forcedAR = outAR$forcedAR)
    return(out_final)
}