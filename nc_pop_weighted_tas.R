#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
  library(lubridate)
})

# Compute population-weighted (and unweighted) zonal mean of `tas` from a NetCDF
# Usage (CLI):
# Rscript nc_pop_weighted_tas.R hunan.nc pop.tif hunan.shp out.csv [time_mode] [var_name] [lon_var] [lat_var] [resample_method]
#  - time_mode: 'avg' (default) | 'first' | 'last' | <number index> | 'all'
#  - var_name: NetCDF variable name for tas (default 'tas')
#  - resample_method: 'near' (default) or 'bilinear'

compute_nc_pop_weighted_tas <- function(nc_path, pop_tif, shp_path,
                                        out_csv = NULL,
                                        time_mode = "avg",
                                        var_name = "tas",
                                        lon_name = "lon",
                                        lat_name = "lat",
                                        resample_method = c("near", "bilinear"),
                                        write_csv = TRUE,
                                        verbose = TRUE) {
  resample_method <- match.arg(resample_method)

  if (!file.exists(nc_path)) stop("NC file not found: ", nc_path)
  if (!file.exists(pop_tif)) stop("pop raster not found: ", pop_tif)
  if (!file.exists(shp_path)) stop("shapefile not found: ", shp_path)

  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

  # read lon/lat
  lon <- tryCatch(ncdf4::ncvar_get(nc, lon_name), error = function(e) NULL)
  if (is.null(lon) && !is.null(nc$dim[[lon_name]]$vals)) lon <- nc$dim[[lon_name]]$vals
  lat <- tryCatch(ncdf4::ncvar_get(nc, lat_name), error = function(e) NULL)
  if (is.null(lat) && !is.null(nc$dim[[lat_name]]$vals)) lat <- nc$dim[[lat_name]]$vals
  if (is.null(lon) || is.null(lat)) stop("Unable to read lon/lat from NetCDF.")
  if (length(dim(lon)) > 1 || length(dim(lat)) > 1) stop("This script supports regular 1D lon/lat grids only.")

  # time
  time_name <- if (!is.null(nc$var$time)) "time" else if (!is.null(nc$dim$time)) "time" else NULL
  time_vals <- if (!is.null(time_name)) tryCatch(ncdf4::ncvar_get(nc, time_name), error = function(e) NULL) else NULL
  time_units <- if (!is.null(time_name)) {
    at <- tryCatch(ncdf4::ncatt_get(nc, time_name), error = function(e) list())
    as.character(at$units %||% "")
  } else ""

  decode_time <- function(vals, units) {
    if (is.null(vals) || !nzchar(units)) return(rep(NA, length(vals)))
    m <- regexec("^(seconds|minutes|hours|days|day) since (.+)$", units, ignore.case = TRUE)
    r <- regmatches(units, m)[[1]]
    if (length(r) < 3) return(rep(NA, length(vals)))
    unit <- tolower(r[2]); origin_str <- r[3]
    origin <- tryCatch(as.POSIXct(origin_str, tz = "UTC"), error = function(e) NA)
    if (is.na(origin)) {
      origin <- tryCatch(lubridate::parse_date_time(origin_str, orders = c("Ymd HMS", "Ymd", "ymd HMS", "ymd"), tz = "UTC"), error = function(e) NA)
    }
    if (is.na(origin)) return(rep(NA, length(vals)))
    sec <- switch(unit, seconds = 1, minutes = 60, hours = 3600, day = 86400, days = 86400, NA)
    if (is.na(sec)) return(rep(NA, length(vals)))
    lubridate::as_datetime(origin) + lubridate::dseconds(vals * sec)
  }
  time_str <- if (!is.null(time_vals)) as.character(decode_time(time_vals, time_units)) else NULL

  # read variable and reorder to [lat, lon, time]
  var <- nc$var[[var_name]]
  if (is.null(var)) stop(sprintf("Variable '%s' not found in NetCDF.", var_name))
  arr <- ncdf4::ncvar_get(nc, var_name)
  dim_names <- sapply(var$dim, function(x) x$name)

  reorder_to <- function(a, from_names, target_names) {
    idx <- match(target_names, from_names)
    idx <- idx[!is.na(idx)]
    if (length(idx) == length(from_names)) {
      aperm(a, idx)
    } else {
      if (all(c("lon","lat") %in% from_names)) {
        idx2 <- match(c("lon","lat"), from_names)
        aperm(a, idx2)
      } else {
        stop("Unexpected variable dimensions; cannot reorder.")
      }
    }
  }
  arr <- reorder_to(arr, dim_names, c("lat","lon","time"))

  nx <- length(lon); ny <- length(lat)
  if (dim(arr)[1] != ny || dim(arr)[2] != nx) stop("Data dimensions do not match lon/lat lengths.")
  nt <- if (length(dim(arr)) >= 3) dim(arr)[3] else 1

  # ensure lon ascending and lat descending (north->south)
  lon_order <- order(lon)
  lon_sorted <- lon[lon_order]
  arr <- arr[, lon_order, , drop = FALSE]
  lat_order_desc <- order(lat, decreasing = TRUE)
  arr <- arr[lat_order_desc, , , drop = FALSE]
  lat_sorted <- lat[lat_order_desc]

  resx <- median(diff(lon_sorted), na.rm = TRUE)
  resy <- abs(median(diff(lat_sorted), na.rm = TRUE))
  xmin <- min(lon_sorted) - resx/2; xmax <- max(lon_sorted) + resx/2
  ymin <- min(lat_sorted) - resy/2; ymax <- max(lat_sorted) + resy/2

  # open pop raster early so we can use its CRS for the NetCDF template
  pop_r <- terra::rast(pop_tif)
  if (terra::nlyr(pop_r) > 1) pop_r <- pop_r[[1]]
  pop_crs <- terra::crs(pop_r, proj = TRUE)
  # Determine whether lon/lat are geographic (degrees)
  is_geog <- all(abs(lon_sorted) <= 360) && all(abs(lat_sorted) <= 90)
  tmpl_crs <- if (is_geog) "EPSG:4326" else if (!is.na(pop_crs) && nzchar(pop_crs)) pop_crs else "EPSG:4326"

  tmpl <- terra::rast(ncols = nx, nrows = ny, ext = terra::ext(xmin, xmax, ymin, ymax), crs = tmpl_crs)
  make_layer <- function(mat_ready) {
    r <- tmpl
    terra::values(r) <- as.vector(t(mat_ready))
    r
  }

  # read shapefile
  vec <- sf::st_read(shp_path, quiet = !verbose)

  # transform vector to pop CRS if possible
  pop_crs <- terra::crs(pop_r, proj = TRUE)
  if (!is.na(pop_crs) && nzchar(pop_crs) && !is.na(sf::st_crs(vec))) {
    vec <- tryCatch(sf::st_transform(vec, pop_crs), error = function(e) {
      warning(sprintf("Failed to transform shapefile to pop CRS: %s. Proceeding without transform.", conditionMessage(e)))
      vec
    })
  }

  # helper for city field
  choose_city_field <- function(df) {
    nms <- names(df)
    candidates <- c("ShiName", "CITY", "CITY_NAME", "NAME", "NAME_1", "NAME_2",
                    "ADMIN_NAME", "ADM2_NAME", "COUNTY", "PREF", "PREFECTURE",
                    "CITY_CN", "NAME_CH", "NAME_ZH")
    nm_upper <- toupper(nms)
    for (cand in candidates) {
      idx <- which(nm_upper == toupper(cand))
      if (length(idx) > 0) return(nms[idx[1]])
    }
    char_cols <- nms[sapply(df, function(x) is.character(x) || is.factor(x))]
    if (length(char_cols) > 0) return(char_cols[1])
    non_geom <- setdiff(nms, attr(df, "sf_column"))
    if (length(non_geom) > 0) return(non_geom[1])
    stop("No suitable city/name field found in shapefile attributes.")
  }
  city_field <- choose_city_field(vec)
  if (verbose) message("Using city field: ", city_field)

  select_time_index <- function(mode, t_str, t_len) {
    if (grepl("^[0-9]+$", mode)) {
      i <- as.integer(mode); if (i < 1 || i > t_len) stop("time index out of range"); return(i)
    }
    m <- tolower(mode)
    if (m == "first") return(1)
    if (m == "last") return(t_len)
    return(NA_integer_)
  }

  v_spat <- terra::vect(vec)

  results <- list()
  process_layer <- function(mat, tlabel = NULL, idx = NA_integer_) {
    r_raw <- make_layer(mat)
    # Diagnostics: report CRS/extents and non-NA counts before resampling
    crs_r_raw <- terra::crs(r_raw, proj = TRUE)
    crs_pop <- terra::crs(pop_r, proj = TRUE)
    ext_r_raw <- terra::ext(r_raw)
    ext_pop <- terra::ext(pop_r)
    ncell_r_raw <- terra::ncell(r_raw)
    ncell_pop <- terra::ncell(pop_r)
    non_na_r_raw <- sum(!is.na(terra::values(r_raw)))
    non_na_pop <- sum(!is.na(terra::values(pop_r)))
    if (verbose) {
      message(sprintf("Layer diagnostics (pre-resample): r_raw CRS=%s pop CRS=%s", crs_r_raw, crs_pop))
      message(sprintf("r_raw extent: xmin=%s xmax=%s ymin=%s ymax=%s", ext_r_raw[1], ext_r_raw[2], ext_r_raw[3], ext_r_raw[4]))
      message(sprintf("pop extent: xmin=%s xmax=%s ymin=%s ymax=%s", ext_pop[1], ext_pop[2], ext_pop[3], ext_pop[4]))
      message(sprintf("ncell r_raw=%d ncell pop=%d nonNA r_raw=%d nonNA pop=%d", ncell_r_raw, ncell_pop, non_na_r_raw, non_na_pop))
    }

    # reproject r_raw to pop CRS if pop CRS is available and r_raw is geographic or different CRS
    if (!is.na(crs_pop) && nzchar(crs_pop)) {
      # only project if CRS differ (and r_raw has a CRS)
      if (is.na(crs_r_raw) || crs_r_raw != crs_pop) {
        if (verbose) message("Projecting NetCDF-derived layer to population CRS before resampling...")
        r_proj <- tryCatch(terra::project(r_raw, crs_pop), error = function(e) {
          warning(sprintf("Project failed: %s; attempting to assign CRS and continue.", conditionMessage(e)))
          terra::crs(r_raw) <- crs_pop
          r_raw
        })
      } else {
        r_proj <- r_raw
      }
    } else {
      r_proj <- r_raw
    }

    tas_rs <- terra::resample(r_proj, pop_r, method = resample_method)
    non_na_rs <- sum(!is.na(terra::values(tas_rs)))
    ext_rs <- terra::ext(tas_rs)
    if (verbose) {
      message(sprintf("Post-resample: ncell_rs=%d nonNA_rs=%d", terra::ncell(tas_rs), non_na_rs))
      message(sprintf("tas_rs extent: xmin=%s xmax=%s ymin=%s ymax=%s", ext_rs[1], ext_rs[2], ext_rs[3], ext_rs[4]))
      # attempt min/max
      mm_raw <- try(terra::minmax(r_raw), silent = TRUE)
      mm_rs <- try(terra::minmax(tas_rs), silent = TRUE)
      if (!inherits(mm_raw, "try-error")) message(sprintf("r_raw min/max: %s", paste(mm_raw, collapse=",")))
      if (!inherits(mm_rs, "try-error")) message(sprintf("tas_rs min/max: %s", paste(mm_rs, collapse=",")))
    }

    mean_tas_df <- terra::extract(tas_rs, v_spat, fun = mean, na.rm = TRUE)
    pop_masked <- terra::mask(pop_r, tas_rs)
    prod_r <- pop_r * tas_rs
    sum_pop_tas <- terra::extract(prod_r, v_spat, fun = sum, na.rm = TRUE)
    sum_pop <- terra::extract(pop_masked, v_spat, fun = sum, na.rm = TRUE)

    df_mean <- as.data.frame(mean_tas_df)
    df_sum_tas <- as.data.frame(sum_pop_tas)
    df_sum_p <- as.data.frame(sum_pop)

    cities <- vec %>% st_drop_geometry() %>% mutate(.row_id = seq_len(n()))

    df <- data.frame(.row_id = df_mean$ID,
                     mean_tas = df_mean[[2]],
                     sum_pop_tas = df_sum_tas[[2]],
                     sum_pop = df_sum_p[[2]],
                     stringsAsFactors = FALSE)
    df <- left_join(cities %>% mutate(.row_id = seq_len(n())), df, by = ".row_id")
    df <- df %>% mutate(band = ifelse(is.null(tlabel), idx, as.character(tlabel)),
                        pw_mean = ifelse(is.na(sum_pop) | sum_pop == 0, NA_real_, sum_pop_tas / sum_pop)) %>%
      select(all_of(city_field), band, mean_tas, sum_pop_tas, sum_pop, pw_mean)
    names(df)[1] <- "city"
    return(df)
  }

  if (tolower(time_mode) == "avg" && nt > 1) {
    m <- apply(arr, MARGIN = c(1,2), FUN = function(x) mean(x, na.rm = TRUE))
    res_df <- process_layer(m, tlabel = "avg", idx = NA_integer_)
    results[[1]] <- res_df
  } else if (nt == 1) {
    m <- arr[,,1, drop = TRUE]
    results[[1]] <- process_layer(m, tlabel = time_str %||% NULL, idx = 1)
  } else if (tolower(time_mode) == "all" && nt > 1) {
    out_list <- list()
    for (ti in seq_len(nt)) {
      m <- arr[,,ti, drop = TRUE]
      out_list[[ti]] <- process_layer(m, tlabel = if (!is.null(time_str)) time_str[ti] else NULL, idx = ti)
    }
    results <- out_list
  } else {
    ti <- select_time_index(time_mode, time_str, nt)
    if (is.na(ti)) ti <- nt
    m <- arr[,,ti, drop = TRUE]
    results[[1]] <- process_layer(m, tlabel = if (!is.null(time_str)) time_str[ti] else NULL, idx = ti)
  }

  res_all <- bind_rows(results) %>% arrange(city, band)

  if (write_csv && !is.null(out_csv)) {
    readr::write_csv(res_all, out_csv)
    if (verbose) message("Wrote results to: ", out_csv)
  }
  return(res_all)
}

# CLI
args <- commandArgs(trailingOnly = TRUE)
if (!interactive() && length(args) >= 3) {
  nc_path <- args[1]
  pop_tif <- args[2]
  shp_path <- args[3]
  out_csv <- if (length(args) >= 4) args[4] else file.path(dirname(shp_path), "nc_pop_weighted_tas.csv")
  time_mode <- if (length(args) >= 5) args[5] else "avg"
  var_name <- if (length(args) >= 6) args[6] else "tas"
  lon_name <- if (length(args) >= 7) args[7] else "lon"
  lat_name <- if (length(args) >= 8) args[8] else "lat"
  resample_method <- if (length(args) >= 9) args[9] else "near"

  compute_nc_pop_weighted_tas(nc_path, pop_tif, shp_path, out_csv = out_csv,
                              time_mode = time_mode, var_name = var_name,
                              lon_name = lon_name, lat_name = lat_name,
                              resample_method = resample_method, write_csv = TRUE, verbose = TRUE)
} else if (interactive()) {
  message("Function compute_nc_pop_weighted_tas() loaded; call it in R with appropriate arguments.")
} else {
  cat("Usage: Rscript nc_pop_weighted_tas.R hunan.nc pop.tif hunan.shp out.csv [time_mode] [var_name] [lon_var] [lat_var] [resample_method]\n")
}
