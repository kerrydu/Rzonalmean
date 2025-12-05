#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
})

#' Compute population-weighted lights for polygons
#'
#' This script resamples the DMSP-like light raster to the population raster grid,
#' then for each polygon in the shapefile computes:
#'  - sum_pop_light = sum(pop * light) (exclude cells where pop or light is NA)
#'  - sum_pop = sum(pop) over cells where light is not NA
#'  - pw_mean = sum_pop_light / sum_pop
#'
#' Usage (CLI):
#' Rscript pop_weighted_light.R pop.tif DMSP-like2020.tif hunan.shp out.csv [resample_method] [band_sel]
#'  - resample_method: 'near' (default) or 'bilinear'
#'  - band_sel: 'all' (default) or a single band index (1-based)
#'

compute_pop_weighted_light <- function(pop_tif, dmsp_tif, shp_path,
                                       out_csv = NULL,
                                       resample_method = c("near", "bilinear"),
                                       band_sel = "all",
                                       write_csv = TRUE,
                                       verbose = TRUE) {
  resample_method <- match.arg(resample_method)

  if (!file.exists(pop_tif)) stop("pop raster not found: ", pop_tif)
  if (!file.exists(dmsp_tif)) stop("dmsp raster not found: ", dmsp_tif)
  if (!file.exists(shp_path)) stop("shapefile not found: ", shp_path)

  pop_r <- terra::rast(pop_tif)
  if (terra::nlyr(pop_r) > 1) {
    warning("pop raster has multiple layers; using first band")
    pop_r <- pop_r[[1]]
  }

  dmsp_r <- terra::rast(dmsp_tif)
  # read shapefile early so we can ensure consistent CRS
  vec <- sf::st_read(shp_path, quiet = !verbose)
  pop_crs <- terra::crs(pop_r, proj = TRUE)
  if (is.na(pop_crs) || !nzchar(pop_crs)) {
    warning("Population raster has no CRS; proceeding but results may be incorrect for reprojection/overlay.")
  } else {
    # attempt to transform shapefile to pop raster CRS; catch errors rather than string-compare
    if (is.na(sf::st_crs(vec))) {
      warning("Shapefile has no CRS; cannot reliably transform to population raster CRS. Proceeding without transforming.")
    } else {
      if (verbose) message("Attempting to transform shapefile to pop raster CRS...")
      vec <- tryCatch(sf::st_transform(vec, pop_crs), error = function(e) {
        warning(sprintf("Transforming shapefile failed: %s. Proceeding without transform.", conditionMessage(e)))
        vec
      })
    }
  }

  # ensure crs for DMSP: if different, reproject DMSP to pop crs
  if (!is.na(pop_crs) && nzchar(pop_crs)) {
    if (terra::crs(dmsp_r, proj = TRUE) != pop_crs) {
      if (verbose) message("Reprojecting DMSP to population CRS...")
      dmsp_r <- terra::project(dmsp_r, pop_crs)
    }
  }

  # Resample DMSP to pop grid
  if (verbose) message("Resampling DMSP to pop raster grid using method: ", resample_method)
  dmsp_rs <- terra::resample(dmsp_r, pop_r, method = resample_method)

  # helper: pick city name field
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
    return(non_geom[1])
  }

  city_field <- choose_city_field(vec)
  if (verbose) message("Using city field: ", city_field)

  # determine bands to process
  nb <- terra::nlyr(dmsp_rs)
  bands <- if (identical(band_sel, "all") || is.null(band_sel)) seq_len(nb) else {
    i <- as.integer(band_sel)
    if (is.na(i) || i < 1 || i > nb) stop("Invalid band_sel")
    i
  }

  v_spat <- terra::vect(vec)

  results <- list()
  for (b in bands) {
    if (verbose) message(sprintf("Processing band %d/%d", b, nb))
    light_r <- dmsp_rs[[b]]

    # compute unweighted mean of light per polygon (ignore NA cells)
    mean_light_df <- terra::extract(light_r, v_spat, fun = mean, na.rm = TRUE)
    # masked pop: set pop to NA where light is NA so sum(pop) excludes those cells
    pop_masked <- terra::mask(pop_r, light_r)

    # product raster pop * light (NA if either NA)
    prod_r <- pop_r * light_r

    # sum(pop * light)
    sum_pop_light <- terra::extract(prod_r, v_spat, fun = sum, na.rm = TRUE)
    # sum(pop) over cells where light is not NA
    sum_pop <- terra::extract(pop_masked, v_spat, fun = sum, na.rm = TRUE)

    # sum_* are returned as data.frames with ID and the value
    df_sum_pl <- as.data.frame(sum_pop_light)
    df_sum_p  <- as.data.frame(sum_pop)
    df_mean_light <- as.data.frame(mean_light_df)

    # extract cities
    cities <- vec %>% st_drop_geometry() %>% mutate(.row_id = seq_len(n()))

    df <- data.frame(.row_id = df_mean_light$ID,
                     mean_light = df_mean_light[[2]],
                     sum_pop_light = df_sum_pl[[2]],
                     sum_pop = df_sum_p[[2]],
                     stringsAsFactors = FALSE)
    df <- left_join(cities %>% mutate(.row_id = seq_len(n())), df, by = ".row_id")
    df <- df %>% mutate(band = b,
                        pw_mean = ifelse(is.na(sum_pop) | sum_pop == 0, NA_real_, sum_pop_light / sum_pop)) %>%
      select(all_of(city_field), band, mean_light, sum_pop_light, sum_pop, pw_mean)
    names(df)[1] <- "city"
    results[[length(results) + 1]] <- df
  }

  res_all <- bind_rows(results) %>% arrange(city, band)

  if (write_csv && !is.null(out_csv)) {
    readr::write_csv(res_all, out_csv)
    if (verbose) message("Wrote results to: ", out_csv)
  }

  return(res_all)
}

# CLI entry
args <- commandArgs(trailingOnly = TRUE)
if (!interactive() && length(args) >= 3) {
  pop_tif <- args[1]
  dmsp_tif <- args[2]
  shp_path <- args[3]
  out_csv <- if (length(args) >= 4) args[4] else file.path(dirname(shp_path), "pop_weighted_light.csv")
  resample_method <- if (length(args) >= 5) args[5] else "near"
  band_sel <- if (length(args) >= 6) args[6] else "all"
  compute_pop_weighted_light(pop_tif, dmsp_tif, shp_path, out_csv = out_csv,
                             resample_method = resample_method, band_sel = band_sel,
                             write_csv = TRUE, verbose = TRUE)
} else if (interactive()) {
  message("Function compute_pop_weighted_light() is available for interactive use.")
} else {
  cat("Usage: Rscript pop_weighted_light.R pop.tif DMSP-like2020.tif hunan.shp out.csv [resample_method] [band_sel]\n")
}
