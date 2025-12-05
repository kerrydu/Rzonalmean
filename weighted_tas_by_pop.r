#Rscript ".\weighted_tas_by_pop.r" ".\hunan10_resampled_to_pop.tif" ".\pop.tif" ".\hunan.shp" ".\weighted_per_time.csv" "" "" "per_time"

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("用法: Rscript weighted_tas_by_pop.r <tas_stack.tif> <pop.tif> <hunan.shp> [out_csv] [city_field] [force_vector_crs] [mode]\n")
  cat("  说明: 读取重采样后的温度栅格(多层)与人口密度栅格，按城市计算加权均值: sum(tas*pop)/sum(pop)。\n")
  cat("  mode: overall(默认，跨时间汇总) | per_time(逐层输出)\n")
  quit(status = 1)
}

tas_path <- args[1]
pop_path <- args[2]
shp_path <- args[3]
out_csv <- ifelse(length(args) >= 4 && nzchar(args[4]), args[4], file.path(dirname(tas_path), "weighted_tas_by_pop.csv"))
city_field_arg <- ifelse(length(args) >= 5 && nzchar(args[5]), args[5], "")
force_vector_crs <- ifelse(length(args) >= 6 && nzchar(args[6]), args[6], "")
mode <- ifelse(length(args) >= 7 && nzchar(args[7]), args[7], "overall")

if (!file.exists(tas_path)) stop(paste("找不到温度栅格:", tas_path))
if (!file.exists(pop_path)) stop(paste("找不到人口栅格:", pop_path))
if (!file.exists(shp_path)) stop(paste("找不到SHP:", shp_path))

# Read rasters and vector
r_tas <- terra::rast(tas_path)
r_pop <- terra::rast(pop_path)
v <- sf::st_read(shp_path, quiet = TRUE)

# Align CRS: force vector CRS if missing
if (nzchar(force_vector_crs) && is.na(sf::st_crs(v))) {
  sf::st_crs(v) <- force_vector_crs
  message(paste("已为矢量设置CRS:", force_vector_crs))
}

# Reproject vector to tas CRS (tas and pop should already align; we ensure match to tas)
r_crs <- terra::crs(r_tas, proj = TRUE)
if (!is.na(sf::st_crs(v))) {
  if (sf::st_crs(v)$wkt != r_crs) {
    v <- sf::st_transform(v, r_crs)
    message("已将矢量重投影至温度栅格CRS")
  }
} else {
  warning("矢量无CRS，且未提供 force_vector_crs，可能导致错位")
}

# Ensure pop matches tas grid
if (terra::crs(r_pop, proj = TRUE) != r_crs) {
  r_pop <- terra::project(r_pop, r_crs, method = "near")
}
if (!terra::compareGeom(r_pop, r_tas, stopOnError = FALSE)) {
  r_pop <- terra::resample(r_pop, r_tas, method = "near")
}

# Choose city field
choose_city_field <- function(df) {
  nms <- names(df)
  candidates <- c("ShiName", "CITY", "CITY_NAME", "NAME", "NAME_1", "NAME_2",
                  "ADMIN_NAME", "ADM2_NAME", "COUNTY", "PREF", "PREFECTURE",
                  "CITY_CN", "NAME_CH", "NAME_ZH")
  if (nzchar(city_field_arg) && city_field_arg %in% nms) return(city_field_arg)
  nm_upper <- toupper(nms)
  for (cand in candidates) {
    idx <- which(nm_upper == toupper(cand))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  char_cols <- nms[sapply(df, function(x) is.character(x) || is.factor(x))]
  if (length(char_cols) > 0) return(char_cols[1])
  non_geom <- setdiff(nms, attr(df, "sf_column"))
  non_geom[1]
}
city_field <- choose_city_field(v)
message(paste("使用城市字段:", city_field))

# Convert vector to SpatVector
v_spat <- terra::vect(v)

# Helper: weighted mean per layer
weighted_extract <- function(r_layer, r_weight, v_spat) {
  # Mask weights by where r_layer is non-NA so denominator excludes TAS NAs
  w_mask <- terra::mask(r_weight, r_layer)
  num <- terra::extract(r_layer * w_mask, v_spat, fun = sum, na.rm = TRUE)
  den <- terra::extract(w_mask, v_spat, fun = sum, na.rm = TRUE)
  out <- data.frame(ID = num$ID,
                    num = num[[ncol(num)]],
                    den = den[[ncol(den)]])
  out$wmean <- ifelse(out$den > 0, out$num / out$den, NA_real_)
  out
}

nl <- terra::nlyr(r_tas)
if (nl == 1 && mode != "per_time") {
  # Single layer simple weighted mean
  res_df <- weighted_extract(r_tas, r_pop, v_spat)
  joined <- v |> sf::st_drop_geometry() |> mutate(.row_id = seq_len(n())) |>
    left_join(res_df |> mutate(.row_id = ID), by = ".row_id")
  out_df <- joined |>
    transmute(city = as.character(.data[[city_field]]), weighted_mean = .data$wmean) |>
    arrange(city)
  readr::write_csv(out_df, out_csv)
  message(paste("已输出加权均值CSV:", out_csv))
  print(utils::head(out_df, 10))
} else if (mode == "per_time") {
  # Per-layer weighted mean
  out_list <- vector("list", nl)
  for (i in seq_len(nl)) {
    r_layer <- r_tas[[i]]
    res_df <- weighted_extract(r_layer, r_pop, v_spat)
    tmp <- v |> sf::st_drop_geometry() |> mutate(.row_id = seq_len(n())) |>
      left_join(res_df |> mutate(.row_id = ID), by = ".row_id") |>
      transmute(city = as.character(.data[[city_field]]), layer = i, weighted_mean = .data$wmean)
    out_list[[i]] <- tmp
  }
  out_all <- dplyr::bind_rows(out_list) |> arrange(city, layer)
  readr::write_csv(out_all, out_csv)
  message(paste("已输出逐层加权均值CSV:", out_csv))
  print(utils::head(out_all, 10))
} else {
  # overall across time: sum numerators and denominators over time with NA-safe masking per layer
  n_poly <- nrow(v)
  acc_num <- numeric(n_poly)
  acc_den <- numeric(n_poly)
  for (i in seq_len(nl)) {
    res_df <- weighted_extract(r_tas[[i]], r_pop, v_spat)
    acc_num[res_df$ID] <- acc_num[res_df$ID] + ifelse(is.finite(res_df$num), res_df$num, 0)
    acc_den[res_df$ID] <- acc_den[res_df$ID] + ifelse(is.finite(res_df$den), res_df$den, 0)
  }
  wmean_overall <- ifelse(acc_den > 0, acc_num / acc_den, NA_real_)
  out_df <- v |> sf::st_drop_geometry() |>
    mutate(weighted_mean = wmean_overall) |>
    transmute(city = as.character(.data[[city_field]]), weighted_mean) |>
    arrange(city)
  readr::write_csv(out_df, out_csv)
  message(paste("已输出跨时间加权均值CSV:", out_csv))
  print(utils::head(out_df, 10))
}
