#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript plot_tif_hunan_maps.r <dmsp_like_2020.tif> <hunan.shp> [out_prefix] [city_field] [band_index] [force_vector_crs]\n")
  cat("  out_prefix: 输出文件前缀(默认与tif同目录: maps)\n")
  cat("  city_field: 可选，城市名称字段(默认自动检测)\n")
  cat("  band_index: 可选，选择多波段tif的第几波段(默认1)\n")
  cat("  force_vector_crs: 可选，当shp无CRS时强制设定，如 'EPSG:4326'\n")
  quit(status = 1)
}

tif_path <- args[1]
shp_path <- args[2]
out_prefix <- ifelse(length(args) >= 3 && nzchar(args[3]), args[3], file.path(dirname(tif_path), "maps"))
city_field_arg <- ifelse(length(args) >= 4 && nzchar(args[4]), args[4], "")
band_index <- ifelse(length(args) >= 5 && nzchar(args[5]), as.integer(args[5]), 1L)
force_vector_crs <- ifelse(length(args) >= 6 && nzchar(args[6]), args[6], "")

if (!file.exists(tif_path)) stop(paste("找不到TIFF文件:", tif_path))
if (!file.exists(shp_path)) stop(paste("找不到SHP文件:", shp_path))

# IO helpers
.dir.create.silent <- function(d) { if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE) }
.dir.create.silent(dirname(out_prefix))

# Read data
r_full <- terra::rast(tif_path)
# Select band if multilayer
if (terra::nlyr(r_full) < band_index || band_index < 1) stop("band_index 超出范围")
r <- r_full[[band_index]]

v <- sf::st_read(shp_path, quiet = TRUE)

# Force vector CRS if requested
if (nzchar(force_vector_crs) && is.na(sf::st_crs(v))) {
  sf::st_crs(v) <- force_vector_crs
  message(paste("已为矢量设置CRS:", force_vector_crs))
}

# Reproject vector to raster CRS if needed
r_crs_wkt <- terra::crs(r, proj = TRUE)
if (!is.na(sf::st_crs(v))) {
  if (sf::st_crs(v)$wkt != r_crs_wkt) {
    v <- sf::st_transform(v, r_crs_wkt)
    message("已将矢量重投影至栅格CRS")
  }
} else {
  warning("矢量无CRS，且未提供 force_vector_crs，可能导致错位")
}

# Overlap check
ext_r <- terra::ext(r)
bbox_v <- sf::st_bbox(v)
has_overlap <- !(bbox_v["xmin"] > ext_r$xmax || bbox_v["xmax"] < ext_r$xmin ||
                 bbox_v["ymin"] > ext_r$ymax || bbox_v["ymax"] < ext_r$ymin)
if (!has_overlap) {
  stop(sprintf("矢量与栅格范围不相交。\nRaster ext: [%.6f, %.6f, %.6f, %.6f]\nVector bbox: [%.6f, %.6f, %.6f, %.6f]",
               ext_r$xmin, ext_r$xmax, ext_r$ymin, ext_r$ymax,
               bbox_v[["xmin"]], bbox_v[["xmax"]], bbox_v[["ymin"]], bbox_v[["ymax"]]))
}

# Pick city field
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

# Map A: crop raster to Hunan bbox, then heatmap + boundaries
# Build extent from vector bbox (already in raster CRS)
vb <- sf::st_bbox(v)
ext_v <- terra::ext(vb["xmin"], vb["xmax"], vb["ymin"], vb["ymax"])
rc <- tryCatch(terra::crop(r, ext_v), error = function(e) {
  warning(paste("裁剪失败，使用原始栅格:", e$message)); r
})
rdf <- terra::as.data.frame(rc, xy = TRUE, cells = FALSE)
val_col <- names(rdf)[3]

pA <- ggplot() +
  geom_raster(data = rdf, aes(x = x, y = y, fill = .data[[val_col]])) +
  scale_fill_viridis_c(option = "C", na.value = NA, name = val_col) +
  geom_sf(data = v, fill = NA, color = "white", linewidth = 0.3) +
  coord_sf(expand = FALSE, datum = sf::st_crs(v)) +
  labs(title = "图A: 栅格热力图 + Hunan边界", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

out_map_a <- paste0(out_prefix, "_mapA_overlay.png")
ggsave(out_map_a, pA, width = 9, height = 7, dpi = 300)
message(paste("已输出:", out_map_a))

# Zonal mean by city
v_spat <- terra::vect(v)
ex <- terra::extract(r, v_spat, fun = mean, na.rm = TRUE)
df_left <- v |> sf::st_drop_geometry()
df_left$.row_id <- seq_len(nrow(df_left))
ex_df <- ex |> as.data.frame()
ex_df$.row_id <- ex_df$ID
joined <- dplyr::left_join(df_left, ex_df, by = ".row_id")

val_name <- names(joined)[ncol(joined)]
# 保持与几何顺序一致的结果用于绘图
res_map <- joined |>
  dplyr::transmute(.row_id = .data$.row_id,
                   city = as.character(.data[[city_field]]),
                   mean_value = .data[[val_name]])
# 排序后的结果用于CSV更易读
res_csv <- res_map |>
  dplyr::select(city, mean_value) |>
  dplyr::arrange(city)

# Save CSV
out_csv <- paste0(out_prefix, "_city_mean.csv")
readr::write_csv(res_csv, out_csv)
message(paste("已输出城市均值CSV:", out_csv))

# Map B: choropleth by mean
v$.row_id <- seq_len(nrow(v))
v_plot <- dplyr::left_join(v, res_map[, c(".row_id", "mean_value")], by = ".row_id")

pB <- ggplot(v_plot) +
  geom_sf(aes(fill = mean_value), color = "grey30", linewidth = 0.2) +
  scale_fill_viridis_c(option = "C", na.value = "#f0f0f0", name = "Mean") +
  labs(title = "图B: 各城市均值(填色图)", fill = "Mean", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

out_map_b <- paste0(out_prefix, "_mapB_city_mean.png")
ggsave(out_map_b, pB, width = 9, height = 7, dpi = 300)
message(paste("已输出:", out_map_b))

# Quick preview head
print(utils::head(res_csv, 10))
