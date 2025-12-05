# # Basic
# Rscript ".\crop_tif_by_shp_bbox.r" ".\pop.tif" ".\hunan.shp"

# # If hunan.shp has no CRS, force one (e.g., EPSG:4326)
# Rscript ".\crop_tif_by_shp_bbox.r" ".\pop.tif" ".\hunan.shp" "EPSG:4326"

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript crop_tif_by_shp_bbox.r <pop.tif> <hunan.shp> [force_vector_crs]\n")
  cat("  说明: 读取pop.tif与hunan.shp，按shp的bbox裁剪并覆盖写回pop.tif\n")
  cat("  force_vector_crs: 当shp无CRS时强制设定，如 'EPSG:4326'\n")
  quit(status = 1)
}

tif_path <- args[1]
shp_path <- args[2]
force_vector_crs <- ifelse(length(args) >= 3 && nzchar(args[3]), args[3], "")

if (!file.exists(tif_path)) stop(paste("找不到TIFF文件:", tif_path))
if (!file.exists(shp_path)) stop(paste("找不到SHP文件:", shp_path))

# Read raster and vector
r <- terra::rast(tif_path)
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

# Crop by bbox
vb <- sf::st_bbox(v)
ext_v <- terra::ext(vb["xmin"], vb["xmax"], vb["ymin"], vb["ymax"])
rc <- tryCatch(terra::crop(r, ext_v), error = function(e) {
  stop(paste("裁剪失败:", e$message))
})

# Overwrite original file (create temp then replace)
out_tmp <- paste0(tif_path, ".tmp.tif")
terra::writeRaster(rc, out_tmp, overwrite = TRUE)
# Replace
file.remove(tif_path)
file.rename(out_tmp, tif_path)
message(paste("已裁剪并覆盖写回:", tif_path))
