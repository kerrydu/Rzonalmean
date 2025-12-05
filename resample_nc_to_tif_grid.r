# # Default var/dim names
# Rscript ".\resample_nc_to_tif_grid.r" ".\hunan10.nc" ".\Pop.tif"

# # Explicit names and output
# Rscript ".\resample_nc_to_tif_grid.r" ".\hunan10.nc" ".\Pop.tif" ".\hunan10_on_pop.tif" "tas" "lon" "lat"

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript resample_nc_to_tif_grid.r <hunan10.nc> <Pop.tif> [out_tif] [var_name] [lon_name] [lat_name]\n")
  cat("  说明: 将NetCDF按Pop.tif的x/y格网进行最近邻重采样，输出为多层GeoTIFF。\n")
  quit(status = 1)
}

nc_path <- args[1]
tif_path <- args[2]
out_tif <- ifelse(length(args) >= 3 && nzchar(args[3]), args[3], file.path(dirname(nc_path), "hunan10_resampled_to_pop.tif"))
var_name <- ifelse(length(args) >= 4 && nzchar(args[4]), args[4], "tas")
lon_name <- ifelse(length(args) >= 5 && nzchar(args[5]), args[5], "lon")
lat_name <- ifelse(length(args) >= 6 && nzchar(args[6]), args[6], "lat")

if (!file.exists(nc_path)) stop(paste("找不到NC文件:", nc_path))
if (!file.exists(tif_path)) stop(paste("找不到TIFF文件:", tif_path))

# Read target raster grid
r_template <- terra::rast(tif_path)

# Open NC
nc <- ncdf4::nc_open(nc_path)
on.exit(ncdf4::nc_close(nc), add = TRUE)

# Read lon/lat
lon <- tryCatch(ncdf4::ncvar_get(nc, lon_name), error = function(e) NULL)
if (is.null(lon) && !is.null(nc$dim[[lon_name]]$vals)) lon <- nc$dim[[lon_name]]$vals
lat <- tryCatch(ncdf4::ncvar_get(nc, lat_name), error = function(e) NULL)
if (is.null(lat) && !is.null(nc$dim[[lat_name]]$vals)) lat <- nc$dim[[lat_name]]$vals
if (is.null(lon) || is.null(lat)) stop("无法读取lon/lat坐标数组。")
if (length(dim(lon)) > 1 || length(dim(lat)) > 1) stop("仅支持规则1D lon/lat 网格。")

# Read variable array
arr <- ncdf4::ncvar_get(nc, var_name)
var <- nc$var[[var_name]]
dim_names <- sapply(var$dim, function(x) x$name)

# Ensure order lon/lat/time if present
reorder_to <- function(a, from_names, target_names) {
  idx <- match(target_names, from_names)
  idx <- idx[!is.na(idx)]
  if (length(idx) == length(from_names)) {
    aperm(a, idx)
  } else {
    # if no time, just lon/lat
    if (all(c(lon_name, lat_name) %in% from_names)) {
      idx2 <- match(c(lon_name, lat_name), from_names)
      aperm(a, idx2)
    } else {
      stop("变量维度非预期，无法重排。")
    }
  }
}
arr <- reorder_to(arr, dim_names, c(lon_name, lat_name, "time"))

nx <- length(lon); ny <- length(lat)
nt <- if (length(dim(arr)) >= 3) dim(arr)[3] else 1

# Sort lon ascending and flip lat to ascending
lon_order <- order(lon)
lon_sorted <- lon[lon_order]
arr <- arr[lon_order, , , drop = FALSE]
lat_desc <- isTRUE(all(diff(lat) < 0))
if (lat_desc) {
  lat <- rev(lat)
  arr <- arr[, ny:1, , drop = FALSE]
}
lat_sorted <- lat

# Build source raster template
resx <- median(diff(lon_sorted), na.rm = TRUE)
resy <- median(diff(lat_sorted), na.rm = TRUE)
xmin <- min(lon_sorted) - resx/2; xmax <- max(lon_sorted) + resx/2
ymin <- min(lat_sorted) - resy/2; ymax <- max(lat_sorted) + resy/2
src_crs <- "EPSG:4326"  # assume geographic lon/lat
src_tmpl <- terra::rast(ncols = nx, nrows = ny, ext = terra::ext(xmin, xmax, ymin, ymax), crs = src_crs)

# Create SpatRaster stack per time slice
make_layer <- function(mat_lon_lat) {
  mat <- t(mat_lon_lat)
  mat <- mat[nrow(mat):1, , drop = FALSE]
  r <- src_tmpl
  terra::values(r) <- as.vector(mat)
  r
}

layers <- vector("list", nt)
for (ti in seq_len(nt)) {
  layers[[ti]] <- make_layer(arr[,,ti, drop = TRUE])
}
src_stack <- terra::rast(layers)

# Reproject source stack to target CRS if needed
src_wkt <- terra::crs(src_stack, proj = TRUE)
tgt_wkt <- terra::crs(r_template, proj = TRUE)
if (!is.na(tgt_wkt) && src_wkt != tgt_wkt) {
  # First project to target CRS using nearest neighbor
  src_stack <- terra::project(src_stack, tgt_wkt, method = "near")
}

# Resample to match target grid (resolution/extent) using nearest neighbor
resampled <- terra::resample(src_stack, r_template, method = "near")

# Write output as multi-layer GeoTIFF
terra::writeRaster(resampled, out_tif, overwrite = TRUE)
message(paste("已写出重采样TIFF:", out_tif, "层数:", terra::nlyr(resampled)))
