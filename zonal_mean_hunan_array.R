#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
  library(sf)
  library(dplyr)
  library(readr)
  library(terra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript zonal_mean_hunan_array.R <hunan.nc> <hunan.shp> [out_csv] [time_mode] [time_sel] [var_name] [lon_var] [lat_var] [force_crs]\n")
  cat("  time_mode: avg(默认) | first | last | <数字索引(1-based)>\n")
  cat("  示例: Rscript zonal_mean_hunan_array.R hunan.nc hunan.shp out.csv avg \"\" tas lon lat EPSG:4326\n")
  quit(status = 1)
}

nc_path <- args[1]
shp_path <- args[2]
out_csv <- ifelse(length(args) >= 3, args[3], file.path(dirname(shp_path), "hunan_light_mean_array.csv"))
time_mode <- ifelse(length(args) >= 4, args[4], "avg")
time_sel  <- ifelse(length(args) >= 5, args[5], "")
var_name  <- ifelse(length(args) >= 6 && nzchar(args[6]), args[6], "tas")
lon_name  <- ifelse(length(args) >= 7 && nzchar(args[7]), args[7], "")
lat_name  <- ifelse(length(args) >= 8 && nzchar(args[8]), args[8], "")
force_crs <- ifelse(length(args) >= 9 && nzchar(args[9]), args[9], "EPSG:4326")

if (!file.exists(nc_path)) stop(paste("找不到NC文件:", nc_path))
if (!file.exists(shp_path)) stop(paste("找不到SHP:", shp_path))

nc <- ncdf4::nc_open(nc_path)
on.exit(ncdf4::nc_close(nc), add = TRUE)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# 识别lon/lat/time名称
is_lon_name <- function(n) grepl("^(lon|longitude|x)$", n, ignore.case = TRUE)
is_lat_name <- function(n) grepl("^(lat|latitude|y)$", n, ignore.case = TRUE)

pick_coord <- function(provided, fallback_fun) {
  if (!is.null(provided) && nzchar(provided)) return(provided)
  fallback_fun()
}

lon_name <- pick_coord(lon_name, function(){
  vnames <- names(nc$var); dnames <- names(nc$dim)
  cand <- c(vnames[is_lon_name(vnames)], dnames[is_lon_name(dnames)])
  if (length(cand) == 0) stop("未找到lon坐标名，请作为第7参数提供。")
  cand[1]
})
lat_name <- pick_coord(lat_name, function(){
  vnames <- names(nc$var); dnames <- names(nc$dim)
  cand <- c(vnames[is_lat_name(vnames)], dnames[is_lat_name(dnames)])
  if (length(cand) == 0) stop("未找到lat坐标名，请作为第8参数提供。")
  cand[1]
})

lon <- tryCatch(ncdf4::ncvar_get(nc, lon_name), error = function(e) NULL)
if (is.null(lon) && !is.null(nc$dim[[lon_name]]$vals)) lon <- nc$dim[[lon_name]]$vals
lat <- tryCatch(ncdf4::ncvar_get(nc, lat_name), error = function(e) NULL)
if (is.null(lat) && !is.null(nc$dim[[lat_name]]$vals)) lat <- nc$dim[[lat_name]]$vals
if (is.null(lon) || is.null(lat)) stop("无法读取lon/lat坐标数组。")
if (length(dim(lon)) > 1 || length(dim(lat)) > 1) stop("检测到二维经纬度数组(曲线网格)，本脚本仅支持规则1D网格。")

# 读取time
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
  unit <- tolower(r[2]); origin <- tryCatch(as.POSIXct(r[3], tz = "UTC"), error = function(e) NA)
  if (is.na(origin)) return(rep(NA, length(vals)))
  sec <- switch(unit, seconds=1, minutes=60, hours=3600, day=86400, days=86400, NA)
  if (is.na(sec)) return(rep(NA, length(vals)))
  as.POSIXct(origin + vals * sec, tz = "UTC", origin = "1970-01-01")
}
time_str <- if (!is.null(time_vals)) as.character(decode_time(time_vals, time_units)) else NULL

# 读取数据变量(期望维序: (lon, lat, time) 或 (lon, lat))
var <- nc$var[[var_name]]
if (is.null(var)) stop(paste0("未找到变量 '", var_name, "'。"))
arr <- ncdf4::ncvar_get(nc, var_name)
dim_names <- sapply(var$dim, function(x) x$name)

# 统一到 dim顺序: lon, lat, time
reorder_to <- function(a, from_names, target_names) {
  idx <- match(target_names, from_names)
  idx <- idx[!is.na(idx)]
  if (length(idx) == length(from_names)) {
    aperm(a, idx)
  } else {
    # 若无time维，仅重排前两维
    if (all(c("lon","lat") %in% from_names)) {
      idx2 <- match(c("lon","lat"), from_names)
      aperm(a, idx2)
    } else {
      stop("变量维度非预期，无法重排。")
    }
  }
}
arr <- reorder_to(arr, dim_names, c("lon","lat","time"))

nx <- length(lon); ny <- length(lat)
if (dim(arr)[1] != nx || dim(arr)[2] != ny) stop("数据与lon/lat尺寸不匹配。")
nt <- if (length(dim(arr)) >= 3) dim(arr)[3] else 1

# 排序lon/lat到递增方向，并相应重排数组
lon_order <- order(lon)
lon_sorted <- lon[lon_order]
arr <- arr[lon_order, , , drop = FALSE]

lat_desc <- isTRUE(all(diff(lat) < 0))
if (lat_desc) {
  lat <- rev(lat)
  arr <- arr[, ny:1, , drop = FALSE]
}
lat_sorted <- lat

resx <- median(diff(lon_sorted), na.rm = TRUE)
resy <- median(diff(lat_sorted), na.rm = TRUE)
xmin <- min(lon_sorted) - resx/2; xmax <- max(lon_sorted) + resx/2
ymin <- min(lat_sorted) - resy/2; ymax <- max(lat_sorted) + resy/2

# 按时间选择
select_time_index <- function(mode, t_str, t_len) {
  if (grepl("^[0-9]+$", mode)) {
    i <- as.integer(mode); if (i < 1 || i > t_len) stop("time索引超范围"); return(i)
  }
  m <- tolower(mode)
  if (m == "first") return(1)
  if (m == "last")  return(t_len)
  if (nzchar(time_sel) && !is.null(t_str)) {
    idx <- which(grepl(time_sel, t_str, ignore.case = TRUE))
    if (length(idx) == 0) stop(paste0("未匹配到 time_sel=", time_sel))
    return(idx[1])
  }
  return(NA_integer_)
}

# 创建SpatRaster模板
tmpl <- terra::rast(ncols = nx, nrows = ny, ext = terra::ext(xmin, xmax, ymin, ymax), crs = force_crs)

make_layer <- function(mat_lon_lat) {
  # 输入矩阵维度: [lon, lat]，需转为 [row(top->bottom), col(left->right)]
  mat <- t(mat_lon_lat)            # 现为 [lat, lon]
  mat <- mat[nrow(mat):1, , drop = FALSE]  # 上->下
  r <- tmpl
  terra::values(r) <- as.vector(mat)
  r
}

vec <- sf::st_read(shp_path, quiet = TRUE)

# 生成目标栅格
if (tolower(time_mode) == "avg" && nt > 1) {
  # 对 time 维求均值
  m <- apply(arr, MARGIN = c(1,2), FUN = function(x) mean(x, na.rm = TRUE))
  r <- make_layer(m)
} else if (nt == 1) {
  r <- make_layer(arr[,,1, drop = TRUE])
} else {
  # 选择单个时间层
  ti <- select_time_index(time_mode, time_str, nt)
  if (is.na(ti)) ti <- nt             # 默认最后一层
  r <- make_layer(arr[,,ti, drop = TRUE])
}

# 对齐CRS并提取
if (!is.na(sf::st_crs(vec))) {
  if (sf::st_crs(vec)$wkt != terra::crs(r, proj = TRUE)) vec <- sf::st_transform(vec, terra::crs(r, proj = TRUE))
}

# 猜测城市字段
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
message(paste("使用城市字段:", city_field))

ext_res <- terra::extract(r, terra::vect(vec), fun = mean, na.rm = TRUE)
joined <- vec %>% st_drop_geometry() %>% mutate(.row_id = seq_len(n())) %>%
  left_join(ext_res %>% as.data.frame() %>% mutate(.row_id = ID), by = ".row_id")

res <- joined %>% select(city = all_of(city_field), mean_light = tidyselect::last_col()) %>%
  mutate(city = as.character(city)) %>% arrange(city)

readr::write_csv(res, out_csv)
message(paste("已写出结果:", out_csv))
print(utils::head(res, 10))
