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
	cat("用法: Rscript nc_zonal_hunan.r <hunan.nc> <hunan.shp> [out_csv] [time_mode] [time_sel] [lon_var] [lat_var] [force_crs]\n")
	cat("  time_mode: avg(默认) | first | last | <数字索引(1-based)> | daily(逐日输出)\n")
	quit(status = 1)
}

nc_path <- args[1]
shp_path <- args[2]
out_csv <- ifelse(length(args) >= 3, args[3], file.path(dirname(shp_path), "hunan_tas_mean.csv"))
time_mode <- ifelse(length(args) >= 4, args[4], "avg")
time_sel  <- ifelse(length(args) >= 5, args[5], "")
lon_name  <- ifelse(length(args) >= 6 && nzchar(args[6]), args[6], "lon")
lat_name  <- ifelse(length(args) >= 7 && nzchar(args[7]), args[7], "lat")
force_crs <- ifelse(length(args) >= 8 && nzchar(args[8]), args[8], "EPSG:4326")

if (!file.exists(nc_path)) stop(paste("找不到NC文件:", nc_path))
if (!file.exists(shp_path)) stop(paste("找不到SHP:", shp_path))

nc <- ncdf4::nc_open(nc_path)
on.exit(ncdf4::nc_close(nc), add = TRUE)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# 读取 lon/lat
lon <- tryCatch(ncdf4::ncvar_get(nc, lon_name), error = function(e) NULL)
if (is.null(lon) && !is.null(nc$dim[[lon_name]]$vals)) lon <- nc$dim[[lon_name]]$vals
lat <- tryCatch(ncdf4::ncvar_get(nc, lat_name), error = function(e) NULL)
if (is.null(lat) && !is.null(nc$dim[[lat_name]]$vals)) lat <- nc$dim[[lat_name]]$vals
if (is.null(lon) || is.null(lat)) stop("无法读取lon/lat坐标数组。")
if (length(dim(lon)) > 1 || length(dim(lat)) > 1) stop("检测到二维经纬度数组(曲线网格)，本脚本仅支持规则1D网格。")

# 读取 time
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

# 读取 tas 变量并确保维序 (lat, lon, time)
var_name <- "tas"
var <- nc$var[[var_name]]
if (is.null(var)) stop("未找到变量 'tas'。")
arr <- ncdf4::ncvar_get(nc, var_name)
dim_names <- sapply(var$dim, function(x) x$name)

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
arr <- reorder_to(arr, dim_names, c("lat","lon","time"))

nx <- length(lon); ny <- length(lat)
if (dim(arr)[1] != ny || dim(arr)[2] != nx) stop("数据与lat/lon尺寸不匹配。")
nt <- if (length(dim(arr)) >= 3) dim(arr)[3] else 1

# 对 lon排序到递增方向并重排数组以匹配空间顺序
lon_order <- order(lon)
lon_sorted <- lon[lon_order]
arr <- arr[, lon_order, , drop = FALSE]

## 将纬度统一为北->南(递减)，经度为西->东(递增)
lat_order_desc <- order(lat, decreasing = TRUE)
arr <- arr[lat_order_desc, , , drop = FALSE]
lat_sorted <- lat[lat_order_desc]
## 当前 arr 维度顺序为 [lat, lon, time]

resx <- median(diff(lon_sorted), na.rm = TRUE)
resy <- abs(median(diff(lat_sorted), na.rm = TRUE))
xmin <- min(lon_sorted) - resx/2; xmax <- max(lon_sorted) + resx/2
ymin <- min(lat_sorted) - resy/2; ymax <- max(lat_sorted) + resy/2

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

tmpl <- terra::rast(ncols = nx, nrows = ny, ext = terra::ext(xmin, xmax, ymin, ymax), crs = force_crs)

make_layer <- function(mat_ready) {
	r <- tmpl
	# Terra expects values in row-major (left→right, top→bottom).
	# R matrices flatten column-wise, so transpose before vectorizing.
	terra::values(r) <- as.vector(t(mat_ready))
	r
}

vec <- sf::st_read(shp_path, quiet = TRUE)

if (tolower(time_mode) == "avg" && nt > 1) {
	# 此时 arr 已为 [lat, lon, time] 且行北->南、列西->东
	m <- apply(arr, MARGIN = c(1,2), FUN = function(x) mean(x, na.rm = TRUE))
	r <- make_layer(m)
} else if (nt == 1) {
	m <- arr[,,1, drop = TRUE]
	r <- make_layer(m)
} else {
	ti <- select_time_index(time_mode, time_str, nt)
	if (is.na(ti)) ti <- nt
	m <- arr[,,ti, drop = TRUE]
	r <- make_layer(m)
}

if (!is.na(sf::st_crs(vec))) {
	if (sf::st_crs(vec)$wkt != terra::crs(r, proj = TRUE)) vec <- sf::st_transform(vec, terra::crs(r, proj = TRUE))
}

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

if (tolower(time_mode) == "daily" && nt > 1) {
	# 逐日输出长表: city, date, tas_mean
	v_spat <- terra::vect(vec)
	out_df <- list()
	for (ti in seq_len(nt)) {
		m <- arr[,,ti, drop = TRUE]
		rti <- make_layer(m)
		ers <- terra::extract(rti, v_spat, fun = mean, na.rm = TRUE)
		dti <- data.frame(city = as.character(vec[[city_field]]),
											date = if (!is.null(time_str)) time_str[ti] else ti,
											tas_mean = ers[[ncol(ers)]])
		out_df[[ti]] <- dti
	}
	res_all <- dplyr::bind_rows(out_df) %>% arrange(city, date)
	readr::write_csv(res_all, out_csv)
	message(paste("已写出逐日结果:", out_csv))
	print(utils::head(res_all, 10))
} else {
	ext_res <- terra::extract(r, terra::vect(vec), fun = mean, na.rm = TRUE)
	joined <- vec %>% st_drop_geometry() %>% mutate(.row_id = seq_len(n())) %>%
		left_join(ext_res %>% as.data.frame() %>% mutate(.row_id = ID), by = ".row_id")

	res <- joined %>% select(city = all_of(city_field), mean_light = tidyselect::last_col()) %>%
		mutate(city = as.character(city)) %>% arrange(city)

	readr::write_csv(res, out_csv)
	message(paste("已写出结果:", out_csv))
	print(utils::head(res, 10))
}

