#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(ncdf4)
	library(sf)
	library(dplyr)
	library(readr)
	library(terra)
	library(lubridate)
})


#' Compute zonal statistics from a NetCDF file for Hunan shapefile (interactive API)
#'
#' This function wraps the script logic so you can call it from an interactive R session.
#' It returns a data.frame (either per-city summary or a long table for all times).
#'
#' @param nc_path path to NetCDF file
#' @param shp_path path to Hunan shapefile (folder or .shp)
#' @param out_csv optional path to write CSV (if NULL, will write to working dir only if write_csv=TRUE)
#' @param time_mode one of 'avg','first','last','all' or a numeric index as string
#' @param var_name variable name in NetCDF (default 'tas')
#' @param lon_name name of longitude variable/dimension (default 'lon')
#' @param lat_name name of latitude variable/dimension (default 'lat')
#' @param force_crs coordinate string for output raster (default 'EPSG:4326')
#' @param write_csv whether to write CSV to `out_csv` when result produced
#' @param verbose print progress messages
#' @return data.frame with results (either summary or long table)
nc_zonal_hunan <- function(nc_path, shp_path, out_csv = NULL, time_mode = "avg",
													 var_name = "tas", lon_name = "lon", lat_name = "lat",
													 force_crs = "EPSG:4326", write_csv = TRUE, verbose = TRUE) {
	if (!file.exists(nc_path)) stop(sprintf("找不到NC文件: %s", nc_path))
	if (!file.exists(shp_path)) stop(sprintf("找不到SHP: %s", shp_path))

	nc <- ncdf4::nc_open(nc_path)
	on.exit(ncdf4::nc_close(nc), add = TRUE)

	`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

	lon <- tryCatch(ncdf4::ncvar_get(nc, lon_name), error = function(e) NULL)
	if (is.null(lon) && !is.null(nc$dim[[lon_name]]$vals)) lon <- nc$dim[[lon_name]]$vals
	lat <- tryCatch(ncdf4::ncvar_get(nc, lat_name), error = function(e) NULL)
	if (is.null(lat) && !is.null(nc$dim[[lat_name]]$vals)) lat <- nc$dim[[lat_name]]$vals
	if (is.null(lon) || is.null(lat)) stop("无法读取lon/lat坐标数组。")
	if (length(dim(lon)) > 1 || length(dim(lat)) > 1) stop("检测到二维经纬度数组(曲线网格)，本脚本仅支持规则1D网格。")

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

	var <- nc$var[[var_name]]
	if (is.null(var)) stop(sprintf("未找到变量 '%s'。", var_name))
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
				stop("变量维度非预期，无法重排。")
			}
		}
	}
	arr <- reorder_to(arr, dim_names, c("lat","lon","time"))

	nx <- length(lon); ny <- length(lat)
	if (dim(arr)[1] != ny || dim(arr)[2] != nx) stop("数据与lat/lon尺寸不匹配。")
	nt <- if (length(dim(arr)) >= 3) dim(arr)[3] else 1

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

	select_time_index <- function(mode, t_str, t_len) {
		if (grepl("^[0-9]+$", mode)) {
			i <- as.integer(mode); if (i < 1 || i > t_len) stop("time索引超范围"); return(i)
		}
		m <- tolower(mode)
		if (m == "first") return(1)
		if (m == "last")  return(t_len)
		return(NA_integer_)
	}

	tmpl <- terra::rast(ncols = nx, nrows = ny, ext = terra::ext(xmin, xmax, ymin, ymax), crs = force_crs)
	make_layer <- function(mat_ready) {
		r <- tmpl
		terra::values(r) <- as.vector(t(mat_ready))
		r
	}

	vec <- sf::st_read(shp_path, quiet = !verbose)

	if (tolower(time_mode) == "avg" && nt > 1) {
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
	if (verbose) message(paste("使用城市字段:", city_field))

	if (tolower(time_mode) == "all" && nt > 1) {
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
		if (write_csv && !is.null(out_csv)) readr::write_csv(res_all, out_csv)
		if (verbose && !is.null(out_csv)) message(paste("已写出所有时间层结果:", out_csv))
		if (verbose) print(utils::head(res_all, 10))
		return(res_all)
	} else {
		ext_res <- terra::extract(r, terra::vect(vec), fun = mean, na.rm = TRUE)
		joined <- vec %>% st_drop_geometry() %>% mutate(.row_id = seq_len(n())) %>%
			left_join(ext_res %>% as.data.frame() %>% mutate(.row_id = ID), by = ".row_id")

		res <- joined %>% select(city = all_of(city_field), mean_light = tidyselect::last_col()) %>%
			mutate(city = as.character(city)) %>% arrange(city)

		if (write_csv && !is.null(out_csv)) readr::write_csv(res, out_csv)
		if (verbose && !is.null(out_csv)) message(paste("已写出结果:", out_csv))
		if (verbose) print(utils::head(res, 10))
		return(res)
	}
}


## CLI: keep backward-compatible script invocation
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2 && !interactive()) {
	nc_path <- args[1]
	shp_path <- args[2]
	out_csv <- ifelse(length(args) >= 3, args[3], file.path(dirname(shp_path), "hunan_tas_mean.csv"))
	time_mode <- ifelse(length(args) >= 4, args[4], "avg")
	var_name  <- ifelse(length(args) >= 5 && nzchar(args[5]), args[5], "tas")
	lon_name  <- ifelse(length(args) >= 6 && nzchar(args[6]), args[6], "lon")
	lat_name  <- ifelse(length(args) >= 7 && nzchar(args[7]), args[7], "lat")
	force_crs <- ifelse(length(args) >= 8 && nzchar(args[8]), args[8], "EPSG:4326")

	nc_zonal_hunan(nc_path, shp_path, out_csv = out_csv, time_mode = time_mode,
								 var_name = var_name, lon_name = lon_name, lat_name = lat_name,
								 force_crs = force_crs, write_csv = TRUE, verbose = TRUE)
} else if (interactive()) {
	message("函数 `nc_zonal_hunan()` 已加载；可在交互式会话中调用它。示例: nc_zonal_hunan('hunan.nc','hunan.shp', out_csv=NULL, time_mode='avg')")
} else {
	cat("用法: Rscript nc_zonal_hunan.r <hunan.nc> <hunan.shp> [out_csv] [time_mode] [var_name] [lon_var] [lat_var] [force_crs]\n")
}

#R -e "source('D:/xwechat_files/wxid_vu6eu977ndi121_3458/msg/migrate/File/2025-03/zonalstatistic_data/nc_zonal_hunan.r'); res <- nc_zonal_hunan("D:/xwechat_files/wxid_vu6eu977ndi121_3458/msg/migrate/File/2025-03/zonalstatistic_data/hunan.nc",'D:/xwechat_files/wxid_vu6eu977ndi121_3458/msg/migrate/File/2025-03/zonalstatistic_data/hunan.shp', out_csv = NULL, time_mode='avg'); print(head(res))"

#Rscript "d:/xwechat_files/.../zonalstatistic_data/nc_zonal_hunan.r" "d:/.../hunan.nc" "d:/.../hunan.shp" "d:/.../result_main.csv" all tas lon lat EPSG:4326