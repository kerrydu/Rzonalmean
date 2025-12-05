#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("用法: Rscript nc_step_by_step.r <file.nc>\n")
  quit(status = 1)
}

nc_path <- args[1]
if (!file.exists(nc_path)) stop(paste("找不到NC文件:", nc_path))

nc <- ncdf4::nc_open(nc_path)
on.exit(ncdf4::nc_close(nc), add = TRUE)

cat("\n== 基本信息 ==\n")
cat("文件:", nc$filename, "\n")
cat("格式:", nc$format, "\n")
cat("维度数量:", length(nc$dim), " 变量数量:", length(nc$var), "\n\n")

name_vec <- function(x) paste(x, collapse = ", ")

get_units <- function(varname) {
  at <- tryCatch(ncdf4::ncatt_get(nc, varname), error = function(e) list())
  u <- at$units
  if (is.null(u) || length(u) == 0) "" else as.character(u)
}

cat("-- 维度信息 --\n")
for (dn in names(nc$dim)) {
  d <- nc$dim[[dn]]
  # 尝试从同名变量或维度属性中获取单位
  units <- ""
  if (!is.null(nc$var[[dn]])) units <- get_units(dn)
  if (units == "" && !is.null(d$units)) units <- as.character(d$units)
  cat("  ", dn, ": len=", d$len, ", units=", units, "\n", sep = "")
}
cat("\n")

cat("-- 变量信息 --\n")
for (vn in names(nc$var)) {
  v <- nc$var[[vn]]
  dims <- v$dim
  dim_names <- sapply(dims, function(x) x$name)
  sizes <- paste(v$size, collapse = ",")
  at <- tryCatch(ncdf4::ncatt_get(nc, vn), error = function(e) list())
  units <- at$units %||% "-"
  longn <- at$long_name %||% at$long %||% "-"
  stdn  <- at$standard_name %||% "-"
  gridm <- at$grid_mapping %||% "-"
  fillv <- at$`_FillValue` %||% at$missing_value %||% "-"
  scale <- at$scale_factor %||% "-"
  offset<- at$add_offset %||% "-"
  cat("  ", vn, ": dims=(", name_vec(dim_names), ") size=(", sizes, ") units=", units,
      " long=", longn, " std=", stdn, " grid_mapping=", gridm,
      " _FillValue=", fillv, " scale=", scale, " offset=", offset, "\n", sep = "")
}
cat("\n")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

is_lon_name <- function(n) grepl("^(lon|longitude|x)$", n, ignore.case = TRUE)
is_lat_name <- function(n) grepl("^(lat|latitude|y)$", n, ignore.case = TRUE)

get_coord_vals <- function(name) {
  # 优先尝试作为变量读取；否则看维度的 vals
  if (!is.null(nc$var[[name]])) {
    vals <- tryCatch(ncdf4::ncvar_get(nc, name), error = function(e) NULL)
    return(vals)
  }
  if (!is.null(nc$dim[[name]]) && !is.null(nc$dim[[name]]$vals)) {
    return(nc$dim[[name]]$vals)
  }
  return(NULL)
}

summarize_axis <- function(vals) {
  if (is.null(vals)) return(list(summary = "(无坐标值)", monotonic = NA))
  if (length(dim(vals)) > 1) return(list(summary = paste0("2D数组 ", paste(dim(vals), collapse = "x")), monotonic = NA))
  m <- suppressWarnings(stats::median(diff(vals), na.rm = TRUE))
  mono <- isTRUE(all(diff(vals) > 0)) || isTRUE(all(diff(vals) < 0))
  headv <- paste(head(vals, 5), collapse = ", ")
  tailv <- paste(tail(vals, 3), collapse = ", ")
  sm <- paste0("len=", length(vals), ", min=", min(vals, na.rm = TRUE), ", max=", max(vals, na.rm = TRUE),
               ", step≈", round(m, 6), ", head=", headv, ", tail=", tailv)
  list(summary = sm, monotonic = mono)
}

cat("-- 坐标变量 --\n")
lon_name <- NULL; lat_name <- NULL

# 1) 首先在变量中找标准名/典型名
for (vn in names(nc$var)) {
  at <- tryCatch(ncdf4::ncatt_get(nc, vn), error = function(e) list())
  stdn <- tolower(as.character(at$standard_name %||% ""))
  if (stdn %in% c("longitude")) lon_name <- lon_name %||% vn
  if (stdn %in% c("latitude"))  lat_name <- lat_name %||% vn
  if (is.null(lon_name) && is_lon_name(vn)) lon_name <- vn
  if (is.null(lat_name) && is_lat_name(vn)) lat_name <- vn
}

# 2) 若变量里没找到，使用同名维度作为坐标轴
if (is.null(lon_name)) {
  dnames <- names(nc$dim)
  cand <- dnames[is_lon_name(dnames)]
  if (length(cand) > 0) lon_name <- cand[1]
}
if (is.null(lat_name)) {
  dnames <- names(nc$dim)
  cand <- dnames[is_lat_name(dnames)]
  if (length(cand) > 0) lat_name <- cand[1]
}

if (is.null(lon_name) && is.null(lat_name)) {
  cat("  未检测到标准 lon/lat 坐标轴 (在变量或维度中)。\n")
} else {
  if (!is.null(lon_name)) {
    vals <- get_coord_vals(lon_name)
    sm <- summarize_axis(vals)
    units <- if (!is.null(nc$var[[lon_name]])) get_units(lon_name) else (nc$dim[[lon_name]]$units %||% "")
    cat("  lon轴:", lon_name, " units=", units, " ", sm$summary,
        if (!is.na(sm$monotonic)) paste0(" monotonic=", sm$monotonic) else "", "\n", sep = "")
  }
  if (!is.null(lat_name)) {
    vals <- get_coord_vals(lat_name)
    sm <- summarize_axis(vals)
    units <- if (!is.null(nc$var[[lat_name]])) get_units(lat_name) else (nc$dim[[lat_name]]$units %||% "")
    cat("  lat轴:", lat_name, " units=", units, " ", sm$summary,
        if (!is.na(sm$monotonic)) paste0(" monotonic=", sm$monotonic) else "", "\n", sep = "")
  }
}
cat("\n")

cat("-- 时间维 --\n")
time_name <- NULL
if (!is.null(nc$var$time)) time_name <- "time"
if (is.null(time_name)) {
  # 在维度中查找名为 time
  if (!is.null(nc$dim$time)) time_name <- "time"
}

decode_time <- function(vals, units) {
  if (is.null(vals) || is.null(units) || !nzchar(units)) return(rep(NA, length(vals)))
  m <- regexec("^(seconds|minutes|hours|days|day) since (.+)$", units, ignore.case = TRUE)
  r <- regmatches(units, m)[[1]]
  if (length(r) < 3) return(rep(NA, length(vals)))
  unit <- tolower(r[2])
  origin <- tryCatch(as.POSIXct(r[3], tz = "UTC"), error = function(e) NA)
  if (is.na(origin)) return(rep(NA, length(vals)))
  sec <- switch(unit,
                seconds = 1,
                minutes = 60,
                hours   = 3600,
                day     = 86400,
                days    = 86400,
                NA)
  if (is.na(sec)) return(rep(NA, length(vals)))
  as.POSIXct(origin + vals * sec, origin = "1970-01-01", tz = "UTC")
}

if (is.null(time_name)) {
  cat("  未检测到 time 变量/维度。\n")
} else {
  tlen <- if (!is.null(nc$var[[time_name]])) nc$var[[time_name]]$len else nc$dim[[time_name]]$len
  tunits <- if (!is.null(nc$var[[time_name]])) get_units(time_name) else (nc$dim[[time_name]]$units %||% "")
  tvals <- tryCatch(ncdf4::ncvar_get(nc, time_name), error = function(e) NULL)
  if (is.null(tvals) && !is.null(nc$dim[[time_name]]$vals)) tvals <- nc$dim[[time_name]]$vals
  cat("  dim ", time_name, ": len=", tlen, ", units=", tunits, "\n", sep = "")
  if (!is.null(tvals)) {
    decoded <- decode_time(tvals, tunits)
    show_n <- min(10, length(tvals))
    cat("  time样本: raw=", paste(head(tvals, show_n), collapse = ", "),
        "; decoded=", paste(utils::head(as.character(decoded), show_n), collapse = ", "), "\n", sep = "")
  }
}
cat("\n")

cat("-- 可能的数据变量(依赖 lon/lat 维) --\n")
ll_dims <- c(lon_name, lat_name)
ll_dims <- ll_dims[!is.null(ll_dims)]
found <- FALSE
if (length(ll_dims) >= 2) {
  for (vn in names(nc$var)) {
    v <- nc$var[[vn]]
    dnames <- sapply(v$dim, function(x) x$name)
    if (all(c(lon_name, lat_name) %in% dnames)) {
      cat("  ", vn, ": dims=(", name_vec(dnames), ") size=(", paste(v$size, collapse = ","), ")\n", sep = "")
      found <- TRUE
    }
  }
}
if (!found) {
  cat("  (未发现或未能识别。注意: 若仅有维度无同名变量, 仍可基于维度长度构建规则网格。)\n")
}
cat("\n完成。\n")
