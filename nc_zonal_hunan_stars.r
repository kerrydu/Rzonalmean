#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stars)
  library(exactextractr)
  library(sf)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript nc_zonal_hunan_stars.r <hunan.nc> <hunan.shp> [out_csv] [time_mode] [time_sel] [var_name] [city_field] [force_vector_crs]\n")
  cat("  time_mode: avg(默认) | first | last | <数字索引(1-based)> | daily(逐日输出)\n")
  quit(status = 1)
}

nc_path <- args[1]
shp_path <- args[2]
out_csv <- ifelse(length(args) >= 3, args[3], file.path(dirname(shp_path), "hunan_tas_mean_stars.csv"))
time_mode <- ifelse(length(args) >= 4, args[4], "avg")
time_sel  <- ifelse(length(args) >= 5, args[5], "")
var_name  <- ifelse(length(args) >= 6 && nzchar(args[6]), args[6], "tas")
city_field_arg <- ifelse(length(args) >= 7 && nzchar(args[7]), args[7], "")
force_vector_crs <- ifelse(length(args) >= 8 && nzchar(args[8]), args[8], "")

if (!file.exists(nc_path)) stop(paste("找不到NC文件:", nc_path))
if (!file.exists(shp_path)) stop(paste("找不到SHP:", shp_path))

# 读取 NetCDF 为 stars 对象（支持规则与曲线网格）
x <- stars::read_ncdf(nc_path, var = var_name)
polys <- sf::st_read(shp_path, quiet = TRUE)
# 若NC无投影信息，假定为地理坐标并设为EPSG:4326
if (is.na(sf::st_crs(x))) {
  sf::st_crs(x) <- sf::st_crs(4326)
}

# 若矢量缺少CRS，按需强制设置（默认不强制）；随后进行重投影以匹配栅格
if (is.na(sf::st_crs(polys))) {
  if (nzchar(force_vector_crs)) {
    sf::st_crs(polys) <- sf::st_crs(force_vector_crs)
    message(paste("矢量CRS缺失，已强制设为:", force_vector_crs))
  } else {
    message("警告: 矢量缺少CRS且未指定 force_vector_crs，假定与栅格一致。")
  }
}
if (!is.na(sf::st_crs(x)) && !is.na(sf::st_crs(polys)) && sf::st_crs(x)$wkt != sf::st_crs(polys)$wkt) {
  polys <- sf::st_transform(polys, crs = sf::st_crs(x))
}

# 运行前检查范围重叠
rb <- sf::st_bbox(x)
vb <- sf::st_bbox(polys)
overlap <- !(rb["xmax"] <= vb["xmin"] || rb["xmin"] >= vb["xmax"] || rb["ymax"] <= vb["ymin"] || rb["ymin"] >= vb["ymax"])
if (!isTRUE(overlap)) {
  stop(paste0("栅格范围与矢量范围不重叠。\n",
              "栅格bbox: xmin=", rb["xmin"], ", xmax=", rb["xmax"], ", ymin=", rb["ymin"], ", ymax=", rb["ymax"], "\n",
              "矢量bbox: xmin=", vb["xmin"], ", xmax=", vb["xmax"], ", ymin=", vb["ymin"], ", ymax=", vb["ymax"], "\n"))
}

# stars dims: names(x) gives attribute name; dimensions via st_dimensions(x)
dims <- stars::st_dimensions(x)
if (!"time" %in% names(dims)) {
  # 无时间维，按单层处理
  r1 <- x
} else {
  # 提取时间字符串
  tvals <- tryCatch(as.character(dims$time$values), error = function(e) NULL)
  nt <- length(dims$time$values)
  pick_time_index <- function(mode, t_str, t_len) {
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

  if (tolower(time_mode) == "avg") {
    # 对时间维进行平均，生成单层 stars
    r1 <- stars::st_apply(x, MARGIN = c("lon","lat"), FUN = function(v) mean(v, na.rm = TRUE))
    names(r1) <- var_name
  } else if (tolower(time_mode) == "daily") {
    # 逐日输出
    out_list <- vector("list", nt)
    # 选择城市字段
    choose_city_field <- function(df) {
      nms <- names(df)
      candidates <- c("ShiName","CITY","CITY_NAME","NAME","NAME_1","NAME_2",
                      "ADMIN_NAME","ADM2_NAME","COUNTY","PREF","PREFECTURE",
                      "CITY_CN","NAME_CH","NAME_ZH")
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
    city_field <- if (nzchar(city_field_arg)) city_field_arg else choose_city_field(polys)
    for (ti in seq_len(nt)) {
      xi <- x[,,,ti, drop = TRUE]
      # 确保CRS一致
      if (!is.na(sf::st_crs(xi)) && !is.na(sf::st_crs(polys)) && sf::st_crs(xi)$wkt != sf::st_crs(polys)$wkt) {
        polys <- sf::st_transform(polys, crs = sf::st_crs(xi))
      }
      vals <- exactextractr::exact_extract(xi, polys, fun = "mean")
      if (!(city_field %in% names(polys))) stop(paste0("指定的城市字段不存在: ", city_field))
      if (length(vals) == nrow(polys)) {
        dti <- data.frame(city = as.character(polys[[city_field]]),
                          date = if (!is.null(tvals)) tvals[ti] else ti,
                          tas_mean = vals)
        out_list[[ti]] <- dti
      } else {
        warning(paste0("第 ", ti, " 个时间层提取结果大小异常: ", length(vals), " vs ", nrow(polys), "。跳过该层。"))
      }
    }
    res_all <- dplyr::bind_rows(out_list)
    if (all(c("city","date") %in% names(res_all))) {
      res_all <- dplyr::arrange(res_all, city, date)
    }
    readr::write_csv(res_all, out_csv)
    message(paste("已写出逐日结果:", out_csv))
    print(utils::head(res_all, 10))
    quit(status = 0)
  } else {
    ti <- pick_time_index(time_mode, tvals, nt)
    if (is.na(ti)) ti <- nt
    r1 <- x[,,,ti, drop = TRUE]
  }
}

# 统一CRS：exactextractr支持 stars 对象，若多边形与栅格CRS不同，重投影多边形
crs_r <- sf::st_crs(r1)
crs_v <- sf::st_crs(polys)
if (!is.na(crs_r) && !is.na(crs_v) && crs_r$wkt != crs_v$wkt) {
  polys <- sf::st_transform(polys, crs = crs_r)
}

# 计算区域均值（单层）
if (!is.na(sf::st_crs(r1)) && !is.na(sf::st_crs(polys)) && sf::st_crs(r1)$wkt != sf::st_crs(polys)$wkt) {
  polys <- sf::st_transform(polys, crs = sf::st_crs(r1))
}
vals <- exactextractr::exact_extract(r1, polys, fun = "mean")

choose_city_field <- function(df) {
  nms <- names(df)
  candidates <- c("ShiName","CITY","CITY_NAME","NAME","NAME_1","NAME_2",
                  "ADMIN_NAME","ADM2_NAME","COUNTY","PREF","PREFECTURE",
                  "CITY_CN","NAME_CH","NAME_ZH")
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
city_field <- if (nzchar(city_field_arg)) city_field_arg else choose_city_field(polys)
if (!(city_field %in% names(polys))) stop(paste0("指定的城市字段不存在: ", city_field))

res <- data.frame(city = as.character(polys[[city_field]]), tas_mean = vals) %>%
  arrange(city)

readr::write_csv(res, out_csv)
message(paste("已写出结果:", out_csv))
print(utils::head(res, 10))
