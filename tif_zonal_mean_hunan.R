#!/usr/bin/env Rscript

# 读取TIF和SHP，计算湖南各市平均灯光亮度
# 依赖：terra, sf, dplyr, readr

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
})

# 简单的命令行参数解析：
# 用法：Rscript zonal_mean_hunan.R dmsp_like_tif_path hunan_shp_path [out_csv] [band_mode]
# band_mode 可选：
#   数字N   -> 选择第N个band
#   "first" -> 使用第1个band
#   "last"  -> 使用最后一个band（默认）
#   "avg"   -> 对所有band求像元平均，得到单层
#   "all"   -> 保留所有band，分别统计每band的区域均值，并给出行平均(mean_light)
# 示例：
# E:\working\R\R-4.5.0\bin\Rscript.exe zonal_mean_hunan.r  DMSP-like2020.tif hunan.shp output.csv avg
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("用法: Rscript zonal_mean_hunan.R dmsp_like_tif_path hunan_shp_path [out_csv]\n")
  cat("示例: Rscript zonal_mean_hunan.R dmsp-like2020.tif hunan.shp hunan_light_mean.csv\n")
  quit(status = 1)
}

tif_path <- args[1]
shp_path <- args[2]
out_csv <- ifelse(length(args) >= 3, args[3],
                  file.path(dirname(shp_path), "hunan_light_mean.csv"))
band_mode <- ifelse(length(args) >= 4, args[4], "last")

# 读取栅格与矢量
if (!file.exists(tif_path)) stop(paste("找不到TIF:", tif_path))
if (!file.exists(shp_path)) stop(paste("找不到SHP:", shp_path))

ras <- terra::rast(tif_path)
vec <- sf::st_read(shp_path, quiet = TRUE)

# 打印TIF元信息
print_raster_meta <- function(r, path_hint = NULL) {
  fn <- tryCatch(terra::sources(r), error = function(e) NULL)
  if (is.null(fn) || length(fn) == 0 || all(is.na(fn))) {
    fn <- if (!is.null(path_hint) && nzchar(path_hint)) path_hint else "(内存栅格)"
  } else {
    fn <- paste(fn[!is.na(fn)], collapse = ", ")
  }
  lyr_names <- names(r)
  sz <- c(nrow = terra::nrow(r), ncol = terra::ncol(r), nlyr = terra::nlyr(r))
  rs <- terra::res(r)
  ex <- terra::ext(r)
  crs_txt <- terra::crs(r, proj = TRUE)
  dtypes <- terra::datatype(r)
  naflags <- tryCatch(terra::NAflag(r), error = function(e) rep(NA, terra::nlyr(r)))
  rng <- tryCatch(terra::global(r, c("min","max"), na.rm = TRUE), error = function(e) NULL)

  cat("\n==== TIF元信息 ====\n")
  cat("文件:", fn, "\n")
  cat("层数:", sz["nlyr"], "\n")
  cat("层名:", paste(lyr_names, collapse = ", "), "\n")
  cat("尺寸:", sz["nrow"], "行 x", sz["ncol"], "列\n")
  cat("分辨率:", paste(rs, collapse = ", "), "\n")
  cat("范围(extent): xmin=", ex$xmin, ", xmax=", ex$xmax, ", ymin=", ex$ymin, ", ymax=", ex$ymax, "\n", sep = "")
  cat("CRS:", ifelse(is.na(crs_txt), "(缺失)", crs_txt), "\n")
  cat("数据类型:", paste(dtypes, collapse = ", "), "\n")
  cat("NA值标记:", paste(naflags, collapse = ", "), "\n")
  if (!is.null(rng)) {
    rng_df <- as.data.frame(rng)
    names(rng_df) <- c("min","max")
    cat("范围(每层):\n")
    for (i in seq_len(nrow(rng_df))) {
      cat("  ", lyr_names[i], ": min=", rng_df$min[i], ", max=", rng_df$max[i], "\n", sep = "")
    }
  }
  cat("===================\n\n")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
print_raster_meta(ras, path_hint = tif_path)

# 根据 band_mode 选择或聚合band
nl <- terra::nlyr(ras)
if (!is.null(band_mode) && nzchar(band_mode)) {
  if (grepl("^[0-9]+$", band_mode)) {
    bi <- as.integer(band_mode)
    if (bi < 1 || bi > nl) stop(paste0("band索引超出范围: 1..", nl))
    ras <- ras[[bi]]
    message(paste("已选择band索引:", bi))
  } else if (tolower(band_mode) == "first") {
    ras <- ras[[1]]
    message("已选择第1个band")
  } else if (tolower(band_mode) == "last") {
    ras <- ras[[nl]]
    message("已选择最后一个band")
  } else if (tolower(band_mode) == "avg") {
    ras <- terra::mean(ras, na.rm = TRUE)
    message("已对所有band求平均(像元层面)")
  } else if (tolower(band_mode) == "all") {
    message("保留所有band以输出分band统计")
  } else {
    warning(paste("未知band_mode:", band_mode, "，保持原始多层栅格"))
  }
}

# 坐标系处理：若两者坐标系不同，重投影矢量到栅格坐标系
crs_r <- terra::crs(ras, proj = TRUE)
crs_v <- sf::st_crs(vec)
if (!is.na(crs_r) && !is.na(crs_v)) {
  if (crs_v$wkt != crs_r) {
    vec <- sf::st_transform(vec, crs = crs_r)
  }
}

# 自动选择城市名称字段
choose_city_field <- function(df) {
  nms <- names(df)
  # 常见字段名候选（不区分大小写）
  candidates <- c("CITY", "CITY_NAME", "NAME", "NAME_1", "NAME_2",
                  "ADMIN_NAME", "ADM2_NAME", "COUNTY", "PREF", "PREFECTURE",
                  "CITY_CN", "NAME_CH", "NAME_ZH")
  nm_upper <- toupper(nms)
  for (cand in candidates) {
    idx <- which(nm_upper == cand)
    if (length(idx) > 0) return(nms[idx[1]])
  }
  # 若没有匹配，优先选择字符型或因子型字段
  char_cols <- nms[sapply(df, function(x) is.character(x) || is.factor(x))]
  if (length(char_cols) > 0) return(char_cols[1])
  # 否则退回第一个非几何字段
  non_geom <- setdiff(nms, attr(df, "sf_column"))
  return(non_geom[1])
}

city_field <- choose_city_field(vec)
if (is.null(city_field) || is.na(city_field)) {
  stop("未能识别城市名称字段，请在脚本中手动指定。")
}

message(paste("使用城市字段:", city_field))

# 计算区域平均值（zonal mean）
# terra::extract 默认按多边形提取，使用 fun=mean, na.rm=TRUE
vec_spat <- terra::vect(vec)  # 转为SpatVector

ext_res <- terra::extract(ras, vec_spat, fun = mean, na.rm = TRUE)

# 提取结果包含ID和提取值列（栅格层名）。
# 将其与城市名称拼接。
joined <- vec %>%
  st_drop_geometry() %>%
  mutate(.row_id = seq_len(n())) %>%
  left_join(ext_res %>% as.data.frame() %>% mutate(.row_id = ID), by = ".row_id")

# 输出逻辑：
# - 若当前ras为单层：mean_light取该层
# - 若为多层（例如 band_mode=all）：mean_light为各band均值的行平均，并同时输出各band列
current_layers <- names(ras)
if (terra::nlyr(ras) == 1) {
  res <- joined %>%
    select(city = all_of(city_field), mean_light = tidyselect::last_col())
} else {
  band_cols <- intersect(current_layers, names(joined))
  if (length(band_cols) == 0) stop("未找到band列，提取结果异常")
  res <- joined %>%
    mutate(mean_light = rowMeans(dplyr::select(., dplyr::all_of(band_cols)), na.rm = TRUE)) %>%
    select(city = all_of(city_field), mean_light, dplyr::all_of(band_cols))
}

# 清理并排序
res <- res %>%
  mutate(city = as.character(city)) %>%
  arrange(city)

# 写出CSV
readr::write_csv(res, out_csv)

message(paste("已写出结果:", out_csv))

# 打印前几行预览
print(utils::head(res, 10))


######################

if (FALSE){
  suppressPackageStartupMessages({
    library(terra)
    library(sf)
    library(dplyr)
    library(readr)
  })
  setwd("D:/xwechat_files/wxid_vu6eu977ndi121_3458/msg/migrate/File/2025-03/zonalstatistic_data")

  ras <- terra::rast("DMSP-like2020.tif")
  vec <- sf::st_read("hunan.shp", quiet = TRUE)
  terra::sources(ras)
  r <- ras
  lyr_names <- names(r)
  sz <- c(nrow = terra::nrow(r), ncol = terra::ncol(r), nlyr = terra::nlyr(r))
  rs <- terra::res(r)
  ex <- terra::ext(r)
  crs_txt <- terra::crs(r, proj = TRUE)
  dtypes <- terra::datatype(r)
  naflags <- tryCatch(terra::NAflag(r), error = function(e) rep(NA, terra::nlyr(r)))
  rng <- tryCatch(terra::global(r, c("min","max"), na.rm = TRUE), error = function(e) NULL)

  cat("\n==== TIF元信息 ====\n")
  cat("层数:", sz["nlyr"], "\n")
  cat("层名:", paste(lyr_names, collapse = ", "), "\n")
  cat("尺寸:", sz["nrow"], "行 x", sz["ncol"], "列\n")
  cat("分辨率:", paste(rs, collapse = ", "), "\n")
  cat("范围(extent): xmin=", ex$xmin, ", xmax=", ex$xmax, ", ymin=", ex$ymin, ", ymax=", ex$ymax, "\n", sep = "")
  cat("CRS:", ifelse(is.na(crs_txt), "(缺失)", crs_txt), "\n")
  cat("数据类型:", paste(dtypes, collapse = ", "), "\n")
  cat("NA值标记:", paste(naflags, collapse = ", "), "\n")

  # 坐标系处理：若两者坐标系不同，重投影矢量到栅格坐标系
  crs_r <- terra::crs(ras, proj = TRUE)
  crs_v <- sf::st_crs(vec)
  if (!is.na(crs_r) && !is.na(crs_v)) {
    if (crs_v$wkt != crs_r) {
      vec <- sf::st_transform(vec, crs = crs_r)
    }
  }

  # 查看vec的字段
  names(vec)

  city_field <- "ShiName"

  # 计算区域平均值（zonal mean）
  # terra::extract 默认按多边形提取，使用 fun=mean, na.rm=TRUE
  vec_spat <- terra::vect(vec)  # 转为SpatVector

  ext_res <- terra::extract(ras, vec_spat, fun = mean, na.rm = TRUE)
  # 提取结果包含ID和提取值列（栅格层名）。
  # 将其与城市名称拼接。
  joined <- vec %>%
    st_drop_geometry() %>%
    mutate(.row_id = seq_len(n())) %>%
    left_join(ext_res %>% as.data.frame() %>% mutate(.row_id = ID), by = ".row_id") 

}