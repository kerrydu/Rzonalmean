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
# 用法：Rscript zonal_mean_hunan.R raster_path hunan_shp_path [out_csv] [band_mode] [nc_var] [time_sel] [force_crs] [lon_var] [lat_var]
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
  cat("用法: Rscript zonal_mean_hunan.R raster_path hunan_shp_path [out_csv] [band_mode] [nc_var] [time_sel] [force_crs] [lon_var] [lat_var]\n")
  cat("示例(TIF): Rscript zonal_mean_hunan.R dmsp-like2020.tif hunan.shp hunan_light_mean.csv last\n")
  cat("示例(NC):  Rscript zonal_mean_hunan.R hunan.nc hunan.shp out.csv avg varname 2020-01 EPSG:4326 lon lat\n")
  quit(status = 1)
}

tif_path <- args[1]
shp_path <- args[2]
out_csv <- ifelse(length(args) >= 3, args[3],
                  file.path(dirname(shp_path), "hunan_light_mean.csv"))
band_mode <- ifelse(length(args) >= 4, args[4], "last")
nc_var <- ifelse(length(args) >= 5, args[5], NA)
time_sel <- ifelse(length(args) >= 6, args[6], NA)
force_crs_arg <- ifelse(length(args) >= 7, args[7], NA)
lon_var_arg <- ifelse(length(args) >= 8, args[8], NA)
lat_var_arg <- ifelse(length(args) >= 9, args[9], NA)

# 读取栅格与矢量
if (!file.exists(tif_path)) stop(paste("找不到TIF:", tif_path))
if (!file.exists(shp_path)) stop(paste("找不到SHP:", shp_path))

if (grepl("\\.nc$", tif_path, ignore.case = TRUE)) {
  s <- tryCatch(terra::sds(tif_path), error = function(e) NULL)
  s_names <- tryCatch(names(s), error = function(e) character(0))
  ras <- NULL
  # 若提供了 nc_var，优先匹配子数据集；若无子数据集，再按层名/索引选取
  if (!is.na(nc_var) && nzchar(nc_var)) {
    if (!is.null(s) && length(s) > 0) {
      # 子数据集模式
      if (grepl("^[0-9]+$", nc_var)) {
        bi <- as.integer(nc_var)
        if (bi < 1 || bi > length(s)) stop(paste0("nc子数据集索引超出范围: 1..", length(s)))
        ras <- terra::rast(s[[bi]])
        message(paste("已选择NetCDF变量(按索引):", s_names[bi]))
      } else {
        match_idx <- which(toupper(s_names) == toupper(nc_var))
        if (length(match_idx) == 0) match_idx <- grep(nc_var, s_names, ignore.case = TRUE)
        if (length(match_idx) == 0) stop(paste0("未找到变量 '", nc_var, "'。可选: ", paste(s_names, collapse = ", ")))
        ras <- terra::rast(s[[match_idx[1]]])
        message(paste("已选择NetCDF变量(按名称):", s_names[match_idx[1]]))
      }
    } else {
      # 无子数据集：直接读取并按层选择（常见于按time维展开的多层）
      ras0 <- terra::rast(tif_path)
      if (grepl("^[0-9]+$", nc_var)) {
        bi <- as.integer(nc_var)
        if (bi < 1 || bi > terra::nlyr(ras0)) stop(paste0("nc层索引超出范围: 1..", terra::nlyr(ras0)))
        ras <- ras0[[bi]]
        message(paste("已选择NetCDF第", bi, "层"))
      } else {
        nm <- names(ras0)
        match_idx <- which(toupper(nm) == toupper(nc_var))
        if (length(match_idx) == 0) match_idx <- grep(nc_var, nm, ignore.case = TRUE)
        if (length(match_idx) == 0) stop(paste0("未在层名中找到 '", nc_var, "'。可选: ", paste(nm, collapse = ", ")))
        ras <- ras0[[match_idx[1]]]
        message(paste("已按层名选择:", nm[match_idx[1]]))
      }
    }
  } else {
    # 未提供 nc_var：优先提示子数据集，否则直接读取为多层（按时间）
    if (!is.null(s) && length(s) > 0) {
      message(paste0("检测到NetCDF子数据集: ", paste(s_names, collapse = ", "),
                     "。默认使用第1个: ", s_names[1],
                     "。如需指定，请提供第5个参数(nc_var): 变量名或索引。"))
      ras <- terra::rast(s[[1]])
    } else {
      ras <- terra::rast(tif_path)
      message(paste0("读取NetCDF为多层栅格 (可能按时间展开)，层数:", terra::nlyr(ras)))
    }
  }
} else {
  ras <- terra::rast(tif_path)
}
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
  tvals <- tryCatch(terra::time(r), error = function(e) NULL)

  cat("\n==== TIF元信息 ====\n")
  cat("文件:", fn, "\n")
  cat("层数:", sz["nlyr"], "\n")
  cat("层名:", paste(lyr_names, collapse = ", "), "\n")
  cat("尺寸:", sz["nrow"], "行 x", sz["ncol"], "列\n")
  cat("分辨率:", paste(rs, collapse = ", "), "\n")
  cat("范围(extent): xmin=", ex$xmin, ", xmax=", ex$xmax, ", ymin=", ex$ymin, ", ymax=", ex$ymax, "\n", sep = "")
  cat("CRS:", if (!is.null(crs_txt) && nzchar(crs_txt)) crs_txt else "(缺失)", "\n")
  cat("数据类型:", paste(dtypes, collapse = ", "), "\n")
  cat("NA值标记:", paste(naflags, collapse = ", "), "\n")
  if (!is.null(tvals)) {
    tv <- tryCatch(as.character(tvals), error = function(e) NULL)
    if (!is.null(tv) && length(tv) > 0) {
      show_n <- min(10, length(tv))
      cat("时间维(前", show_n, "): ", paste(tv[seq_len(show_n)], collapse = ", "), if (length(tv) > show_n) " ..." else "", "\n", sep = "")
    }
  }
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

# 若存在时间维，与图层一一对应，则用时间字符串重命名图层，便于输出阅读
tv_all <- tryCatch(terra::time(ras), error = function(e) NULL)
if (!is.null(tv_all)) {
  tv_char <- tryCatch(as.character(tv_all), error = function(e) NULL)
  if (!is.null(tv_char) && length(tv_char) == terra::nlyr(ras)) {
    nm <- make.unique(tv_char)
    names(ras) <- nm
  }
}

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

# 若提供 time_sel（仅对含时间维的多层有意义），则优先按时间字符串子串匹配选择层
if (!is.na(time_sel) && nzchar(time_sel) && terra::nlyr(ras) > 1) {
  tv <- tryCatch(as.character(terra::time(ras)), error = function(e) NULL)
  if (!is.null(tv) && length(tv) == terra::nlyr(ras)) {
    idx <- which(grepl(time_sel, tv, ignore.case = TRUE))
    if (length(idx) == 0) {
      stop(paste0("未在时间维中匹配到 '", time_sel, "'。可选(示例前10): ", paste(head(tv, 10), collapse = ", "), if (length(tv) > 10) " ..." else ""))
    }
    sel <- idx[1]
    ras <- ras[[sel]]
    message(paste0("按时间匹配选择图层: ", tv[sel], " (索引 ", sel, ")"))
  }
}

# 坐标系处理：若两者坐标系不同，重投影矢量到栅格坐标系；若栅格缺失CRS，默认设为 WGS84(EPSG:4326)
crs_r_wkt <- tryCatch(terra::crs(ras, proj = TRUE), error = function(e) "")
crs_v_obj <- sf::st_crs(vec)
r_has_crs <- !is.null(crs_r_wkt) && nzchar(crs_r_wkt)
v_has_crs <- !is.na(crs_v_obj)

if (!r_has_crs) {
  terra::crs(ras) <- "EPSG:4326"
  crs_r_wkt <- tryCatch(terra::crs(ras, proj = TRUE), error = function(e) "")
  r_has_crs <- nzchar(crs_r_wkt)
  message("栅格缺失CRS，已自动设为 EPSG:4326 (WGS84)")
}

if (r_has_crs && v_has_crs) {
  if (crs_v_obj$wkt != crs_r_wkt) {
    vec <- sf::st_transform(vec, crs = crs_r_wkt)
  }
} else if (!v_has_crs) {
  message("警告: 矢量缺少CRS，无法重投影。假定两者已在同一坐标系。")
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

# 额外健壮性检查：
# 1) 检查当前栅格是否存在非NA数值
non_na <- tryCatch(terra::global(ras, 'count', na.rm = TRUE), error = function(e) NULL)
if (!is.null(non_na)) {
  nn <- tryCatch(as.numeric(non_na[1,1]), error = function(e) NA_real_)
  if (is.na(nn) || nn == 0) {
    stop("当前选择的图层无有效像元值(均为NA)。请检查 nc_var/time_sel 或 band_mode (例如 avg)，以及数据的缺失值标记。")
  }
}

# 2) 检查范围是否重叠
rx <- terra::ext(ras)
vb <- sf::st_bbox(vec)
overlap <- !(rx$xmax <= vb["xmin"] || rx$xmin >= vb["xmax"] || rx$ymax <= vb["ymin"] || rx$ymin >= vb["ymax"])
if (!isTRUE(overlap)) {
  # 若是NetCDF，尝试使用lon/lat变量自动重建地理参考
  attempted_fix <- FALSE
  if (grepl("\\.nc$", tif_path, ignore.case = TRUE)) {
    message("栅格与矢量不重叠，尝试用 NetCDF lon/lat 重建地理参考...")
    georef_from_nc <- function(r, nc_path, lon_var = NULL, lat_var = NULL, force_crs = NULL) {
      if (!requireNamespace("ncdf4", quietly = TRUE)) {
        stop("需要包 ncdf4 以从NetCDF读取经纬度变量，请先安装：install.packages('ncdf4')")
      }
      nc <- ncdf4::nc_open(nc_path)
      on.exit(ncdf4::nc_close(nc), add = TRUE)
      vnames <- names(nc$var)
      pick_name <- function(cands, vnames) {
        cands <- cands[!is.na(cands) & nzchar(cands)]
        cands[cands %in% vnames][1]
      }
      lon_name <- pick_name(c(lon_var, "lon", "longitude", "x", "Lon", "Longitude", "X"), vnames)
      lat_name <- pick_name(c(lat_var, "lat", "latitude", "y", "Lat", "Latitude", "Y"), vnames)
      if (is.na(lon_name) || is.na(lat_name)) {
        stop(paste0("未找到 lon/lat 变量。变量列表: ", paste(vnames, collapse = ", "),
                    "。可通过命令行提供第8/9参数分别指定 lon_var/lat_var。"))
      }
      lon <- ncdf4::ncvar_get(nc, lon_name)
      lat <- ncdf4::ncvar_get(nc, lat_name)
      if (length(dim(lon)) > 1 || length(dim(lat)) > 1) {
        stop("检测到二维经纬度数组(曲线/不规则网格)，当前脚本仅支持规则1D网格。")
      }
      nx <- terra::ncol(r); ny <- terra::nrow(r)
      if (length(lon) != nx || length(lat) != ny) {
        stop(paste0("lon/lat 长度与栅格尺寸不匹配: lon=", length(lon), " vs ncol=", nx,
                    "; lat=", length(lat), " vs nrow=", ny))
      }
      mono_lon <- isTRUE(all(diff(lon) > 0)) || isTRUE(all(diff(lon) < 0))
      mono_lat <- isTRUE(all(diff(lat) > 0)) || isTRUE(all(diff(lat) < 0))
      if (!mono_lon) message("警告: lon 非单调，将依据最小/最大与步长设置范围，列顺序可能与空间顺序不一致。")
      if (!mono_lat) message("警告: lat 非单调，将依据最小/最大与步长设置范围，行顺序可能与空间顺序不一致。")
      # 翻转以保证经纬度递增
      if (isTRUE(all(diff(lat) < 0))) {
        r <- terra::flip(r, direction = "vertical")
        lat <- rev(lat)
      }
      if (isTRUE(all(diff(lon) < 0))) {
        r <- terra::flip(r, direction = "horizontal")
        lon <- rev(lon)
      }
      resx <- stats::median(diff(lon), na.rm = TRUE)
      resy <- stats::median(diff(lat), na.rm = TRUE)
      xmin <- min(lon) - resx/2
      xmax <- max(lon) + resx/2
      ymin <- min(lat) - resy/2
      ymax <- max(lat) + resy/2
      terra::ext(r) <- terra::ext(xmin, xmax, ymin, ymax)
      terra::crs(r) <- if (!is.null(force_crs) && nzchar(force_crs)) force_crs else "EPSG:4326"
      message(sprintf("已根据 lon/lat 设置范围: [%.6f, %.6f] x [%.6f, %.6f], CRS=%s",
                      xmin, xmax, ymin, ymax,
                      if (!is.null(force_crs) && nzchar(force_crs)) force_crs else "EPSG:4326"))
      r
    }

    r_try <- try(
      georef_from_nc(
        ras, tif_path,
        lon_var = ifelse(is.na(lon_var_arg), "", lon_var_arg),
        lat_var = ifelse(is.na(lat_var_arg), "", lat_var_arg),
        force_crs = ifelse(is.na(force_crs_arg), "", force_crs_arg)
      ),
      silent = TRUE
    )
    if (!inherits(r_try, "try-error")) {
      ras <- r_try
      print_raster_meta(ras, path_hint = tif_path)
      # 重算CRS对齐
      crs_r_wkt <- tryCatch(terra::crs(ras, proj = TRUE), error = function(e) "")
      r_has_crs <- nzchar(crs_r_wkt)
      if (r_has_crs && !is.na(sf::st_crs(vec))) {
        if (sf::st_crs(vec)$wkt != crs_r_wkt) vec <- sf::st_transform(vec, crs = crs_r_wkt)
      }
      rx <- terra::ext(ras)
      vb <- sf::st_bbox(vec)
      overlap <- !(rx$xmax <= vb["xmin"] || rx$xmin >= vb["xmax"] || rx$ymax <= vb["ymin"] || rx$ymin >= vb["ymax"])
      attempted_fix <- TRUE
    } else {
      message("lon/lat 重建失败: ", as.character(r_try))
    }
  }

  if (!isTRUE(overlap)) {
    msg <- paste0(
      "栅格范围与矢量范围不重叠。请确认CRS与地理参考是否正确。\n",
      "栅格范围: xmin=", rx$xmin, ", xmax=", rx$xmax, ", ymin=", rx$ymin, ", ymax=", rx$ymax, "\n",
      "矢量范围: xmin=", vb["xmin"], ", xmax=", vb["xmax"], ", ymin=", vb["ymin"], ", ymax=", vb["ymax"], "\n",
      if (grepl("\\.nc$", tif_path, ignore.case = TRUE))
        paste0("提示: 该NC可能为规则网格但缺少空间参考。可尝试提供第7-9参数 force_crs/lon_var/lat_var (例如: EPSG:4326 lon lat)。\n",
               if (!attempted_fix) "本次已自动尝试重建失败，详见上方日志。\n" else "")
      else
        "提示: 请核对数据CRS，必要时重投影到一致坐标系后再试。",
      ""
    )
    stop(msg)
  }
}

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
