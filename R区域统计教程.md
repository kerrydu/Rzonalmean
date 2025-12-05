# R语言栅格数据区域统计教程

## 目录
1. [概述](#概述)
2. [R栅格数据处理包介绍](#r栅格数据处理包介绍)
3. [基本工作流程](#基本工作流程)
4. [任务1：计算Zonal Mean](#任务1计算zonal-mean)
5. [任务2：计算人口加权Zonal Mean](#任务2计算人口加权zonal-mean)
6. [完整示例代码](#完整示例代码)

---

## 概述

本教程介绍如何使用R语言处理栅格数据（TIF/NetCDF格式），并计算基于矢量边界（Shapefile）的区域统计。主要内容包括：
- 基础的区域平均值（Zonal Mean）计算
- 人口加权的区域平均值计算

**适用场景**：
- 气候数据（温度、降水等）的区域统计
- 夜间灯光数据的城市级汇总
- 任何需要将栅格数据按行政区划汇总的场景

---

## R栅格数据处理包介绍

### 1. **terra** 包（推荐）

`terra` 是 `raster` 包的现代化替代品，性能更优，语法更简洁。

**主要功能**：
- 读取和处理TIF、NetCDF等栅格数据
- 栅格运算（裁剪、重采样、代数运算）
- 区域统计（`extract()`、`zonal()`）
- 支持大数据量的内存优化处理

**安装**：
```r
install.packages("terra")
```

**基本用法**：
```r
library(terra)

# 读取栅格
r <- rast("data.tif")           # 单波段或多波段TIF
r_nc <- rast("climate.nc")       # NetCDF文件

# 查看信息
print(r)                         # 打印元数据
nlyr(r)                          # 波段数
res(r)                           # 分辨率
crs(r)                           # 坐标系统
ext(r)                           # 范围

# 基本操作
r_subset <- r[[1:3]]             # 选择前3个波段
r_mean <- mean(r)                # 所有波段求平均
r_crop <- crop(r, extent_obj)    # 裁剪
r_resamp <- resample(r, r_target, method="bilinear")  # 重采样
```

---

### 2. **sf** 包

`sf`（Simple Features）是现代化的矢量数据处理包，替代了`sp`包。

**主要功能**：
- 读取Shapefile、GeoJSON等矢量数据
- 空间操作（投影转换、缓冲区、叠加分析）
- 与dplyr完美集成

**基本用法**：
```r
library(sf)

# 读取矢量
vec <- st_read("regions.shp", quiet = TRUE)

# 投影转换
vec_reproj <- st_transform(vec, crs = st_crs(4326))  # 转为WGS84
vec_reproj <- st_transform(vec, crs = crs(raster))   # 匹配栅格CRS

# 查看和操作
st_crs(vec)                      # 查看坐标系
st_bbox(vec)                     # 查看边界框
plot(st_geometry(vec))           # 绘制几何
```

---

### 3. **stars** 包（处理多维数组）

`stars` 专门处理时空数组数据，特别适合NetCDF多时间步数据。

**主要功能**：
- 原生支持NetCDF的时间维度
- 支持曲线坐标网格（curvilinear grid）
- 与`sf`紧密集成

**基本用法**：
```r
library(stars)

# 读取NetCDF
nc <- read_ncdf("climate.nc", var = "tas")

# 查看维度
st_dimensions(nc)

# 时间维操作
nc_mean <- st_apply(nc, MARGIN = c("lon","lat"), FUN = mean, na.rm = TRUE)

# 与sf结合进行区域统计（需配合exactextractr）
library(exactextractr)
zonal_stats <- exact_extract(nc, polygons, fun = "mean")
```

---

### 4. **ncdf4** 包（底层NetCDF操作）

提供对NetCDF文件的底层访问，适合复杂数据结构的处理。

**基本用法**：
```r
library(ncdf4)

# 打开文件
nc <- nc_open("climate.nc")

# 读取变量
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
data <- ncvar_get(nc, "tas")   # 读取温度数据

# 读取属性
units <- ncatt_get(nc, "time", "units")

# 关闭文件
nc_close(nc)
```

---

### 5. **exactextractr** 包（精确区域统计）

提供比`terra::extract()`更精确的区域统计方法，考虑栅格像元的部分重叠。

**基本用法**：
```r
library(exactextractr)

# 精确提取
result <- exact_extract(raster, polygons, fun = "mean")
result_weighted <- exact_extract(raster, polygons, fun = "mean", weights = pop_raster)
```

---

### 6. 辅助包

- **dplyr**：数据框操作和转换
- **readr**：CSV读写
- **lubridate**：时间日期处理

---

## 基本工作流程

### 标准工作流程

```
1. 数据准备
   ├── 读取栅格数据（TIF/NetCDF）
   ├── 读取矢量边界（Shapefile）
   └── 检查并对齐坐标系统

2. 数据预处理
   ├── 投影转换（统一CRS）
   ├── 空间裁剪（可选）
   ├── 时间维度处理（NetCDF）
   └── 栅格重采样（分辨率对齐）

3. 区域统计
   ├── 简单统计：extract() / zonal()
   └── 加权统计：带权重的区域运算

4. 结果输出
   ├── 整理为数据框
   ├── 保存为CSV
   └── 可视化（可选）
```

### 关键注意事项

1. **坐标系统对齐**：确保栅格和矢量使用相同的CRS
2. **范围检查**：验证栅格和矢量是否有空间重叠
3. **内存管理**：大文件处理时使用分块处理或`terra`的内存优化
4. **NA值处理**：统计时注意处理缺失值（`na.rm = TRUE`）

---

## 任务1：计算Zonal Mean

### 方法A：使用terra处理TIF数据

#### 完整流程

```r
#!/usr/bin/env Rscript

library(terra)
library(sf)
library(dplyr)
library(readr)

# ========== 1. 读取数据 ==========
# 读取栅格（TIF格式）
raster <- rast("DMSP-like2020.tif")

# 读取矢量边界
regions <- st_read("hunan.shp", quiet = TRUE)

# ========== 2. 检查和对齐坐标系统 ==========
# 查看栅格CRS
print(paste("栅格CRS:", crs(raster, proj = TRUE)))

# 查看矢量CRS
print(paste("矢量CRS:", st_crs(regions)$wkt))

# 如果矢量缺少CRS，手动设置（假设为WGS84）
if (is.na(st_crs(regions))) {
  st_crs(regions) <- "EPSG:4326"
  message("矢量CRS缺失，已设置为EPSG:4326")
}

# 将矢量投影转换到栅格CRS
if (st_crs(regions)$wkt != crs(raster, proj = TRUE)) {
  regions <- st_transform(regions, crs = crs(raster, proj = TRUE))
  message("已将矢量重投影至栅格CRS")
}

# ========== 3. 检查空间重叠 ==========
raster_bbox <- as.vector(ext(raster))  # c(xmin, xmax, ymin, ymax)
vector_bbox <- st_bbox(regions)

overlap <- !(raster_bbox[2] < vector_bbox["xmin"] || 
             raster_bbox[1] > vector_bbox["xmax"] ||
             raster_bbox[4] < vector_bbox["ymin"] || 
             raster_bbox[3] > vector_bbox["ymax"])

if (!overlap) {
  stop("栅格与矢量边界不重叠，请检查数据！")
}

# ========== 4. 波段选择（可选）==========
# 如果栅格有多个波段，可以选择：
# - 单个波段：raster <- raster[[1]]
# - 多波段平均：raster <- mean(raster)
# - 保留所有波段：不修改

# 示例：使用最后一个波段
if (nlyr(raster) > 1) {
  raster <- raster[[nlyr(raster)]]
  message(paste("使用最后一个波段"))
}

# ========== 5. 区域统计：提取每个区域的平均值 ==========
# 方法1：使用terra::extract（推荐）
zonal_values <- extract(raster, regions, fun = mean, na.rm = TRUE, ID = TRUE)

# 将结果合并到原矢量属性表
# zonal_values$ID 对应 regions 的行号
regions$mean_value <- zonal_values[, 2]  # 第2列是统计值

# ========== 6. 整理输出数据框 ==========
# 选择城市名称字段（根据实际shapefile调整）
city_field <- "NAME"  # 可能是 "CITY", "ShiName", "NAME_1" 等

result <- data.frame(
  city = regions[[city_field]],
  mean_value = regions$mean_value
)

# ========== 7. 保存结果 ==========
write_csv(result, "output_zonal_mean.csv")
print(result)

message("区域统计完成！结果已保存到 output_zonal_mean.csv")
```

#### 关键函数说明

**`terra::extract()`**：
```r
extract(
  x,               # 栅格对象
  y,               # 矢量对象（sf或SpatVector）
  fun = mean,      # 统计函数：mean, sum, min, max, sd等
  na.rm = TRUE,    # 是否忽略NA值
  ID = TRUE,       # 是否返回区域ID
  weights = NULL   # 可选：权重栅格（用于加权统计）
)
```

---

### 方法B：使用terra处理NetCDF数据

```r
library(terra)
library(sf)
library(dplyr)

# ========== 1. 读取NetCDF ==========
nc_raster <- rast("climate.nc")

# 查看NetCDF信息
print(nc_raster)
names(nc_raster)      # 查看变量名
nlyr(nc_raster)       # 查看层数（时间步数）

# 如果nc文件包含多个变量，指定变量名
# nc_raster <- rast("climate.nc", subds = "tas")  # 读取tas变量

# ========== 2. 处理时间维度 ==========
# 选项1：所有时间步求平均
raster_mean <- mean(nc_raster, na.rm = TRUE)

# 选项2：选择特定时间步
raster_t1 <- nc_raster[[1]]        # 第一个时间步
raster_t_last <- nc_raster[[nlyr(nc_raster)]]  # 最后一个时间步

# 选项3：选择时间范围求平均
raster_subset <- nc_raster[[1:10]]  # 前10个时间步
raster_subset_mean <- mean(raster_subset, na.rm = TRUE)

# ========== 3. 读取矢量并对齐 ==========
regions <- st_read("hunan.shp", quiet = TRUE)

# NC文件通常为EPSG:4326，检查并设置
if (is.na(crs(nc_raster))) {
  crs(nc_raster) <- "EPSG:4326"
}

# 投影转换
if (!is.na(st_crs(regions)) && st_crs(regions)$wkt != crs(nc_raster, proj = TRUE)) {
  regions <- st_transform(regions, crs = crs(nc_raster, proj = TRUE))
}

# ========== 4. 区域统计 ==========
# 使用处理后的单层栅格
zonal_result <- extract(raster_mean, regions, fun = mean, na.rm = TRUE, ID = TRUE)

# 合并结果
regions$tas_mean <- zonal_result[, 2]

# 输出
result <- data.frame(
  city = regions$NAME,
  temperature_mean = regions$tas_mean
)

write_csv(result, "nc_zonal_mean.csv")
```

---

### 方法C：使用stars + exactextractr（NetCDF多时间步）

适合需要逐时间步统计或精确计算的场景。

```r
library(stars)
library(exactextractr)
library(sf)
library(dplyr)

# ========== 1. 读取NetCDF为stars对象 ==========
nc_stars <- read_ncdf("hunan.nc", var = "tas")

# 查看维度
st_dimensions(nc_stars)

# 设置CRS（如果缺失）
if (is.na(st_crs(nc_stars))) {
  st_crs(nc_stars) <- st_crs(4326)
}

# ========== 2. 读取矢量 ==========
regions <- st_read("hunan.shp", quiet = TRUE)

# 投影对齐
if (!is.na(st_crs(regions)) && st_crs(regions)$wkt != st_crs(nc_stars)$wkt) {
  regions <- st_transform(regions, crs = st_crs(nc_stars))
}

# ========== 3. 时间维度处理 ==========
# 选项1：对时间维求平均
nc_mean <- st_apply(nc_stars, MARGIN = c("lon", "lat"), FUN = mean, na.rm = TRUE)
names(nc_mean) <- "tas"

# 选项2：保留时间维，逐时间步统计（见下文）

# ========== 4. 区域统计 ==========
# 使用exactextractr进行精确提取
zonal_result <- exact_extract(nc_mean, regions, fun = "mean")

# 添加到矢量属性
regions$tas_mean <- zonal_result

# 输出
result <- data.frame(
  city = regions$NAME,
  temperature_mean = regions$tas_mean
)

write_csv(result, "stars_zonal_mean.csv")
```

#### 逐时间步统计（daily模式）

```r
# 获取时间维度信息
time_dim <- st_dimensions(nc_stars)$time
n_times <- length(time_dim$values)

# 逐时间步提取
result_list <- vector("list", n_times)

for (i in 1:n_times) {
  # 提取第i个时间步
  nc_slice <- nc_stars[,,,i]  # 根据维度顺序调整
  
  # 区域统计
  zonal <- exact_extract(nc_slice, regions, fun = "mean")
  
  # 保存结果
  result_list[[i]] <- data.frame(
    city = regions$NAME,
    time = as.character(time_dim$values[i]),
    tas_mean = zonal
  )
}

# 合并所有时间步
result_all <- bind_rows(result_list)

write_csv(result_all, "daily_zonal_mean.csv")
```

---

## 任务2：计算人口加权Zonal Mean

### 原理说明

人口加权平均考虑了每个栅格单元的人口分布，公式为：

$$
\text{Weighted Mean} = \frac{\sum_{i=1}^{n} (V_i \times P_i)}{\sum_{i=1}^{n} P_i}
$$

其中：
- $V_i$ = 栅格单元的值（如温度）
- $P_i$ = 栅格单元的人口
- $n$ = 区域内的栅格单元数

**应用场景**：
- 计算城市人口加权平均温度（关注人口密集区）
- 人口暴露度分析（污染、灾害等）
- 任何需要考虑人口分布的区域统计

---

### 完整实现流程

#### 步骤1：准备和对齐数据

```r
library(terra)
library(sf)
library(dplyr)
library(readr)

# ========== 1. 读取数据 ==========
# 温度数据（可以是TIF或NetCDF）
tas_raster <- rast("temperature.tif")  # 或 rast("climate.nc")

# 人口密度栅格
pop_raster <- rast("population.tif")

# 行政区划矢量
regions <- st_read("hunan.shp", quiet = TRUE)

# ========== 2. 检查栅格对齐 ==========
# 检查CRS
if (crs(tas_raster, proj = TRUE) != crs(pop_raster, proj = TRUE)) {
  stop("温度栅格和人口栅格的CRS不一致！")
}

# 检查分辨率和范围
if (!compareGeom(tas_raster, pop_raster, stopOnError = FALSE)) {
  warning("温度栅格和人口栅格的分辨率或范围不一致，将进行重采样")
  
  # 将人口栅格重采样到温度栅格
  pop_raster <- resample(pop_raster, tas_raster, method = "bilinear")
  message("人口栅格已重采样至温度栅格的网格")
}

# ========== 3. 矢量投影对齐 ==========
if (is.na(st_crs(regions))) {
  st_crs(regions) <- "EPSG:4326"
}

if (st_crs(regions)$wkt != crs(tas_raster, proj = TRUE)) {
  regions <- st_transform(regions, crs = crs(tas_raster, proj = TRUE))
}
```

#### 步骤2：计算加权和非加权统计

```r
# ========== 4. 提取区域内的栅格值 ==========
# 提取温度值
tas_values <- extract(tas_raster, regions, ID = TRUE, na.rm = FALSE)

# 提取对应的人口值
pop_values <- extract(pop_raster, regions, ID = TRUE, na.rm = FALSE)

# ========== 5. 计算加权平均 ==========
# 合并数据
combined <- data.frame(
  ID = tas_values$ID,
  tas = tas_values[, 2],    # 温度值
  pop = pop_values[, 2]     # 人口值
)

# 移除NA值
combined <- combined[!is.na(combined$tas) & !is.na(combined$pop) & combined$pop > 0, ]

# 按区域分组计算加权平均
weighted_result <- combined %>%
  group_by(ID) %>%
  summarise(
    tas_weighted_mean = sum(tas * pop) / sum(pop),  # 加权平均
    tas_simple_mean = mean(tas),                     # 简单平均
    total_population = sum(pop),                     # 区域总人口
    n_pixels = n()                                   # 像元数
  )

# ========== 6. 合并城市名称 ==========
result <- data.frame(
  city = regions$NAME,
  tas_weighted_mean = weighted_result$tas_weighted_mean,
  tas_simple_mean = weighted_result$tas_simple_mean,
  total_population = weighted_result$total_population,
  n_pixels = weighted_result$n_pixels
)

# ========== 7. 保存结果 ==========
write_csv(result, "population_weighted_zonal_mean.csv")
print(result)
```

---

### 方法2：使用terra的weights参数（更简洁）

`terra::extract()` 支持直接传入权重栅格：

```r
library(terra)
library(sf)

# 读取数据
tas <- rast("temperature.tif")
pop <- rast("population.tif")
regions <- st_read("hunan.shp", quiet = TRUE)

# 对齐数据（同上）
pop <- resample(pop, tas, method = "bilinear")
regions <- st_transform(regions, crs = crs(tas, proj = TRUE))

# ========== 直接使用weights参数 ==========
# 加权提取
weighted <- extract(
  tas, 
  regions, 
  fun = function(x, w) sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE),
  weights = pop,
  ID = TRUE
)

# 简单平均（对比用）
simple <- extract(tas, regions, fun = mean, na.rm = TRUE, ID = TRUE)

# 整理结果
result <- data.frame(
  city = regions$NAME,
  tas_weighted = weighted[, 2],
  tas_simple = simple[, 2]
)

write_csv(result, "weighted_result.csv")
```

---

### 方法3：处理NetCDF多时间步的人口加权

```r
library(ncdf4)
library(terra)
library(sf)
library(dplyr)

# ========== 1. 读取NetCDF ==========
nc <- nc_open("climate.nc")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
tas_array <- ncvar_get(nc, "tas")  # 维度: [lon, lat, time]
nc_close(nc)

# ========== 2. 准备人口栅格 ==========
pop <- rast("population.tif")
regions <- st_read("hunan.shp", quiet = TRUE)

# ========== 3. 逐时间步处理 ==========
n_times <- length(time)
result_list <- vector("list", n_times)

for (t in 1:n_times) {
  # 提取第t个时间步的数据
  tas_slice <- tas_array[,,t]
  
  # 转换为raster（terra格式）
  tas_raster <- rast(
    vals = as.vector(t(tas_slice[, ncol(tas_slice):1])),  # 转置并翻转
    nrows = length(lat),
    ncols = length(lon),
    extent = ext(min(lon), max(lon), min(lat), max(lat)),
    crs = "EPSG:4326"
  )
  
  # 重采样人口栅格到当前温度栅格
  pop_aligned <- resample(pop, tas_raster, method = "bilinear")
  
  # 投影矢量
  regions_proj <- st_transform(regions, crs = crs(tas_raster, proj = TRUE))
  
  # 加权提取
  weighted <- extract(
    tas_raster,
    regions_proj,
    fun = function(x, w) sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE),
    weights = pop_aligned,
    ID = TRUE
  )
  
  # 保存结果
  result_list[[t]] <- data.frame(
    city = regions_proj$NAME,
    time_index = t,
    tas_weighted = weighted[, 2]
  )
}

# 合并所有时间步
result_all <- bind_rows(result_list)
write_csv(result_all, "nc_pop_weighted_all_times.csv")

# 计算每个城市的时间平均
result_summary <- result_all %>%
  group_by(city) %>%
  summarise(
    tas_weighted_mean = mean(tas_weighted, na.rm = TRUE),
    tas_weighted_min = min(tas_weighted, na.rm = TRUE),
    tas_weighted_max = max(tas_weighted, na.rm = TRUE)
  )

write_csv(result_summary, "nc_pop_weighted_summary.csv")
```

---

## 完整示例代码

### 示例1：TIF文件的简单区域统计

```r
#!/usr/bin/env Rscript
# 文件名: tif_zonal_simple.R
# 用法: Rscript tif_zonal_simple.R input.tif regions.shp output.csv

library(terra)
library(sf)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("用法: Rscript tif_zonal_simple.R <input.tif> <regions.shp> <output.csv>")
}

tif_path <- args[1]
shp_path <- args[2]
out_csv <- args[3]

# 读取数据
raster <- rast(tif_path)
regions <- st_read(shp_path, quiet = TRUE)

# CRS对齐
if (is.na(st_crs(regions))) st_crs(regions) <- "EPSG:4326"
regions <- st_transform(regions, crs = crs(raster, proj = TRUE))

# 区域统计
zonal <- extract(raster, regions, fun = mean, na.rm = TRUE, ID = TRUE)

# 输出
result <- data.frame(
  region_id = 1:nrow(regions),
  region_name = regions[[1]],  # 假设第一列是名称
  mean_value = zonal[, 2]
)

write_csv(result, out_csv)
cat("完成！结果已保存到", out_csv, "\n")
```

---

### 示例2：NetCDF人口加权统计（封装函数）

```r
#!/usr/bin/env Rscript
# 文件名: nc_pop_weighted.R

library(terra)
library(sf)
library(dplyr)
library(readr)

#' 计算NetCDF的人口加权区域统计
#'
#' @param nc_path NetCDF文件路径
#' @param pop_path 人口栅格路径
#' @param shp_path Shapefile路径
#' @param out_csv 输出CSV路径
#' @param var_name NetCDF变量名（默认"tas"）
#' @param time_mode 时间处理模式："avg"平均, "first"首个, "last"最后, "all"全部
#' @return data.frame 结果数据框
compute_pop_weighted <- function(nc_path, pop_path, shp_path, out_csv,
                                 var_name = "tas", time_mode = "avg") {
  
  # 1. 读取数据
  nc <- rast(nc_path)
  pop <- rast(pop_path)
  regions <- st_read(shp_path, quiet = TRUE)
  
  # 2. 处理时间维度
  if (time_mode == "avg") {
    nc <- mean(nc, na.rm = TRUE)
  } else if (time_mode == "first") {
    nc <- nc[[1]]
  } else if (time_mode == "last") {
    nc <- nc[[nlyr(nc)]]
  } else if (time_mode != "all") {
    stop("time_mode 必须是 'avg', 'first', 'last' 或 'all'")
  }
  
  # 3. 对齐栅格
  if (!compareGeom(nc, pop, stopOnError = FALSE)) {
    pop <- resample(pop, nc, method = "bilinear")
  }
  
  # 4. 对齐矢量
  if (is.na(st_crs(regions))) st_crs(regions) <- "EPSG:4326"
  regions <- st_transform(regions, crs = crs(nc, proj = TRUE))
  
  # 5. 加权提取
  if (time_mode == "all" && nlyr(nc) > 1) {
    # 逐层处理
    result_list <- lapply(1:nlyr(nc), function(i) {
      weighted <- extract(
        nc[[i]], regions,
        fun = function(x, w) sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE),
        weights = pop, ID = TRUE
      )
      data.frame(
        city = regions[[1]],
        time_index = i,
        value_weighted = weighted[, 2]
      )
    })
    result <- bind_rows(result_list)
  } else {
    # 单层处理
    weighted <- extract(
      nc, regions,
      fun = function(x, w) sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE),
      weights = pop, ID = TRUE
    )
    simple <- extract(nc, regions, fun = mean, na.rm = TRUE, ID = TRUE)
    
    result <- data.frame(
      city = regions[[1]],
      value_weighted = weighted[, 2],
      value_simple = simple[, 2]
    )
  }
  
  # 6. 保存结果
  write_csv(result, out_csv)
  message(sprintf("结果已保存到 %s", out_csv))
  
  return(result)
}

# 命令行调用
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) {
    cat("用法: Rscript nc_pop_weighted.R <nc_file> <pop_file> <shp_file> <out_csv> [time_mode]\n")
    quit(status = 1)
  }
  
  time_mode <- ifelse(length(args) >= 5, args[5], "avg")
  
  compute_pop_weighted(
    nc_path = args[1],
    pop_path = args[2],
    shp_path = args[3],
    out_csv = args[4],
    time_mode = time_mode
  )
}
```

**使用示例**：
```bash
# 计算时间平均的人口加权值
Rscript nc_pop_weighted.R climate.nc population.tif hunan.shp result.csv avg

# 计算所有时间步的人口加权值
Rscript nc_pop_weighted.R climate.nc population.tif hunan.shp result_all.csv all
```

---

## 常见问题和解决方案

### 1. CRS不匹配

**问题**：栅格和矢量坐标系不一致导致结果错误

**解决**：
```r
# 检查CRS
print(crs(raster, proj = TRUE))
print(st_crs(vector)$wkt)

# 统一CRS：将矢量投影到栅格
vector <- st_transform(vector, crs = crs(raster, proj = TRUE))
```

### 2. 栅格分辨率不一致

**问题**：温度栅格和人口栅格分辨率不同

**解决**：
```r
# 重采样到目标分辨率
pop_resampled <- resample(pop, target_raster, method = "bilinear")
```

### 3. 内存不足

**问题**：大文件导致内存溢出

**解决**：
```r
# 方法1：裁剪到研究区域
raster_crop <- crop(raster, regions)

# 方法2：降低分辨率
raster_agg <- aggregate(raster, fact = 4)  # 分辨率降低4倍

# 方法3：分块处理（terra自动优化）
terraOptions(memfrac = 0.6)  # 使用60%内存
```

### 4. NetCDF时间解析错误

**问题**：无法正确解析时间维度

**解决**：
```r
library(ncdf4)
nc <- nc_open("file.nc")

# 读取时间单位
time_units <- ncatt_get(nc, "time", "units")$value
# 例如: "days since 1850-01-01 00:00:00"

# 手动解析
library(lubridate)
time_origin <- as.POSIXct("1850-01-01", tz = "UTC")
time_values <- ncvar_get(nc, "time")
dates <- time_origin + days(time_values)

nc_close(nc)
```

### 5. 区域没有有效像元

**问题**：某些区域提取结果为NA

**解决**：
```r
# 检查重叠
raster_bbox <- st_bbox(st_as_sf(as.polygons(ext(raster))))
vector_bbox <- st_bbox(vector)

# 确保有重叠
st_intersects(raster_bbox, vector_bbox, sparse = FALSE)

# 可视化检查
plot(raster)
plot(st_geometry(vector), add = TRUE, border = "red")
```

---

## 性能优化建议

1. **选择合适的包**：
   - 小数据：`terra` 足够
   - 大数据/多维：`stars` + `exactextractr`
   - 需要底层控制：`ncdf4`

2. **并行处理**：
```r
library(future.apply)
plan(multisession, workers = 4)

result <- future_lapply(1:n_times, function(t) {
  # 处理第t个时间步
})
```

3. **避免重复投影**：
```r
# 预先统一CRS，避免循环中重复投影
regions <- st_transform(regions, target_crs)
```

4. **使用适当的统计方法**：
```r
# 快速但不精确
extract(raster, vector, fun = mean)

# 精确但较慢
exact_extract(raster, vector, fun = "mean")
```

---

## 参考资源

- **terra包文档**: https://rspatial.github.io/terra/
- **sf包文档**: https://r-spatial.github.io/sf/
- **stars包文档**: https://r-spatial.github.io/stars/
- **exactextractr文档**: https://isciences.gitlab.io/exactextractr/

---

## 总结

本教程涵盖了使用R进行栅格区域统计的完整流程：

1. **基础工具**：`terra`处理栅格，`sf`处理矢量，`stars`处理多维数组
2. **标准流程**：数据读取 → CRS对齐 → 区域统计 → 结果输出
3. **简单统计**：使用`extract()`计算区域平均值
4. **加权统计**：结合人口数据计算加权平均，更真实反映人类活动区域的情况

选择合适的方法取决于：
- 数据格式（TIF vs NetCDF）
- 时间维度（单时刻 vs 多时间步）
- 精度要求（快速 vs 精确）
- 数据规模（小文件 vs 大数据）

建议从简单案例开始实践，逐步掌握各种工具和技巧。
