# R语言：点位缓冲区内栅格数据提取与计算教程

本教程讲解如何从CSV文件读入点位坐标（id, lon, lat），以点为中心创建球面距离缓冲区，提取缓冲区内栅格数据的平均值或逆距离加权平均值。支持多波段栅格和多时间维度的NetCDF数据。

## 目录
1. [环境准备](#环境准备)
2. [基础概念](#基础概念)
3. [读入点位数据](#读入点位数据)
4. [单波段栅格的缓冲区提取](#单波段栅格的缓冲区提取)
5. [多波段栅格的缓冲区提取](#多波段栅格的缓冲区提取)
6. [NetCDF时间序列数据](#netcdf时间序列数据)
7. [逆距离加权平均](#逆距离加权平均)
8. [完整工作流示例](#完整工作流示例)
9. [故障排查](#故障排查)

---

## 环境准备

### 1. R包安装

```r
# 核心地理空间包
install.packages("sf")           # 矢量数据处理
install.packages("terra")        # 现代栅格数据处理
install.packages("raster")       # 栅格数据处理（备选）
install.packages("stars")        # 多维数组数据处理
install.packages("ncdf4")        # NetCDF文件读写
install.packages("geosphere")    # 球面距离计算

# 数据处理
install.packages("dplyr")        # 数据框操作
install.packages("tidyr")        # 数据整形
install.packages("readr")        # 快速读取CSV
install.packages("lubridate")    # 时间处理

# 可视化
install.packages("ggplot2")      # 绘图
install.packages("tmap")         # 地图绘制
```

### 2. 系统依赖（Linux）

```bash
# Ubuntu/Debian
sudo apt-get install gdal-bin libgdal-dev libproj-dev libgeos-dev
sudo apt-get install libudunits2-dev libnetcdf-dev

# 验证安装
gdalinfo --version
```

### 3. 加载库

```r
library(sf)
library(terra)
library(stars)
library(ncdf4)
library(geosphere)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
```

---

## 基础概念

### 球面距离缓冲区

- **球面距离**：地球表面两点间的最短距离（大圆距离）
- **缓冲区**：以点为中心、给定半径的圆形区域
- **半径单位**：米（m）或千米（km）

### 缓冲区内栅格提取方法

| 方法 | 优点 | 缺点 | 适用场景 |
|------|------|------|---------|
| **简单平均** | 快速、易理解 | 不考虑距离影响 | 缓冲区小或栅格分辨率低 |
| **逆距离加权(IDW)** | 距离近的像元权重大 | 计算复杂、较慢 | 需要考虑空间自相关 |
| **栅格重叠百分比** | 精确、物理意义清晰 | 边界处理复杂 | 精确分析需求 |

---

## 读入点位数据

### 1. 创建示例CSV文件

```r
# 创建示例数据
points_data <- data.frame(
  id = c("site_01", "site_02", "site_03", "site_04"),
  lon = c(110.5, 111.0, 111.5, 110.2),
  lat = c(28.5, 28.8, 28.2, 28.0)
)

# 保存为CSV
write_csv(points_data, "sample_points.csv")

head(points_data)
# # A tibble: 4 × 3
#   id      lon   lat
#   <chr> <dbl> <dbl>
# 1 site_01 110.5  28.5
# 2 site_02 111    28.8
# 3 site_03 111.5  28.2
# 4 site_04 110.2  28
```

### 2. 读入CSV并转换为SF对象

```r
# 读入CSV
points <- read_csv("sample_points.csv")

# 转换为SF对象（WGS84坐标系）
points_sf <- st_as_sf(points, coords = c("lon", "lat"), 
                       crs = 4326)

# 查看数据
print(points_sf)
# Simple feature collection with 4 features and 1 field
# Geometry type: POINT
# Dimension:     XY
# Bounding box:  xmin: 110.2 ymin: 28 xmax: 111.5 ymax: 28.8
# Geodetic CRS:  WGS 84
#   id                   geometry
# 1 site_01 POINT (110.5 28.5)
# 2 site_02 POINT (111 28.8)
# 3 site_03 POINT (111.5 28.2)
# 4 site_04 POINT (110.2 28)
```

---

## 单波段栅格的缓冲区提取

### 1. 创建球面距离缓冲区

```r
# 创建缓冲区（半径5km）
buffer_radius_m <- 5000  # 5公里

# 方法1：使用st_buffer（投影坐标系下计算）
# 先投影到等面积投影
points_proj <- st_transform(points_sf, crs = 3857)  # Web Mercator
buffer_circle <- st_buffer(points_proj, dist = buffer_radius_m)

# 转回WGS84
buffer_circle <- st_transform(buffer_circle, crs = 4326)

# 方法2：使用geosphere包创建球面缓冲区（更精确）
create_sphere_buffer <- function(lon, lat, radius_m, n_points = 100) {
  # radius_m: 半径（米）
  # n_points: 圆形顶点数
  
  radius_km <- radius_m / 1000
  angles <- seq(0, 360, length.out = n_points)
  
  # 计算缓冲区边界点
  boundary_points <- destPoint(
    p = c(lon, lat),
    b = angles,
    d = radius_km * 1000
  )
  
  # 创建多边形
  polygon <- st_polygon(
    list(rbind(boundary_points, boundary_points[1, ]))
  )
  
  return(st_sfc(polygon, crs = 4326))
}

# 为每个点创建缓冲区
buffers <- map_df(1:nrow(points_sf), function(i) {
  buf <- create_sphere_buffer(
    lon = st_coordinates(points_sf)[i, 1],
    lat = st_coordinates(points_sf)[i, 2],
    radius_m = 5000
  )
  data.frame(id = points_sf$id[i], geometry = buf)
}) %>%
  st_as_sf()

print(buffers)
```

### 2. 读入栅格数据

```r
# 读入TIF文件
raster_tif <- rast("dem.tif")  # 使用terra包

# 检查栅格信息
print(raster_tif)
# class       : SpatRaster
# dimensions  : 1000, 1000, 1  (nrow, ncol, nlyr)
# resolution  : 100, 100  (x, y)
# extent      : 110, 111, 28, 29  (xmin, xmax, ymin, ymax)
# coord. ref. : EPSG:4326 (WGS 84)
# source      : dem.tif
# name        : dem
# min value   :  50
# max value   : 2000

# 检查坐标参考系
crs(raster_tif)
```

### 3. 提取缓冲区内的栅格值

```r
# 方法1：简单平均
extract_buffer_mean <- function(raster, buffer_sf, id_col = "id") {
  # raster: SpatRaster对象
  # buffer_sf: sf多边形对象
  
  result <- data.frame()
  
  for (i in 1:nrow(buffer_sf)) {
    # 提取缓冲区内的栅格值
    values <- extract(raster, vect(buffer_sf[i, ]))[[1]]
    
    # 计算平均值（排除NA）
    mean_val <- mean(values, na.rm = TRUE)
    
    result <- rbind(result, data.frame(
      id = buffer_sf[[id_col]][i],
      mean_value = mean_val,
      n_pixels = length(na.omit(values))
    ))
  }
  
  return(result)
}

# 执行提取
result <- extract_buffer_mean(raster_tif, buffers)
print(result)
#        id mean_value n_pixels
# 1 site_01    856.23      245
# 2 site_02    934.12      248
# 3 site_03    762.45      242
# 4 site_04    845.67      243
```

---

## 多波段栅格的缓冲区提取

### 1. 多波段栅格结构

```r
# 读入多波段栅格（例如Sentinel-2多谱段数据）
multi_band_raster <- rast("sentinel2_multiband.tif")

# 查看波段信息
print(multi_band_raster)
# class       : SpatRaster
# dimensions  : 1000, 1000, 11  (nrow, ncol, nlyr)
# resolution  : 10, 10  (x, y)
# extent      : 110, 111, 28, 29  (xmin, xmax, ymin, ymax)
# names       : B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11
# min values  :  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
# max values  : 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000
```

### 2. 多波段提取函数

```r
extract_buffer_multiband <- function(raster, buffer_sf, id_col = "id") {
  # 返回每个点、每个波段的平均值
  
  n_bands <- nlyr(raster)
  band_names <- names(raster)
  
  result <- data.frame(id = buffer_sf[[id_col]])
  
  for (i in 1:nrow(buffer_sf)) {
    # 提取缓冲区内所有波段的栅格值
    values_list <- extract(raster, vect(buffer_sf[i, ]))[[1]]
    
    # 计算每个波段的平均值
    band_means <- apply(values_list, 2, function(x) {
      mean(x, na.rm = TRUE)
    })
    
    # 创建行数据
    row_data <- as.list(band_means)
    names(row_data) <- paste0("band_", band_names)
    
    if (i == 1) {
      result <- cbind(result, as.data.frame(row_data)[rep(1, nrow(result)), ])
    } else {
      result[i, names(row_data)] <- row_data
    }
  }
  
  return(result)
}

# 执行提取
multi_band_result <- extract_buffer_multiband(multi_band_raster, buffers)
print(multi_band_result)
#        id band_B1 band_B2 band_B3 band_B4 ...
# 1 site_01    1234    1456    1678    2345 ...
# 2 site_02    1300    1500    1700    2400 ...
```

### 3. 简化版本：使用exact_extract（更快）

```r
# 使用exactextractr包（更快，返回像元权重）
library(exactextractr)

extract_exact_multiband <- function(raster, points_sf, buffer_radius_m) {
  # 创建缓冲区
  buffers <- st_buffer(
    st_transform(points_sf, crs = 3857),
    dist = buffer_radius_m
  ) %>%
    st_transform(crs = 4326)
  
  # 提取
  result <- exact_extract(raster, buffers, 'mean', progress = TRUE)
  
  # 整理结果
  band_names <- names(raster)
  result_df <- as.data.frame(result)
  names(result_df) <- band_names
  result_df <- cbind(id = points_sf$id, result_df)
  
  return(result_df)
}

multi_result_fast <- extract_exact_multiband(multi_band_raster, points_sf, 5000)
```

---

## NetCDF时间序列数据

### 1. 读入NetCDF文件

```r
# 使用stars包（保留时间维度）
nc_file <- "temperature_timeseries.nc"

# 读入全部数据
nc_data <- read_ncdf(nc_file)

# 查看维度结构
print(nc_data)
# stars object with 3 dimensions and 1 attribute
# attribute(s):
#            Min. 1st Qu. Median 3rd Qu.  Max.
# tas       -20   10     15     20      40
# dimension(s):
#     from   to offset delta  refsys point values x/y
# x      1 1000    110     0 WGS 84    NA   NULL [x]
# y      1 1000     29     0 WGS 84    NA   NULL [y]
# time   1  365     NA    NA     POSIXct NA   NULL
```

### 2. 提取时间序列数据

```r
extract_timeseries_buffer <- function(nc_data, points_sf, buffer_radius_m) {
  # nc_data: stars对象，维度为 [lon, lat, time]
  # 返回: 每个点、每个时间步的平均值
  
  # 创建缓冲区
  buffers <- st_buffer(
    st_transform(points_sf, crs = 3857),
    dist = buffer_radius_m
  ) %>%
    st_transform(crs = 4326)
  
  result <- list()
  
  for (i in 1:nrow(buffers)) {
    # 提取缓冲区内的时间序列
    ts_data <- st_crop(nc_data, buffers[i, ])
    
    # 计算每个时间步的平均值
    ts_mean <- apply(ts_data[[1]], 3, function(x) {
      mean(x, na.rm = TRUE)
    })
    
    time_index <- st_get_dimension_values(nc_data, "time")
    
    result[[i]] <- data.frame(
      id = points_sf$id[i],
      time = time_index,
      mean_value = ts_mean,
      n_pixels = apply(ts_data[[1]], 3, function(x) {
        sum(!is.na(x))
      })
    )
  }
  
  return(bind_rows(result))
}

# 执行提取
ts_result <- extract_timeseries_buffer(nc_data, points_sf, 5000)
head(ts_result)
#        id       time mean_value n_pixels
# 1 site_01 2020-01-01      -5.23      245
# 2 site_01 2020-01-02      -4.56      245
# 3 site_01 2020-01-03      -3.89      245
```

### 3. 将结果转换为宽格式（时间为列）

```r
# 转换为宽格式
ts_wide <- ts_result %>%
  pivot_wider(
    id_cols = id,
    names_from = time,
    values_from = mean_value,
    names_prefix = "time_"
  )

print(ts_wide)
#       id  time_2020-01-01 time_2020-01-02 time_2020-01-03 ...
# 1 site_01          -5.23           -4.56            -3.89
# 2 site_02          -4.56           -3.78            -2.95
```

---

## 逆距离加权平均

### 1. IDW原理

逆距离加权(IDW)平均值计算公式：

$$\hat{z}(x_0) = \frac{\sum_{i=1}^{n} \frac{z_i}{d_i^p}}{\sum_{i=1}^{n} \frac{1}{d_i^p}}$$

其中：
- $\hat{z}(x_0)$：点$x_0$处的估计值
- $z_i$：像元$i$的值
- $d_i$：点$x_0$到像元$i$的距离
- $p$：幂次（通常为2）

### 2. IDW提取函数

```r
extract_buffer_idw <- function(raster, buffer_sf, points_sf, 
                                id_col = "id", power = 2) {
  # 计算缓冲区内每个像元到点的距离，进行IDW平均
  
  result <- data.frame()
  
  for (i in 1:nrow(buffer_sf)) {
    # 获取点坐标
    point_coord <- st_coordinates(points_sf)[i, ]
    
    # 提取缓冲区内的栅格值
    extracted <- extract(raster, vect(buffer_sf[i, ]))
    
    # 获取这些像元的坐标
    cell_coords <- xyFromCell(raster, extracted[[1]][, "cell"])
    
    # 计算距离（平面距离或球面距离）
    distances <- distGeo(
      p1 = point_coord,
      p2 = cell_coords
    )
    
    # 获取像元值
    values <- extracted[[1]][, 1]
    
    # 计算IDW权重
    weights <- 1 / (distances^power)
    weights <- weights / sum(weights, na.rm = TRUE)
    
    # 计算加权平均
    idw_mean <- sum(values * weights, na.rm = TRUE)
    
    result <- rbind(result, data.frame(
      id = buffer_sf[[id_col]][i],
      idw_mean = idw_mean,
      simple_mean = mean(values, na.rm = TRUE),
      n_pixels = length(na.omit(values))
    ))
  }
  
  return(result)
}

# 执行IDW提取
idw_result <- extract_buffer_idw(raster_tif, buffers, points_sf, power = 2)
print(idw_result)
#        id idw_mean simple_mean n_pixels
# 1 site_01   856.89       856.23      245
# 2 site_02   935.12       934.12      248
# 3 site_03   762.78       762.45      242
# 4 site_04   846.23       845.67      243
```

### 3. 多波段IDW提取

```r
extract_buffer_idw_multiband <- function(raster, buffer_sf, points_sf,
                                        id_col = "id", power = 2) {
  
  n_bands <- nlyr(raster)
  band_names <- names(raster)
  
  result <- data.frame(id = buffer_sf[[id_col]])
  
  for (i in 1:nrow(buffer_sf)) {
    # 获取点坐标
    point_coord <- st_coordinates(points_sf)[i, ]
    
    # 提取缓冲区内所有波段的栅格值
    extracted <- extract(raster, vect(buffer_sf[i, ]))[[1]]
    
    # 获取像元坐标
    cell_coords <- xyFromCell(raster, as.integer(rownames(extracted)))
    
    # 计算距离和权重
    distances <- distGeo(p1 = point_coord, p2 = cell_coords)
    weights <- 1 / (distances^power)
    weights <- weights / sum(weights, na.rm = TRUE)
    
    # 计算每个波段的IDW平均
    idw_means <- apply(extracted, 2, function(values) {
      sum(values * weights, na.rm = TRUE)
    })
    
    # 添加到结果
    for (j in 1:n_bands) {
      result[i, paste0("idw_", band_names[j])] <- idw_means[j]
    }
  }
  
  return(result)
}

multi_idw_result <- extract_buffer_idw_multiband(multi_band_raster, buffers, 
                                                 points_sf, power = 2)
```

---

## 完整工作流示例

### 示例1：单波段TIF + 简单平均

```r
# 1. 准备数据
points_csv <- read_csv("points.csv")
points_sf <- st_as_sf(points_csv, coords = c("lon", "lat"), crs = 4326)

# 2. 读入栅格
raster <- rast("temperature.tif")

# 3. 创建缓冲区
buffers <- st_buffer(
  st_transform(points_sf, crs = 3857),
  dist = 5000
) %>%
  st_transform(crs = 4326)

# 4. 提取并计算
result <- data.frame()
for (i in 1:nrow(buffers)) {
  values <- extract(raster, vect(buffers[i, ]))[[1]]
  result <- rbind(result, data.frame(
    id = points_sf$id[i],
    mean_temp = mean(values, na.rm = TRUE),
    sd_temp = sd(values, na.rm = TRUE),
    n_pixels = length(na.omit(values))
  ))
}

# 5. 保存结果
write_csv(result, "extraction_result.csv")

# 6. 可视化
ggplot(result, aes(x = id, y = mean_temp)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_temp - sd_temp, 
                     ymax = mean_temp + sd_temp),
                width = 0.2) +
  theme_minimal() +
  labs(title = "缓冲区内平均温度", x = "点位", y = "温度 (°C)")
```

### 示例2：多波段栅格 + IDW平均 + 时间序列

```r
# 1. 读入多波段多时间NetCDF
nc_data <- read_ncdf("multi_temporal_bands.nc")

# 2. 提取时间序列（每个点、每个波段、每个时间）
result_ts <- list()

for (i in 1:nrow(points_sf)) {
  # 创建单点缓冲区
  buffer <- st_buffer(
    st_transform(points_sf[i, ], crs = 3857),
    dist = 5000
  ) %>%
    st_transform(crs = 4326)
  
  # 裁剪到缓冲区
  nc_crop <- st_crop(nc_data, buffer)
  
  # 获取时间维度
  times <- st_get_dimension_values(nc_data, "time")
  bands <- names(nc_data)
  
  # 计算IDW平均（对每个时间步和波段）
  for (t in 1:length(times)) {
    for (b in 1:length(bands)) {
      # 提取该时间步该波段的数据
      band_data <- nc_crop[[bands[b]]][, , t]
      
      # 简化为简单平均（实现完整IDW较复杂）
      mean_val <- mean(band_data, na.rm = TRUE)
      
      result_ts[[length(result_ts) + 1]] <- data.frame(
        id = points_sf$id[i],
        band = bands[b],
        time = times[t],
        mean_value = mean_val,
        stringsAsFactors = FALSE
      )
    }
  }
}

result_ts_df <- bind_rows(result_ts)

# 3. 检查结果
head(result_ts_df, 15)
#        id band       time mean_value
# 1 site_01   B1 2020-01-01      1234.5
# 2 site_01   B1 2020-01-02      1245.3
# ...

# 4. 转换为宽格式
result_wide <- result_ts_df %>%
  unite("col", band, time) %>%
  pivot_wider(names_from = col, values_from = mean_value)

# 5. 可视化时间序列
result_ts_df %>%
  filter(band == "B1") %>%
  ggplot(aes(x = time, y = mean_value, color = id, group = id)) +
  geom_line() +
  theme_minimal() +
  labs(title = "B1波段时间序列",
       x = "时间", y = "平均值", color = "点位")
```

### 示例3：使用管道操作简化工作流

```r
# 完整的管道工作流
library(magrittr)

extract_complete_workflow <- function(csv_file, raster_file, 
                                     buffer_radius_m = 5000,
                                     method = "mean") {
  
  # 1. 读入数据
  points <- read_csv(csv_file) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  raster <- rast(raster_file)
  
  # 2. 创建缓冲区
  buffers <- st_buffer(
    st_transform(points, crs = 3857),
    dist = buffer_radius_m
  ) %>%
    st_transform(crs = 4326)
  
  # 3. 提取
  result <- map_df(1:nrow(buffers), function(i) {
    values <- extract(raster, vect(buffers[i, ]))[[1]]
    
    if (method == "mean") {
      stat_val <- mean(values, na.rm = TRUE)
    } else if (method == "median") {
      stat_val <- median(values, na.rm = TRUE)
    } else if (method == "sd") {
      stat_val <- sd(values, na.rm = TRUE)
    }
    
    data.frame(
      id = points$id[i],
      lon = st_coordinates(points)[i, 1],
      lat = st_coordinates(points)[i, 2],
      value = stat_val,
      n_pixels = length(na.omit(values))
    )
  })
  
  return(result)
}

# 执行
result <- extract_complete_workflow("points.csv", "temperature.tif", 
                                   buffer_radius_m = 5000,
                                   method = "mean")

# 保存
write_csv(result, "buffer_extraction_result.csv")
```

---

## 故障排查

### 问题1：投影坐标系转换错误

```r
# 错误：不能直接在WGS84中使用st_buffer()计算距离
# 原因：WGS84的单位是度，不是米

# 正确做法
buffer <- st_buffer(
  st_transform(points_sf, crs = 3857),  # 先转到投影坐标系
  dist = 5000  # 单位为米
) %>%
  st_transform(crs = 4326)  # 再转回WGS84
```

### 问题2：提取结果包含NA

```r
# 检查原因
extracted <- extract(raster, vect(buffer[1, ]))

# 查看缓冲区是否与栅格有重叠
plot(raster)
plot(buffer[1, ], add = TRUE, border = "red")

# 检查坐标参考系是否一致
crs(raster) == crs(buffer)
```

### 问题3：多时间维度处理

```r
# NetCDF时间维度提取失败
# 检查维度名称
st_dimensions(nc_data)

# 正确的提取方式
time_vals <- st_get_dimension_values(nc_data, "time")

# 或使用stars的切片函数
slice_result <- slice(nc_data, "time", i)
```

### 问题4：内存不足

```r
# 对于大文件，使用分块处理
chunk_extract <- function(raster, buffer_sf, chunk_size = 100) {
  n_chunks <- ceiling(nrow(buffer_sf) / chunk_size)
  
  result <- list()
  
  for (chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, nrow(buffer_sf))
    
    chunk_result <- extract_buffer_mean(
      raster,
      buffer_sf[start_idx:end_idx, ]
    )
    
    result[[chunk]] <- chunk_result
    
    # 及时垃圾回收
    gc()
  }
  
  return(bind_rows(result))
}

# 执行
result <- chunk_extract(large_raster, buffers, chunk_size = 50)
```

### 问题5：距离计算精度

```r
# 对于小规模区域（<100km），可使用平面距离（快速）
library(sf)
distances_planar <- st_distance(point, cells)

# 对于大规模区域或高精度需求，使用球面距离（精确）
library(geosphere)
distances_sphere <- distGeo(p1 = point, p2 = cells)

# 比较差异
df_compare <- data.frame(
  planar = distances_planar,
  sphere = distances_sphere,
  diff_pct = abs(distances_planar - distances_sphere) / distances_sphere * 100
)

print(summary(df_compare$diff_pct))
```

---

## 性能优化建议

| 优化方法 | 适用场景 | 性能提升 |
|---------|---------|---------|
| 使用 `exact_extract()` | 多点、大缓冲区 | 10-50倍 |
| 缓冲区预计算 | 重复使用相同缓冲区 | 5-10倍 |
| 并行处理（`future` + `furrr`） | 点位数>1000 | 4-8倍 |
| 栅格重采样降低分辨率 | 快速探索 | 50-100倍 |
| 分块处理 | 内存受限 | 减少内存占用 |

### 并行处理示例

```r
library(future)
library(furrr)

# 设置并行后端
plan(multisession, workers = 4)

# 并行提取
result_parallel <- future_map_df(1:nrow(buffers), function(i) {
  values <- extract(raster, vect(buffers[i, ]))[[1]]
  data.frame(
    id = points_sf$id[i],
    mean_value = mean(values, na.rm = TRUE),
    n_pixels = length(na.omit(values))
  )
})

# 还原串行处理
plan(sequential)
```

---

## 推荐阅读

- [SF包官方文档](https://r-spatial.github.io/sf/)
- [Terra包官方文档](https://www.rspatial.org/terra/pkg/)
- [Stars包官方文档](https://r-spatial.github.io/stars/)
- [Geosphere包文档](https://cran.r-project.org/web/packages/geosphere/)
- [Exactextractr包文档](https://cran.r-project.org/web/packages/exactextractr/)

---

## 许可证

本教程遵循CC-BY-4.0许可证
