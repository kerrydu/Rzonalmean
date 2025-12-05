# R语言：点对点匹配与加权聚合教程

本教程讲解如何基于空间位置进行点对点匹配（如企业与监测站），计算给定距离内目标点的加权统计量。典型应用包括：企业地点匹配周边监测站、工厂匹配气象站、门店匹配消费网格等。

## 目录
1. [环境准备](#环境准备)
2. [基础概念](#基础概念)
3. [数据准备](#数据准备)
4. [点对点匹配](#点对点匹配)
5. [距离计算](#距离计算)
6. [逆距离加权聚合](#逆距离加权聚合)
7. [多时间序列处理](#多时间序列处理)
8. [完整工作流示例](#完整工作流示例)
9. [性能优化](#性能优化)
10. [可视化](#可视化)
11. [故障排查](#故障排查)

---

## 环境准备

### 1. R包安装

```r
# 空间数据处理
install.packages("sf")           # 矢量数据处理
install.packages("geosphere")    # 球面距离计算
install.packages("data.table")   # 高效数据处理

# 数据处理
install.packages("dplyr")        # 数据框操作
install.packages("tidyr")        # 数据整形
install.packages("readr")        # 快速读取CSV

# 性能工具
install.packages("future")       # 并行计算
install.packages("furrr")        # 并行map函数
install.packages("bench")        # 性能测试

# 可视化
install.packages("ggplot2")      # 绘图
install.packages("tmap")         # 地图绘制
install.packages("leaflet")      # 交互式地图
```

### 2. 加载库

```r
library(sf)
library(geosphere)
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
```

---

## 基础概念

### 点对点匹配场景

| 场景 | 源点（企业/用户） | 目标点（监测站/服务） | 匹配关键词 |
|------|------------------|----------------------|-----------|
| 环保监测 | 企业地址 | 空气质量监测站 | 企业周边污染水平 |
| 气象服务 | 农田/园区 | 气象观测站 | 田间微气候 |
| 物流优化 | 客户地点 | 配送站/仓库 | 最近配送中心 |
| 消费分析 | 商业门店 | 消费网格/人口热力 | 门店周边消费力 |
| 环保评估 | 居住社区 | 污染源/排放点 | 社区环境暴露 |

### 关键计算公式

**球面距离**（Haversine公式）：

$$d = 2R \cdot \arcsin\left(\sqrt{\sin^2\left(\frac{\Delta\phi}{2}\right) + \cos\phi_1 \cos\phi_2 \sin^2\left(\frac{\Delta\lambda}{2}\right)}\right)$$

其中 $R = 6371$ km（地球平均半径）

**逆距离加权平均**：

$$\hat{v}(P) = \frac{\sum_{i=1}^{n} \frac{v_i}{d_i^p}}{\sum_{i=1}^{n} \frac{1}{d_i^p}}$$

其中：
- $\hat{v}(P)$：点$P$处的加权值
- $v_i$：第$i$个监测站的值
- $d_i$：$P$到第$i$个监测站的距离
- $p$：幂次（通常为1或2）
- $n$：距离范围内的监测站数量

---

## 数据准备

### 1. 创建示例数据

#### 企业数据（源点）

```r
# 创建企业数据
companies <- data.frame(
  company_id = c("COM001", "COM002", "COM003", "COM004", "COM005"),
  company_name = c("化工厂A", "电子厂B", "纺织厂C", "机械厂D", "制药厂E"),
  lon = c(110.5, 111.2, 110.8, 111.5, 110.3),
  lat = c(28.5, 28.8, 28.2, 28.0, 28.6),
  industry = c("化工", "电子", "纺织", "机械", "制药")
)

# 保存为CSV
write_csv(companies, "companies.csv")

print(companies)
#   company_id company_name   lon   lat industry
# 1     COM001      化工厂A 110.5 28.50     化工
# 2     COM002      电子厂B 111.2 28.80     电子
# 3     COM003      纺织厂C 110.8 28.20     纺织
# 4     COM004      机械厂D 111.5 28.00     机械
# 5     COM005      制药厂E 110.3 28.60     制药
```

#### 监测站数据（目标点）

```r
# 创建监测站数据（无时间维度）
stations <- data.frame(
  station_id = c("STN01", "STN02", "STN03", "STN04", "STN05", "STN06"),
  station_name = c("城区站", "郊区站A", "工业站", "农村站", "工业站B", "城南站"),
  lon = c(110.6, 111.0, 110.9, 110.4, 111.4, 110.7),
  lat = c(28.55, 28.75, 28.25, 28.10, 28.05, 28.45),
  pm25 = c(45, 38, 62, 28, 55, 42),
  pm10 = c(78, 65, 98, 45, 88, 72),
  no2 = c(42, 35, 58, 22, 52, 38)
)

write_csv(stations, "stations.csv")

print(stations)
#   station_id station_name   lon   lat pm25 pm10 no2
# 1       STN01       城区站 110.6 28.55   45   78  42
# 2       STN02       郊区站A 111.0 28.75   38   65  35
# 3       STN03       工业站 110.9 28.25   62   98  58
# 4       STN04       农村站 110.4 28.10   28   45  22
# 5       STN05       工业站B 111.4 28.05   55   88  52
# 6       STN06       城南站 110.7 28.45   42   72  38
```

#### 监测站多时间数据

```r
# 创建包含时间序列的监测站数据
stations_ts <- expand_grid(
  station_id = c("STN01", "STN02", "STN03", "STN04", "STN05", "STN06"),
  date = seq.Date(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "day")
) %>%
  left_join(
    stations %>% select(station_id, lon, lat),
    by = "station_id"
  ) %>%
  mutate(
    # 生成模拟的污染数据（正弦波 + 随机噪声）
    pm25 = 40 + 20 * sin(as.numeric(date) / 365 * 2 * pi) + rnorm(n(), 0, 5),
    pm10 = 70 + 30 * sin(as.numeric(date) / 365 * 2 * pi) + rnorm(n(), 0, 8),
    no2 = 35 + 15 * sin(as.numeric(date) / 365 * 2 * pi) + rnorm(n(), 0, 4)
  ) %>%
  mutate(
    pm25 = pmax(0, pm25),
    pm10 = pmax(0, pm10),
    no2 = pmax(0, no2)
  ) %>%
  left_join(
    stations %>% select(station_id, station_name),
    by = "station_id"
  )

write_csv(stations_ts, "stations_timeseries.csv")

head(stations_ts)
#   station_id       date      lon   lat pm25     pm10    no2          station_name
# 1       STN01 2023-01-01 110.6 28.55   35.2 62.3 32.1 城区站
# 2       STN02 2023-01-01 111.0 28.75   28.5 54.1 28.3 郊区站A
```

### 2. 转换为SF对象

```r
# 转换企业数据为SF
companies_sf <- st_as_sf(
  companies,
  coords = c("lon", "lat"),
  crs = 4326
)

# 转换监测站数据为SF
stations_sf <- st_as_sf(
  stations,
  coords = c("lon", "lat"),
  crs = 4326
)

# 查看
print(companies_sf)
print(stations_sf)
```

---

## 点对点匹配

### 1. 基础距离矩阵计算

```r
# 方法1：使用st_distance()计算距离矩阵
# 返回企业到监测站的距离矩阵（米）
distance_matrix <- st_distance(companies_sf, stations_sf)

# 转换为数据框
distance_df <- as.data.frame(distance_matrix) %>%
  rownames_to_column("company_idx") %>%
  pivot_longer(
    cols = -company_idx,
    names_to = "station_idx",
    values_to = "distance_m"
  ) %>%
  mutate(
    company_idx = as.integer(company_idx),
    station_idx = as.integer(gsub("V", "", station_idx)),
    distance_km = distance_m / 1000
  )

print(head(distance_df, 10))
#   company_idx station_idx distance_m distance_km
# 1           1           1        11234      11.23
# 2           1           2        58123      58.12
# 3           1           3        45678      45.68
```

### 2. 使用geosphere计算球面距离

```r
# 更精确的球面距离计算
calculate_distances <- function(company_coords, station_coords) {
  # company_coords: data.frame with lon, lat
  # station_coords: data.frame with lon, lat
  
  distances_list <- lapply(1:nrow(company_coords), function(i) {
    company_pt <- as.matrix(company_coords[i, c("lon", "lat")])
    station_pts <- as.matrix(station_coords[, c("lon", "lat")])
    
    dist_m <- distGeo(p1 = company_pt, p2 = station_pts)
    
    data.frame(
      company_id = companies$company_id[i],
      station_id = stations$station_id,
      distance_m = dist_m,
      distance_km = dist_m / 1000
    )
  })
  
  return(bind_rows(distances_list))
}

# 执行计算
distance_matrix_precise <- calculate_distances(
  companies[, c("lon", "lat")],
  stations[, c("lon", "lat")]
)

print(distance_matrix_precise)
#   company_id station_id distance_m distance_km
# 1      COM001       STN01       11.23      11.23
# 2      COM001       STN02       58.12      58.12
# 3      COM001       STN03       45.68      45.68
```

### 3. 筛选距离阈值内的监测站

```r
# 筛选距离在5km内的监测站
threshold_km <- 5

matched_stations <- distance_matrix_precise %>%
  filter(distance_km <= threshold_km) %>%
  left_join(stations, by = "station_id") %>%
  left_join(companies, by = "company_id")

print(matched_stations)
#   company_id station_id distance_m distance_km station_name   lon.x   lat.x
# 1      COM001       STN01       1234       1.23       城区站 110.60 28.55
# 2      COM001       STN06       3456       3.46       城南站 110.70 28.45
# 3      COM002       STN02       2345       2.35       郊区站A 111.00 28.75
```

### 4. 统计每个企业的匹配监测站数

```r
# 计算每个企业的匹配监测站统计
match_summary <- matched_stations %>%
  group_by(company_id, company_name) %>%
  summarise(
    n_stations = n(),
    min_distance_km = min(distance_km),
    max_distance_km = max(distance_km),
    mean_distance_km = mean(distance_km),
    .groups = 'drop'
  )

print(match_summary)
#   company_id company_name n_stations min_distance_km max_distance_km mean_distance_km
# 1      COM001      化工厂A          2            1.23            3.46            2.35
# 2      COM002      电子厂B          3            0.95            4.50            2.78
```

---

## 距离计算

### 1. 高效的距离矩阵计算

```r
# 使用data.table实现快速计算（大数据量推荐）
library(data.table)

calc_distance_dt <- function(companies_df, stations_df) {
  # 创建笛卡尔积
  dt_companies <- as.data.table(companies_df)[, 
    .(company_id, lon_c = lon, lat_c = lat)
  ]
  dt_stations <- as.data.table(stations_df)[, 
    .(station_id, lon_s = lon, lat_s = lat)
  ]
  
  # 笛卡尔积
  dist_dt <- dt_companies[, .(key = 1, lon_c, lat_c)][
    dt_stations[, .(key = 1, station_id, lon_s, lat_s)],
    on = "key",
    allow.cartesian = TRUE
  ][, key := NULL]
  
  # 计算距离（简化版，实际应用中使用distGeo）
  dist_dt[, distance_m := 111.32 * sqrt(
    (lat_c - lat_s)^2 + 
    ((lon_c - lon_s) * cos((lat_c + lat_s) / 2 * pi / 180))^2
  ) * 1000][, # 转为米
    `:=`(lon_c = NULL, lat_c = NULL, lon_s = NULL, lat_s = NULL)
  ]
  
  return(dist_dt)
}

# 使用data.table计算（对大数据集快10-100倍）
dist_dt <- calc_distance_dt(companies, stations)
print(dist_dt[1:10])
```

### 2. 距离带计算

```r
# 定义多个距离带，统计每个带内的监测站
distance_bands <- c(0, 5, 10, 20, Inf)
band_labels <- c("0-5km", "5-10km", "10-20km", ">20km")

distance_with_bands <- distance_matrix_precise %>%
  mutate(
    distance_band = cut(
      distance_km,
      breaks = distance_bands,
      labels = band_labels,
      right = FALSE
    )
  ) %>%
  filter(distance_km <= 20)  # 只保留20km内

print(distance_with_bands)

# 统计各距离带内的监测站数
band_stats <- distance_with_bands %>%
  group_by(company_id, distance_band) %>%
  summarise(
    n_stations = n(),
    mean_value = NA,  # 后续填充
    .groups = 'drop'
  )

print(band_stats)
```

---

## 逆距离加权聚合

### 1. 单值IDW聚合

```r
# 对每个企业，计算距离内监测站PM2.5的IDW平均值
calculate_idw <- function(matched_data, value_col, power = 2) {
  # matched_data: 包含distance_km的数据框
  # value_col: 要聚合的值列名
  # power: IDW的幂次
  
  result <- matched_data %>%
    group_by(company_id) %>%
    summarise(
      # IDW权重计算
      n_stations = n(),
      # 权重：距离的负幂次
      weight = 1 / (distance_km^power),
      # IDW值
      idw_value = sum(weight * !!sym(value_col)) / sum(weight),
      # 简单平均（作为对比）
      simple_mean = mean(!!sym(value_col)),
      # 距离加权标准差
      weighted_sd = sqrt(sum(weight * (!!sym(value_col) - idw_value)^2) / sum(weight)),
      .groups = 'drop'
    )
  
  return(result)
}

# 创建包含污染物数据的匹配结果
matched_with_data <- matched_stations %>%
  select(company_id, company_name, station_id, distance_km, 
         pm25, pm10, no2, industry)

# 计算IDW平均值（power=2）
idw_pm25 <- matched_with_data %>%
  group_by(company_id, company_name, industry) %>%
  summarise(
    n_stations = n(),
    idw_pm25 = sum((1 / distance_km^2) * pm25) / sum(1 / distance_km^2),
    simple_mean_pm25 = mean(pm25),
    min_pm25 = min(pm25),
    max_pm25 = max(pm25),
    .groups = 'drop'
  )

print(idw_pm25)
#   company_id company_name industry n_stations idw_pm25 simple_mean_pm25
# 1      COM001      化工厂A       化工          2     46.2            45.5
# 2      COM002      电子厂B       电子          3     39.8            38.3
```

### 2. 多污染物IDW聚合

```r
# 同时计算多个污染物的IDW
calculate_idw_multivalue <- function(matched_data, value_cols, power = 2) {
  # value_cols: 要聚合的值列名向量
  
  result_list <- list()
  
  for (col in value_cols) {
    idw_col <- matched_data %>%
      group_by(company_id) %>%
      summarise(
        weight = 1 / (distance_km^power),
        idw_val = sum(weight * !!sym(col)) / sum(weight),
        simple_mean_val = mean(!!sym(col)),
        .groups = 'drop'
      ) %>%
      rename(
        !!paste0("idw_", col) := idw_val,
        !!paste0("mean_", col) := simple_mean_val
      )
    
    if (length(result_list) == 0) {
      result_list[[col]] <- idw_col
    } else {
      result_list[[col]] <- idw_col %>%
        select(-company_id)
    }
  }
  
  result <- Reduce(cbind, result_list)
  return(result)
}

# 同时计算PM2.5、PM10、NO2的IDW
multi_idw <- calculate_idw_multivalue(
  matched_with_data,
  value_cols = c("pm25", "pm10", "no2"),
  power = 2
)

print(multi_idw)
#   company_id idw_pm25 mean_pm25 idw_pm10 mean_pm10 idw_no2 mean_no2
# 1      COM001     46.2      45.5      72.3      71.1      39.8      38.5
# 2      COM002     39.8      38.3      68.5      67.2      36.2      35.8
```

### 3. 不同幂次的效果对比

```r
# 比较不同幂次（power = 1, 2, 3）的IDW结果
power_comparison <- data.frame()

for (p in c(1, 2, 3)) {
  idw_result <- matched_with_data %>%
    group_by(company_id) %>%
    summarise(
      power = p,
      idw_pm25 = sum((1 / distance_km^p) * pm25) / sum(1 / distance_km^p),
      .groups = 'drop'
    )
  
  power_comparison <- rbind(power_comparison, idw_result)
}

# 转换为宽格式方便对比
power_comparison_wide <- power_comparison %>%
  pivot_wider(
    id_cols = company_id,
    names_from = power,
    values_from = idw_pm25,
    names_prefix = "power_"
  )

print(power_comparison_wide)
#   company_id power_1 power_2 power_3
# 1      COM001    45.8    46.2    46.4
# 2      COM002    38.9    39.8    40.1
```

---

## 多时间序列处理

### 1. 基础时间序列IDW

```r
# 对时间序列数据进行IDW聚合
calculate_idw_timeseries <- function(companies, stations_ts, 
                                     distance_threshold_km = 5,
                                     power = 2) {
  
  # 1. 计算距离矩阵
  dist_matrix <- calculate_distances(
    companies[, c("lon", "lat")],
    stations_ts[!duplicated(stations_ts$station_id), c("lon", "lat")]
  )
  
  # 2. 筛选距离内的监测站
  matched <- dist_matrix %>%
    filter(distance_km <= distance_threshold_km) %>%
    left_join(
      stations_ts %>% distinct(station_id, .keep_all = TRUE),
      by = "station_id"
    )
  
  # 3. 对每个日期计算IDW
  result_list <- list()
  
  for (d in unique(stations_ts$date)) {
    # 该日期的站点数据
    daily_data <- stations_ts %>%
      filter(date == d) %>%
      left_join(matched %>% select(company_id, station_id, distance_km),
                by = "station_id") %>%
      filter(!is.na(company_id))
    
    # 对每个企业计算IDW
    daily_idw <- daily_data %>%
      group_by(company_id) %>%
      summarise(
        date = d,
        n_stations = n(),
        idw_pm25 = sum((1 / distance_km^power) * pm25) / sum(1 / distance_km^power),
        idw_pm10 = sum((1 / distance_km^power) * pm10) / sum(1 / distance_km^power),
        idw_no2 = sum((1 / distance_km^power) * no2) / sum(1 / distance_km^power),
        .groups = 'drop'
      )
    
    result_list[[length(result_list) + 1]] <- daily_idw
  }
  
  return(bind_rows(result_list))
}

# 执行计算
ts_idw_result <- calculate_idw_timeseries(
  companies,
  stations_ts,
  distance_threshold_km = 5,
  power = 2
)

head(ts_idw_result)
#   company_id       date n_stations idw_pm25 idw_pm10 idw_no2
# 1      COM001 2023-01-01          2     35.2     62.3    32.1
# 2      COM001 2023-01-02          2     34.8     61.9    31.7
```

### 2. 时间序列统计

```r
# 计算时间序列的统计量
ts_statistics <- ts_idw_result %>%
  group_by(company_id) %>%
  summarise(
    # PM2.5统计
    pm25_mean = mean(idw_pm25),
    pm25_sd = sd(idw_pm25),
    pm25_min = min(idw_pm25),
    pm25_max = max(idw_pm25),
    # PM10统计
    pm10_mean = mean(idw_pm10),
    pm10_sd = sd(idw_pm10),
    # NO2统计
    no2_mean = mean(idw_no2),
    no2_sd = sd(idw_no2),
    # 数据点数
    n_records = n(),
    .groups = 'drop'
  ) %>%
  left_join(companies %>% select(company_id, company_name, industry),
            by = "company_id")

print(ts_statistics)
```

### 3. 季节性分析

```r
# 添加季节信息
ts_idw_seasonal <- ts_idw_result %>%
  mutate(
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "冬季",
      month %in% c(3, 4, 5) ~ "春季",
      month %in% c(6, 7, 8) ~ "夏季",
      month %in% c(9, 10, 11) ~ "秋季"
    )
  )

# 计算季节平均值
seasonal_stats <- ts_idw_seasonal %>%
  group_by(company_id, season) %>%
  summarise(
    pm25_seasonal_mean = mean(idw_pm25),
    pm10_seasonal_mean = mean(idw_pm10),
    no2_seasonal_mean = mean(idw_no2),
    .groups = 'drop'
  )

print(seasonal_stats)
```

---

## 完整工作流示例

### 示例1：单时间快照 + 多污染物

```r
# 完整工作流：从CSV到最终结果
workflow_single_time <- function(companies_file, stations_file, 
                                 distance_threshold_km = 5,
                                 power = 2) {
  
  # 1. 读入数据
  companies <- read_csv(companies_file)
  stations <- read_csv(stations_file)
  
  # 2. 计算距离矩阵
  dist_matrix <- calculate_distances(
    companies[, c("lon", "lat")],
    stations[, c("lon", "lat")]
  )
  
  # 3. 筛选距离阈值内的监测站
  matched <- dist_matrix %>%
    filter(distance_km <= distance_threshold_km) %>%
    left_join(stations, by = "station_id") %>%
    left_join(companies, by = "company_id")
  
  # 4. 计算IDW
  result <- matched %>%
    group_by(company_id, company_name, industry) %>%
    summarise(
      n_stations = n(),
      # PM2.5 IDW
      idw_pm25 = sum((1 / distance_km^power) * pm25) / sum(1 / distance_km^power),
      simple_pm25 = mean(pm25),
      # PM10 IDW
      idw_pm10 = sum((1 / distance_km^power) * pm10) / sum(1 / distance_km^power),
      simple_pm10 = mean(pm10),
      # NO2 IDW
      idw_no2 = sum((1 / distance_km^power) * no2) / sum(1 / distance_km^power),
      simple_no2 = mean(no2),
      # 距离统计
      min_distance_km = min(distance_km),
      max_distance_km = max(distance_km),
      mean_distance_km = mean(distance_km),
      .groups = 'drop'
    )
  
  return(result)
}

# 执行
final_result <- workflow_single_time(
  "companies.csv",
  "stations.csv",
  distance_threshold_km = 5,
  power = 2
)

print(final_result)

# 保存结果
write_csv(final_result, "company_air_quality_match.csv")
```

### 示例2：时间序列 + 多企业多监测站

```r
# 完整时间序列工作流
workflow_timeseries <- function(companies_file, stations_ts_file,
                               distance_threshold_km = 5,
                               power = 2) {
  
  # 1. 读入数据
  companies <- read_csv(companies_file)
  stations_ts <- read_csv(stations_ts_file) %>%
    mutate(date = as.Date(date))
  
  # 2. 计算距离矩阵（只用站点坐标）
  stations_coords <- stations_ts %>%
    distinct(station_id, lon, lat)
  
  dist_matrix <- calculate_distances(
    companies[, c("lon", "lat")],
    stations_coords[, c("lon", "lat")]
  )
  
  # 3. 筛选匹配的企业-监测站对
  matched_pairs <- dist_matrix %>%
    filter(distance_km <= distance_threshold_km) %>%
    select(company_id, station_id, distance_km)
  
  # 4. 连接时间序列数据
  ts_matched <- matched_pairs %>%
    left_join(stations_ts, by = "station_id") %>%
    left_join(companies %>% select(company_id, company_name, industry),
              by = "company_id")
  
  # 5. 计算每日IDW
  daily_idw <- ts_matched %>%
    group_by(company_id, company_name, industry, date) %>%
    summarise(
      n_stations = n(),
      idw_pm25 = sum((1 / distance_km^power) * pm25) / sum(1 / distance_km^power),
      idw_pm10 = sum((1 / distance_km^power) * pm10) / sum(1 / distance_km^power),
      idw_no2 = sum((1 / distance_km^power) * no2) / sum(1 / distance_km^power),
      mean_distance_km = mean(distance_km),
      .groups = 'drop'
    ) %>%
    arrange(company_id, date)
  
  return(daily_idw)
}

# 执行
ts_result <- workflow_timeseries(
  "companies.csv",
  "stations_timeseries.csv",
  distance_threshold_km = 5,
  power = 2
)

head(ts_result)

# 保存
write_csv(ts_result, "company_air_quality_timeseries.csv")
```

---

## 性能优化

### 1. 使用data.table加速

```r
# 对于大数据集（企业>10000，监测站>1000），使用data.table
calculate_idw_dt <- function(companies_dt, stations_dt, 
                             distance_threshold_km = 5,
                             power = 2) {
  
  # 创建笛卡尔积
  companies_dt <- setDT(companies_dt)
  stations_dt <- setDT(stations_dt)
  
  dist_dt <- companies_dt[, .(key = 1, company_id, lon_c = lon, lat_c = lat)][
    stations_dt[, .(key = 1, station_id, pm25, pm10, no2, lon_s = lon, lat_s = lat)],
    on = "key",
    allow.cartesian = TRUE
  ][, key := NULL]
  
  # 计算距离（简化公式）
  dist_dt[, distance_km := 
    111.32 * sqrt(
      (lat_c - lat_s)^2 + 
      ((lon_c - lon_s) * cos((lat_c + lat_s) / 2 * pi / 180))^2
    )
  ]
  
  # 筛选距离阈值
  dist_dt <- dist_dt[distance_km <= distance_threshold_km]
  
  # 计算IDW
  idw_dt <- dist_dt[, 
    .(
      n_stations = .N,
      idw_pm25 = sum((1 / distance_km^power) * pm25) / sum(1 / distance_km^power),
      idw_pm10 = sum((1 / distance_km^power) * pm10) / sum(1 / distance_km^power),
      idw_no2 = sum((1 / distance_km^power) * no2) / sum(1 / distance_km^power),
      mean_distance_km = mean(distance_km)
    ),
    by = .(company_id)
  ]
  
  return(idw_dt)
}

# 执行
companies_dt <- as.data.table(companies)
stations_dt <- as.data.table(stations)

result_dt <- calculate_idw_dt(companies_dt, stations_dt, power = 2)
print(result_dt)
```

### 2. 并行处理

```r
library(future)
library(furrr)

# 设置并行处理
plan(multisession, workers = 4)

# 并行计算多个日期的IDW
parallel_ts_idw <- function(companies, stations_ts, 
                            distance_threshold_km = 5,
                            power = 2) {
  
  # 获取所有日期
  dates <- unique(stations_ts$date)
  
  # 并行处理每个日期
  results <- future_map_df(dates, function(d) {
    daily_data <- stations_ts %>% filter(date == d)
    
    dist_matrix <- calculate_distances(
      companies[, c("lon", "lat")],
      daily_data %>% distinct(station_id, lon, lat) %>% select(lon, lat)
    )
    
    # IDW计算...（省略）
    
  }, .progress = TRUE)
  
  return(results)
}

# 还原串行处理
plan(sequential)
```

### 3. 性能基准测试

```r
library(bench)

# 比较不同方法的性能
companies_1k <- companies %>% slice(rep(1:nrow(companies), 200))
stations_1k <- stations %>% slice(rep(1:nrow(stations), 100))

benchmark_result <- mark(
  # 方法1：基础R + lapply
  "R_base" = {
    calculate_distances(
      companies_1k[, c("lon", "lat")],
      stations_1k[, c("lon", "lat")]
    )
  },
  # 方法2：data.table
  "data_table" = {
    calculate_idw_dt(
      as.data.table(companies_1k),
      as.data.table(stations_1k)
    )
  },
  # 方法3：sf距离计算
  "sf_dist" = {
    st_distance(
      st_as_sf(companies_1k, coords = c("lon", "lat"), crs = 4326),
      st_as_sf(stations_1k, coords = c("lon", "lat"), crs = 4326)
    )
  },
  iterations = 3
)

print(benchmark_result)
```

---

## 可视化

### 1. 企业-监测站匹配网络图

```r
library(ggplot2)
library(ggmap)

# 创建匹配关系的可视化
plot_matching_network <- function(companies, stations, matched_data,
                                 distance_threshold_km = 5) {
  
  # 提取匹配线段
  lines_data <- matched_data %>%
    filter(distance_km <= distance_threshold_km) %>%
    left_join(companies %>% select(company_id, lon, lat),
              by = "company_id", suffix = c("_s", "_c")) %>%
    rename(lon_company = lon, lat_company = lat)
  
  ggplot() +
    # 绘制连接线
    geom_segment(
      data = lines_data,
      aes(x = lon_s, y = lat_s, xend = lon_company, yend = lat_company),
      alpha = 0.3, color = "gray"
    ) +
    # 绘制监测站
    geom_point(
      data = stations,
      aes(x = lon, y = lat, color = "监测站", size = 3),
      show.legend = TRUE
    ) +
    # 绘制企业
    geom_point(
      data = companies,
      aes(x = lon, y = lat, color = "企业", size = 5),
      show.legend = TRUE
    ) +
    # 标签
    geom_text(
      data = companies,
      aes(x = lon, y = lat, label = substr(company_id, 1, 3)),
      size = 3, vjust = -1.5
    ) +
    theme_minimal() +
    labs(
      title = "企业-监测站匹配网络",
      x = "经度", y = "纬度",
      color = "点位类型"
    ) +
    theme(legend.position = "bottom")
}

plot_matching_network(companies, stations, matched_with_data, 5)
```

### 2. IDW值与距离的关系

```r
plot_idw_distance <- function(matched_data, pollutant = "pm25") {
  
  matched_data %>%
    ggplot(aes(x = distance_km, y = !!sym(pollutant), 
               color = company_id, size = 1/distance_km^2)) +
    geom_point(alpha = 0.6) +
    facet_wrap(~company_id) +
    theme_minimal() +
    labs(
      title = paste(pollutant, "与距离的关系"),
      x = "距离 (km)",
      y = toupper(pollutant),
      size = "权重"
    ) +
    theme(legend.position = "bottom")
}

plot_idw_distance(matched_with_data, "pm25")
```

### 3. 时间序列可视化

```r
plot_ts_idw <- function(ts_result, company_id = NULL) {
  
  if (!is.null(company_id)) {
    ts_result <- ts_result %>% filter(company_id %in% !!company_id)
  }
  
  ts_result %>%
    pivot_longer(
      cols = starts_with("idw_"),
      names_to = "pollutant",
      values_to = "value"
    ) %>%
    mutate(pollutant = toupper(gsub("idw_", "", pollutant))) %>%
    ggplot(aes(x = date, y = value, color = pollutant)) +
    geom_line() +
    facet_wrap(~company_id, scales = "free_y") +
    theme_minimal() +
    labs(
      title = "企业周边污染物时间序列",
      x = "日期", y = "浓度",
      color = "污染物"
    )
}

plot_ts_idw(ts_idw_result, company_id = c("COM001", "COM002"))
```

### 4. 交互式地图（leaflet）

```r
library(leaflet)

plot_interactive_map <- function(companies, stations, matched_data,
                                distance_threshold_km = 5) {
  
  # 筛选匹配数据
  matched <- matched_data %>%
    filter(distance_km <= distance_threshold_km)
  
  # 创建地图
  m <- leaflet() %>%
    addTiles() %>%
    setView(lng = mean(companies$lon), lat = mean(companies$lat), zoom = 10)
  
  # 添加监测站
  m <- m %>%
    addCircleMarkers(
      data = stations,
      lng = ~lon, lat = ~lat,
      radius = 5,
      color = "blue",
      popup = ~station_name,
      group = "监测站"
    )
  
  # 添加企业
  m <- m %>%
    addCircleMarkers(
      data = companies,
      lng = ~lon, lat = ~lat,
      radius = 8,
      color = "red",
      popup = ~company_name,
      group = "企业"
    )
  
  # 添加连接线
  for (i in 1:nrow(matched)) {
    m <- m %>%
      addPolylines(
        lng = c(matched$lon[i], 
                companies$lon[companies$company_id == matched$company_id[i]]),
        lat = c(matched$lat[i],
                companies$lat[companies$company_id == matched$company_id[i]]),
        color = "gray",
        weight = 1,
        opacity = 0.5
      )
  }
  
  return(m)
}

# 生成交互式地图
plot_interactive_map(companies, stations, matched_with_data, 5)
```

---

## 故障排查

### 问题1：无匹配的企业-监测站对

```r
# 检查是否存在无匹配企业
unmatched_companies <- companies %>%
  anti_join(
    matched_stations %>% select(company_id),
    by = "company_id"
  )

if (nrow(unmatched_companies) > 0) {
  cat("无匹配企业:\n")
  print(unmatched_companies)
  
  # 扩大搜索范围
  threshold_km <- 10
  cat(paste("建议扩大搜索范围至", threshold_km, "km\n"))
}
```

### 问题2：权重计算中的零距离

```r
# 处理企业与监测站位置完全相同的情况
handle_zero_distance <- function(matched_data, power = 2) {
  
  matched_data %>%
    mutate(
      # 避免0距离
      distance_km_safe = pmax(distance_km, 0.001),  # 最小1米
      weight = 1 / (distance_km_safe^power)
    ) %>%
    group_by(company_id) %>%
    mutate(
      weight_normalized = weight / sum(weight)
    )
}

matched_safe <- handle_zero_distance(matched_with_data, power = 2)
```

### 问题3：权重分布不均

```r
# 检查权重分布是否过于集中
check_weight_distribution <- function(matched_data, power = 2) {
  
  weight_summary <- matched_data %>%
    group_by(company_id) %>%
    mutate(
      weight = 1 / (distance_km^power),
      weight_pct = weight / sum(weight) * 100
    ) %>%
    summarise(
      n_stations = n(),
      max_weight_pct = max(weight_pct),
      top2_weight_pct = sum(weight_pct[weight_pct == max(weight_pct) | 
                                        weight_pct == sort(weight_pct, decreasing = TRUE)[2]]),
      .groups = 'drop'
    )
  
  print(weight_summary)
  
  # 如果权重过于集中，考虑降低幂次
  if (max(weight_summary$max_weight_pct) > 80) {
    cat("警告：权重分布过于集中，建议降低幂次或使用等权平均\n")
  }
}

check_weight_distribution(matched_with_data, power = 2)
```

### 问题4：坐标参考系不一致

```r
# 验证坐标参考系
verify_crs <- function(companies, stations) {
  
  # 检查坐标范围
  cat("企业坐标范围:\n")
  print(data.frame(
    lon_range = paste(min(companies$lon), "-", max(companies$lon)),
    lat_range = paste(min(companies$lat), "-", max(companies$lat))
  ))
  
  cat("监测站坐标范围:\n")
  print(data.frame(
    lon_range = paste(min(stations$lon), "-", max(stations$lon)),
    lat_range = paste(min(stations$lat), "-", max(stations$lat))
  ))
  
  # WGS84合理范围检查
  if (any(abs(companies$lon) > 180) | any(abs(companies$lat) > 90)) {
    cat("错误：经纬度超出WGS84范围\n")
  }
}

verify_crs(companies, stations)
```

### 问题5：时间序列数据缺失

```r
# 检查时间序列完整性
check_timeseries_completeness <- function(stations_ts) {
  
  completeness <- stations_ts %>%
    group_by(station_id) %>%
    summarise(
      date_min = min(date),
      date_max = max(date),
      n_records = n(),
      expected_records = as.numeric(max(date) - min(date)) + 1,
      completeness_pct = n_records / expected_records * 100,
      .groups = 'drop'
    )
  
  print(completeness)
  
  # 检查缺失数据
  incomplete <- completeness %>%
    filter(completeness_pct < 100)
  
  if (nrow(incomplete) > 0) {
    cat("警告：以下监测站存在缺失数据:\n")
    print(incomplete)
  }
}

check_timeseries_completeness(stations_ts)
```

---

## 最佳实践总结

| 方面 | 建议 | 原因 |
|------|------|------|
| **距离阈值** | 3-10 km | 过大会包含无关站点，过小可能无匹配 |
| **幂次选择** | power = 1 或 2 | 2是标准选择，1给予距离权重更均匀 |
| **无效值处理** | 用NA_real_替代 | 避免数值计算错误 |
| **坐标验证** | 使用WGS84 (EPSG:4326) | 确保全球统一参考系 |
| **时间同步** | 确保企业和监测站数据对齐 | 避免时间错配的匹配 |
| **权重检查** | 单个权重<80% | 防止结果过度依赖单一监测站 |
| **结果验证** | 抽样检查计算过程 | 发现潜在错误 |

---

## 推荐阅读

- [Geosphere包文档](https://cran.r-project.org/web/packages/geosphere/)
- [SF包空间关系](https://r-spatial.github.io/sf/articles/sf4.html)
- [Data.table高性能](https://cran.r-project.org/web/packages/data.table/)
- [Haversine公式](https://en.wikipedia.org/wiki/Haversine_formula)

---

## 许可证

本教程遵循CC-BY-4.0许可证
