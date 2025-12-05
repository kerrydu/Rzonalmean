#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
读取TIF和SHP，计算湖南各市平均灯光亮度（Python版）
功能与 R 脚本等价：
- 打印TIF元信息（文件、层数、尺寸、分辨率、范围、CRS、dtype、nodata、每层min/max）
- 自动识别城市字段（常见字段名；否则选择第一个字符/因子列；再否则第一个非几何列）
- 支持 band_mode: N, first, last(默认), avg, all
- 计算分市平均灯光亮度，写出 CSV
依赖：rasterio, geopandas, shapely, numpy, pandas, rasterstats
"""

import argparse
import os
from typing import List, Optional

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio import features
import geopandas as gpd
from rasterstats import zonal_stats


def print_raster_meta(src: rasterio.io.DatasetReader, path_hint: Optional[str] = None) -> None:
    try:
        sources = getattr(src, 'files', None) or [src.name]
    except Exception:
        sources = [src.name]
    fn = ", ".join([s for s in sources if s]) if sources else (path_hint or "(内存栅格)")

    count = src.count
    width, height = src.width, src.height
    res = src.res
    bounds = src.bounds
    crs = src.crs
    dtypes = [src.dtypes[i - 1] for i in range(1, count + 1)]
    nodata = src.nodata

    print("\n==== TIF元信息 ====")
    print(f"文件: {fn}")
    print(f"层数: {count}")
    print(f"尺寸: {height} 行 x {width} 列")
    print(f"分辨率: {res}")
    print(f"范围(extent): xmin={bounds.left}, xmax={bounds.right}, ymin={bounds.bottom}, ymax={bounds.top}")
    print(f"CRS: {crs or '(缺失)'}")
    print(f"数据类型: {', '.join(map(str, dtypes))}")
    print(f"NA值标记: {nodata}")

    # 每层min/max（可能较慢）
    print("范围(每层):")
    for i in range(1, count + 1):
        arr = src.read(i, masked=True)
        # 把masked值设为nan以便计算
        a = np.where(arr.mask, np.nan, arr.data)
        mn = float(np.nanmin(a)) if np.isnan(a).sum() < a.size else np.nan
        mx = float(np.nanmax(a)) if np.isnan(a).sum() < a.size else np.nan
        print(f"  band_{i}: min={mn}, max={mx}")
    print("===================\n")


def choose_city_field(gdf: gpd.GeoDataFrame) -> str:
    nms = list(gdf.columns)
    # 常见字段名候选（不区分大小写）
    candidates = [
        "CITY", "CITY_NAME", "NAME", "NAME_1", "NAME_2",
        "ADMIN_NAME", "ADM2_NAME", "COUNTY", "PREF", "PREFECTURE",
        "CITY_CN", "NAME_CH", "NAME_ZH", "ShiName"
    ]
    upper = [c.upper() for c in nms]
    for cand in candidates:
        if cand in upper:
            return nms[upper.index(cand)]
    # 字符型优先
    for nm in nms:
        if nm == gdf.geometry.name:
            continue
        vals = gdf[nm]
        if pd.api.types.is_string_dtype(vals) or pd.api.types.is_object_dtype(vals):
            return nm
    # 退回第一个非几何列
    for nm in nms:
        if nm != gdf.geometry.name:
            return nm
    raise ValueError("未能识别城市名称字段")


def band_selection(src: rasterio.io.DatasetReader, mode: str) -> np.ndarray:
    """根据band_mode返回二维数组；若mode='all'返回三维数组(count, h, w)"""
    count = src.count
    if mode is None or mode == "":
        mode = "last"
    m = mode.lower()
    if m.isdigit():
        bi = int(m)
        if bi < 1 or bi > count:
            raise ValueError(f"band索引超出范围: 1..{count}")
        arr = src.read(bi, masked=True)
        return np.where(arr.mask, np.nan, arr.data)
    elif m == "first":
        arr = src.read(1, masked=True)
        return np.where(arr.mask, np.nan, arr.data)
    elif m == "last":
        arr = src.read(count, masked=True)
        return np.where(arr.mask, np.nan, arr.data)
    elif m == "avg":
        stack = src.read(np.arange(1, count + 1), masked=True)  # (count, h, w)
        data = np.where(stack.mask, np.nan, stack.data)
        return np.nanmean(data, axis=0)
    elif m == "all":
        stack = src.read(np.arange(1, count + 1), masked=True)
        return np.where(stack.mask, np.nan, stack.data)  # (count, h, w)
    else:
        raise ValueError(f"未知band_mode: {mode}")


def compute_zonal_means(tif_path: str, shp_path: str, out_csv: str, band_mode: str = "last", city_field: Optional[str] = None) -> pd.DataFrame:
    if not os.path.exists(tif_path):
        raise FileNotFoundError(f"找不到TIF: {tif_path}")
    if not os.path.exists(shp_path):
        raise FileNotFoundError(f"找不到SHP: {shp_path}")

    gdf = gpd.read_file(shp_path)
    with rasterio.open(tif_path) as src:
        print_raster_meta(src, path_hint=tif_path)
        # 坐标系统一
        if gdf.crs is not None and src.crs is not None and gdf.crs != src.crs:
            gdf = gdf.to_crs(src.crs)
        # 选择城市字段
        if not city_field:
            city_field = choose_city_field(gdf)
        print(f"使用城市字段: {city_field}")

        # 选择或聚合band
        arr = band_selection(src, band_mode)
        transform = src.transform
        nodata = src.nodata

        # zonal stats
        if arr.ndim == 2:  # 单层
            zs = zonal_stats(
                gdf.geometry,
                arr,
                affine=transform,
                nodata=nodata,
                stats=["mean"],
                all_touched=False
            )
            means = [z["mean"] if z and ("mean" in z) else np.nan for z in zs]
            df = pd.DataFrame({"city": gdf[city_field].astype(str), "mean_light": means})
        else:  # 多层（all）
            band_means = {}
            for i in range(arr.shape[0]):
                zs = zonal_stats(
                    gdf.geometry,
                    arr[i, :, :],
                    affine=transform,
                    nodata=nodata,
                    stats=["mean"],
                    all_touched=False
                )
                band_means[f"band_{i+1}"] = [z["mean"] if z and ("mean" in z) else np.nan for z in zs]
            df = pd.DataFrame({"city": gdf[city_field].astype(str)})
            for k, v in band_means.items():
                df[k] = v
            df["mean_light"] = df[[c for c in df.columns if c.startswith("band_")]].mean(axis=1, skipna=True)
        df = df.sort_values("city").reset_index(drop=True)
    # 写出
    out_dir = os.path.dirname(out_csv) or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(out_csv, index=False, encoding='utf-8-sig')
    print(f"已写出结果: {out_csv}")
    print(df.head(10))
    return df


def main():
    parser = argparse.ArgumentParser(description="计算湖南各市平均灯光亮度（Python）")
    parser.add_argument("tif_path", help="GeoTIFF路径，例如 DMSP-like2020.tif")
    parser.add_argument("shp_path", help="湖南行政区划Shapefile路径，例如 hunan.shp")
    parser.add_argument("out_csv", nargs="?", default=None, help="输出CSV路径，默认与SHP同目录下的 hunan_light_mean.csv")
    parser.add_argument("band_mode", nargs="?", default="last", help="band模式: 数字N/first/last/avg/all，默认last")
    parser.add_argument("city_field", nargs="?", default=None, help="可选：显式指定城市字段名")
    args = parser.parse_args()

    out_csv = args.out_csv or os.path.join(os.path.dirname(args.shp_path), "hunan_light_mean.csv")
    compute_zonal_means(args.tif_path, args.shp_path, out_csv, band_mode=args.band_mode, city_field=args.city_field)


if __name__ == "__main__":
    main()
