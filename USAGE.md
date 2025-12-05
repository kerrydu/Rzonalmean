# Hunan Zonal Mean – Usage Guide

This folder contains scripts to compute per-city mean nightlight brightness for Hunan using a GeoTIFF (e.g., `DMSP-like2020.tif`) and a shapefile (`hunan.shp`). Three options are provided: R, Python, and Google Earth Engine (GEE).

- Raster: `DMSP-like2020.tif` (single or multi-band)
- Vector: `hunan.shp` (+ .dbf/.prj/.shx)
- Output: CSV with columns like `city` and `mean_light` (plus per-band means for some modes)

## R Script (`zonal_mean_hunan.R`)

- Prints TIFF metadata (bands, size, resolution, extent, CRS, dtype, NA, per-band min/max).
- Auto-detects a city name field; if incorrect, edit the script to fix the column name.
- Reprojects shapefile to raster CRS when needed.

Install packages (once):

```powershell
R -e "install.packages(c('terra','sf','dplyr','readr'), repos='https://cloud.r-project.org')"
```

Run examples (from this folder):

```powershell
# Default: use last band; output defaults to hunan_light_mean.csv next to the shapefile
Rscript ".\zonal_mean_hunan.R" ".\DMSP-like2020.tif" ".\hunan.shp"

# Explicit output path
Rscript ".\zonal_mean_hunan.R" ".\DMSP-like2020.tif" ".\hunan.shp" ".\output.csv"

# Band modes (optional 4th arg):
#   number N -> Nth band
#   first    -> band 1
#   last     -> last band (default)
#   avg      -> average all bands to a single layer
Rscript ".\zonal_mean_hunan.R" ".\DMSP-like2020.tif" ".\hunan.shp" ".\output.csv" avg

### NetCDF (.nc)

When using a NetCDF raster (e.g., `hunan.nc`):

- 5th argument `nc_var` (optional): select variable/subdataset by name or 1-based index. If omitted, the script uses the first subdataset (if any) or reads the file as a multi-layer time series.
- 6th argument `time_sel` (optional, for time-series): substring to match a specific time layer (e.g., `2020-01`, `2020-01-15`). The script tries to match against the time strings and selects the first match.

```powershell
# Use first variable by default (prints available subdatasets)
Rscript ".\zonal_mean_hunan.R" ".\hunan.nc" ".\hunan.shp" ".\output.csv" last

# Specify variable name (recommended if multiple variables exist)
Rscript ".\zonal_mean_hunan.R" ".\hunan.nc" ".\hunan.shp" ".\output.csv" avg varname

# Specify variable by index (1-based)
Rscript ".\zonal_mean_hunan.R" ".\hunan.nc" ".\hunan.shp" ".\output.csv" first 2

# Select a specific time layer by substring (e.g., January 2020)
Rscript ".\zonal_mean_hunan.R" ".\hunan.nc" ".\hunan.shp" ".\output.csv" last varname 2020-01
```
```

Notes:
- If the script picked the wrong city column, open the script and set `city_field` explicitly after the shapefile is read.
- Metadata prints right after reading the raster.

## Python Script (`zonal_mean_hunan.py`)

- Prints GeoTIFF metadata as above.
- Auto-detects the city field (includes common names like `ShiName`, `NAME`, etc.); you can also pass it as an argument.
- Band modes: `N`, `first`, `last` (default), `avg` (pixel-wise average across bands to one layer), `all` (keep all bands and output per-band means plus overall mean).
- Handles CRS alignment.

Recommended: create a clean conda-forge environment (Windows, avoids GDAL DLL issues):

```powershell
conda config --set channel_priority strict
conda create -n geozones311 -c conda-forge -y python=3.11 gdal rasterio fiona geopandas rasterstats pandas numpy
conda activate geozones311
python -c "import rasterio, geopandas; print('OK', rasterio.__version__)"
```

Run examples (from this folder):

```powershell
# Default (last band), output defaults to hunan_light_mean.csv next to the shapefile
python ".\zonal_mean_hunan.py" ".\DMSP-like2020.tif" ".\hunan.shp"

# Explicit output path
python ".\zonal_mean_hunan.py" ".\DMSP-like2020.tif" ".\hunan.shp" ".\output.csv"

# Use band 1
python ".\zonal_mean_hunan.py" ".\DMSP-like2020.tif" ".\hunan.shp" ".\output.csv" 1

# Average all bands to one layer
python ".\zonal_mean_hunan.py" ".\DMSP-like2020.tif" ".\hunan.shp" ".\output.csv" avg

# Keep all bands and output per-band means + mean_light
python ".\zonal_mean_hunan.py" ".\DMSP-like2020.tif" ".\hunan.shp" ".\output.csv" all

# Explicit city field (e.g., ShiName) as 5th arg
python ".\zonal_mean_hunan.py" ".\DMSP-like2020.tif" ".\hunan.shp" ".\output.csv" last ShiName
```

Troubleshooting GDAL/rasterio DLLs on Windows:
- Prefer a fresh `conda-forge` environment with `channel_priority=strict`.
- Avoid mixing `pip` and `conda` for `rasterio/gdal`.
- If imports fail with `_vsiopener` or missing DLLs, reinstall `gdal` and `rasterio` from `conda-forge` in the same environment.

## GEE Script (`zonal_mean_hunan_gee.js`)

- Runs in Google Earth Engine Code Editor: https://code.earthengine.google.com/
- Upload assets:
  - Image upload: `dmsp-like2020.tif` → Image asset ID (e.g., `users/you/dmsp_like2020`).
  - Table upload: `hunan.shp` zipped with `.dbf/.prj/.shx` → FeatureCollection asset ID (e.g., `users/you/hunan`).
- Set at top of script:
  - `ASSET_IMAGE`, `ASSET_HUNAN`
  - Optional: `CITY_FIELD`, `BAND_MODE` (`N`, `first`, `last`, `avg`, `all`)
  - `EXPORT_TO_DRIVE`, `EXPORT_FILE_NAME`, `EXPORT_FOLDER`

Run:
- Paste the script into the Code Editor, set asset IDs, click Run.
- Preview prints in the Console; an export task appears in Tasks if `EXPORT_TO_DRIVE=true`.

## Band Modes Summary
- `N`: use Nth band (1-based)
- `first`: band 1
- `last`: last band (R/Python default)
- `avg`: average across all bands to a single-band image, then zonal mean
- `all` (Python/GEE): keep all bands, output per-band means plus `mean_light` as row average

## Output
- CSV with at least: `city`, `mean_light`
- When `all`: additional columns per band (e.g., `band_1`, `band_2`, ...)

## Notes
- CRS is aligned by reprojecting the vector to raster CRS (R/Python); GEE handles projection internally.
- NA/nodata is ignored in means.
- Polygons: default uses pixel center inclusion (`all_touched=False` in Python). Adjust if needed.
