conda list gdal
conda list rasterio
conda list proj
conda list fiona
conda list libgdal
conda config --show channel_priority
where gdal*.dll
where proj*.dll
python -c "import rasterio, sys; print('rasterio file:', rasterio.__file__); import fiona; print('fiona file:', fiona.__file__)"