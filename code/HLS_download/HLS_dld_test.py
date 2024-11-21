"""
_summary_
"""

# Load packages
import os
from datetime import datetime
import numpy as np
import pandas as pd
import geopandas as gp
from skimage import io
import matplotlib.pyplot as plt
from osgeo import gdal
import xarray as xr
import rioxarray as rxr
import hvplot.xarray
import hvplot.pandas
import earthaccess

# Configure Earthaccess
auth = earthaccess.Auth()
auth.login(strategy = "netrc")
print(auth.authenticated)

# Load AOIs
aois = gp.read_file("../data/feature_layers/aois_test.geojson")

bbox = tuple(list(aois.total_bounds))

# Define time window
temporal = ("2020-06-01T00:00:00", "2020-06-30T23:59:59")

# Get results
max_img = 100 # nr. of maximum returned images
results = earthaccess.search_data(
    short_name=['HLSL30','HLSS30'],
    bounding_box=bbox,
    temporal=temporal,
    count=max_img
)

pd.json_normalize(results).head(5)

# retrieve URLS
print("Retrieving granules...")
hls_results_urls = [granule.data_links() for granule in results]

browse_urls = [granule.dataviz_links()[0] for granule in results] # 0 retrieves only the https links

# GDAL configurations used to successfully access LP DAAC Cloud Assets via vsicurl 
gdal.SetConfigOption('GDAL_HTTP_COOKIEFILE','~/cookies.txt')
gdal.SetConfigOption('GDAL_HTTP_COOKIEJAR', '~/cookies.txt')
gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN','EMPTY_DIR')
gdal.SetConfigOption('CPL_VSIL_CURL_ALLOWED_EXTENSIONS','TIF')
gdal.SetConfigOption('GDAL_HTTP_UNSAFESSL', 'YES')
gdal.SetConfigOption('GDAL_HTTP_MAX_RETRY', '10')
gdal.SetConfigOption('GDAL_HTTP_RETRY_DELAY', '0.5')

# List of returned granules from URLs
h = hls_results_urls[10]

evi_band_links = []

# Define which HLS product is being accessed
if h[0].split('/')[4] == 'HLSS30.020':
    evi_bands = ['B8A', 'B04', 'B02', 'Fmask'] # NIR RED BLUE for S30
else:
    evi_bands = ['B05', 'B04', 'B02', 'Fmask'] # NIR RED BLUE for L30

# Subset the assets in the item down to only the desired bands
for a in h: 
    if any(b in a for b in evi_bands):
        evi_band_links.append(a)

# Use vsicurl to load the data directly into memory (be patient, may take a few seconds)
chunk_size = dict(band=1, x=512, y=512) # Tiles have 1 band and are divided into 512x512 pixel chunks
for e in evi_band_links:
    print(e)
    # Open and build datasets
    if e.rsplit('.', 2)[-2] == evi_bands[0]:      # NIR index
        nir = rxr.open_rasterio(e, chunks=chunk_size, masked=True).squeeze('band', drop=True)
        nir.attrs['scale_factor'] = 0.0001        # hard coded the scale_factor attribute 
    elif e.rsplit('.', 2)[-2] == evi_bands[1]:    # red index
        red = rxr.open_rasterio(e, chunks=chunk_size, masked=True).squeeze('band', drop=True)
        red.attrs['scale_factor'] = 0.0001        # hard coded the scale_factor attribute
    elif e.rsplit('.', 2)[-2] == evi_bands[2]:    # blue index
        blue = rxr.open_rasterio(e, chunks=chunk_size, masked=True).squeeze('band', drop=True)
        blue.attrs['scale_factor'] = 0.0001       # hard coded the scale_factor attribute
    elif e.rsplit('.', 2)[-2] == evi_bands[3]:
        fmask = rxr.open_rasterio(e, chunks=chunk_size, masked=True).squeeze('band', drop=True)
        # No Scaling for Fmask

print("The COGs have been loaded into memory!")

nir.spatial_ref.crs_wkt