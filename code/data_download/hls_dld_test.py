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
hvplot.extension('bokeh')
from bokeh.plotting import show

# Configure Earthaccess
auth = earthaccess.Auth()
auth.login(strategy = "netrc")
print(auth.authenticated)

# Load AOIs
aois = gp.read_file("./data/feature_layers/roi.geojson")

bbox = tuple(list(aois.total_bounds))

# Define time window
temporal = ("2020-05-01T00:00:00", "2020-08-31T23:59:59")

# Get results
MAX_IMG = 100 # nr. of maximum returned images
results = earthaccess.search_data(
    short_name=['HLSL30','HLSS30'],
    bounding_box=bbox,
    temporal=temporal,
    count=MAX_IMG
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
# Tiles have 1 band and are divided into 512x512 pixel chunks
chunk_size = dict(band=1, x=512, y=512)
for e in evi_band_links:
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

# reproject AOIs to HLS CRS
aoiUTM = aois.to_crs(nir.spatial_ref.crs_wkt)

# Crop data with aois and include any pixels touched by the polygon
nir_cropped = nir.rio.clip(aoiUTM.geometry.values, aoiUTM.crs, all_touched=True)

# Define function to scale 
def scaling(band):
    scale_factor = band.attrs['scale_factor'] 
    band_out = band.copy()
    band_out.data = band.data*scale_factor
    band_out.attrs['scale_factor'] = 1
    return(band_out)

nir_cropped_scaled = scaling(nir_cropped)

# Red
red_cropped = red.rio.clip(aoiUTM.geometry.values, aoiUTM.crs, all_touched=True)
red_cropped_scaled = scaling(red_cropped)
# Blue
blue_cropped = blue.rio.clip(aoiUTM.geometry.values, aoiUTM.crs, all_touched=True)
blue_cropped_scaled = scaling(blue_cropped)
# Fmask
fmask_cropped = fmask.rio.clip(aoiUTM.geometry.values, aoiUTM.crs, all_touched=True)

print('Data is loaded and scaled!')

def calc_evi(red, blue, nir):
    # Create EVI xarray.DataArray that has the same coordinates and metadata
    evi = red.copy()
    # Calculate the EVI
    evi_data = 2.5 * ((nir.data - red.data) / (nir.data + 6.0 * red.data - 7.5 * blue.data + 1.0))
    # Replace the Red xarray.DataArray data with the new EVI data
    evi.data = evi_data
    # exclude the inf values
    evi = xr.where(evi != np.inf, evi, np.nan, keep_attrs=True)
    # change the long_name in the attributes
    evi.attrs['long_name'] = 'EVI'
    evi.attrs['scale_factor'] = 1
    return evi

# Generate EVI array
evi_cropped = calc_evi(red_cropped_scaled, blue_cropped_scaled, nir_cropped_scaled)

# Quick plot of EVI in aoi
img = evi_cropped.hvplot.image(cmap='YlGn', 
                               frame_width= 800, 
                               fontscale=1.6, 
                               crs='EPSG:32610', 
                               tiles='ESRI'
                               ).opts(title=f'HLS-derived EVI, {evi_cropped.SENSING_TIME}')

show(hvplot.render(img))

# define bits to mask out
"""
7-6 = aerosols
5 = water
4 = snow/ice
3 = cloud shadow
2 = cloud/shadow adjacent
1 = cloud
0 = cirrus
"""
bit_nums = [1,2,3,4,5]

# Define bit masking function
def create_quality_mask(quality_data, bit_nums: list = [1, 2, 3, 4, 5]):
    """
    Uses the Fmask layer and bit numbers to create a binary mask of good pixels.
    By default, bits 1-5 are used.
    """
    mask_array = np.zeros((quality_data.shape[0], quality_data.shape[1]))
    # Remove/Mask Fill Values and Convert to Integer
    quality_data = np.nan_to_num(quality_data, 0).astype(np.int8)
    for bit in bit_nums:
        # Create a Single Binary Mask Layer
        mask_temp = np.array(quality_data) & 1 << bit > 0
        mask_array = np.logical_or(mask_array, mask_temp)
    return mask_array

# Apply mask and filter cropped image
mask_layer = create_quality_mask(fmask_cropped.data, bit_nums)
evi_cropped_qf = evi_cropped.where(~mask_layer)

# Export to COG tiffs
original_name = evi_band_links[0].split('/')[-1]

# Generate output name from the original filename
out_name = f"{original_name.split('v2.0')[0]}v2.0_EVI_cropped.tif"

out_folder = './data/hls/'
evi_cropped.rio.to_raster(raster_path=f'{out_folder}{out_name}', driver='COG')

# delete unused files
del evi_cropped, out_folder, out_name, red_cropped, blue_cropped, nir_cropped, red_cropped_scaled, blue_cropped_scaled, nir_cropped_scaled

