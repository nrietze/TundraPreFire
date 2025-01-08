# %%
"""
_summary_
"""

# Load packages
import os
import numpy as np
import geopandas as gp
from osgeo import gdal
import xarray as xr
import rioxarray as rxr
import earthaccess

# GDAL configurations used to successfully access LP DAAC Cloud Assets via vsicurl
gdal.SetConfigOption('GDAL_HTTP_COOKIEFILE','~/cookies.txt')
gdal.SetConfigOption('GDAL_HTTP_COOKIEJAR', '~/cookies.txt')
gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN','EMPTY_DIR')
gdal.SetConfigOption('CPL_VSIL_CURL_ALLOWED_EXTENSIONS','TIF')
gdal.SetConfigOption('GDAL_HTTP_UNSAFESSL', 'YES')
gdal.SetConfigOption('GDAL_HTTP_MAX_RETRY', '10')
gdal.SetConfigOption('GDAL_HTTP_RETRY_DELAY', '0.5')

# Configure Earthaccess
auth = earthaccess.Auth()
auth.login(strategy = "netrc")
if auth.authenticated:
    print("Successfully authenticated Earthaccess.")
# %% User configurations

# Define bands/indices to download
BAND_INDEX = ["GEMI"]

# Define time search window
START_DATE = "2020-09-10T00:00:00"
# START_DATE = "2019-09-12T00:00:00"
END_DATE = "2020-09-12T23:59:59"
# END_DATE = "2019-09-13T23:59:59"

# nr. of maximum returned images
MAX_IMG = 10000

# define chunk size for data loading
chunk_size = dict(band=1, x=512, y=512)

# Load AOIs
aois = gp.read_file("data/feature_layers/roi.geojson")
bbox = tuple(list(aois.total_bounds))
print("AOI loaded.")

# Define time window
temporal = (START_DATE, END_DATE)

# %% Define user functions

def get_matching_url(url_list, pattern):
    """Function to extract the link of a desired spectral band.

    Args:
        url_list (list): List with URLs of all filtered images.
        pattern (str): HLS product band name for the desired band. E.g., "B05" (NIR)

    Returns:
        str: the URL to that band
    """
    return next((url for url in url_list if pattern in url), None)

# Define function to scale 
def scaling(band):
    scale_factor = band.attrs['scale_factor'] 
    band_out = band.copy()
    band_out.data = band.data*scale_factor
    band_out.attrs['scale_factor'] = 1
    return(band_out)

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

# Define raster loading function
def load_rasters(index_band_links, band_name,
                 band_dict = {"Fmask" : "Fmask"}, 
                 region = None, 
                 chunk_size = dict(band=1, x=512, y=512)):
    """Function to load COG rasters of selected bands into memory.

    Args:
        index_band_links (list): List of URLs
        band_name (str): Spectral band abbreviaiton. E.g., SWIR1, NIR, RED.
        band_dict (dict): Dictionary with name pairs of band_name-type names and the HLS band names ("B05")
        region (gp): Region to crop imagery with. Will be reprojected to data if not same.
        chunk_size (dict, optional): Settings for handling data chunking. Defaults to dict(band=1, x=512, y=512).

    Returns:
        xr.DataArray: The raster of the desired band raster.
    """
    
    hls_band_name = band_dict[band_name]
    
    # Get URL of matching spectral band
    lnk = get_matching_url(index_band_links, hls_band_name)

    # Open the raster
    # Use vsicurl to load the data directly into memory (be patient, may take a few seconds)
    # Tiles have 1 band and are divided into 512x512 pixel chunks
    raster = rxr.open_rasterio(lnk, chunks=chunk_size, masked=True).squeeze('band', drop=True)
    
    if region is not None:
        # reproject AOIs to HLS CRS
        regionUTM = region.to_crs(raster.spatial_ref.crs_wkt)
        
        raster = raster.rio.clip(regionUTM.geometry.values, regionUTM.crs, all_touched=True)
    
    # print("The COGs have been loaded into memory!")

    if hls_band_name != "Fmask":
        raster.attrs['scale_factor'] = 0.0001  # hard coded scale factor
        raster = scaling(raster)
    
    # print('Data is loaded and scaled!')

    # Return the raster
    return raster

# Calculate spectral index
def calc_index(urls, INDEX_NAME, bit_nums, region = None): 
    # Band name pairs
    if urls[0].split('/')[4] == 'HLSS30.020':
        HLS_BAND_DICT = {
            "BLUE" : "B02",
            "GREEN" : "B03",
            "RED" : "B04",
            "NIR" : "B8A",
            "SWIR1" : "B11",
            "SWIR2" : "B12",
            "Fmask" : "Fmask"
        }
    else:
        HLS_BAND_DICT = {
            "BLUE" : "B02",
            "GREEN" : "B03",
            "RED" : "B04",
            "NIR" : "B05",
            "SWIR1" : "B06",
            "SWIR2" : "B07",
            "TIR1" : "B10",
            "TIR2" : "B11",
            "Fmask" : "Fmask"
        }
    
    if INDEX_NAME == "NDMI":
        # load spectral bands needed for NDMI
        nir = load_rasters(urls, "NIR",
                           band_dict = HLS_BAND_DICT, region = region,
                           chunk_size = chunk_size)
        
        swir1 = load_rasters(urls, "SWIR1",
                             band_dict = HLS_BAND_DICT, region = region,
                             chunk_size = chunk_size)
        
        spectral_index = nir.copy()
        spectral_index_data = (nir - swir1) / (nir + swir1)
        
        SCALE_FACTOR = 1

    elif INDEX_NAME == "NDVI":
        # load spectral bands needed for NDMI
        nir = load_rasters(urls, "NIR",
                           band_dict = HLS_BAND_DICT, region = region,
                           chunk_size = chunk_size)
        
        red = load_rasters(urls, "RED",
                           band_dict = HLS_BAND_DICT, region = region,
                           chunk_size = chunk_size)
        
        spectral_index = nir.copy()
        spectral_index_data = (nir - red) / (nir + red)
        
        SCALE_FACTOR = 1
        
    elif INDEX_NAME == "NBR":
        # load spectral bands needed for NDMI
        nir = load_rasters(urls, "NIR",
                           band_dict = HLS_BAND_DICT, region = region,
                           chunk_size = chunk_size)
        
        swir2 = load_rasters(urls, "SWIR2",
                           band_dict = HLS_BAND_DICT, region = region,
                           chunk_size = chunk_size)
        
        spectral_index = nir.copy()
        spectral_index_data = (nir - swir2) / (nir + swir2)
        
        SCALE_FACTOR = 1
        
    elif INDEX_NAME == "GEMI":
        # load spectral bands needed for NDMI
        nir = load_rasters(urls, "NIR",
                           band_dict = HLS_BAND_DICT, region = region,
                           chunk_size = chunk_size)
        
        red = load_rasters(urls, "RED",
                           band_dict = HLS_BAND_DICT, region = region,
                           chunk_size = chunk_size)
        
        spectral_index = nir.copy()
        
        term1 = (2 * (nir**2 - red**2) + 1.5 * nir + 0.5 * red) / (nir + red + 0.5)
        spectral_index_data = term1 * (1 - 0.25 * term1) - ((red - 0.125) / (1 - red))
        
        SCALE_FACTOR = 1
        
    # Replace the dummy xarray.DataArray data with the new spectral index data
    spectral_index.data = spectral_index_data
    # exclude the inf values
    spectral_index = xr.where(spectral_index != np.inf, spectral_index, np.nan, keep_attrs=True)
    # change the long_name in the attributes
    spectral_index.attrs['long_name'] = INDEX_NAME
    spectral_index.attrs['scale_factor'] = SCALE_FACTOR
    
    # get Fmask
    fmask = load_rasters(urls, "Fmask",
                         region = region,
                         band_dict = HLS_BAND_DICT,
                         chunk_size = dict(band=1, x=512, y=512))
    
    # Apply mask and filter spectral_index image
    mask_layer = create_quality_mask(fmask.data, bit_nums)
    spectral_index_qf = spectral_index.where(~mask_layer)
    
    return spectral_index_qf

# %% Data search
# Get results
results = earthaccess.search_data(
    short_name=['HLSL30','HLSS30'],
    bounding_box=bbox,
    temporal=temporal,
    count=MAX_IMG
)

# retrieve URLS
print("Retrieving granules...")
hls_results_urls = [granule.data_links() for granule in results]

print(f"{len(hls_results_urls)} granules found.")

browse_urls = [granule.dataviz_links()[0] for granule in results] # 0 retrieves only the https links

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
bit_nums = [0,1,2,3,4,5]

# Generate EVI array
# h = hls_results_urls[10]
# ndmi = calc_index(h, "NDMI", region = aois, bit_nums = bit_nums)

# %% automate and run over all found granules
# output directory
out_folder = 'data/raster/hls/'

for j, h in enumerate(hls_results_urls):
    original_name = h[0].split('/')[-1]

    print(f"Processing: {original_name.split('v2.0')[0]}", end = "\n")
    
    for INDEX_NAME in BAND_INDEX:
        # Generate output name from the original filename
        out_name = f"{original_name.split('v2.0')[0]}v2.0_{INDEX_NAME}_cropped.tif"
        
        out_path = f'{out_folder}{out_name}'
        
        # Check if file already exists in output directory, if yes--skip that file and move to the next observation
        if os.path.exists(out_path):
            print(f"{out_name} has already been processed and is available in this directory, moving to next file.")
            continue
        
        # Calculate spectral index
        spectral_index = calc_index(h, INDEX_NAME, 
                                    region = aois, 
                                    bit_nums = bit_nums)
        
        print(f"{INDEX_NAME} index calculated.", end = "\n")
        
        # Remove any observations that are entirely fill value
        if np.nansum(spectral_index.data) == 0.0:
            print(f"File: {h[0].split('/')[-1].rsplit('.', 1)[0]} was entirely fill values and will not be exported.")
            continue
        
        # Export to COG tiffs
        spectral_index.rio.to_raster(raster_path = out_path, driver = 'COG')
    
    # Load and export scene's fmask
    fmask = load_rasters(h, "Fmask",
                            region = aois,
                            chunk_size = dict(band=1, x=512, y=512))
    
    out_name = f"{original_name.split('v2.0')[0]}v2.0_Fmask_cropped.tif"
        
    fmask.rio.to_raster(raster_path = f'{out_folder}{out_name}', driver = 'COG')
    
    print(f"Processed file {j+1} of {len(hls_results_urls)}")
