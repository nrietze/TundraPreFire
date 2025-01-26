"""
Script to download and preprocess HLS data.
- HLS data search
- Cloud  masking using Fmask provided with imagery
- Water masking done with dynamic NDWI thresholding (Otsu's method)
- Automated implementation using joblib.Parallel.

Author: Nils Rietze (nils.rietze@uzh.ch), github.com/nrietze

Date: 15.01.2025 (last update of this docstring)
"""

# Load packages
import os
import re
from glob import glob
import time
import datetime
import numpy as np
import dask.array as da
from joblib import Parallel, delayed
import multiprocessing

import geopandas as gp
from osgeo import gdal
import xarray as xr
from xrspatial import multispectral
import rioxarray as rxr
import earthaccess

from skimage.filters import threshold_multiotsu, threshold_otsu

# %% Define user functions
# Function to gather a list of lists with tiff files
def get_hls_tif_list(main_dir):
    hls_granule_tiffs = []
    
    for root, dirs, files in os.walk(main_dir):
        
        if not dirs:
            tif_files = glob(os.path.join(root, '*.tif'))
            
            hls_granule_tiffs.append(tif_files)

    return hls_granule_tiffs

def get_doy(date_str):
    date_obj = datetime.datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S")
    return date_obj.timetuple().tm_yday

def search_files_by_doy_range(file_lists, start_date, end_date):
    """
    Function to extract HLS granules that fall in between two dates, defined by their doy. This is only applied if the script is executed to process NBR.
    """
    matching_lists = []
    
    start_doy = get_doy(start_date)  
    end_doy = get_doy(end_date)  
    
    year = datetime.datetime.strptime(start_date,
                                      "%Y-%m-%dT%H:%M:%S").year
    
    start_doy_str = f"{start_doy:03d}"
    end_doy_str = f"{end_doy:03d}"
    
    for file_list in file_lists:
        for file in file_list:
            # Extract the DOY from the file path
            match = re.search(r'(\d{4})(\d{3})', file)
            if match:
                file_year = match.group(1)  
                file_doy = match.group(2)  
                if file_year == str(year) and start_doy_str <= file_doy <= end_doy_str:  # Check if file_doy is in the range
                    matching_lists.append(file_list)
                    break 
        
    return matching_lists
    
def get_matching_url(url_list, pattern):
    """Function to extract the link of a desired spectral band.

    Args:
        url_list (list): List with files of all filtered images.
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
def create_quality_mask(quality_data, 
                        bit_nums: list = [1, 2, 3, 4]):
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
    """
    Function to load COG rasters of selected bands into memory.

    Args:
        index_band_links (list): List of files
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
def calc_index(files, 
               index_name,
               bit_nums,
               out_folder,
               region = None): 
    
    # Band name pairs
    if 'HLS.S30' in files[0]:
        hls_band_dict = {
            "COASTAL-AEROSOL": "B01",
            "BLUE": "B02",
            "GREEN": "B03",
            "RED": "B04",
            "RED-EDGE1": "B05",
            "RED-EDGE2": "B06",
            "RED-EDGE3": "B07",
            "NIR-Broad": "B08",
            "NIR1": "B8A",
            "WATER-VAPOR": "B09",
            "CIRRUS": "B10",
            "SWIR1": "B11",
            "SWIR2": "B12",
            "FMASK": "Fmask",
            "VZA": "VZA",
            "VAA": "VAA",
            "SZA": "SZA",
            "SAA": "SAA",
        }
    else:
        hls_band_dict = {
            "COASTAL-AEROSOL": "B01",
            "BLUE": "B02",
            "GREEN": "B03",
            "RED": "B04",
            "NIR1": "B05",
            "SWIR1": "B06",
            "SWIR2": "B07",
            "CIRRUS": "B09",
            "TIR1": "B10",
            "TIR2": "B11",
            "FMASK": "Fmask",
            "VZA": "VZA",
            "VAA": "VAA",
            "SZA": "SZA",
            "SAA": "SAA",
        }
    
    print(f"Calculating: {index_name}")
    
    if index_name == "NDMI":
        # load spectral bands needed for NDMI
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict, region = region,
                           chunk_size = chunk_size)
        
        swir1 = load_rasters(files, "SWIR1",
                             band_dict = hls_band_dict, region = region,
                             chunk_size = chunk_size)
        spectral_index = multispectral.ndmi(nir, swir1)
        
        # Exclude data outside valid value range
        spectral_index = spectral_index.where(
            (spectral_index >= -1) & (spectral_index <= 1),np.nan)
        
        SCALE_FACTOR = 1

    elif index_name == "NDVI":
        # load spectral bands needed for NDVI
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict, region = region,
                           chunk_size = chunk_size)
        
        red = load_rasters(files, "RED",
                           band_dict = hls_band_dict, region = region,
                           chunk_size = chunk_size)
        
        spectral_index = multispectral.ndvi(nir, red)
        
        # Exclude data outside valid value range
        spectral_index = spectral_index.where(
            (spectral_index >= -1) & (spectral_index <= 1),np.nan)
        
        SCALE_FACTOR = 1

    elif index_name == "NBR":
        # load spectral bands needed for NBR
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict, region = region,
                           chunk_size = chunk_size)
        
        swir2 = load_rasters(files, "SWIR2",
                           band_dict = hls_band_dict, region = region,
                           chunk_size = chunk_size)
        
        # spectral_index = nir.copy()
        spectral_index = multispectral.nbr(nir, swir2)
        
        # Exclude data outside valid value range
        spectral_index = spectral_index.where(
            (spectral_index >= -1) & (spectral_index <= 1),np.nan)
        
        SCALE_FACTOR = 1
        
    elif index_name == "GEMI":
        # load spectral bands needed for GEMI
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict, region = region,
                           chunk_size = chunk_size)
        
        red = load_rasters(files, "RED",
                           band_dict = hls_band_dict, 
                           region = region,
                           chunk_size = chunk_size)
        
        spectral_index = nir.copy()
        
        term1 = (2 * (nir**2 - red**2) + 1.5 * nir + 0.5 * red) / (nir + red + 0.5)
        spectral_index_data = term1 * (1 - 0.25 * term1) - ((red - 0.125) / (1 - red))
        
        # Replace the dummy xarray.DataArray data with the new spectral index data
        spectral_index.data = spectral_index_data
        
        SCALE_FACTOR = 1
        
    # change the long_name in the attributes
    spectral_index.attrs['long_name'] = index_name
    spectral_index.attrs['scale_factor'] = SCALE_FACTOR

    # ---- MASKING ----
    # Compute Water mask
    tile_name = files[0].split('/')[-1]
    filename = f"{tile_name.split('v2.0')[0]}v2.0_watermask.tif"
    wm_path = f"{out_folder}{filename}"
    
    # Check if that tile's water mask exists
    if os.path.exists(wm_path):
        print("Loading water mask.")
        water_mask = rxr.open_rasterio(wm_path,
                                       chunks=chunk_size,
                                       masked=True
                                       ).squeeze('band', drop=True
                                                 ).astype(bool)
    else:
        print("Computing NDWI water mask.")
        # load spectral bands needed for NDWI
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict,
                           region = region,
                           chunk_size = chunk_size)
        
        green = load_rasters(files, "GREEN",
                             band_dict = hls_band_dict, region = region,
                             chunk_size = chunk_size)
        
        ndwi = nir.copy()
        ndwi.data = (green - nir) / (green + nir)
        
        # Exclude data outside valid value range
        ndwi = ndwi.where((ndwi >= -1) & (ndwi <= 1), np.nan)
        
        ndwi.attrs['long_name'] = "NDWI"
        
        # Export NDWI as COG tiff
        ndwi_filename = f"{tile_name.split('v2.0')[0]}v2.0_NDWI.tif"
        ndwi_path = f"{out_folder}{ndwi_filename}"
        ndwi.rio.to_raster(raster_path = ndwi_path, driver = 'COG')
    
        # Apply Otsu's thresholding to NDWI
        hist_ndwi = da.histogram(ndwi, bins=2, range=[-1, 1])
        # otsu_thresh = threshold_otsu(hist_ndwi[1])[1]
        otsu_thresh = threshold_multiotsu(hist_ndwi[1])[1]
        
        # Apply threshold to mask out water
        water_mask = ndwi > otsu_thresh
        
        # Export water mask as COG tiff
        water_mask.astype(float).rio.to_raster(
            raster_path = wm_path, driver = 'COG')
    
    # Fmask
    tile_name = files[0].split('/')[-1]
    filename = f"{tile_name.split('v2.0')[0]}v2.0_Fmask.tif"
    fmask_path = f"{out_folder}{filename}"
    
    # Check if that tile's Fmask tiff exists
    if os.path.exists(fmask_path):
        print("Loading existing Fmask.")
        fmask = rxr.open_rasterio(fmask_path,
                                  chunks=chunk_size,
                                  masked=True).squeeze('band', drop=True)
    
    else:
        print("Downloading Fmask.")
        # Load tile's Fmask
        fmask = load_rasters(files, "Fmask",region=region,chunk_size=chunk_size)

        # Export scene's Fmask
        fmask.rio.to_raster(raster_path = fmask_path, driver = 'COG')
    
    # Convert Fmask to quality mask
    mask_layer = create_quality_mask(fmask.data, bit_nums)
    
    # Apply mask and filter spectral_index image
    merged_mask = np.logical_and(mask_layer,water_mask)
    spectral_index_qf = spectral_index.where(~merged_mask)
    
    # exclude the inf values
    spectral_index_qf = xr.where(spectral_index_qf != np.inf, 
                                 spectral_index_qf, np.nan, 
                                 keep_attrs=True)
    
    return spectral_index_qf

def joblib_hls_preprocessing(files:list, 
                             band_index:str,
                             bit_nums:list,
                             out_folder: str,
                             REMASK_DATA:bool):
    start_time = time.time()
    
    original_name = files[0].split('/')[-1]

    print(f"Processing: {original_name.split('v2.0')[0]}", end = "\n")
    
    for index_name in band_index:
        # Generate output name from the original filename
        out_name = f"{original_name.split('v2.0')[0]}v2.0_{index_name}.tif"
        
        out_path = f'{out_folder}{out_name}'
        
        # Check if file already exists in output directory, if yes--skip that file and move to the next observation
        if os.path.exists(out_path):
            
            if OVERWRITE_DATA:
                print(f"{out_name} exists but will be overwritten.")
                
                # Load existing data
                spectral_index = calc_index(files, index_name,
                                            out_folder = out_folder,
                                            bit_nums = bit_nums)
                print(f"{index_name} index calculated.", end = "\n")
                
            else:
                print(f"{out_name} has already been processed and is available in this directory, moving to next file.")
                continue
        
        # "else" necessary to handle remasking of existing tiles in previous "if"
        else: 
            # Calculate spectral index
            spectral_index = calc_index(files, index_name, 
                                        out_folder = out_folder,
                                        bit_nums = bit_nums)
            
            print(f"{index_name} index calculated.", end = "\n")
        
        # Remove any observations that are entirely fill value
        if np.nansum(spectral_index.data) == 0.0:
            print(f"File: {files[0].split('/')[-1].rsplit('.', 1)[0]} was entirely fill values and will not be exported.")
            continue
        
        # Export to COG tiffs
        spectral_index.rio.to_raster(raster_path = out_path, driver = 'COG')
    
    print("--- %s seconds ---" % (time.time() - start_time))


# %% Execute pararellized processing
if __name__ == "__main__":
    # GDAL configurations used to successfully access LP DAAC Cloud Assets via vsicurl
    gdal.SetConfigOption('GDAL_HTTP_COOKIEFILE','~/cookies.txt')
    gdal.SetConfigOption('GDAL_HTTP_COOKIEJAR', '~/cookies.txt')
    gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN','EMPTY_DIR')
    gdal.SetConfigOption('CPL_VSIL_CURL_ALLOWED_EXTENSIONS','TIF')
    gdal.SetConfigOption('GDAL_HTTP_UNSAFESSL', 'YES')
    gdal.SetConfigOption('GDAL_HTTP_MAX_RETRY', '10')
    gdal.SetConfigOption('GDAL_HTTP_RETRY_DELAY', '0.5')

    # Load AOIs
    # aois = gp.read_file("data/feature_layers/roi.geojson")

    # fire_perimeters_2020 = gp.read_file("data/feature_layers/fire_atlas/final_viirs2020.gpkg")

    # # Intersect aoi & perimeter
    # fire_within_roi = gp.overlay(fire_perimeters_2020, aois,
    #                             how='intersection')

    # convert end day of fire to date
    # fire_within_roi['end_date'] = fire_within_roi.apply(
    #     lambda row: datetime.date(int(row['ted_year']), 
    #                             int(row['ted_month']), 
    #                             int(row['ted_day'])), axis=1)

    # END_OF_FIRE = fire_within_roi['end_date'].max()

    # bbox = tuple(list(fire_within_roi.total_bounds))
    bbox = (144.4001779629389, 71.24588582272926,
            146.16765743760948, 71.6168891179392)
    print("AOI loaded.")

    # Load MGRS tile centroids to find msot suitable HLS tile
    # mgrs_tile_centroids = gp.read_file("data/feature_layers/MGRS_centroids.geojson")

    # Find best UTM tile for the feature by identifying closest tile centroid to feature centroid
    # sindex = mgrs_tile_centroids.geometry.sindex.nearest(aois.geometry.centroid)

    # # Query tile in all MGRS tiles
    # nearest_mgrs_tile_centroid = mgrs_tile_centroids.iloc[sindex[1],:]

    # OPTIMAL_TILE_NAME = nearest_mgrs_tile_centroid.Name.item()
    OPTIMAL_TILE_NAME = "54WXE"

    # Define bands/indices to process
    band_index = ["NBR"]

    # Overwrite existing tiles?
    OVERWRITE_DATA = True
    
    # define chunk size for data loading
    chunk_size = dict(band=1, x=3600, y=3600)

    # Set original data paths
    HLS_PARENT_PATH = "/home/nrietz/scratch/raster/hls/"

    hls_granules_paths = get_hls_tif_list(HLS_PARENT_PATH)

    if "NBR" in band_index:
        # Define time search window
        # START_DATE = "2020-05-01T00:00:00" #full growing season
        START_DATE = "2020-09-10T00:00:00"
        # START_DATE = "2019-09-12T00:00:00"

        # END_DATE = "2020-10-15T23:59:59" #full growing season
        END_DATE = "2020-09-12T23:59:59"
        # END_DATE = "2019-09-13T23:59:59"

        hls_granules_paths = search_files_by_doy_range(hls_granules_paths, START_DATE, END_DATE)

    print(f"{len(hls_granules_paths)} granules found to process.")

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
    bit_nums = [0,1,2,3,4]

    # output directory
    out_folder = '/data/nrietz/raster/hls/'
    N_CORES = -1 # -1 = all are used, -2 all but one

    multiprocessing.set_start_method('spawn')
    
    Parallel(n_jobs=N_CORES, backend='loky')(
        delayed(joblib_hls_preprocessing)(files,band_index,bit_nums,out_folder,OVERWRITE_DATA) for files in hls_granules_paths)