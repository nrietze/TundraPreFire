"""
Script to download and preprocess HLS data.
- HLS data search
- Cloud  masking using Fmask provided with imagery
- Water masking done with dynamic NDWI thresholding (Otsu's method)
- Automated implementation using joblib.Parallel.

Author: Nils Rietze (nils.rietze@uzh.ch), github.com/nrietze

Date: 15.01.2025 (last update of this docstring)
"""

#%% Load packages
import os
import sys
import ast
import re
from glob import glob
import time
import datetime
import multiprocessing
import platform
import numpy as np
import dask.array as da
from joblib import Parallel, delayed
import pandas as pd
from dateutil.relativedelta import relativedelta
from osgeo import gdal
import xarray as xr
from xrspatial import multispectral
import rioxarray as rxr
from pyproj import CRS

from skimage.filters import threshold_multiotsu

# %% Define user functions
# Function to gather a list of lists with tiff files
def get_hls_tif_list(main_dir):
    hls_granule_tiffs = []
    
    for root, dirs, files in os.walk(main_dir):
        
        if not dirs:
            tif_files = [f for f in glob(os.path.join(root, '*.tif')) if "checkpoint" not in f]
            
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
def create_quality_mask(quality_data: np.array,
                        bit_nums: list):
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
def load_rasters(index_band_links: list,
                 band_name: str,
                 band_dict: dict,
                 chunk_size: dict,
                 region = None):
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
        region_utm = region.to_crs(raster.spatial_ref.crs_wkt)
        
        raster = raster.rio.clip(region_utm.geometry.values, region_utm.crs, all_touched=True)
    
    # print("The COGs have been loaded into memory!")

    if hls_band_name != "Fmask":
        raster.attrs['scale_factor'] = 0.0001  # hard coded scale factor
        raster = scaling(raster)
    
    # print('Data is loaded and scaled!')

    # Return the raster
    return raster

# Calculate spectral index
def calc_index(files: list,
               index_name: str,
               bit_nums: list,
               out_folder: str,
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
    
    tile_name = os.path.basename(files[0])
    
    # Extract UTM EPSG code for reprojection
    file_utm_zone = tile_name.split(".")[2]
    utm_epsg_code = 32600 + int(file_utm_zone[1:3])
    
    print(f"Calculating: {index_name}")
    
    if index_name == "NDMI":
        # load spectral bands needed for NDMI
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict,
                           region = region,
                           chunk_size = chunk_size)
        
        swir1 = load_rasters(files, "SWIR1",
                             band_dict = hls_band_dict,
                             region = region,
                             chunk_size = chunk_size)
        with np.errstate(divide='ignore'): #to ignore divide by zero
            spectral_index = multispectral.ndmi(nir, swir1)
        
        # Exclude data outside valid value range
        spectral_index = spectral_index.where(
            (spectral_index >= -1) & (spectral_index <= 1),np.nan)
        
        SCALE_FACTOR = 1

    elif index_name == "NDVI":
        # load spectral bands needed for NDVI
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict,
                           region = region,
                           chunk_size = chunk_size)
        
        red = load_rasters(files, "RED",
                           band_dict = hls_band_dict,
                           region = region,
                           chunk_size = chunk_size)
        
        with np.errstate(divide='ignore'): #to ignore divide by zero
            spectral_index = multispectral.ndvi(nir, red)
        
        # Exclude data outside valid value range
        spectral_index = spectral_index.where(
            (spectral_index >= -1) & (spectral_index <= 1),np.nan)
        
        SCALE_FACTOR = 1
        
    elif index_name == "NMDI":
        # load spectral bands needed for NMDI
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict,
                           region = region,
                           chunk_size = chunk_size)
        
        swir1 = load_rasters(files, "SWIR1",
                             band_dict = hls_band_dict,
                             region = region,
                             chunk_size = chunk_size)
        
        swir2 = load_rasters(files, "SWIR2",
                             band_dict = hls_band_dict,
                             region = region,
                             chunk_size = chunk_size)
        
        spectral_index = nir.copy()
        
        with np.errstate(divide='ignore'): #to ignore divide by zero
            spectral_index_data = (nir - (swir1 - swir2)) / (nir + (swir1 - swir2))
        
        # Replace the dummy xarray.DataArray data with the new spectral index data
        spectral_index.data = spectral_index_data
        
        # Exclude data outside valid value range
        spectral_index = spectral_index.where(
            (spectral_index >= -1) & (spectral_index <= 1),np.nan)
        
        SCALE_FACTOR = 1

    elif index_name == "NBR":
        # load spectral bands needed for NBR
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict,
                           region = region,
                           chunk_size = chunk_size)
        
        swir2 = load_rasters(files, "SWIR2",
                             band_dict = hls_band_dict,
                             region = region,
                             chunk_size = chunk_size)
        
        # spectral_index = nir.copy()
        with np.errstate(divide='ignore'): #to ignore divide by zero
            spectral_index = multispectral.nbr(nir, swir2)
        
        # Exclude data outside valid value range
        spectral_index = spectral_index.where(
            (spectral_index >= -1) & (spectral_index <= 1),np.nan)
        
        SCALE_FACTOR = 1
        
    elif index_name == "GEMI":
        # load spectral bands needed for GEMI
        nir = load_rasters(files, "NIR1",
                           band_dict = hls_band_dict,
                           region = region,
                           chunk_size = chunk_size)
        
        red = load_rasters(files, "RED",
                           band_dict = hls_band_dict,
                           region = region,
                           chunk_size = chunk_size)
        
        spectral_index = nir.copy()
        
        with np.errstate(divide='ignore'): #to ignore divide by zero
            term1 = (2 * (nir**2 - red**2) + 1.5 * nir + 0.5 * red) / (nir + red + 0.5)
            spectral_index_data = term1 * (1 - 0.25 * term1) - ((red - 0.125) / (1 - red))
        
        # Replace the dummy xarray.DataArray data with the new spectral index data
        spectral_index.data = spectral_index_data
        
        # Exclude data outside valid value range
        spectral_index = spectral_index.where(
            (spectral_index >= -1) & (spectral_index <= 1),np.nan)
        
        SCALE_FACTOR = 1
        
    # change the long_name in the attributes
    spectral_index.attrs['long_name'] = index_name
    spectral_index.attrs['scale_factor'] = SCALE_FACTOR
    
    # ---- MASKING ----
    # Compute Water mask
    filename = f"{tile_name.split('v2.0')[0]}v2.0_watermask.tif"
    wm_path = os.path.join(out_folder,filename)
    
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
                             band_dict = hls_band_dict,
                             region = region,
                             chunk_size = chunk_size)
        
        ndwi = nir.copy()
        with np.errstate(divide='ignore'): #to ignore divide by zero
            ndwi.data = (green - nir) / (green + nir)
        
        # Exclude data outside valid value range
        # ndwi = ndwi.where((ndwi >= -1) & (ndwi <= 1), np.nan)
        
        ndwi.attrs['long_name'] = "NDWI"
        
        # Reproject if not EPSG conform
        if ndwi.rio.crs is not CRS.from_epsg(utm_epsg_code):
            with np.errstate(divide='ignore'):
                ndwi_rprj = ndwi.rio.reproject(f"EPSG:{utm_epsg_code}",
                                            chunks=chunk_size)
            
        # Export NDWI as COG tiff
        ndwi_filename = f"{tile_name.split('v2.0')[0]}v2.0_NDWI.tif"
        ndwi_path = os.path.join(out_folder,ndwi_filename)
        ndwi_rprj.rio.to_raster(raster_path = ndwi_path, driver = 'COG')
    
        # Apply Otsu's thresholding to NDWI
        hist_ndwi = da.histogram(ndwi, bins=2, range=[-1, 1])
        # otsu_thresh = threshold_otsu(hist_ndwi[1])[1]
        otsu_thresh = threshold_multiotsu(hist_ndwi[1])[1]
        
        # Apply threshold to mask out water
        water_mask = ndwi > otsu_thresh
        
        # Reproject if not EPSG conform
        if water_mask.rio.crs is not CRS.from_epsg(utm_epsg_code):
            water_mask_rprj = water_mask.astype(int).rio.reproject(f"EPSG:{utm_epsg_code}")
            
        # Export water mask as COG tiff
        water_mask_rprj.rio.to_raster(raster_path = wm_path, driver = 'COG')
    
    # Fmask
    filename = f"{tile_name.split('v2.0')[0]}v2.0_Fmask.tif"
    fmask_path = os.path.join(out_folder,filename)
    
    # Check if that tile's Fmask tiff exists
    if os.path.exists(fmask_path):
        print("Loading existing Fmask.")
        fmask = rxr.open_rasterio(fmask_path,
                                  chunks=chunk_size,
                                  masked=True).squeeze('band', drop=True)
    
    else:
        print("Loading Fmask.")
        # Load tile's Fmask
        fmask = load_rasters(files, "FMASK",
                             band_dict = hls_band_dict,
                             region=region,
                             chunk_size=chunk_size)

        # Reproject if not EPSG conform
        if fmask.rio.crs is not CRS.from_epsg(utm_epsg_code):
            fmask_rprj = fmask.rio.reproject(f"EPSG:{utm_epsg_code}")
        
        # Export scene's Fmask
        fmask_rprj.rio.to_raster(raster_path = fmask_path, driver = 'COG')
    
    # Convert Fmask to quality mask
    mask_layer = create_quality_mask(fmask.data, bit_nums)
    
    # Apply mask and filter spectral_index image
    merged_mask = np.logical_or(mask_layer,water_mask.data)
    spectral_index_qf = spectral_index.where(~merged_mask)
    
    # exclude the inf values
    spectral_index_qf = xr.where(spectral_index_qf != np.inf,
                                 spectral_index_qf, np.nan,
                                 keep_attrs=True)
    
    return spectral_index_qf

def joblib_hls_preprocessing(files: list,
                             band_index: str,
                             bit_nums: list,
                             out_folder: str,
                             OVERWRITE_DATA = False,
                             skip_source = None):
    start_time = time.time()
    
    original_name = os.path.basename(files[0])

    print(f"Processing: {original_name.split('v2.0')[0]}", end = "\n")
    
    for index_name in band_index:
        # Generate output name from the original filename
        out_name = f"{original_name.split('v2.0')[0]}v2.0_{index_name}.tif"
        
        out_path = os.path.join(out_folder,out_name)
        
        if skip_source and skip_source in files[0]:
            print("not reprocessing Sentinel data, only for Landsat..")
            continue
        
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

        # Extract UTM EPSG code for reprojection
        file_utm_zone = original_name.split(".")[2]
        utm_epsg_code = 32600 + int(file_utm_zone[1:3])

        # reproject if not epsg conform
        if spectral_index.rio.crs is not CRS.from_epsg(utm_epsg_code):
            spectral_index = spectral_index.rio.reproject(f"EPSG:{utm_epsg_code}")
        
        # Export to COG tiffs
        spectral_index.rio.to_raster(raster_path = out_path, 
                                     driver = 'COG')
    
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

    if platform.system() == "Windows":
        DATA_FOLDER = 'data/' # on local machine
        HLS_PARENT_PATH = "data/raster/hls/test"
    else:
        DATA_FOLDER = '~/data/' # on sciencecluster
        HLS_PARENT_PATH = "/home/nrietz/scratch/raster/hls/" # Set original data paths
    
    # Load Processing look-up-table to match UTM tiles to fire perimeter IDs
    processing_lut = pd.read_csv(
        os.path.join(DATA_FOLDER,"tables/processing_LUT.csv"),
        index_col=0)
    
    # Create date objects
    processing_lut['tst_date'] = processing_lut.apply(
        lambda row: datetime.date(row['tst_year'], row['tst_month'], row['tst_day']), axis=1)
    processing_lut['ted_date'] = processing_lut.apply(
        lambda row: datetime.date(row['ted_year'], row['ted_month'], row['ted_day']), axis=1)

    # Define bands/indices to process
    band_index = ast.literal_eval(sys.argv[2])
    
    print("Band metrcis parsed from slurm:",band_index)
    
    # band_index = ["NDVI","NDMI"]

    # Overwrite existing tiles?
    OVERWRITE_DATA = False
    
    # define chunk size for data loading
    chunk_size = dict(band=1, x=3600, y=3600)

    hls_granules_paths = get_hls_tif_list(HLS_PARENT_PATH)

    # Load UTM tiles
    UTM_TILE_FILE = sys.argv[1]
    
    with open(UTM_TILE_FILE, "r") as file:
        UTM_TILE_LIST = [line.strip() for line in file] 
        
    # UTM_TILE_LIST = ["54WXE","53WMU"] # for testing
    
    # Execute preprocessing for each UTM tile
    # ----
    for UTM_TILE_NAME in UTM_TILE_LIST:
        print("Preprocessing HLS data for UTM tile: ", UTM_TILE_NAME)
        
        if UTM_TILE_NAME:
            hls_granules_paths = [
                sublist for sublist in hls_granules_paths if sublist and UTM_TILE_NAME in sublist[0]
                ]
        
        # Execute preprocessing for each year
        # ----
        
        # Filter fire perimeters that use this UTM tile
        fire_perimeters_in_utm = processing_lut.loc[processing_lut.opt_UTM_tile == UTM_TILE_NAME]
        
        # extract years of burning (one UTM tile can have multiple years)
        utm_years = fire_perimeters_in_utm.tst_year.unique()
        
        # Loop through each year
        for year in utm_years:
            
            # For burn severity raster processing -
            # Create time range from LUT to subset HLS granules temporally
            if any(pattern in band_index for pattern in ["NBR","GEMI"]):
                # 15 d maximum offset from end of fire to search HLS imagery
                MAX_TIMEDELTA = 31
                
                # Get latest end date of fire in that year
                latest_tile_fire_end = max(fire_perimeters_in_utm.ted_date)
                
                # Define time search window
                # START_DATE = "2020-05-01T00:00:00" # full growing season
                START_DATE_POST = latest_tile_fire_end.strftime("%Y-%m-%dT%H:%M:%S")
                START_DATE_PRE = (latest_tile_fire_end -
                                  pd.Timedelta(days=MAX_TIMEDELTA) -
                                  relativedelta(years=1)).strftime("%Y-%m-%dT%H:%M:%S")

                # END_DATE = "2020-10-15T23:59:59" # full growing season
                END_DATE_POST = (latest_tile_fire_end +
                                 pd.Timedelta(days=MAX_TIMEDELTA)).strftime("%Y-%m-%dT%H:%M:%S")
                END_DATE_PRE = (latest_tile_fire_end +
                                pd.Timedelta(days=MAX_TIMEDELTA)-
                                relativedelta(years=1)).strftime("%Y-%m-%dT%H:%M:%S")

                hls_granules_paths_post = search_files_by_doy_range(hls_granules_paths,
                                                                    START_DATE_POST, END_DATE_POST)
                hls_granules_paths_pre = search_files_by_doy_range(hls_granules_paths,
                                                                   START_DATE_PRE, END_DATE_PRE)
                
                hls_granules_paths = hls_granules_paths_pre + hls_granules_paths_post

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
            # OUT_FOLDER = '/data/nrietz/raster/hls/'
            OUT_FOLDER = '/scratch/nrietz/raster/hls/processed/'
            
            num_workers = (
                int(os.environ.get("SLURM_CPUS_PER_TASK"))
                if os.environ.get("SLURM_CPUS_PER_TASK")
                else os.cpu_count()
            )
            # num_workers = -1 # -1 = all are used, -2 all but one

            multiprocessing.set_start_method('spawn')
            
            Parallel(n_jobs=num_workers, backend='loky')(
                delayed(joblib_hls_preprocessing)(files,band_index,bit_nums,OUT_FOLDER,OVERWRITE_DATA,skip_source=None) for files in hls_granules_paths)