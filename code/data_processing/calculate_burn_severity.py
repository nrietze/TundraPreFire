# %%
import os
import multiprocessing
import numpy as np
import geopandas as gpd
from osgeo import gdal
import xarray as xr
import rioxarray as rxr
import rasterio as rio
from rasterio.mask import mask
from shapely.geometry import mapping
import datetime
from joblib import Parallel, delayed

os.chdir("code/data_processing")
import HLS_preprocessing as hlsPrep
os.chdir("../..")

# %% 0. Configure functions
def get_raster_metadata(filepath):
    """Extract metadata for cloud and spatial coverage from a raster file."""
    with rio.open(filepath) as src:
        # Extract raster metadata attributes lazily
        tags = src.tags()
        cloud_coverage = float(tags.get("cloud_coverage", 100))  # default to 100 if missing
        spatial_coverage = float(tags.get("spatial_coverage", 0))  # default to 0 if missing
        return cloud_coverage, spatial_coverage

def replace_string_in_filename(hls_filename, string):
    """Replace the second last part of the filename with a string."""
    out_string = f"{hls_filename.split('v2.0')[0]}v2.0_{string}.tif"
    
    return out_string

def calculate_clear_pixel_percentage(files:list,
                                     hls_processed_path:str,
                                     polygon:gpd.GeoDataFrame,
                                     bit_nums = [0,1,2,3,4]):
    """Calculate percentage of non-NaN, non-water pixels in the polygon."""
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
    
    image = hlsPrep.load_rasters(files, "RED",
                                 band_dict = hls_band_dict,
                                 region=polygon,
                                 chunk_size=dict(band=1, x=3600, y=3600))
    
    
    # sample_file = [file for file in files if "B03" in file][0]
    
    # # Load the sample raster
    # with rio.open(sample_file) as src:
    #     polygon_geom = polygon.to_crs(src.crs).geometry.iloc[0]
        
    #     # Crop to the polygon area
    #     image, _ = mask(src, [polygon_geom], crop=True)
    #     nodata_value = src.nodata if src.nodata is not None else np.nan
    
    # # Load corresponding Fmask
    # fmask_file = sample_file.replace('B03', 'Fmask')
    # with rio.open(fmask_file) as fmask_src:
    #     fmask, _ = mask(fmask_src, [polygon_geom], crop=True)
    
    fmask = hlsPrep.load_rasters(files, "FMASK",
                                 band_dict = hls_band_dict,
                                 region=polygon,
                                 chunk_size=dict(band=1, x=3600, y=3600))
    
    
    # Load corresponding water mask
    water_mask_file = os.path.join(hls_processed_path,
                                    os.path.basename(
                                        replace_string_in_filename(files[0], "watermask")))
    
    with rio.open(water_mask_file) as watermask_src:
        polygon_geom = polygon.to_crs(watermask_src.crs).geometry.iloc[0]
        water_mask, _ = mask(watermask_src, [mapping(polygon_geom)], crop=True)

    # Convert Fmask to binary
    fmask_valid = hlsPrep.create_quality_mask(fmask.data, bit_nums)
    
    # Flatten arrays for easier pixel counting
    image_flat = image.data.flatten()
    fmask_valid_flat = fmask_valid.flatten()
    water_mask_flat = water_mask.flatten()

    # Mask water and NaN values
    valid_pixels = np.isfinite(image_flat) 
    non_water_pixels = water_mask_flat == 0  # Pixels that are not water

    # Combine valid and non-water conditions
    clear_pixels = valid_pixels & non_water_pixels & (fmask_valid_flat == 1)

    # Compute percentage of clear pixels
    percentage_clear = (np.sum(clear_pixels) / np.sum(valid_pixels)) * 100 if np.sum(valid_pixels) > 0 else 0
    
    return percentage_clear

def filter_granules(granules:list,
                    polygon:gpd.GeoDataFrame,
                    hls_processed_path:str,
                    UTM_TILEID: str,
                    MAX_CLOUD_COVER: int,
                    MIN_TILE_COVER: int,
                    MIN_VALID_COVER: int):
    """
    Function that filters HLS granules on cloud and tile coverage, and the percent coverage of valid data over the fire perimeter. 

    Args:
        granules (list): _description_
        polygon (gpd.GeoDataFrame): _description_
        UTM_TILEID (str): _description_
        MAX_CLOUD_COVER (int): _description_
        MIN_TILE_COVER (int): _description_
        MIN_VALID_COVER (int): _description_

    Returns:
        _type_: _description_
    """
    filtered_granules = []
    
    # Filter granules first by UTM tile ID
    utm_filtered_granules = [
        sublist for sublist in granules if UTM_TILEID in sublist[0]
    ]
    
    for granule_files in utm_filtered_granules:

        # Extract metadata without loading the raster
        cloud_coverage, spatial_coverage = get_raster_metadata(granule_files[0])
        
        clear_percentage = calculate_clear_pixel_percentage(granule_files, 
                                                            hls_processed_path,
                                                            polygon)
        
        # Apply filtering conditions
        if cloud_coverage <= MAX_CLOUD_COVER and \
            spatial_coverage >= MIN_TILE_COVER and \
                clear_percentage >= MIN_VALID_COVER:
            filtered_granules.append(granule_files)

    return filtered_granules

def time_index_from_filenames(files):
    '''
    Helper function to create a pandas DatetimeIndex
    '''
    return [datetime.datetime.strptime(f.split('.')[-4], '%Y%jT%H%M%S') for f in files]

# %% 1. Load data
PATH_VIIRS_PERIMETERS = "../data/feature_layers/fire_atlas/"

# Load VIIRS perimeters in Siberian tundra
FN_VIIRS_CAVM_PERIMETERS = os.path.join(PATH_VIIRS_PERIMETERS,
                                        "viirs_perimeters_in_cavm_e113.gpkg")

if os.path.exists(FN_VIIRS_CAVM_PERIMETERS):
    print("Loading VIIRS fire perimeters in CAVM zone.")
    fire_within_cavm = gpd.read_file(FN_VIIRS_CAVM_PERIMETERS)
else:
    print("Please prepare and filter the VIIRS fire perimeters using\n \"0_preprocess_ancillary_data.py\" ")

HLS_PARENT_PATH = "/home/nrietz/scratch/raster/hls/"
PROCESSED_HLS_DIR = '/home/nrietz/data/raster/hls/'

TEST_ID = 14211 # fire ID for part of the large fire scar

# %% 2. Find HLS imagery pre- and post-fire

# Select individual fire perimeter from ID
perimeter = fire_within_cavm.loc[fire_within_cavm.fireid==TEST_ID]

# get the end of burn date
end_date = datetime.datetime(perimeter.ted_year.item(),
                             perimeter.ted_month.item(),
                             perimeter.ted_day.item())

season_end = datetime.datetime(perimeter.ted_year.item(),10,31) # end of season Oct 31

# retrieve optimal UTM tile ID
utm_tileid = perimeter.opt_UTM_tile.unique().item()
utm_tileid = "54WXE"

# Get list of HLS tiles for this fire perimeter
hls_granules_paths = hlsPrep.get_hls_tif_list(HLS_PARENT_PATH)

# Filter for post-fire period
hls_granules_paths = hlsPrep.search_files_by_doy_range(
    hls_granules_paths,
    end_date.strftime("%Y-%m-%dT%H:%M:%S"), 
    season_end.strftime("%Y-%m-%dT%H:%M:%S"))

# Filter tiles for < 20 % cloud cover, > 20 % tile coverage, and > 20% of valid data
MIN_VALID_COVERAGE = 80
MIN_TILE_COVERAGE = 20
MAX_CLOUD_COVERAGE = 80

print(f"Filtering out HLS granules with more than {MAX_CLOUD_COVERAGE} % cloud cover\
     \n and less than {MIN_TILE_COVERAGE} % tile coverage.")

filtered_hls_granules = filter_granules(hls_granules_paths,
                                        polygon = perimeter,
                                        hls_processed_path=PROCESSED_HLS_DIR,
                                        UTM_TILEID=utm_tileid,
                                        MAX_CLOUD_COVER=MAX_CLOUD_COVERAGE,MIN_TILE_COVER=MIN_TILE_COVERAGE,
                                        MIN_VALID_COVER=MIN_VALID_COVERAGE)

print(f"Found {len(filtered_hls_granules)} good granules post-fire in UTM tile {utm_tileid}.")

# %% 3. Calculate NBR for the filtered granules
# Create list of unprocessed HLS granules for which no NBR file exists
granules_without_nbr = [
    granule for granule in filtered_hls_granules
    if not os.path.exists(replace_string_in_filename(granule[0],'NBR'))  
]

# Calculate NBR and GEMI for the granules without NBR data
print("Start calculation of NBR and GEMI for HLS granules missing that data.")

N_CORES = -1 # -1 = all are used, -2 all but one
OVERWRITE_DATA = False
OUT_FOLDER = '/data/nrietz/raster/hls/'
bit_nums = [0,1,2,3,4]

multiprocessing.set_start_method('spawn')

Parallel(n_jobs=N_CORES, backend='loky')(
    delayed(hlsPrep.joblib_hls_preprocessing)(files,
                                            ["GEMI","NBR"],
                                            bit_nums,
                                            OUT_FOLDER,OVERWRITE_DATA) for files in granules_without_nbr)
print("NBR and GEMI processing complete.")

# %% 4. Stack NBR imagery
# List NBR COGs
nbr_files = [
    os.path.abspath(os.path.join(PROCESSED_HLS_DIR, o)) 
    for o in os.listdir(PROCESSED_HLS_DIR)  # list dir is not recursive, so dNBR subfolder will not be searched :)
    if o.endswith('NBR.tif') and utm_tileid in o # find files that are in the tile and are NBR
    ]

print(f"There are {len(nbr_files)} NBR files for this feature.")

time = xr.Variable('time', time_index_from_filenames(nbr_files))

# define chunk size for data loading
chunk_size = dict(band=1, x=3600, y=3600)

# create image stack of post-fire NBR
nbr_ts = xr.concat([rxr.open_rasterio(f, mask_and_scale=True, 
                                      chunks=chunk_size).squeeze('band', drop=True) for f in nbr_files], dim=time)
nbr_ts.name = 'NBR'