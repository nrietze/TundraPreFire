# %%
import os
import multiprocessing
import geopandas as gpd
from osgeo import gdal
import xarray as xr
import rioxarray as rxr
import rasterio as rio
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

def filter_granules(granules, 
                    UTM_TILEID: str,
                    MAX_CLOUD_COVER: int,
                    MIN_TILE_COVER: int):
    """Filter granules based on cloud and spatial coverage criteria."""
    filtered_granules = []
    
    # Filter granules first by UTM tile ID
    utm_filtered_granules = [
        sublist for sublist in granules if UTM_TILEID in sublist[0]
    ]
    
    for granule_files in utm_filtered_granules:
        sample_file = granule_files[0]  

        # Extract metadata without loading the raster
        cloud_coverage, spatial_coverage = get_raster_metadata(sample_file)
        
        # Apply filtering conditions
        if cloud_coverage <= MAX_CLOUD_COVER and spatial_coverage >= MIN_TILE_COVER:
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

# Filter tiles for < 20 % cloud cover and > 20 % tile coverage
MIN_TILE_COVERAGE = 20
MAX_CLOUD_COVERAGE = 80

print(f"Filtering out HLS granules with more than {MAX_CLOUD_COVERAGE} % cloud cover \n and less than {MIN_TILE_COVERAGE} % tile coverage.")

filtered_hls_granules = filter_granules(hls_granules_paths,
                                        UTM_TILEID=utm_tileid,
                                        MAX_CLOUD_COVER=MAX_CLOUD_COVERAGE,MIN_TILE_COVER=MIN_TILE_COVERAGE)

print(f"Found {len(filtered_hls_granules)} good granules post-fire in UTM tile {utm_tileid}.")

# %% 3. Calculate NBR for the filtered granules
N_CORES = -1 # -1 = all are used, -2 all but one
OVERWRITE_DATA = False
OUT_FOLDER = '/data/nrietz/raster/hls/'
bit_nums = [0,1,2,3,4]

multiprocessing.set_start_method('spawn')

Parallel(n_jobs=N_CORES, backend='loky')(
    delayed(hlsPrep.joblib_hls_preprocessing)(files,
                                              ["GEMI","NBR"],
                                              bit_nums,
                                              OUT_FOLDER,OVERWRITE_DATA) for files in filtered_hls_granules)


# %% 4. Stack NBR imagery
PROCESSED_HLS_DIR = '/home/nrietz/data/raster/hls/'

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