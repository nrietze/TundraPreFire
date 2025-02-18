import os
import subprocess
import platform
from osgeo import gdal
import geopandas as gpd
import pandas as pd
import shlex
import tarfile
from glob import glob
import multiprocessing
from joblib import Parallel, delayed
import re

#%% 0. configure functions
def find_optimal_dem_tile(dem_tile_centroids: gpd.GeoDataFrame,
                          fire_polygon: gpd.GeoDataFrame) -> str:
    """
    Find best ArcticDEM mosaic tile for the feature by identifying closest tile centroid to feature centroid

    Args:
        dem_tile_centroids (gpd.GeoDataFrame): feature layer with all ArcticDEM tile centroids. 
        fire_polygon (gpd.GeoDataFrame): One final fire polygon (from the VIIRS fire atlas) to which the nearest UTM tile will be searched

    Returns:
        str: Name of the optimal UTM tile for that fire polygon.
    """
    # get fire polygon's centroid
    fire_centroid = fire_polygon.geometry.centroid
    
    # get index of the nearest UTM centroid to the fire centroid
    sindex = dem_tile_centroids.geometry.sindex.nearest(fire_centroid)
    
    # select that nearest tile
    nearest_mgrs_tile_centroid = dem_tile_centroids.iloc[sindex[1],:] 

    # reportt the tilename
    return nearest_mgrs_tile_centroid.Name.item()

def download_and_extract_arcticdem(url, output_dir="data/raster/arcticDEM"):
    """Function to download ArcticDEM mosaic tiles with wget

    Args:
        tile_index (str: string of the supertile index.
        output_dir (str, optional): Output directory where DEM data is dummped. Defaults to "data/raster/arcticDEM".
    """
    # wget command
    wget_command = f"wget -r -N -nH -np -R index.html* --cut-dirs=8 -P {shlex.quote(output_dir)} {shlex.quote(url)}"
    
    # Run the wget command
    try:
        subprocess.run(wget_command, shell=True, check=True)  # shell=True is necessary here because of the wildcard
        print(f"Successfully downloaded files from {url} to {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading files from {url}: {e}")
        return 
    
    # Extract tar.gz files
    for filename in os.listdir(output_dir):
        if filename.endswith(".tar.gz"):
            filepath = os.path.join(output_dir, filename)
            print(f"Extracting: {filepath}")
            try:
                with tarfile.open(filepath, 'r:gz') as tar:
                    tar.extractall(path=output_dir)
                print(f"Successfully extracted {filepath} to {output_dir}")

                #Optionally remove the tar.gz file after extraction to save space
                os.remove(filepath)
                print(f"Removed {filepath}")

            except tarfile.ReadError as e:
                 print(f"Error extracting {filepath}: {e}")

# Function to extract and rename the s3 url in the ArcticDEM LUT
def transform_url(url):
    match = re.search(r'external/(.+)', url)
    return f"https://{match.group(1).replace('.json', '_dem.tif')}" if match else None

#%% 1. Load data
if platform.system() == "Windows":
    DATA_FOLDER = 'data/' # on local machine
else:
    DATA_FOLDER = '~/data/' # on sciencecluster

DOWNLOAD_DATA = True

# Load fire perimeters
fire_perimeters = gpd.read_file(
    os.path.join(DATA_FOLDER,
                 "feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")
    )

fname_dem_lut = os.path.join(DATA_FOLDER,"tables/ArcticDEM_LUT.csv")

if not os.path.exists(fname_dem_lut):
    print("Matching ArcticDEM tiles with fire perimeters")
  
    # ArcticDEM mosaic tile index
    arctic_dem_tile_index = gpd.read_file(
        os.path.join(DATA_FOLDER,
                     "feature_layers/ArcticDEM_Mosaic_Index_v4_1.gpkg")
        )

    # reproject to ArcticDEM tile crs (Polar Stereographic)
    fire_perimeters_stereo = fire_perimeters.to_crs(arctic_dem_tile_index.crs)
       
    # perform spatial intersect of centroids with CAVM
    dem_tiles_intersect = gpd.overlay(arctic_dem_tile_index,fire_perimeters_stereo,
                                    how='intersection')

    # Selecting only neccessary columns for export to LUT csv
    dem_tile_df = dem_tiles_intersect[["dem_id","tile","supertile",
                                       "s3url","fileurl","fireid","UniqueID"]]

    # Export dataframe as csv
    dem_tile_df.to_csv(fname_dem_lut)

if DOWNLOAD_DATA:
    # Read tile list file
    dem_tile_df = pd.read_csv(fname_dem_lut,index_col=0)

    s3_urls = dem_tile_df["s3url"].apply(transform_url)
    
    # Download ArcticDEM rasters
    num_workers = (
            int(os.environ.get("SLURM_CPUS_PER_TASK"))
            if os.environ.get("SLURM_CPUS_PER_TASK")
            else os.cpu_count()
        )

    multiprocessing.set_start_method('spawn')

    Parallel(n_jobs=num_workers, backend='loky')(
        delayed(download_and_extract_arcticdem)(url,
                                                output_dir = "/scratch/nrietz/raster/arcticDEM") for url in s3_urls)