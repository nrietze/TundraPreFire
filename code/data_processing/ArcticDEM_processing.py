import os
import numpy as np
import geopandas as gpd
import subprocess
import shlex
import tarfile

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

def download_and_extract_arcticdem(tile_index, output_dir="data/raster/arcticDEM"):
    """Function to download ArcticDEM mosaic tiles with wget

    Args:
        tile_index (str: string of the supertile index.
        output_dir (str, optional): Output directory where DEM data is dummped. Defaults to "data/raster/arcticDEM".
    """
    url = f'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v4.1/2m/{tile_index}/'

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

#%% 1. Load data
DATA_FOLDER = '/data/nrietz/' # on sciencecluster
# DATA_FOLDER = 'data/' # on local machine

# ArcticDEM mosaic tile index
index_features = gpd.read_file(os.path.join(DATA_FOLDER,
                                            "feature_layers/ArcticDEM_Mosaic_Index_v4_1_gpkg.gpkg"))

# Load fire perimeters
fire_perimeters = gpd.read_file(os.path.join(DATA_FOLDER,
                                             "feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"))

# reproject to ArcticDEM tile crs
fire_perimeters_stereo = fire_perimeters.to_crs(index_features.crs)
    
# perform spatial intersect of centroids with CAVM
dem_tiles_intersect = gpd.overlay(index_features,fire_perimeters_stereo,how='intersection')

# get tile index list for download
print("Finding DEM tiles overlapping with fire perimeters.")
tile_indices = dem_tiles_intersect.supertile.unique()

# Write tile index list to file
fname_tile_list = os.path.join(DATA_FOLDER,"tables/ArcticDEM_tileid_file.txt")
with open(fname_tile_list,"w") as outfile:
  outfile.write('\n'.join(str(i) for i in tile_indices))

# Read tile list file
with open(fname_tile_list, "r") as tile_file:
    # lines = tile_file.readlines()
    tile_indices_fromfile = [line.replace("\n","") for line in tile_file]
  
for tile_index in tile_indices_fromfile[:2]:
  print("Downloading ArcticDEM data for tile:", tile_index)
  download_and_extract_arcticdem(tile_index,
                                 output_dir = "~/scratch/raster/arcticDEM")
