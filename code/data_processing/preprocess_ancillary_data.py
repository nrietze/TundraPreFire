# %%
import os
from glob import glob

import geopandas as gpd
import numpy as np
import pandas as pd

# %% 0. Configure stuff

def find_optimal_utm_tile(mgrs_tile_centroids: gpd.GeoDataFrame,
                          fire_polygon: gpd.GeoDataFrame) -> str:
    """
    Find best UTM tile for the feature by identifying closest tile centroid to feature centroid

    Args:
        mgrs_tile_centroids (gpd.GeoDataFrame): feature layer with all UTM tile centroids. The UTM tile glossary can be downloaded from https://hls.gsfc.nasa.gov/wp-content/uploads/2016/03/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml
        fire_polygon (gpd.GeoDataFrame): One final fire polygon (from the VIIRS fire atlas) to which the nearest UTM tile will be searched

    Returns:
        str: Name of the optimal UTM tile for that fire polygon.
    """
    # get fire polygon's centroid
    fire_centroid = fire_polygon.geometry.centroid
    
    # get index of the nearest UTM centroid to the fire centroid
    sindex = mgrs_tile_centroids.geometry.sindex.nearest(fire_centroid)
    
    # select that nearest tile
    nearest_mgrs_tile_centroid = mgrs_tile_centroids.iloc[sindex[1],:] 

    # reportt the tilename
    return nearest_mgrs_tile_centroid.Name.item()

# %% 1. Load data
DATA_FOLDER = '/data/nrietz/'

# Load features (fire perimeters and ROIs)
fire_perimeters_2020 = gpd.read_file(os.path.join(DATA_FOLDER,
                                                  "feature_layers/fire_atlas/final_viirs2020.gpkg"))

# test region of interest
roi = gpd.read_file(os.path.join(DATA_FOLDER,
                                 "feature_layers/roi.geojson"))

# Load MGRS tile centroids to find msot suitable HLS tile
mgrs_tile_centroids = gpd.read_file(os.path.join(DATA_FOLDER,
                                                 "feature_layers/MGRS_centroids.geojson"))

# Intersect
fire_within_roi = gpd.overlay(fire_perimeters_2020, roi, how='intersection')
fire_within_roi.head()

# CAVM outline
cavm_outline = gpd.read_file(os.path.join(DATA_FOLDER,
                                          "feature_layers/cavm_outline.gpkg"))

# %% 2. Extract all VIIRS fire events in CAVM and east of 113° E
fname_perims_in_cavm = os.path.join(DATA_FOLDER,
                                    "feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")

# Apply spatial intersection but check if feature layer exists
if os.path.exists(fname_perims_in_cavm):
    print("File of fire perimeters in CAVM exists, loading.")
    merged_fire_perimeters = gpd.read_file(fname_perims_in_cavm)
else:
    print("File of fire perimeters in CAVM doesn't exist, processing.")
    # List all final VIIRS perimeters
    gpkg_files = glob(os.path.join(DATA_FOLDER,
                                   "feature_layers/fire_atlas/final_viirs*.gpkg"))

    fire_perimeters = []
    
    # Grab all VIIRS perimeters and add to list
    for gpkg_file in gpkg_files:
        gdf = gpd.read_file(gpkg_file)
        fire_perimeters.append(gdf)
    
    # concat to one dataframe
    merged_fire_perimeters = gpd.GeoDataFrame(pd.concat(fire_perimeters,
                                                        ignore_index=True))
    
    # add UniqueID column (some FireIDs are duplicates despite large distances)
    merged_fire_perimeters["UniqueID"] = range(len(merged_fire_perimeters))
   
    # extract fires east of 113°N (too many flase detections towards W Siberia with industry)
    fires_east = merged_fire_perimeters[
        merged_fire_perimeters.geometry.apply(lambda x: x.centroid.x > 113)]
    fires_east['centroid'] = fires_east.geometry.centroid
    
    # Extract centroids for intersection (faster and avoids holes in polygon intersections)
    fire_centroids = gpd.GeoDataFrame(fires_east[['UniqueID', 'centroid']],
                                      geometry='centroid')
    fire_centroids = fire_centroids.to_crs(cavm_outline.crs)
    
    # perform spatial intersect of centroids with CAVM
    centroids_within_cavm = gpd.overlay(fire_centroids, 
                                        cavm_outline, 
                                        how='intersection')

    # select perimeters in CAVM zone
    fire_within_cavm = fires_east[fires_east['UniqueID'].
                                  isin(centroids_within_cavm.UniqueID.values)]
    
    
    # Assign optimal UTM tile name to the fire perimeter
    fire_within_cavm["opt_UTM_tile"] = fire_within_cavm.apply(lambda row :
        find_optimal_utm_tile(mgrs_tile_centroids, row), axis=1)
    
    # export filtered fire perimeters to file ( exclude "centroids" column from export)
    fire_within_cavm.drop(columns=["centroid"]).to_file(fname_perims_in_cavm)
    
    # Export uniwue tile names as a txt file for HLS dowloading
    unique_utm_tile_ids = fire_within_cavm.opt_UTM_tile.unique()
    np.savetxt("code/data_processing/tileid_file.txt",
               unique_utm_tile_ids,fmt="%s")
# %%
