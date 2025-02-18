# %%
import os
from glob import glob
from tqdm import tqdm
from osgeo import gdal
import rasterio as rio
from shapely.geometry import Point
import geopandas as gpd
import platform
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

if platform.system() == "Windows":
    DATA_FOLDER = 'data/' # on local machine
else:
    DATA_FOLDER = '~/data/' # on sciencecluster

# Load features (fire perimeters and ROIs)
fire_perimeters_2020 = gpd.read_file(os.path.join(DATA_FOLDER,
                                                  "feature_layers/fire_atlas/final_viirs2020.gpkg"))

# Load MGRS tile centroids to find msot suitable HLS tile
mgrs_tile_centroids = gpd.read_file(os.path.join(DATA_FOLDER,
                                                 "feature_layers/MGRS_centroids.geojson"))
# CAVM outline
cavm_outline = gpd.read_file(os.path.join(DATA_FOLDER,
                                          "feature_layers/cavm_outline.gpkg"))

# %% 2. Extract all VIIRS fire events in CAVM and east of 113° E
fname_perims_in_cavm = os.path.join(DATA_FOLDER,
                                    "feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")

# List all final VIIRS perimeters
gpkg_files = glob(os.path.join(DATA_FOLDER,
                                "feature_layers/fire_atlas/final_viirs*.gpkg"))

fire_perimeters = []

# Grab all VIIRS perimeters and add to list (only perimeters since 2016, e.g., 4:)
for gpkg_file in gpkg_files[4:]:
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
centroids_in_cavm = gpd.overlay(fire_centroids, 
                                    cavm_outline, 
                                    how='intersection')

# select perimeters in CAVM zone
fire_perims_in_cavm = fires_east[fires_east['UniqueID'].
                                isin(centroids_in_cavm.UniqueID.values)]


# Assign optimal UTM tile name to the fire perimeter
fire_perims_in_cavm["opt_UTM_tile"] = fire_perims_in_cavm.apply(lambda row :
    find_optimal_utm_tile(mgrs_tile_centroids, row), axis=1)

# Reset index
fire_perims_in_cavm = fire_perims_in_cavm.reset_index(drop=True)
    
# Create final format of look-up-table
final_lut = fire_perims_in_cavm.drop(columns = ['mergid', 'n_pixels', 'farea', 'fperim', 
                                             'duration','lcc_final', 'geometry',
                                             'centroid'])
final_lut[["descals_file"]] = "n"

# %% 2. Add Descals burned area tiles to final look-up-table
tiff_list = glob(os.path.join(DATA_FOLDER,"raster/burned_area_descals/*.tif"))

# Get centroids of descals tiles and format to geopandas data series
coordinates = [rio.open(file).lnglat() for file in tiff_list]

points = [Point(x, y) for x, y in coordinates]

descals_tile_centroids = gpd.GeoSeries(points)
descals_tile_centroids = gpd.GeoDataFrame(
    geometry=points,
    data={
        'filename': tiff_list,
        'longitude': [p.x for p in points],
        'latitude': [p.y for p in points]
    }
)

# Match Descals burned area tiles and ArcticDEM files to fire perimeters
for row in tqdm(fire_perims_in_cavm.iterrows(), total = len(fire_perims_in_cavm)):
    i = row[0]
    vect = row[1]

    # Calculate fire perimeter centroid
    fire_centroid = vect.geometry.centroid
    
    # Find nearest Descals tile centroid to fire perimeter (in EPSG:4326)
    sindex = descals_tile_centroids.geometry.sindex.nearest(fire_centroid)
    
     # select that nearest tile
    nearest_descals_tile_centroid = descals_tile_centroids.iloc[sindex[1,:]]
    
    final_lut.loc[final_lut.fireid == vect.fireid, 'descals_file'] = os.path.basename(nearest_descals_tile_centroid.filename.item())
    
# Export dataframe as csv
final_lut.to_csv(os.path.join(DATA_FOLDER,"tables/processing_LUT.csv"))   