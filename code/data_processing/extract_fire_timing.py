import os
from glob import glob

import geopandas as gpd
import pandas as pd
import datetime
import numpy as np

from osgeo import gdal
import xarray as xr
import rioxarray as rxr
import rasterio as rio


# %% 0. Configure functions
def extract_datetime(PATH):
    FN = os.path.basename(PATH)
    datestr = FN.split(".")[0]
    return datetime.datetime.strptime(datestr, "%Y%m%d%p")

def spatial_join_multiple_files(random_points: gpd.GeoDataFrame,
                                df_viirs_filtered: pd.DataFrame,
                                perimeter_id: str) -> gpd.GeoDataFrame:
    """
    Perform spatial joins between random points and polygons from filepaths,
    makes sure that each point is joined only once.
    
    Parameters:
    - random_points: GeoDataFrame of points
    - df_viirs_filtered: DataFrame with a 'filepath' column containing polygon files
    - perimeter_id: ID used to filter each polygon layer
    
    Returns:
    - GeoDataFrame containing the points with polygon information attached
    """
    # Create column to track if assignment was successful
    random_points['joined'] = False
    
    result_list = []

    for filepath in df_viirs_filtered['filepath']:
        # Read sub-daily perimeter file
        sd_perimeters = gpd.read_file(filepath)
        
        # Filter by perimeter_id
        sd_perimeter_filtered = sd_perimeters[sd_perimeters.fireid == perimeter_id]
        
        if sd_perimeter_filtered.empty:
            continue
        
        # Select only points that haven't been joined yet
        points_to_join = random_points.loc[~random_points['joined']]

        # Apply spatial join
        point_in_poly = gpd.sjoin(points_to_join, sd_perimeter_filtered, predicate='within')
        
        # Mark points as joined
        if not point_in_poly.empty:
            random_points.loc[point_in_poly.index, 'joined'] = True

        # Append to results
        result_list.append(point_in_poly)

    # Combine results
    final_result = pd.concat(result_list, ignore_index=True) if result_list else gpd.GeoDataFrame()
    
    # Drop 'joined' column
    final_result = final_result.drop(columns=['joined'], errors='ignore')
    
    return final_result

#%% 1. Load data
PATH_VIIRS_PERIMETERS = "../data/feature_layers/fire_atlas/"

# Load VIIRS perimeters in Siberian tundra
FN_VIIRS_CAVM_PERIMETERS = os.path.join(PATH_VIIRS_PERIMETERS,"viirs_perimeters_in_cavm_e113.gpkg")

if os.path.exists(FN_VIIRS_CAVM_PERIMETERS):
    print("CAVM Fire perimeter file exists, loading.")
    merged_fire_perimeters = gpd.read_file(FN_VIIRS_CAVM_PERIMETERS)
else:
    print("Please extract VIIRS fire perimeters in CAVM extent first.")
    pass

TEST_ID = 14211 # fire ID for part of the large fire scar

# Load the random points file for the currently selected fire sacr
random_points = gpd.read_file("../data/feature_layers/random_points.shp")

# %%
# Load list of sub-daily perimeters and assign datetime
FLIST_VIIRS_SUBDAILY = glob(os.path.join(PATH_VIIRS_PERIMETERS,
                                         "sub_daily/*/Snapshot/*M.gpkg"))
df_viirs_sub_daily = pd.DataFrame(data={"filepath":FLIST_VIIRS_SUBDAILY})
df_viirs_sub_daily['datetime'] = [extract_datetime(FN) for FN in FLIST_VIIRS_SUBDAILY]

df_viirs_sub_daily = df_viirs_sub_daily.sort_values("datetime")

perimeter = merged_fire_perimeters.loc[merged_fire_perimeters.fireid==TEST_ID]
start_date = datetime.datetime(perimeter.tst_year.item(),
                               perimeter.tst_month.item(),
                               perimeter.tst_day.item())

end_date = datetime.datetime(perimeter.ted_year.item(),
                             perimeter.ted_month.item(),
                             perimeter.ted_day.item())

# filter sub-daily perimeters that are within the fire's start and end dates
cond = np.logical_and(df_viirs_sub_daily.datetime >= start_date,
                      df_viirs_sub_daily.datetime <= end_date)
df_viirs_filtered = df_viirs_sub_daily[cond]
# %%
# Load sample perimeter file for CRS info
sample_perimeters = gpd.read_file(df_viirs_filtered.iloc[0,0])

# reproject random points to the sub-daily perimeter CRS
random_points = random_points.to_crs(sample_perimeters.crs)

# apply spatial join
random_points_sjoin = spatial_join_multiple_files(
    random_points,
    df_viirs_filtered,
    TEST_ID)