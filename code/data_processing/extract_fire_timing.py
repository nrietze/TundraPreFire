import os
import re
from glob import glob
import platform
import geopandas as gpd
import pandas as pd
import datetime
import numpy as np
from tqdm import tqdm
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

    print("Performing spatial join of sample points and sub-daily perimeters.")
    for filepath in tqdm(df_viirs_filtered['filepath']):
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
if platform.system() == "Windows":
    DATA_FOLDER = 'data/' # on local machine
else:
    DATA_FOLDER = '/home/nrietz/data/' # on sciencecluster

PATH_VIIRS_PERIMETERS = os.path.join(DATA_FOLDER,"feature_layers/fire_atlas/")

# Load VIIRS perimeters in Siberian tundra
FN_VIIRS_CAVM_PERIMETERS = os.path.join(PATH_VIIRS_PERIMETERS,
                                        "viirs_perimeters_in_cavm_e113.gpkg")

if os.path.exists(FN_VIIRS_CAVM_PERIMETERS):
    print("CAVM Fire perimeter file exists, loading.")
    fire_perims_in_cavm = gpd.read_file(FN_VIIRS_CAVM_PERIMETERS)
else:
    print("Please prepare and filter the VIIRS fire perimeters using\n \"0_preprocess_ancillary_data.py\" ")

gpkg_list = glob(os.path.join(DATA_FOLDER,"feature_layers","*_sample_points.gpkg"))

#%% 2. Apply burn date extraction to all fire perimeters
for path in gpkg_list:
    fname = os.path.basename(path)
    fname_out = re.sub(r"\\*.gpkg$", r"_burn_date.gpkg", path)
    
    FIRE_ID = int(fname.split("_")[0])

    print(f"Extracting time of burn for fire: {FIRE_ID}")
    
    # Load the random points file for the currently selected fire scar
    random_points = gpd.read_file(path)
    
    # Load list of sub-daily perimeters and assign datetime
    FLIST_VIIRS_SUBDAILY = glob(os.path.join(PATH_VIIRS_PERIMETERS,
                                             "sub_daily/*/Snapshot/*M.gpkg"))
    df_viirs_sub_daily = pd.DataFrame(data={"filepath":FLIST_VIIRS_SUBDAILY})
    df_viirs_sub_daily['datetime'] = [extract_datetime(FN) for FN in FLIST_VIIRS_SUBDAILY]
    
    df_viirs_sub_daily = df_viirs_sub_daily.sort_values("datetime")
    
    perimeter = fire_perims_in_cavm.loc[fire_perims_in_cavm.fireid==FIRE_ID]
    
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
    
    # Load a fire perimeter file for CRS info
    sample_perimeters = gpd.read_file(df_viirs_filtered.iloc[0,0])
    
    # reproject random points to the sub-daily perimeter CRS
    random_points = random_points.to_crs(sample_perimeters.crs)
    
    # apply spatial join
    random_points_sjoin = spatial_join_multiple_files(
        random_points,
        df_viirs_filtered,
        FIRE_ID)
    
    # Create start DOY column
    random_points_sjoin['burn_date'] = pd.to_datetime(
        random_points_sjoin[['ted_year', 'ted_month', 'ted_day']].
        rename(columns={'ted_year': 'year', 
                        'ted_month': 'month', 
                        'ted_day': 'day'}))
    random_points_sjoin['burn_doy'] = random_points_sjoin['burn_date'].dt.dayofyear
    
    # export data
    random_points_sjoin.to_file(fname_out)
