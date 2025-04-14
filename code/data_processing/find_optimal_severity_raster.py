# %%
import os
import sys
import platform
import datetime
import multiprocessing
from glob import glob
from tqdm import tqdm
import numpy as np
import pandas as pd
import geopandas as gpd
from osgeo import gdal
import xarray as xr
import rioxarray as rxr
import rasterio as rio
from rasterio.mask import mask
from shapely.geometry import mapping
from dateutil.relativedelta import relativedelta
from joblib import Parallel, delayed

# %% 0. Configure functions
def calculate_clear_pixel_percentage(data,
                                     polygon:gpd.GeoDataFrame):
    """Calculate percentage of non-NaN, non-water pixels in the polygon.

    Args:
        data: Can be a filename string or xr.DataArray of the raster to be looked up.
        polygon (gpd.GeoDataFrame): The polygon of the fire perimeter that is used to calculate the pixel percentage.

    Returns:
        int: Valid pixel percentage within the polygon.
    """
    # Apply function to existing raster with rasterio
    if isinstance(data, str):
        # Load the sample raster
        with rio.open(data) as src:
            geometry = gpd.GeoSeries([polygon.geometry.iloc[0]],crs = 4326)
            polygon_geom = geometry.to_crs(src.crs).geometry.iloc[0]
        
            # Crop to the polygon area
            image, _ = mask(src, [polygon_geom], crop=True, nodata = np.nan)

    # Apply function to memory loaded xarray
    elif isinstance(data, xr.DataArray):
        # Reproject polygon to match xarray CRS
        crs = data.rio.crs
        geometry = gpd.GeoSeries([polygon.geometry], crs=4326).to_crs(crs)
        
        # Clip the xarray to the polygon
        clipped = data.rio.clip(geometry, all_touched=True)
        image = clipped.values
    
    # Flatten array for easier pixel counting
    image_flat = image.flatten()

    non_nan_count = np.count_nonzero(~np.isnan(image_flat))
    
    # Calculate the total number of elements
    total_count = image_flat.size
    
    # Calculate the percentage of non-NaN values
    percentage_clear = (non_nan_count / total_count) * 100 if np.sum(non_nan_count) > 0 else 0
    
    return percentage_clear


def calculate_optimality_statistics(data, polygon:gpd.GeoDataFrame, stat='mean'):
    """Calculate summary statistics of the optimality raster non-NaN pixels in the polygon.

    Args:
        data: Filename string of the raster to be looked up.
        polygon (gpd.GeoDataFrame): The polygon of the fire perimeter that is used to calculate the pixel percentage.

    Returns:
        int: Valid pixel percentage within the polygon.
    """
    # Load the sample raster
    with rio.open(data) as src:
        geometry = gpd.GeoSeries([polygon.geometry.iloc[0]],crs = 4326)
        polygon_geom = geometry.to_crs(src.crs).geometry.iloc[0]
    
        # Crop to the polygon area
        image, _ = mask(src, [polygon_geom], crop=True, nodata = np.nan)

    # Flatten array for easier pixel counting
    image_flat = image.flatten()

    stat_value = np.nanmean(image_flat) if stat == 'mean' else np.nanmedian(image_flat)
    
    return stat_value

# %% Execute pararellized processing
if __name__ == "__main__":
    print("Loading data...")

    OVERWRTIE_DATA = True
    
    if platform.system() == "Windows":
        DATA_FOLDER = 'data/' # on local machine
        HLS_FOLDER = './data/raster/hls'
    else:
        DATA_FOLDER = '/home/nrietz/data/' # on sciencecluster
        HLS_FOLDER = "/home/nrietz/scratch/raster/hls"
    
    # Load Processing look-up-table to match UTM tiles to fire perimeter IDs
    processing_lut = pd.read_csv(
        os.path.join(DATA_FOLDER,"tables/processing_LUT.csv"),
        index_col=0)
    
    # Create date objects
    processing_lut['tst_date'] = processing_lut.apply(
        lambda row: datetime.date(row['tst_year'], row['tst_month'], row['tst_day']), axis=1)
    processing_lut['ted_date'] = processing_lut.apply(
        lambda row: datetime.date(row['ted_year'], row['ted_month'], row['ted_day']), axis=1)

    # Load VIIRS perimeters in Siberian tundra
    FN_VIIRS_CAVM_PERIMETERS = os.path.join(DATA_FOLDER,
                                            "feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")
    fire_polygons = gpd.read_file(FN_VIIRS_CAVM_PERIMETERS)
    
    # Only process HLS data for fires starting in 2017
    TEST_ID = [14211,14664,10792,17548]
    fire_polygons = fire_polygons.loc[(fire_polygons.tst_year>=2017) & (np.isin(fire_polygons.fireid,TEST_ID))]
    
    # Create output tables
    output_lut = pd.DataFrame(columns=[
        "fireid",
        "fname_severity_raster",
        "fname_optimality_raster",
        "date_severity_raster",
        "pct_clear_pixel",
        "mean_optimality",
        "severity_index"
    ])
    output_table = output_lut.copy()
    
    for i in fire_polygons.index:
        perimeter = fire_polygons.loc[[i]]
        
        # Set Time
        YEAR_START = perimeter.tst_year.item()
        YEAR_END = perimeter.ted_year.item()
        MONTH_START = 5
        MONTH_END = 10
        
        # Get attributes of this perimeter
        FIREID = perimeter.fireid.item()
        fire_perimeter_attrs = processing_lut.loc[processing_lut.fireid == FIREID]
        UTM_TILE_NAME = fire_perimeter_attrs.opt_UTM_tile.item()

        for severity_index in ["dNBR","dGEMI"]:
            
            # Get list of severity rasters for this fire perimeter
            burn_severity_files = glob(
                os.path.join(HLS_FOLDER,"severity_rasters",f"{severity_index}_{UTM_TILE_NAME}*.tif")
                )
            severity_dates = [os.path.basename(f).split('.')[0][-10:] for f in burn_severity_files]
            
            df_sev = pd.DataFrame({
                "fname_severity_raster":burn_severity_files,
                "date_severity_raster":severity_dates
            })
            
            # and optimality rasters
            optimality_files = glob(
                os.path.join(HLS_FOLDER,"optimality_rasters",f"*{UTM_TILE_NAME}*.tif")
                )
            optimality_dates = [os.path.basename(f).split('.')[0][-10:] for f in optimality_files]
            
            df_opt = pd.DataFrame({
                "fname_optimality_raster":optimality_files,
                "date_severity_raster":optimality_dates
            })
            
            df_temp = pd.merge(df_sev,df_opt, on="date_severity_raster")
            df_temp["fireid"] = FIREID
            df_temp["severity_index"] = severity_index
            
            if not burn_severity_files:
                print(f"No {severity_index} rasters found for fire {FIREID}.")
                continue
            
            print(f"Searching optimal {severity_index} raster for fire {FIREID}.")
            
            # da_burn_severity = ConcatRasters(burn_severity_files, severity_index)
            
            # Calculate clear pixel percentage for all burn severity files of this scar
            df_temp['pct_clear_pixel'] = [calculate_clear_pixel_percentage(fname, polygon = perimeter) for fname in tqdm(df_temp.fname_severity_raster)]
            
            # Calculate optimality summary statistics for all files of this scar
            df_temp['mean_optimality'] = [calculate_optimality_statistics(fname, polygon = perimeter, stat = "mean") for fname in tqdm(df_temp.fname_optimality_raster)]
            
            # Find optimal pair of optimality & clear pixel percentage (normalizing both to range [0, 1])
            norm_opt = (df_temp["mean_optimality"] - df_temp["mean_optimality"].min()) / (df_temp["mean_optimality"].max() - df_temp["mean_optimality"].min())
            norm_pct = (df_temp["pct_clear_pixel"] - df_temp["pct_clear_pixel"].min()) / (df_temp["pct_clear_pixel"].max() - df_temp["pct_clear_pixel"].min())

            # Get index of best pair with wieghted average
            score = 0.2 * norm_opt + 0.8 * norm_pct
            
            if all(np.isnan(score)):
                continue
            best_idx = np.argmax(score)
            
            df_best = df_temp.iloc[best_idx,:].to_frame().T 
            output_lut = pd.concat([df_best,output_lut])
            
            # output_lut.loc[output_lut["fireid"] == FIREID, "severity_index"] = severity_index
            # output_lut.loc[output_lut["fireid"] == FIREID, "pct_clear_pixel"] = df_temp['pct_clear_pixel'][best_idx]
            # output_lut.loc[output_lut["fireid"] == FIREID, "mean_optimality"] = df_temp['mean_optimality'][best_idx]
            # output_lut.loc[output_lut["fireid"] == FIREID, "severity_raster_filename"] = df_temp['fname_severity_raster'][best_idx]
            # output_lut.loc[output_lut["fireid"] == FIREID, "date_severity_raster"] = df_temp['date'][best_idx]
            
            output_table = pd.concat([df_temp,output_table])
    
    # Write to ouput tables
    output_lut.to_csv(os.path.join(DATA_FOLDER,"tables","optimality_LUT.csv"),sep=";")        
    output_table.to_csv(os.path.join(DATA_FOLDER,"tables","optimality_overview_table.csv"),sep=";")        