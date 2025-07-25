# %%
import os
import sys
import platform
import datetime
from glob import glob
from tqdm import tqdm
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
from osgeo import gdal
import xarray as xr
import rioxarray as rxr
import rasterio as rio
from rasterio.mask import mask

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
            geometry = gpd.GeoSeries([polygon.geometry.iloc[0]],crs = 3413)
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
        geometry = gpd.GeoSeries([polygon.geometry.iloc[0]],crs = 3413)
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
    fire_polygons = fire_polygons.loc[(fire_polygons.tst_year>=2017)]
    
    USE_TOP_25 = True
    if USE_TOP_25:
        fire_polygons = fire_polygons.sort_values("farea",ascending=False).head(25)
    
        print(f"Using top {len(fire_polygons)} largest fire polygons. \n")

    FALSE_FIRES_ID = [23633,21461,15231,15970,17473,13223,
                     14071,12145,10168,24037,13712]
    
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

    print(f"Processing {len(fire_polygons)} polygons.")
    
    for i,_ in enumerate(fire_polygons.index):
        perimeter = fire_polygons.iloc[[i]]
        
        # Set Time
        YEAR_START = perimeter.tst_year.item()
        YEAR_END = perimeter.ted_year.item()
        MONTH_START = 5
        MONTH_END = 10
        
        # Get attributes of this perimeter
        FIREID = perimeter.fireid.item()
        print(f"Processing fire: {FIREID} ({YEAR_START})...")
        
        if np.isin(FIREID,FALSE_FIRES_ID):
            print("False positive fire, skipping.")
            continue
        
        fire_perimeter_attrs = processing_lut.loc[processing_lut.fireid == FIREID]
        
        if len(fire_perimeter_attrs)==0:
            print(f"Fire not in processing lookup table.")
            continue
        
        UTM_TILE_NAME = fire_perimeter_attrs.opt_UTM_tile.item()

        for severity_index in ["dNBR","dGEMI"]:

            # Get list of severity rasters for this fire perimeter
            burn_severity_files = glob(
                os.path.join(HLS_FOLDER,"severity_rasters",f"{severity_index}_{UTM_TILE_NAME}_{YEAR_END}*.tif")
                )

            # skip perimeters where no rasters exist yet
            if not burn_severity_files:
                print(f"No {severity_index} rasters exist for fire {FIREID}.")
                continue
            
            severity_datestrings = [os.path.basename(f).split('.')[0][-10:] for f in burn_severity_files]
            severity_dates = [datetime.datetime.strptime(s,"%Y-%m-%d") for s in severity_datestrings]
            
            df_sev = pd.DataFrame({
                "fname_severity_raster":burn_severity_files,
                "datestring_raster":severity_datestrings,
                "date_raster":severity_dates
            })

            df_sev["doy"] = df_sev.date_raster.dt.day_of_year
            
            # and optimality rasters
            optimality_files = glob(
                os.path.join(HLS_FOLDER,"optimality_rasters",f"*{UTM_TILE_NAME}*.tif")
                )
            optimality_datestrings = [os.path.basename(f).split('.')[0][-10:] for f in optimality_files]
            
            df_opt = pd.DataFrame({
                "fname_optimality_raster":optimality_files,
                "datestring_raster":optimality_datestrings
            })
            
            df_temp = pd.merge(df_sev,df_opt, on="datestring_raster")
            df_temp["fireid"] = FIREID
            df_temp["severity_index"] = severity_index
            
            if not burn_severity_files:
                print(f"No {severity_index} rasters found.")
                continue
            
            # Calculate clear pixel percentage for all burn severity files of this scar
            df_temp['pct_clear_pixel'] = [calculate_clear_pixel_percentage(fname, polygon = perimeter) for fname in tqdm(df_temp.fname_severity_raster)]
            
            # Calculate optimality summary statistics for all files of this scar
            df_temp['mean_optimality'] = [calculate_optimality_statistics(fname, polygon = perimeter, stat = "mean") for fname in tqdm(df_temp.fname_optimality_raster)]
            
            
            # Filter out rows for burn severity rasters before the end of burn for each fire event (not >= bc. some fires were still burning on last day)
            df_temp = df_temp[df_temp['doy'] > fire_perimeter_attrs.ted_date.item().timetuple().tm_yday]
            
            # Find optimal combination of optimality, clear pixel percentage, and doy
            score = df_temp["mean_optimality"] * df_temp["pct_clear_pixel"] * (366 - df_temp["doy"])/366
            
            if all(np.isnan(score)):
                print("Nan optimality score...")
                continue
                
            best_idx = np.argmax(score)
            
            df_best = df_temp.iloc[best_idx,:].to_frame().T
            output_lut = pd.concat([df_best,output_lut])
            
            df_temp["score"] = score
            
            output_table = pd.concat([df_temp,output_table])

    # Write to ouput tables
    output_lut.to_csv(os.path.join(DATA_FOLDER,"tables","optimality_LUT.csv"),sep=";")      
    output_table.to_csv(os.path.join(DATA_FOLDER,"tables","optimality_overview_table.csv"),sep=";")        
