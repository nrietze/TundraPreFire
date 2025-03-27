# %%
import os
import sys
import platform
import multiprocessing
from glob import glob
import numpy as np
import pandas as pd
import geopandas as gpd
from osgeo import gdal
import xarray as xr
import rioxarray as rxr
import rasterio as rio
from rasterio.mask import mask
from shapely.geometry import mapping
import datetime
from dateutil.relativedelta import relativedelta
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
            geometry = gpd.GeoSeries([polygon.geometry],crs = 4326)
            polygon_geom = geometry.to_crs(src.crs).geometry.iloc[0]
        
            # Crop to the polygon area
            image, _ = mask(src, [polygon_geom], crop=True)

        # Load corresponding water mask
        water_mask_file = os.path.join(os.path.dirname(data),
                                        os.path.basename(
                                            replace_string_in_filename(data, "watermask")))
        
        with rio.open(water_mask_file) as watermask_src:
            eometry = gpd.GeoSeries([polygon.geometry],crs = 4326)
            polygon_geom = geometry.to_crs(watermask_src.crs).geometry.iloc[0]
            water_mask, _ = mask(watermask_src, [mapping(polygon_geom)], crop=True)
        
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

def time_index_from_filenames(files):
    '''
    Helper function to create a pandas DatetimeIndex
    '''
    return [datetime.datetime.strptime(f.split('.')[-4], '%Y%jT%H%M%S') for f in files]

def calculate_severity_metrics(gemi_prefire_composite,
                               gemi_postfire_composite,
                               nbr_prefire_composite, nbr_postfire_composite,
                               date_postfire,
                               index_names: list,
                               utm_tileid: str,
                               polygon:gpd.GeoDataFrame,
                               MIN_VALID_PERCENTAGE:int,
                               OUT_DIR=None):
    
    for index_name in index_names:
        print(f"Calculating {index_name}.")
        if index_name == "dNBR":
            # Calculate dNBR (pre - post fire)
            final_index = nbr_prefire_composite - nbr_postfire_composite
            
        elif index_name == "dGEMI":
            # Calculate dGEMI (pre - post fire)
            final_index = gemi_prefire_composite - gemi_postfire_composite
            
        elif index_name == "RdNBR":
            # Calculate RdNBR (pre - post fire)
            dnbr = nbr_prefire_composite - nbr_postfire_composite
            final_index = dnbr / np.sqrt(abs(nbr_prefire_composite))
            
        elif index_name == "RBR":
            # Calculate relativized burn ratio (pre - post fire)
            dnbr = nbr_prefire_composite - nbr_postfire_composite
            final_index = dnbr / (nbr_prefire_composite + 1.001)

        # Exclude infinite values
        final_index = xr.where(final_index != np.inf, final_index, np.nan, keep_attrs=True)
    
        # Change the long_name attribute
        final_index.attrs['long_name'] = index_name
        
        # Check if output raster has > 80% valid pixels in fire perimeter
        clear_percentage = calculate_clear_pixel_percentage(final_index,polygon=polygon)
        
        print(f"Clear pixel percentage in current perimeter: {clear_percentage}")
        
        if clear_percentage > MIN_VALID_PERCENTAGE and OUT_DIR is not None:
            out_name = os.path.join(OUT_DIR, f"{index_name}_{utm_tileid}_{date_postfire}.tif")
    
            # Export as Cloud Optimized GeoTIFF
            final_index.rio.to_raster(raster_path=out_name, driver="COG")
        else:
            print(f"Severity raster doesn't have > {MIN_VALID_PERCENTAGE} % valid pixels in perimter. Skipping export.\n")
            continue
    return None

def joblib_fct_calculate_severity(PROCESSED_HLS_DIR:str, 
                                  utm_tileid:str,
                                  processing_lut:pd.DataFrame,
                                  fire_polygons:gpd.GeoDataFrame,
                                  OUT_DIR = None):

    # List of severity indices that are calculated
    index_names = ["dNBR","dGEMI","RdNBR","RBR"]

    # extract pattern to search for in existing HLS files
    search_index_names = ["NBR","GEMI"]

    # Filter fire perimeters that use this UTM tile
    fires_in_utm = processing_lut.loc[processing_lut.opt_UTM_tile == utm_tileid]
    
    fire_polygons_in_utm = fire_polygons.loc[np.isin(fire_polygons.fireid, fires_in_utm.fireid)]
    
    # Loop through each year
    for (_,perimeter) in fire_polygons_in_utm.iterrows():
        fire_perimeter_attrs = processing_lut.loc[processing_lut.fireid == perimeter.fireid]
        
        # extract years of burning
        year = perimeter.tst_year

        # Only process HLS data for fires starting in 2017
        if year < 2017:
            continue
        
        print(f"processing data for fire perimeter {perimeter.fireid} in {year} (UTM tile {utm_tileid}). \n")
        
        # Create list of NBR and GEMI files for that tile and year
        nbr_files = glob(os.path.join(PROCESSED_HLS_DIR,f"*{utm_tileid}*NBR.tif"))
        gemi_files = glob(os.path.join(PROCESSED_HLS_DIR,f"*{utm_tileid}*GEMI.tif"))

        # Only use scenes with <= 20 % clear pixels in the perimeter to process further
        good_nbr_files = []
        good_gemi_files = []

        MIN_VALID_PERCENTAGE = 20
        
        for file in nbr_files:
            clear_percentage = calculate_clear_pixel_percentage(file,
                                                                polygon = perimeter)
            
            if clear_percentage > MIN_VALID_PERCENTAGE:
                good_nbr_files.append(file)

        for file in gemi_files:
            clear_percentage = calculate_clear_pixel_percentage(file,
                                                                polygon = perimeter)
            
            if clear_percentage > MIN_VALID_PERCENTAGE:
                good_gemi_files.append(file)

        if len(good_nbr_files) == 0:
            print(f"Found no good NBR files with >{MIN_VALID_PERCENTAGE}% valid pixels for fire perimeter: {perimeter.fireid} (UTM tile {utm_tileid}).")
            continue

        print(f"Found {len(good_nbr_files)} NBR files with >{MIN_VALID_PERCENTAGE}% valid pixels for fire perimeter: {perimeter.fireid}.")
        
        time = xr.Variable('time', time_index_from_filenames(good_nbr_files))
        
        # define chunk size for data loading
        chunk_size = dict(band=1, x=3600, y=3600)
    
        # create image stack of NBR/GEMI scenes for that tile
        nbr_ts = xr.concat([rxr.open_rasterio(f, mask_and_scale=True,
                                                chunks=chunk_size
                                                ).squeeze('band',drop=True) for f in good_nbr_files], dim=time)
        nbr_ts.name = "NBR"

        time = xr.Variable('time', time_index_from_filenames(good_gemi_files))
        gemi_ts = xr.concat([rxr.open_rasterio(f, mask_and_scale=True,
                                                chunks=chunk_size
                                                ).squeeze('band',drop=True) for f in good_gemi_files], dim=time)
        gemi_ts.name = "GEMI"

        gemi_ts = gemi_ts.sortby("time")
        nbr_ts = nbr_ts.sortby("time")
        
        # get date from earliest fire end in that tile, from then on start matching NBR pairs until Oct 31
        MAX_TIMEDELTA = 7
        earliest_tile_fire_end = min(fire_perimeter_attrs.ted_date)
        first_search_date_post = earliest_tile_fire_end.strftime("%Y-%m-%d")
        first_search_date_pre = (earliest_tile_fire_end - pd.Timedelta(days=MAX_TIMEDELTA) -
                                 relativedelta(years=1)).strftime("%Y-%m-%d")
        last_search_date_post = f"{year}-10-31"
        last_search_date_pre = f"{year-1}-10-31"
        
        # Resample to daily resolution and compute max composite
        nbr_ts.coords['date'] = nbr_ts.time.dt.floor('1D') # create date coord. for grouping
        nbr_daily = nbr_ts.groupby("date").max()
        gemi_ts.coords['date'] = gemi_ts.time.dt.floor('1D')
        gemi_daily = gemi_ts.groupby("date").max()
        
        # Select post-fire time series
        postfire_nbr = nbr_daily.sel(date=slice(first_search_date_post,
                                                last_search_date_post))

        # Iterate through post-fire dates
        for postfire_date in postfire_nbr["date"].values:
            postfire_date = pd.Timestamp(postfire_date)

            # Retrieve post-fire composite for both GEMI and NBR
            gemi_postfire_composite = gemi_daily.sel(date=postfire_date.strftime("%Y-%m-%d"))
            nbr_postfire_composite = nbr_daily.sel(date=postfire_date.strftime("%Y-%m-%d"))
            
            # Define pre-fire search window (±7 days of same day last year)
            prefire_start = postfire_date - relativedelta(years=1) - pd.Timedelta(days=MAX_TIMEDELTA)
            prefire_end = postfire_date - relativedelta(years=1) + pd.Timedelta(days=MAX_TIMEDELTA)
        
            # Select pre-fire composites
            prefire_gemi = gemi_daily.sel(date=slice(prefire_start, prefire_end))
            prefire_nbr = nbr_daily.sel(date=slice(prefire_start, prefire_end))
        
            if prefire_gemi.date.size > 0 and prefire_nbr.date.size > 0:  # Ensure pre-fire data exists
                # Median composite over available days
                gemi_prefire_composite = prefire_gemi.mean(dim="date")
                nbr_prefire_composite = prefire_nbr.mean(dim="date")

                calculate_severity_metrics(
                    gemi_prefire_composite, gemi_postfire_composite,
                    nbr_prefire_composite, nbr_postfire_composite,
                    date_postfire=postfire_date.strftime("%Y-%m-%d"),
                    index_names=index_names, utm_tileid=utm_tileid,
                    polygon=perimeter,
                    MIN_VALID_PERCENTAGE=1,
                    OUT_DIR=OUT_DIR
                )
                
            else:
                print(f"No pre-fire match found for {postfire_date}")

        return None

# %% Execute pararellized processing
if __name__ == "__main__":
    print("Loading data...")
    
    if platform.system() == "Windows":
        DATA_FOLDER = 'data/' # on local machine
        PROCESSED_HLS_DIR = "data/raster/hls/processed"
        OUT_FOLDER = './data/raster/hls/severity_rasters'
    else:
        DATA_FOLDER = '~/data/' # on sciencecluster
        PROCESSED_HLS_DIR = "/home/nrietz/scratch/raster/hls/processed/" # Set original data paths
        OUT_FOLDER = '/data/nrietz/raster/hls/severity_rasters'
    
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
    fire_polygons = fire_polygons.loc[fire_polygons.tst_year>=2017]
    
    # Load UTM tiles from argument ingested through SLURM bash
    UTM_TILE_FILE = sys.argv[1]
    
    with open(UTM_TILE_FILE, "r") as file:
        UTM_TILE_LIST = [line.strip() for line in file]

    # 3. Set up parallelization to calculate burn severity rasters
    num_workers = (
                    int(os.environ.get("SLURM_CPUS_PER_TASK"))
                    if os.environ.get("SLURM_CPUS_PER_TASK")
                    else os.cpu_count()
                )
    
    print("Starting parallelized calculation...")
    
    multiprocessing.set_start_method('spawn')
    
    Parallel(n_jobs=num_workers, backend='loky')(
        delayed(joblib_fct_calculate_severity)(PROCESSED_HLS_DIR,
                                               UTM_TILE_NAME,
                                               processing_lut,
                                               fire_polygons,
                                               OUT_DIR = OUT_FOLDER) for UTM_TILE_NAME in UTM_TILE_LIST)
