# %%
"""
This sciprt downloads Landsat-8 & 9 Collection-2 Level-2 imagery from Google Earth Engine.


Author: nils.rietze@uzh.ch
Date: 27.11.2024
"""
import ee
import geemap
import datetime
import os
from tqdm import tqdm
import platform
from glob import glob
from zipfile import ZipFile
import geopandas as gpd
import pandas as pd

ee.Initialize(project = "ee-nrietze")

#%% 0. Configure parameters
EXPORT_TO_DRIVE = False # False --> download locally
MAX_CLOUD_COVER = 80
PIXEL_RESOLUTION = 30
RASTER_BAND = 'TIR'

if platform.system() == "Windows":
    DATA_FOLDER = 'data/' # on local machine
    OUT_FOLDER = "data/raster/landsat8"
else:
    DATA_FOLDER = '~/data/' # on sciencecluster
    OUT_FOLDER = '/data/nrietz/raster/landsat8'

# Load VIIRS perimeters in Siberian tundra
FN_VIIRS_CAVM_PERIMETERS = os.path.join(DATA_FOLDER,
                                        "feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")
fire_polygons = gpd.read_file(FN_VIIRS_CAVM_PERIMETERS)
fire_polygons = fire_polygons.loc[fire_polygons.tst_year>=2017]

# Load Processing look-up-table to match UTM tiles to fire perimeter IDs
processing_lut = pd.read_csv(
    os.path.join(DATA_FOLDER,"tables/processing_LUT.csv"),
    index_col=0)

# Create date objects
processing_lut['tst_date'] = processing_lut.apply(
    lambda row: datetime.date(row['tst_year'], row['tst_month'], row['tst_day']), axis=1)
processing_lut['ted_date'] = processing_lut.apply(
    lambda row: datetime.date(row['ted_year'], row['ted_month'], row['ted_day']), axis=1)

# %% Define user functions
# Algorithm constants
NDVI_V = 0.6 # vegetation
NDVI_S = 0.2 # soil

# emissivity values
EPSILON_V = 0.985
EPSILON_S = 0.97
EPSILON_W = 0.99

# Constants
cs_l8 = [0.04019, 0.02916, 1.01523,
         -0.38333, -1.50294, 0.20324,
         0.00918, 1.36072, -0.27514]
cs_l7 = [0.06518, 0.00683, 1.02717,
         -0.53003, -1.25866, 0.10490,
         -0.01965, 1.36947, -0.24310]
cs_l5 = [0.07518, -0.00492, 1.03189,
         -0.59600, -1.22554, 0.08104,
         -0.02767, 1.43740, -0.25844]

# Rename and scale bands
def fun_bands(img):
    """
    Function to rename spectral bands and scale them acc. to User guide. 

    Args:
        img (ee.Image): ee Image object

    Returns:
        img: ee.Image with renamed and scaled C2-L2 bands (for Landsat-8 & 9)
    """
    bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
    thermal_band = ['ST_B10']
    new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
    new_thermal_bands = ['TIR']
    vnirswir = img.select(bands)\
        .multiply(0.0000275).add(-0.2)\
            .rename(new_bands)
    tir = img.select(thermal_band)\
        .multiply(0.00341802).add(149.0).subtract(273.15)\
            .rename(new_thermal_bands)
    return vnirswir.addBands(tir).copyProperties(img, ['system:time_start'])

# Mask OLI-TIRS images
def fun_mask_ls_sr(img):
    """
    7-6 = aerosols
    5 = water
    4 = snow/ice
    3 = cloud shadow
    2 = cloud/shadow adjacent
    1 = cloud
    0 = cirrus
    """
    # bits used in HLS filtering: [0,1,2,3,4]
    
    # Create a bitmask to filter out bits 0 to 4
    bitmask = (ee.Number(2).pow(0)        # Cirrus (bit 0)
               .add(ee.Number(2).pow(1))  # Cloud (bit 1)
               .add(ee.Number(2).pow(2))  # Cloud/Shadow Adjacent (bit 2)
               .add(ee.Number(2).pow(3))  # Cloud Shadow (bit 3)
               .add(ee.Number(2).pow(4))) # Snow/Ice (bit 4)
    
    qa = img.select('QA_PIXEL')
    
    mask = qa.bitwiseAnd(bitmask).eq(0)
    
    return img.updateMask(mask)

# calculate NDVI
def fun_ndvi(img):
    ndvi = img.normalizedDifference(['NIR', 'R']).rename('NDVI')
    return img.addBands(ndvi)

# Fraction Vegetation Cover (FVC)
def fun_fvc(img):
    fvc = img.expression(
        '((NDVI-NDVI_S)/(NDVI_V-NDVI_S))**2',
        {
            'NDVI': img.select('NDVI'),
            'NDVI_S': NDVI_S,
            'NDVI_V': NDVI_V
        }
    ).rename('FVC')
    return img.addBands(fvc)

def fun_date(img):
    return ee.Date(ee.Image(img).date().format("YYYY-MM-dd"))


def fun_getdates(imgCol):
    return ee.List(imgCol.toList(imgCol.size()).map(fun_date))


def fun_mosaic(date, newList):
    # cast list & date
    newList = ee.List(newList)
    date = ee.Date(date)

    # filter img-collection
    filtered = ee.ImageCollection(subCol.filterDate(date, date.advance(1, 'day')))

    # check duplicate
    img_previous = ee.Image(newList.get(-1))
    img_previous_datestring = img_previous.date().format("YYYY-MM-dd")
    img_previous_millis = ee.Number(ee.Date(img_previous_datestring).millis())

    img_new_datestring = filtered.select(RASTER_BAND).first().date().format("YYYY-MM-dd")
    img_new_date = ee.Date(img_new_datestring).millis()
    img_new_millis = ee.Number(ee.Date(img_new_datestring).millis())

    date_diff = img_previous_millis.subtract(img_new_millis)

    # mosaic
    img_mosaic = ee.Algorithms.If(
        date_diff.neq(0),
        filtered.select(RASTER_BAND).mosaic().set('system:time_start', img_new_date),
        ee.ImageCollection(subCol.filterDate(pseudodate, pseudodate.advance(1, 'day')))
    )

    tester = ee.Algorithms.If(date_diff.neq(0), ee.Number(1), ee.Number(0))

    return ee.Algorithms.If(tester, newList.add(img_mosaic), newList)


def fun_timeband(img):
    time = ee.Image(img.metadata('system:time_start', 'TIME').divide(86400000))
    timeband = time.updateMask(img.select(RASTER_BAND).mask())
    return img.addBands(timeband)

# %% Run operations on GEE cloud

# Buffer all fire perimeters by 2km to avoid clipping data from perimeter edge
for i in fire_polygons.index:
    perimeter = fire_polygons.loc[[i]]
    
    # Set Time
    YEAR_START = perimeter.tst_year.item()
    YEAR_END = perimeter.ted_year.item()
    MONTH_START = 5
    MONTH_END = 10
    
    # Only process HLS data for fires starting in 2017
    if YEAR_START < 2017:
        continue
    
    # Get attributes of this perimeter
    FIREID = perimeter.fireid.item()
    fire_perimeter_attrs = processing_lut.loc[processing_lut.fireid == FIREID]
    TILE_NAME = fire_perimeter_attrs.opt_UTM_tile.item()
    
    # Set output EPSG for this perimeter by extracting UTM EPSG code for reprojection
    utm_epsg_code = 32600 + int(TILE_NAME[:2])
    EPSG = f'EPSG:{utm_epsg_code}'
    
    FILENAME_PREFIX = f"{FIREID}_{TILE_NAME}_"
    
    buffered_perimeter = perimeter.to_crs(utm_epsg_code)
    buffered_perimeter.geometry = buffered_perimeter.geometry.buffer(2000)
    
    bbox = list(buffered_perimeter.to_crs("EPSG:4326").total_bounds)
    # bbox = geemap.gdf_bounds(buffered_perimeter.to_crs("EPSG:4326"))
    select_roi = ee.Geometry.Rectangle(bbox)
    
    print(f"Processing perimeter nr. {FIREID} ({YEAR_START},{TILE_NAME})")
    
    # Fetch and mask Landsat 8 OLI-TIRS C2-L2 imagery
    imgCol_L8_SR = (
        ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
        .filterBounds(select_roi)
        .filter(ee.Filter.calendarRange(YEAR_START,YEAR_END,'year'))
        .filter(ee.Filter.calendarRange(MONTH_START,MONTH_END,'month'))
        .filter(ee.Filter.lt('CLOUD_COVER_LAND', MAX_CLOUD_COVER))
        .map(fun_mask_ls_sr))
    print("Searched & filtered data on GEE...")
    
    # Rename bands
    imgCol_L8_SR = imgCol_L8_SR.map(fun_bands)

    # Calculate parameters and indices
    # imgCol_L8_SR = imgCol_L8_SR.map(fun_ndvi)
    # imgCol_L8_SR = imgCol_L8_SR.map(fun_fvc)

    # Format time in Image Collection
    pseudodate = ee.Date('1960-01-01')
    subCol = ee.ImageCollection(imgCol_L8_SR.select(RASTER_BAND))
    dates = fun_getdates(subCol)
    ini_date = ee.Date(dates.get(0))
    ini_merge = subCol.filterDate(ini_date, ini_date.advance(1, 'day'))
    ini_merge = ini_merge.select(RASTER_BAND).mosaic().set('system:time_start',
                                                           ini_date.millis())
    ini = ee.List([ini_merge])
    
    # Mosaic imgCollection
    imgCol_mosaic = ee.ImageCollection(ee.List(dates.iterate(fun_mosaic, ini)))
    imgCol_mosaic = imgCol_mosaic.map(fun_timeband)

    # Download entire image collection
    print("Exporting LSt data...")
    
    # construct list of filenames with image date
    datelist = dates.getInfo()
    FILENAMES = [f"{FILENAME_PREFIX}{datetime.datetime.fromtimestamp(d['value'] / 1000).strftime('%Y-%m-%d')}" for d in datelist]
    
    if EXPORT_TO_DRIVE:
        geemap.ee_export_image_collection_to_drive(imgCol_mosaic,
                                                folder = "L8_C2L2_LST",
                                                fileNamePrefix = FILENAME_PREFIX,
                                                scale = PIXEL_RESOLUTION,
                                                maxPixels = 1e13,
                                                region = select_roi, #['coordinates'][0],
                                                crs = EPSG)
    else:
        geemap.ee_export_image_collection(imgCol_mosaic,
                                          out_dir=OUT_FOLDER,
                                          filenames = FILENAMES,
                                          scale = PIXEL_RESOLUTION,
                                          region = select_roi,
                                          crs = EPSG)
        
    print("done.")

#%% [OPTIONAL] Generate temporal statistics & Export
lookup_metrics = {
    'mean': ee.Reducer.mean(),
    'min': ee.Reducer.min(),
    'max': ee.Reducer.max(),
    'std': ee.Reducer.stdDev(),
    'median': ee.Reducer.median(),
    'ts': ee.Reducer.sensSlope()
}

select_metrics = ['nobs']
percentiles = []
ROI_FILENAME = 'KYT_TEST'

for metric in select_metrics:
    if metric == 'ts':
        temp = imgCol_mosaic.select(['TIME', RASTER_BAND]).reduce(ee.Reducer.sensSlope())
        temp = temp.select('slope')
        temp = temp.multiply(365)
        temp = temp.multiply(100000000).int32()
    elif metric == 'nobs':
        temp = imgCol_mosaic.select(RASTER_BAND).count()
        temp = temp.int16()
    else:
        if metric == 'percentile':
            temp = imgCol_mosaic.select(RASTER_BAND).reduce(ee.Reducer.percentile(percentiles))
        else:
            reducer = lookup_metrics[metric]
            temp = imgCol_mosaic.select(RASTER_BAND).reduce(reducer)
        # if RASTER_BAND == 'TIR':
        #     temp = temp.multiply(100).int16()
        # else:
        #     temp = temp.multiply(10000).int16()

    # Export to Drive
    filename = RASTER_BAND + '_'+ \
        ROI_FILENAME + '_GEE_' + \
            str(YEAR_START) + '-' + \
                str(YEAR_END)+ '_' + \
                    str(MONTH_START) + '-' + \
                        str(MONTH_END) + '_' + metric
    
    task = ee.batch.Export.image.toDrive(
        image = temp, 
        description = filename,
        folder = "Tundra_PreFire",
        scale = PIXEL_RESOLUTION,
        maxPixels = 1e13,
        region = select_roi['coordinates'][0],
        crs = EPSG,
        formatOptions={
            'cloudOptimized': True
            }
        )
    process = ee.batch.Task.start(task)