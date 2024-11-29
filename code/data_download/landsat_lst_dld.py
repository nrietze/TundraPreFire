# %%
"""
This sciprt downloads Landsat-8 & 9 Colelction-2 Level-2 imagery from Google Earth Engine.


Author: nils.rietze@uzh.ch
Date: 27.11.2024
"""
import ee
import geopandas as gp
import geemap
import os

ee.Initialize()

#%% 0. Configure parameters

select_parameters = ['TIR']
select_metrics = ['nobs']
percentiles = []

# Download images or imageCollection?
DOWNLOAD_IMG = False
DOWNLOAD_IMGCOL = True

# Time
YEAR_START = 2020
YEAR_END = 2020
MONTH_START = 5
MONTH_END = 8

# Space
aois = gp.read_file("./data/feature_layers/roi.geojson")

bbox = list(aois.total_bounds)
select_roi = ee.Geometry.Rectangle(bbox)
MAX_CLOUD_COVER = 60
EPSG = 'EPSG:32655'
PIXEL_RESOLUTION = 30
ROI_FILENAME = 'KYT_TEST'

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

# %% Define user functions

lookup_metrics = {
    'mean': ee.Reducer.mean(),
    'min': ee.Reducer.min(),
    'max': ee.Reducer.max(),
    'std': ee.Reducer.stdDev(),
    'median': ee.Reducer.median(),
    'ts': ee.Reducer.sensSlope()
}

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
    cloudShadowBitMask = ee.Number(2).pow(3).int()
    cloudsBitMask = ee.Number(2).pow(5).int()
    snowBitMask = ee.Number(2).pow(4).int()
    qa = img.select('QA_PIXEL')
    mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(
        qa.bitwiseAnd(cloudsBitMask).eq(0)).And(
            qa.bitwiseAnd(snowBitMask).eq(0)
            )
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

    img_new_datestring = filtered.select(parameter).first().date().format("YYYY-MM-dd")
    img_new_date = ee.Date(img_new_datestring).millis()
    img_new_millis = ee.Number(ee.Date(img_new_datestring).millis())

    date_diff = img_previous_millis.subtract(img_new_millis)

    # mosaic
    img_mosaic = ee.Algorithms.If(
        date_diff.neq(0),
        filtered.select(parameter).mosaic().set('system:time_start', img_new_date),
        ee.ImageCollection(subCol.filterDate(pseudodate, pseudodate.advance(1, 'day')))
    )

    tester = ee.Algorithms.If(date_diff.neq(0), ee.Number(1), ee.Number(0))

    return ee.Algorithms.If(tester, newList.add(img_mosaic), newList)


def fun_timeband(img):
    time = ee.Image(img.metadata('system:time_start', 'TIME').divide(86400000))
    timeband = time.updateMask(img.select(parameter).mask())
    return img.addBands(timeband)

# %% Run operations on GEE cloud

# Fetch and mask Landsat 8 OLI-TIRS C2-L2 imagery
imgCol_L8_SR = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')\
    .filterBounds(select_roi)\
    .filter(ee.Filter.calendarRange(YEAR_START,YEAR_END,'year'))\
    .filter(ee.Filter.calendarRange(MONTH_START,MONTH_END,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER_LAND', MAX_CLOUD_COVER))\
    .map(fun_mask_ls_sr)

# Rename bands
imgCol_L8_SR = imgCol_L8_SR.map(fun_bands)

# Calculate parameters and indices
imgCol_L8_SR = imgCol_L8_SR.map(fun_ndvi)
imgCol_L8_SR = imgCol_L8_SR.map(fun_fvc)


#%% Generate temporal statistics & Export
for parameter in select_parameters:

    # Mosaic imgCollection
    pseudodate = ee.Date('1960-01-01')
    subCol = ee.ImageCollection(imgCol_L8_SR.select(parameter))
    dates = fun_getdates(subCol)
    ini_date = ee.Date(dates.get(0))
    ini_merge = subCol.filterDate(ini_date, ini_date.advance(1, 'day'))
    ini_merge = ini_merge.select(parameter).mosaic().set('system:time_start', ini_date.millis())
    ini = ee.List([ini_merge])
    imgCol_mosaic = ee.ImageCollection(ee.List(dates.iterate(fun_mosaic, ini)))
    imgCol_mosaic = imgCol_mosaic.map(fun_timeband)

    for metric in select_metrics:
        if metric == 'ts':
            temp = imgCol_mosaic.select(['TIME', parameter]).reduce(ee.Reducer.sensSlope())
            temp = temp.select('slope')
            temp = temp.multiply(365)
            temp = temp.multiply(100000000).int32()
        elif metric == 'nobs':
            temp = imgCol_mosaic.select(parameter).count()
            temp = temp.int16()
        else:
            if metric == 'percentile':
                temp = imgCol_mosaic.select(parameter).reduce(ee.Reducer.percentile(percentiles))
            else:
                reducer = lookup_metrics[metric]
                temp = imgCol_mosaic.select(parameter).reduce(reducer)
            # if parameter == 'TIR':
            #     temp = temp.multiply(100).int16()
            # else:
            #     temp = temp.multiply(10000).int16()

        if DOWNLOAD_IMG:
            # Export to Drive
            filename = parameter + '_'+ \
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
        
    # Download entire image collection
    if DOWNLOAD_IMGCOL:
        geemap.ee_export_image_collection_to_drive(imgCol_mosaic,
                                                   folder = "Tundra_PreFire/LST/",
                                                   scale = PIXEL_RESOLUTION,
                                                   region = select_roi['coordinates'][0],
                                                   crs = EPSG)
