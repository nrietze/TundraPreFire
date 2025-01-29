# -*- coding: utf-8 -*-
"""
===============================================================================
HLS Processing and Exporting Reformatted Data (HLS_PER)

This module contains functions to conduct subsetting and quality filtering of 
search results.
-------------------------------------------------------------------------------
Authors: Cole Krehbiel, Mahsa Jami, and Erik Bolch
Last Updated: 2024-09-18
===============================================================================
"""

import os
import sys
import logging

from osgeo import gdal
import numpy as np
from datetime import datetime as dt
import xarray as xr
import rioxarray as rxr
from xrspatial import multispectral
import dask.distributed
import dask.array as da
from skimage.filters import threshold_multiotsu, threshold_otsu


def create_output_name(url, band_dict):
    """
    Uses HLS default naming scheme to generate an output name with common band names.
    This allows for easier stacking of bands from both collections.
    """
    # Get Necessary Strings
    prod = url.split("/")[4].split(".")[0]
    asset = url.split("/")[-1].split(".")[-2]
    # Hard-coded one off for Fmask name incase it is not in the band_dict but is needed for masking
    if asset == "Fmask":
        output_name = f"{'.'.join(url.split('/')[-1].split('.')[:-2])}_Fmask.tif"
    else:
        for key, value in band_dict[prod].items():
            if value == asset:
                output_name = (
                    f"{'.'.join(url.split('/')[-1].split('.')[:-2])}_{key}.tif"
                )
    return output_name


def open_hls(url, roi=None, scale=True, chunk_size=dict(band=1, x=512, y=512)):
    """
    Generic Function to open an HLS COG and clip to ROI. For consistent scaling, this must be done manually.
    Some HLS Landsat scenes have the metadata in the wrong location.
    """
    # Open using rioxarray
    da_hls = rxr.open_rasterio(url, chunks=chunk_size, mask_and_scale=False).squeeze(
        "band", drop=True
    )

    # Reproject ROI and Clip if ROI is provided
    if roi is not None:
        roi = roi.to_crs(da_hls.spatial_ref.crs_wkt)
        da_hls = da_hls.rio.clip(roi.geometry.values, roi.crs, all_touched=True)

    # Apply Scale Factor if desired for non-quality layer
    if scale and "Fmask" not in url:
        # Mask Fill Values
        da_hls = xr.where(da_hls == -9999, np.nan, da_hls)
        # Scale Data
        da_hls = da_hls * 0.0001
        # Remove Scale Factor After Scaling - Prevents Double Scaling
        da_hls.attrs["scale_factor"] = 1.0

    # Add Scale Factor to Attributes Manually - This will overwrite/add if the data is missing.
    if not scale and "Fmask" not in url:
        da_hls.attrs["scale_factor"] = 0.0001

    return da_hls


def create_quality_mask(quality_data, bit_nums: list = [0, 1, 2, 3, 4, 5]):
    """
    Uses the Fmask layer and bit numbers to create a binary mask of good pixels.
    By default, bits 0-5 are used.
    """
    mask_array = np.zeros((quality_data.shape[0], quality_data.shape[1]))
    # Remove/Mask Fill Values and Convert to Integer
    quality_data = np.nan_to_num(quality_data, 0).astype(np.int8)
    for bit in bit_nums:
        # Create a Single Binary Mask Layer
        mask_temp = np.array(quality_data) & 1 << bit > 0
        mask_array = np.logical_or(mask_array, mask_temp)
    return mask_array


def process_granule(
    granule_urls,
    roi,
    quality_filter,
    scale,
    output_dir,
    band_dict,
    bit_nums=[0, 1, 2, 3, 4],
    chunk_size=dict(band=1, x=512, y=512),
):
    """
    Processes a list of HLS asset urls for a single granule.
    """

    # Setup Logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s:%(asctime)s ||| %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    # Check if all Outputs Exist for a Granule
    if not all(
        os.path.isfile(f"{output_dir}/{create_output_name(url, band_dict)}")
        for url in granule_urls
    ):

        # First Handle Quality Layer
        if quality_filter:
            # Generate Quality Layer URL
            split_asset = granule_urls[0].split("/")[-1].split(".")
            split_asset[-2] = "Fmask"
            quality_url = (
                f"{'/'.join(granule_urls[0].split('/')[:-1])}/{'.'.join(split_asset)}"
            )

            # Check if File exists in Output Directory
            output_name = create_output_name(quality_url, band_dict)
            output_file = f"{output_dir}/{output_name}"

            # Open Quality Layer
            qa_da = open_hls(quality_url, roi, scale, chunk_size)

            # Check if quality asset is already processed
            if not os.path.isfile(output_file):
                # Write Output
                qa_da.rio.to_raster(raster_path=output_file, driver="COG")
            else:
                logging.info(
                    f"Existing file {output_name} found in {output_dir}. Skipping."
                )

            # Remove Quality Layer from Granule Asset List if Present
            granule_urls = [asset for asset in granule_urls if asset != quality_url]

            # Create Quality Mask
            qa_mask = create_quality_mask(qa_da, bit_nums=bit_nums)

            
            # Compute Water mask 
            tile_name = granule_urls[0].split('/')[-1]
            filename = f"{tile_name.split('v2.0')[0]}v2.0_watermask.tif"
            wm_path = f"{output_dir}{filename}"
            
            # Check if that tile's water mask exists
            if os.path.exists(wm_path):
                logging.info("Loading water mask.")
                water_mask = rxr.open_rasterio(wm_path,
                                            chunks=chunk_size,
                                            masked=True
                                            ).squeeze('band', drop=True
                                                        ).astype(bool)
            else:
                logging.info("Computing NDWI water mask.")
                if os.path.basename(quality_url).split('.')[1] == "S30":
                    nir_url = quality_url.replace("Fmask","B8A")
                else: 
                    nir_url = quality_url.replace("Fmask","B05")
                green_url = quality_url.replace("Fmask","B03")
                
                # load spectral bands needed for NDWI
                nir = open_hls(nir_url, roi, scale, chunk_size)
                green = open_hls(green_url, roi, scale, chunk_size)
                
                ndwi = nir.copy()
                ndwi.data = (green - nir) / (green + nir)
                
                # Exclude data outside valid value range
                ndwi = ndwi.where((ndwi < -1) & (ndwi > 1), np.nan)
                
                ndwi.attrs['long_name'] = "NDWI"
                
                # Export NDWI as COG tiff
                ndwi_filename = filename.replace("watermask","NDWI")
                ndwi_path = f"{output_dir}{ndwi_filename}"
                ndwi.rio.to_raster(raster_path = ndwi_path, driver = 'COG')
            
                # Apply Otsu's thresholding to NDWI 
                hist_ndwi = da.histogram(ndwi, bins=2, range=[-1, 1])
                # otsu_thresh = threshold_otsu(hist_ndwi[1])[1]
                otsu_thresh = threshold_multiotsu(hist_ndwi[1])[1]
                
                # Apply threshold to mask out water
                water_mask = ndwi > otsu_thresh
                
                # Export water mask as COG tiff
                water_mask.astype(int).rio.to_raster(
                    raster_path = wm_path, driver = 'COG')
            
        # Process Remaining Assets
        for url in granule_urls:
            # Check if File exists in Output Directory
            output_name = create_output_name(url, band_dict)
            output_file = f"{output_dir}/{output_name}"

            # Check if scene is already processed
            if not os.path.isfile(output_file):
                # Open Asset
                da_hls = open_hls(url, roi, scale, chunk_size)

                # Apply Quality Mask if Desired
                if quality_filter:
                    merged_mask = np.logical_and(qa_da,water_mask)
                    da_hls = da_hls.where(~merged_mask)

                # Write Output
                da_hls.rio.to_raster(raster_path=output_file, driver="COG")
            else:
                logging.info(
                    f"Existing file {output_name} found in {output_dir}. Skipping."
                )
    else:
        logging.info(
            f"All assets related to {granule_urls[0].split('/')[-1]} are already processed, skipping."
        )


def build_hls_xarray_timeseries(
    hls_cog_list, mask_and_scale=True, chunk_size=dict(band=1, x=512, y=512)
):
    """
    Builds a single band timeseries using xarray for a list of HLS COGs. Dependent on file naming convention.
    Works on SuPERScript named files. Files need common naming bands corresponding HLSS and HLSL bands,
    e.g. HLSL30 Band 5 (NIR1) and HLSS30 Band 8A (NIR1)
    """
    # Define Band(s)
    bands = [filename.split(".")[6] for filename in hls_cog_list]

    # Make sure all files in list are the same band
    if not all(band == bands[0] for band in bands):
        raise ValueError("All listed files must be of the same band.")

    band_name = bands[0]

    # Create Time Variable
    try:
        time_list = [
            dt.strptime(filename.split(".")[3], "%Y%jT%H%M%S")
            for filename in hls_cog_list
        ]
    except ValueError:
        print("A COG does not have a valid date string in the filename.")

    time = xr.Variable("time", time_list)

    timeseries_da = xr.concat(
        [
            rxr.open_rasterio(
                filename, mask_and_scale=mask_and_scale, chunks=chunk_size
            ).squeeze("band", drop=True)
            for filename in hls_cog_list
        ],
        dim=time,
    )
    timeseries_da.name = band_name

    return timeseries_da


def create_timeseries_dataset(hls_file_dir, output_type, output_dir=None):
    """
    Creates an xarray dataset timeseries from a directory of HLS COGs.
    Writes to a netcdf output. Currently only works for HLS SuPER outputs.
    """

    # Setup Logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s:%(asctime)s ||| %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    # List Files in Directory
    all_files = [file for file in os.listdir(hls_file_dir) if file.endswith(".tif")]

    # Create Dictionary of Files by Band
    file_dict = {}
    for file in all_files:
        tile = file.split(".")[2]
        band = file.split(".")[6]
        full_path = os.path.join(hls_file_dir, file)
        if tile not in file_dict:
            file_dict[tile] = {}
        if band not in file_dict[tile]:
            file_dict[tile][band] = []
        file_dict[tile][band].append(full_path)

    # logging.info(f"{file_dict}")

    # Check that all bands within each tile have the same number of observations
    for tile, bands in file_dict.items():
        q_obs = {band: len(files) for band, files in bands.items()}
        if not all(q == list(q_obs.values())[0] for q in q_obs.values()):
            logging.info(
                f"Not all bands in {tile} have the same number of observations."
            )
            logging.info(f"{q_obs}")

    # Loop through each tile and build timeseries output

    for tile, bands in file_dict.items():
        dataset = xr.Dataset()

        timeseries_dict = {
            band: dask.delayed(build_hls_xarray_timeseries)(files)
            for band, files in bands.items()
        }
        timeseries_dict = dask.compute(timeseries_dict)[0]
        dataset = xr.Dataset(timeseries_dict)

        # Set up CF-Compliant Coordinate Attributes
        dataset.attrs["Conventions"] = "CF-1.6"
        dataset.attrs["title"] = "HLS SuPER Timeseries Dataset"
        dataset.attrs["institution"] = "LP DAAC"

        dataset.x.attrs["axis"] = "X"
        dataset.x.attrs["standard_name"] = "projection_x_coordinate"
        dataset.x.attrs["long_name"] = "x-coordinate in projected coordinate system"
        dataset.x.attrs["units"] = "m"

        dataset.y.attrs["axis"] = "Y"
        dataset.y.attrs["standard_name"] = "projection_y_coordinate"
        dataset.y.attrs["long_name"] = "y-coordinate in projected coordinate system"
        dataset.y.attrs["units"] = "m"

        dataset.time.attrs["axis"] = "Z"
        dataset.time.attrs["standard_name"] = "time"
        dataset.time.attrs["long_name"] = "time"

        # Get first and last date
        first_date = (
            dataset.time.data[0].astype("M8[ms]").astype(dt).strftime("%Y-%m-%d")
        )
        final_date = (
            dataset.time.data[-1].astype("M8[ms]").astype(dt).strftime("%Y-%m-%d")
        )

        # Write Outputs
        # if output_type == "NC4":
        output_path = os.path.join(
            output_dir, f"HLS.{tile}.{first_date}.{final_date}.subset.nc"
        )
        dataset.to_netcdf(output_path)
        # elif output_type == "ZARR":
        #     output_path = os.path.join(output_dir, "hls_timeseries_dataset.zarr")
        #     dataset.to_zarr(output_path)
        logging.info(f"Output saved to {output_path}")
