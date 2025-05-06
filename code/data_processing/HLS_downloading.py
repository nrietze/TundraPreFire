import pandas as pd
import subprocess
import sys
import platform
import os
import numpy as np

# Set paths
if platform.system() == "Windows":
    DATA_FOLDER = 'data/' # on local machine
    OUTPUT_DIR = "data/raster/hls/test"
else:
    DATA_FOLDER = '~/data/' # on sciencecluster
    OUTPUT_DIR = "/home/nrietz/scratch/raster/hls/" # Set original data paths

# Load Processing look-up-table to match UTM tiles to fire perimeter IDs
lut_df = pd.read_csv(
    os.path.join(DATA_FOLDER,"tables/processing_LUT.csv"),
    index_col=0)

tileid_file = "code/data_processing/tileid_file.txt"
download_script = "code/data_processing/getHLS.sh"

# Fixed date range for HLS bulk download
start_date = "05-01"
end_date = "10-31"

# Load UTM tiles
UTM_TILE_FILE = sys.argv[1]

with open(UTM_TILE_FILE, "r") as file:
    UTM_TILE_LIST = [line.strip() for line in file] 

# Iterate through tiles and find corresponding year
for TILE_ID in UTM_TILE_LIST:
    matching_row = lut_df[lut_df["opt_UTM_tile"] == TILE_ID]
    if not matching_row.empty:
        years = np.unique(matching_row["tst_year"].values)

        for year in years:
            # Format dynamic date range
            start_date_dynamic = f"{year}-{start_date}"
            end_date_dynamic = f"{year}-{end_date}"
    
            print(f"Processing tile: {TILE_ID} for {start_date_dynamic} to {end_date_dynamic}")

            # Create a temporary file with the single tile ID
            TASK_ID = os.environ.get("SLURM_ARRAY_TASK_ID")
            
            temp_tile_file = f"tmp/temp_tile{TASK_ID}.txt"
            with open(temp_tile_file, "w") as temp_file:
                temp_file.write(TILE_ID + "\n")
            
            # Execute bash command
            subprocess.run(["bash", download_script, temp_tile_file, start_date_dynamic, end_date_dynamic, OUTPUT_DIR])
            
            # Remove the temporary file
            os.remove(temp_tile_file)
    else:
        print(f"No matching year found for tile: {TILE_ID}")
