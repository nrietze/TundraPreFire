#!/bin/bash

### Here are the SBATCH parameters that you should always consider:
# SBATCH --time=0-00:20:00   ## days-hours:minutes:seconds
# SBATCH --mem 16000M         ## 16GB ram (hardware ratio is < 4GB/core)
# SBATCH --ntasks=1          ## Not strictly necessary because default is 1
# SBATCH --cpus-per-task=32   ## ask for 32 cpus (Use greater than 1 for parallelized jobs)
# SBATCH --output=--output=home/$USER/scratch/logs/job_%j.out

module load mamba
source activate lpdaac_vitals

echo 'HLS processing starting.'

## python code/data_download/hls_downloader.py
## python code/data_download/HLS_SuPER.py -roi '144.4001779629389, 71.24588582272926,146.16765743760948, 71.6168891179392' -dir /home/nrietz/data/raster/hls/ -start 2020-05-01 -end 2020-10-31 -bands RED,GREEN,BLUE,NIR1,SWIR1,SWIR2,FMASK -qf True
python code/data_processing/HLS_preprocessing.py

echo 'finished'