#!/bin/bash

### Here are the SBATCH parameters that you should always consider:
#SBATCH --time=0-00:05:00   ## days-hours:minutes:seconds
#SBATCH --mem 3000M         ## 3GB ram (hardware ratio is < 4GB/core)
#SBATCH --ntasks=1          ## Not strictly necessary because default is 1
#SBATCH --cpus-per-task=1   ## Use greater than 1 for parallelized jobs

module load mamba
source activate lpdaac_vitals

echo 'HLS download starting.'
date        ## Prints the system date

#python code/data_download/hls_downloader.py
python code/data_download/HLS_SuPER.py -roi '144.4001779629389, 71.24588582272926,146.16765743760948, 71.6168891179392' -dir /home/nrietz/data/raster/hls/ -start 2020-05-01 -end 2020-10-31 -bands RED,GREEN,BLUE,NIR1,SWIR1,SWIR2,FMASK -qf True


--output=home/$USER/scratch/logs/job_%j.out

echo 'finished'