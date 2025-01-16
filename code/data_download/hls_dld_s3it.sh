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

python code/data_download/hls_downloader.py

--output=home/$USER/scratch/logs/job_%j.out

echo 'finished'