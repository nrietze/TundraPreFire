# Bottom-Up Factors Best Predict Tundra Burn Severity Shortly Before Ignition
This repository contains the code and tabular data to generate the output for Rietze et al. (in prep.): Bottom-Up Factors Best Predict Tundra Burn Severity Shortly Before Ignition.

Last update to this readme: 19 November 2025.

- [1. Repository structure](#1-repo-structure)
- [2. Required data and software](#2-required-data-and-software)
- [3. Contact](#3-contact)
- [4. Acknowledgements](#4-acknowledgements)
- [5. Citation](#5-citation)
- [6. License](#6-license)

## 1. Repository structure
Here is the structure of this repo, files have been excluded from this tree for readability.

```bash
├───code
│   ├───data_processing
│   ├───exploratory
│   └───figures_and_tables
├───data
│   ├───feature_layers
│   │   ├───1pct
│   │   ├───fire_atlas
│   │   │   ├───ignitions
│   │   │   └───sub_daily
│   │   │       ├───2020
│   │   │       │   ├───NFP
│   │   │       │   └───Snapshot
│   │   │       └───2022
│   │   │           ├───NFP
│   │   │           └───Snapshot
│   │   ├───sample_points_filtered
│   │   └───Sentinel2_MGRS_tiles.gdb
│   ├───models
│   ├───raster
│   │   ├───arcticDEM
│   │   ├───burned_area_descals
│   │   ├───hls
│   │   │   ├───optimality_rasters
│   │   │   ├───processed
│   │   │   ├───severity_rasters
│   │   ├───landsat8
│   │   │   └───lst
│   └───tables
│       ├───model_dataframes
│       │   └───1pct
│       └───sampled_data
│           └───1pct
├───figures
│   ├───bad_filtering
│   ├───burn_severity_maps
│   ├───dnbr_correction
│   ├───Figure_workflow
│   ├───FRP_vs_dNBR
│   └───manuscript
└───tmp
```

- The scripts in `code` are used to download and process HLS imagery. A visual overview of the processing steps is shown below:
![Processing workflow for the project.](https://github.com/nrietze/TundraPreFire/blob/figures/Figure_workflow/code_workflow.png?raw=true)

- `data_processing` contains the code to download HLS imagery from LPDAAC and Landsat LST, preprocess into spectral index rasters, and calculating burn severity rasters. It also contains R scripts to run the data extraction.
- `exploratory` contains code for the spline fitting and statistical analysis.
- The scripts in `figures_and_tables` are used to generate the main and supplementary figures as well as the supporting tables.
- The folder `data` contains the data generated throughout the analysis and should contain data that can be downloaded from the external sources (e.g., ArcticDEM, Arctic-Boreal fire atlas), and Zenodo.
- The folder `figures` contains the figures that are produced in the correspoinding scripts.
- The folder `tables` contains the tables that are produced as intermediate dataframes in the data sampling steps during the data extraction.

[\[back to content\]](#content)

## 2. Required data and software
Raster data is automatically downloaded:
- HLS from LPDAAC
- ArcticDEM from the PGC server
- Landsat-8 LST from Google Earth Engine

Feature layers delineating sub-daily (12h) and final fire perimeters can be downloaded from the Arctic-Boreal fire atlas:
  Scholten, R., Chen, Y., Veraverbeke, S., & Randerson, J. (2024). Arctic-boreal fire atlas: 12-hourly perimeters of individual fires in the Arctic-boreal domain from 2012 to 2023 [Dataset]. PANGAEA. [https://doi.org/10.1594/PANGAEA.967653](https://doi.org/10.1594/PANGAEA.967653)

The data pre-processing and data analysis was using python version 3.10.14 and R 4.2.2 (2022-10-31 ucrt). Newer versions of these software packages will likely work, but have not been tested.

The remote processing on the S3IT servers of the University of Zurich was carried out on the same versions of R and python but running on a Linux system. The batch processing scripts were executed using SLURM [https://en.wikipedia.org/wiki/Slurm_Workload_Manager](https://en.wikipedia.org/wiki/Slurm_Workload_Manager) and not developed for other workload managing software.

[\[back to content\]](#content)

## 3. Contact
Code development and maintenance: Nils Rietze (nils.rietze [at] pm.me)

[\[back to content\]](#content)

## 4. Acknowledgements

From the manuscript:
N.R. was supported through the TRISHNA Science and Electronics Contribution (T‐ SEC), ESA PRODEX Trishna T‐SEC project (PEA C4000133711). We would like to thank David S. Schimel and Charles E. Miller for constructive discussions regarding the interpretation of our results and Max van Gerrevink for his help regarding the optimality calculation.*

[\[back to content\]](#content)

## 5. Citation
When citing elements in this repository, please cite as:

N. Rietze, G. Schaepman-Strub, K. R. Miner, and J.J. Assmann (in prep.). Bottom-Up Factors Best Predict Tundra Burn Severity Shortly Before Ignition.

[\[back to content\]](#content)

## 6. License
The scripts in this repository (*.R files) are licensed under the MIT license (see [license text](https://github.com/nrietze/TundraFires/blob/main/LICENSE)).<br>

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />The remaining content in this repo is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

[\[back to content\]](#content)
