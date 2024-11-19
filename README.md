# Fine-scale burn patterns in Siberian tundra fires - Code and Data Repository
This repository contains the code and tabular data to generate the output for Rietze et al. (in prep.): Pre-fire Vegetation Conditions and Topography Shape Burn Mosaics of Siberian Tundra Fire Scars.

Last update to this readme: 11 November 2024.

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
│   │
│   │  ZOIB_model.R
│   │
│   ├───classification
│   │   │  prep.R
│   │   │  predict_burned_area.R
│   │   │  predict_water.py
│   │
│   └───figures_and_tables
│          Table_1.R
│          figure_1.R
│          figure_2.R
│          figure_3.R
│    
├───data
│   └───geodata
│       ├───feature_layers
│       └───raster
│
├───figures  
│
└───tables
```

- The scripts in `classification` are used to prepare the PlanetScope imagery and run the random forest classification.
  - `prep.R` is used to crop, rename and prepare raster files for the image classification.
  - `predict_burned_area.R` is used to execute the classification of burned areas, performing validation and predicting the burned area maps.
- The scripts in `figures_and_tables` are used to generate the main and supplementary figures as well as the supporting tables and are named appropriately.
- The script `ZOIB_model.R` contains the code for the zero-one inflated beta regression. 
- The folder `data` is empty and should contain the data that can be downloaded from Zenodo (see link on top).
- The folder `figures` contains the figures that are produced in the correspoinding scripts.
- The folder `tables` contains the tables that are produced in the correspoinding scripts.

[\[back to content\]](#content)

## 2. Required data and software
The necessary data to run the code is publicly available under [![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.12650945-blue)](https://doi.org/10.5281/zenodo.12650945).

The data pre-processing and data analysis was using R 4.2.2 (2022-10-31 ucrt). Newer versions of these software packages will likely work, but have not been tested.

Code development and processing were carried out in Windows 10 (64 bit), but execution should (in theory) be platform independent.

[\[back to content\]](#content)

## 3. Contact
Code development and maintenance: Nils Rietze ([nils.rietze@uzh.ch](nils.rietze@uzh.ch))

[\[back to content\]](#content)

## 4. Acknowledgements

From the manuscript:
*N.R. was supported through the TRISHNA Science and Electronics Contribution (T-SEC), ESA PRODEX Trishna T-SEC project (PEA C4000133711). Field work and vegetation sample processing were conducted in the scope of State Assignment of the Ministry of Science and Higher Education of the Russian Federation (Project  АААА-А21-121012190038-0), using the equipment of the Centre for collective use of Federal Research Centre «Yakut Scientific Centre» (grant no. 13.TsKP.21.0016). We would like to thank Planet Labs for free access to PlanetScope imagery. We would like to thank Tim Gyger for helpful discussions regarding our statistical analysis. The authors declare no competing interests.*

[\[back to content\]](#content)

## 5. Citation
When citing elements in this repository, please cite as:

Rietze, N., Heim, R., Troeva, E., Schaepman-Strub, G., Assmann, J. J. (in prep.). Pre-fire Vegetation Conditions and Topography Shape Burn Mosaics of Siberian Tundra Fire Scars. 

[\[back to content\]](#content)

## 6. License
The scripts in this repository (*.R files) are licensed under the MIT license (see [license text](https://github.com/nrietze/SiberiaFires/blob/main/LICENSE)).<br>

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />The remaining content in this repo is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

[\[back to content\]](#content)
