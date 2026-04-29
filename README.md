# Project Title

Measuring the Impact of Urban Morphology on Residential Energy Consumption at the Neighborhood-level in the Netherlands 

## Project overview

This repository contains the workflow and scripts used to investigate the relationship between urban morphology and residential energy consumption at the neighborhood level in the Netherlands.

## Workflow

### Step 1. Energy data processing (`main_energy_updated.ipynb`)

Process and clean residential energy consumption data.

See folder `energy`

**Input:** Raw residential energy consumption and demographic data from CBS   
**Output:** Cleaned energy indicators, building characteristics, and selected demographic data for neighborhood-level analysis

---

### Step 2. Landscape metrics calculation (`main_landscape_updated.ipynb`)

Calculate landscape metrics for neighborhood green spaces.

See folder `landscape`

**Input:** - Reclassified land-use map based on LGN land-use data (water and built-up)
           - Reclassified vegetation map based on RIVM vegetation data 
           
**Output:** landscape metrics dataset

---

### Step 3. Statistical analysis 

**Data merging and preliminary analysis:**

`main_analysis_updated.ipynb`, merged energy and landscape data together, and processing correlation analysis and VIF and OLS regression


**Mainstream of analysis:**

`main_SDM_updated.R`, estimates the Spatial Durbin Model (SDM) using a K-nearest neighbor spatial weights matrix with `k = 6`

`main_comparisons_updated.R`, compares OLS, SAR, SEM, and SDM models to identify the most appropriate spatial model

`main_matrix_AIC_updated.R`, conducts block-wise model comparison based on variable groups


**Robustness checks:**

`KNN4+8_updated.R`, estimates SDM models using KNN spatial weights with `k = 4` and `k = 8`

`QUEEN_updated.R`, estimates the SDM using a queen contiguity spatial weights matrix


**Figure production:**

`main_figure_updated.ipynb`, processes model outputs from R and produces figures in Python

---

### Step 4. Supplementary analysis 

- `main_addition_updated.ipynb`: performs additional Spearman correlation analysis with land surface temperature (LST), corresponding to Supplementary Material S14
- 
- `2022_winter_LST_updated.ipynb`: generates the 2022 winter LST map used in the supplementary analysis

---

## How to run

Each folder contains a dedicated README file with detailed instructions for running the corresponding scripts.

## Data availability

Due to data licensing restrictions, some raw datasets are not included in this repository. The scripts are provided to document the workflow and support reproducibility where data access is available.

## Software requirements

The analysis uses both Python and R.

Main Python packages include:

- pandas
- geopandas
- rasterio
- numpy
- scipy
- statsmodels
- matplotlib

Main R packages include:

- sf
- spdep
- spatialreg
- dplyr
- readr

## Repository structure

```text
Repository
├── energy/
├── landscape/
├── analysis/
└── supplementary_analysis/
