# Project Title

## Project overview

This project investigates the relationship between urban morphology and residential energy consumption.

## Workflow

### Step 1. Energy processing (`main_energy_updated.ipynb`)

Process and clean residential energy consumption data.

**Input:** raw energy dataset
**Output:** cleaned energy indicators

### Step 2. Landscape metrics calculation (`main_landscape_updated.ipynb`)

Calculate landscape metrics for neighborhood green spaces.

**Input:** reclassified land-use map (main land-use types from LGN and main vegetation map from RIVM)
**Output:** landscape metrics dataset

### Step 3. Statistical analysis 

**Model selection:**

`main_comparisonss.R`, compare OLS, SAR/SEM, SDM to select the most proper model

**Group priorities:**

`main_matrix_AIC.R`, block-wise model comparisons

**Mainstream of analysis:**

`main_analysis_updated.ipynb`, merged energy and landscape data together, and processing correlation analysis and VIF and OLS regression

`main_SDM_updated.R`, run SDM (KNN=6)in R

`main_matrix_AIC_updated.R`, run block-wised model comparisons

`main_figure_updated.ipynb`, transfer the results from R to process in Python for figures

**Robustness checks:**

`KNN4+8_updated.R`, run SDM (KNN=4/8) in R as comparisons

`QUEEN_updated.R`, run SDM (QUEEN) in R as a comparison

### Step 4. Supplementary analysis 

Additional analysis of LST spearman correlation `main_addition_updated.ipynb`，see supplementary material S14

The 2022 winter LST map is generated from code: `2022_winter_LST_updated.ipynb`


## How to run

Every folder has their own README as instruction
