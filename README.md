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

`main_comparisonss.R`, compare OLS，SAR/SEM, SDM to select the most proper model

**Group priorities:**

`main_matrix_AIC.R`, block-wise model comparisons

**Mainstream of analysis:**

`main_analysis_updated.ipynb`, merged energy and landscape data together, and processing correlation analysis and VIF and OLS regression

`0311_main_SDM.R`, run SDM (KNN=6)in R

`main_figure_updated.ipynb`, transfer the results from R to process in Python for figures

**Robustness checks:**

`1007_KNN4+8.R`, run SDM (KNN=4/8) in R as comparisons

`1228_QUEEN.R`, run SDM (QUEEN) in R as a comparison

### Step 4. Supplementary analysis (`main_addition_updated.py`)

Additional analysis of LST spearman correlation，see supplementary material S14

## Data structure

## Requirements

## How to run

1. Run energy.py
2. Run landscape.py
3. Run analysis.py
4. Run figure.py
5. Run supplementary.py
