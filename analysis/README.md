# Analysis Workflow

This folder contains the scripts for statistical analysis, spatial modelling, robustness checks, and figure generation.

The workflow combines Python and R scripts. Python is mainly used for data preparation, OLS, correlation, VIF diagnostics, and figure generation. R is used for Spatial Durbin Model (SDM) estimation, model selection, robustness checks, and variable-group comparison.

## Scripts

```text
main_analysis_updated.ipynb
main_SDM_updated.R
main_comparisons_updated.R
KNN4+8_updated.R
Queen_updated.R
main_matrix_AIC_updated.R
main_figure_updated.ipynb
```

## Workflow

1. `main_analysis_updated.ipynb`  
   Merges energy/building and landscape data, cleans the modelling dataset, filters the shapefile, runs summary statistics, Spearman correlation, VIF diagnostics, and OLS regression.

2. `main_SDM_updated.R`  
   Estimates the main Spatial Durbin Model using `KNN = 6`, calculates direct, indirect, and total impacts, saves model diagnostics and residuals.

3. `main_comparisons_updated.R`  
   Compares OLS, SAR, SEM, and SDM using Moran’s I, LM tests, LR tests, log-likelihood, AIC, and BIC.

4. `KNN4+8_updated.R`  
   Tests robustness using alternative spatial weights: `KNN = 4` and `KNN = 8`.
   
5. `Queen_updated.R`  
   Tests robustness using alternative spatial weights: `QUEEN`.
   
6. `main_matrix_AIC_updated.R`  
   Compares the full SDM with models excluding demographic, building, or landscape variable groups.

7. `main_figure_updated.ipynb`  
   Generates urbanity-level violin plots and SDM impact plots.

## Input Files

Main input files:

```text
output/neighborhoods_energy_building_characteristics_merged_more.csv
output/form_merged_green_cleaned.csv
data/main_buurten selection.shp
data/kwb-2022.xlsx
```

Main modelling inputs generated during the workflow:

```text
model/form_energy_final_cleaned_renamed.csv
shp/10580_filtered_neighborhood_boundariesx.shp
output/sdm_impacts_with_p.xlsx
```

## Main Outputs

```text
model/form_energy_final_cleaned_renamed.csv
output/sdm_impacts_with_p.xlsx
models/
tables/
logs/
figures/
```

Important figure outputs:

```text
figures/violin_urbanity_*.png
figures/sdm_impacts_SocioEconomic_Building.png
figures/sdm_impacts_LandscapeMetrics.png
```

## Dependencies

Python:

```text
pandas
numpy
geopandas
matplotlib
seaborn
statsmodels
scikit-learn
openpyxl
```

R:

```text
sf
spdep
spatialreg
dplyr
readr
tidyr
ggplot2
writexl
car
```

## Run Order

```text
1. main_analysis_updated.ipynb
2. main_SDM_updated.R
3. main_comparisons_updated.R
4. KNN4+8_updated.R and Queen_updated.R
5. main_matrix_AIC_updated.R
6. main_figure_updated.ipynb
```
