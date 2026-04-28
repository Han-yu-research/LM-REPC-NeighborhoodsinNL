# Analysis Workflow

This folder contains the scripts used for statistical analysis, spatial modelling, robustness checks, and figure generation.

The analysis workflow combines Python and R scripts.

## Scripts

```text
analysis_OLS.py
sdm_analysis.R
sdm_model_selection.R
sdm_knn_robustness.R
sdm_dropoff_analysis.R
figure.py
```

---

## Workflow overview

The analysis workflow consists of six steps.

---

## Step 1. Preliminary statistical analysis (`analysis_OLS.py`)

This Python script prepares the final modelling dataset and performs preliminary statistical analysis.

Main tasks:

* merge processed energy/building and landscape datasets
* select modelling variables
* calculate summary statistics
* clean missing values
* filter the neighborhood shapefile
* rename variables
* perform Spearman correlation analysis
* perform VIF diagnostics
* estimate OLS regression

Main output:

```text
model/form_energy_final_cleaned_renamed.csv
```

This file is used as the input for all SDM analyses.

---

## Step 2. Main Spatial Durbin Model (`sdm_analysis.R`)

This R script estimates the main Spatial Durbin Model.

Main tasks:

* read the final modelling dataset
* read the filtered shapefile
* construct KNN spatial weights (`KNN = 6`)
* standardize explanatory variables
* log-transform the dependent variable
* estimate the SDM
* calculate direct, indirect, and total impacts
* export coefficients
* export impacts with significance
* export model diagnostics
* save residuals
* generate residual violin plots

Main outputs:

```text
output/sdm_impacts_with_p.xlsx
data/knn6_after_clean.gal
data/02_listw_knn6.rds
data/03_model_input_knn6.rds
models/
tables/
logs/
```

The impact file is used for figure generation.

---

## Step 3. SDM model selection (`sdm_model_selection.R`)

This script evaluates whether SDM is the most appropriate spatial specification.

Main tasks:

* Moran’s I test on OLS residuals
* LM tests
* robust LM tests
* SARMA test
* estimate SAR and SEM models
* compare OLS, SAR, SEM, and SDM
* LR simplification tests
* Wald simplification tests

Purpose:

To justify the use of SDM over alternative spatial models.

Main outputs:

```text
models/LM_tests_OLS.csv
models/LR_tests_SDM.csv
models/logLik_models.csv
models/model_comparison.csv
```

---

## Step 4. KNN robustness checks (`sdm_knn_robustness.R`)

This script tests whether the SDM results are robust to alternative spatial weight specifications.

Main model:

```text
KNN = 6
```

Robustness models:

```text
KNN = 4
KNN = 8
```

Main tasks:

* reconstruct spatial weights
* fit alternative SDM models
* compare model fit

Comparison metrics:

* log-likelihood
* AIC
* BIC
* rho
* pseudo R²
* adjusted pseudo R²

Main outputs:

```text
data/knn4_after_clean.gal
data/knn8_after_clean.gal
models/sdm_knn4.rds
models/sdm_knn8.rds
tables/sdm_knn_results_compare.csv
```

---

## Step 5. Variable-group drop-off analysis (`sdm_dropoff_analysis.R`)

This script evaluates the contribution of variable groups.

Compared models:

* Full model
* Drop demographic variables
* Drop building variables
* Drop landscape variables

Comparison metrics:

* log-likelihood
* AIC
* BIC
* ΔAIC
* ΔBIC

Purpose:

To assess the relative explanatory contribution of each variable group.

Main outputs:

```text
tables/{run_id}/sdm_matrix_bigcomp_metrics.csv
tables/{run_id}/sdm_matrix_bigcomp_dropoff_compare.csv
tables/{run_id}/sdm_dropoff_delta_barplot_styled_eigen_all_500grid.png
models/{run_id}/
logs/{run_id}/
data/{run_id}/
```

---

## Step 6. Figure generation (`figure.py`)

This Python script generates figures for interpretation and presentation.

Main tasks:

* violin plots across urbanity levels
* SDM impact plots

Input files:

```text
model/form_energy_final_cleaned_renamed.csv
output/sdm_impacts_with_p.xlsx
```

Output:

```text
figures/
```

---

## Input files

```text
output/neighborhoods_energy_building_characteristics_merged_more.csv
output/form_merged_green_cleaned.csv
data/main_buurten selection.shp
data/kwb-2022.xlsx
```

---

## Output folders

```text
data/
output/
weights/
src/
shp/
model/
models/
tables/
logs/
figures/
```

---

## Dependencies

### Python

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

### R

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

---

## Workflow order

Run the scripts in the following order:

```text
1. energy.py
2. landscape.py
3. FRAGSTATS processing
4. analysis_OLS.py
5. sdm_analysis.R
6. sdm_model_selection.R
7. sdm_knn_robustness.R
8. sdm_dropoff_analysis.R
9. figure.py
```

---

## How to run

### Python analysis

```bash
python analysis/analysis_OLS.py
```

### Main SDM

```bash
Rscript analysis/sdm_analysis.R
```

### SDM model selection

```bash
Rscript analysis/sdm_model_selection.R
```

### KNN robustness checks

```bash
Rscript analysis/sdm_knn_robustness.R
```

### Variable-group drop-off analysis

```bash
Rscript analysis/sdm_dropoff_analysis.R
```

### Figure generation

```bash
python analysis/figure.py
```
