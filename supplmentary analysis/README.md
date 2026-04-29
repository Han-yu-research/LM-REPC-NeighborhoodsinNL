# Supplementary Analysis

This folder contains the additional Land Surface Temperature (LST) analysis conducted for the reviewer response.

The analysis explores the relationship between residential energy consumption per capita (REPC), neighborhood-level LST, and selected landscape metrics.

## Workflow

The script performs the following steps:

1. Merge the cleaned modelling dataset with the filtered neighborhood shapefile
2. Check possible shapefile field-name truncation caused by the 10-character shapefile limit
3. Calculate neighborhood-level median LST from the full LST raster
4. Calculate Spearman correlation between REPC and LST
5. Repeat the REPC-LST correlation by urbanity level
6. Calculate median LST for selected land-cover classes
7. Calculate Spearman correlations between class-specific LST and related landscape metrics

## Input Files

Required input files:

### 1. Cleaned modelling dataset

Path: model/form_energy_final_cleaned.csv

---

### 2. Neighborhood boundary shapefile

Path: shp/10580_filtered_neighborhood_boundariesx.shp

---

### 3. Winter LST raster data

The full LST raster was generated using code:

```text
data/NL_LST_2022_winter_m_3.tif
```

Description: These raster files are used to calculate neighborhood-level median LST values through zonal statistics.

## Main Outputs

```text
shp/merged_cleaned_shapefile.shp
model/Spearman_REPC_LST_by_ste_mvs_and_all.csv
model/Median_REPC_LST_by_ste_mvs.csv
model/Spearman_neighborhood_LST_landscape_metrics.csv
```

## Main Analyses

The supplementary analysis includes two parts:

1. **REPC and full-neighborhood LST**

   Spearman correlations are calculated between REPC and neighborhood-level median LST for all neighborhoods and by urbanity level.

2. **Neighborhood-level LST and landscape metrics**

   Median LST is calculated for selected land-cover classes:


## Dependencies

```text
pandas
numpy
geopandas
rasterio
rasterstats
```

## Run

```bash
python additional_LST_analysis.py
```
