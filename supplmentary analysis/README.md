## Additional LST Analysis (`additional_LST_analysis.py`)

This script was added as an additional analysis in response to reviewer comments. It examines the relationship between Land Surface Temperature (LST), residential energy consumption per capita (REPC), and selected landscape metrics.

---

## Purpose

The analysis aims to explore whether LST is associated with:

* residential energy consumption per capita
* landscape metrics for specific land-cover types
* urbanity-level differences

The analysis is conducted for all neighborhoods and separately by `ste_mvs` urbanity level.

---

## Input files

```text
model/form_energy_final_cleaned.csv
shp/10580_filtered_neighborhood_boundariesx.shp
data/NL_LST_2022_winter_m_3.tif
data/Extract_NL_water3.tif
data/Extract_NL_built3.tif
data/Extract_NL_Dense_Short3.tif
data/Extract_NL_Dense_Tall3.tif
```

---

## Main workflow

### 1. Merge cleaned data with shapefile

The script merges the cleaned modelling dataset with the filtered neighborhood shapefile using `buurtcode`.

**Output:**

```text
shp/merged_cleaned_shapefile.shp
```

---

### 2. Check shapefile field-name truncation

Because shapefiles only support field names up to 10 characters, some variables may be shortened after saving.

The script reports possible mappings between original CSV field names and truncated shapefile field names.

Examples:

```text
ec_inten_percapita → ec_inten_p
PLAND_Dense_Short → PLAND_Dens
PLAND_Dense_Tall  → PLAND_De_1
AI_Dense_Short    → AI_Dense_S
AI_Dense_Tall     → AI_Dense_T
```

---

### 3. Full-neighborhood LST and REPC correlation

The script calculates median LST for each neighborhood using the full LST raster.

Then it calculates Spearman correlations between:

```text
REPC
LST_median_2022
```

The analysis is conducted for:

* all neighborhoods
* each `ste_mvs` urbanity group

**Outputs:**

```text
model/Spearman_REPC_LST_by_ste_mvs_and_all.csv
model/Median_REPC_LST_by_ste_mvs.csv
```

---

### 4. Land-cover-specific LST analysis

The script also calculates median LST for selected land-cover types:

```text
water
built
dense_short
dense_tall
```

Each land-cover type is linked to corresponding landscape metrics.

### Water

```text
LST raster: Extract_NL_water3.tif
Primary metric: AI_water
Secondary metric: IJI_water
```

### Built-up areas

```text
LST raster: Extract_NL_built3.tif
Primary metric: AI_built
Secondary metric: IJI_built
```

### Dense short vegetation

```text
LST raster: Extract_NL_Dense_Short3.tif
Primary metric: PLAND_Dens
Secondary metric: AI_Dense_S
```

### Dense tall vegetation

```text
LST raster: Extract_NL_Dense_Tall3.tif
Primary metric: PLAND_De_1
Secondary metric: AI_Dense_T
```

---

### 5. Spearman correlation by urbanity level

For each land-cover type, the script calculates Spearman correlations between land-cover-specific LST and the corresponding landscape metrics.

The correlations are calculated for:

* all neighborhoods
* each `ste_mvs` group

**Outputs:**

```text
model/LST_landcover_metric_spearman_summary.csv
model/LST_landcover_metric_spearman_grouped.csv
```

---

## Output files

```text
shp/merged_cleaned_shapefile.shp
model/Spearman_REPC_LST_by_ste_mvs_and_all.csv
model/Median_REPC_LST_by_ste_mvs.csv
model/LST_landcover_metric_spearman_summary.csv
model/LST_landcover_metric_spearman_grouped.csv
```

---

## Dependencies

Required Python packages:

```text
pandas
numpy
geopandas
rasterio
rasterstats
```

---

## How to run

From the project root directory:

```bash
python analysis/additional_LST_analysis.py
```

Or, if the script is stored directly in the project root:

```bash
python additional_LST_analysis.py
```

---

## Notes

This analysis is designed as an additional exploratory analysis for reviewer response.

The script uses median LST within each neighborhood or land-cover mask because median values are less sensitive to extreme raster values.

Spearman correlation is used because it does not assume a linear relationship and is suitable for assessing monotonic associations.

For land-cover-specific analyses, only valid observations with available LST values and positive landscape metric values are included.
