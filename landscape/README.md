# Landscape Module

This module preprocesses land-cover and green-space raster data and prepares landscape metrics at the neighborhood level.

The workflow clips raster datasets by neighborhood boundaries, generates Fragstats batch files, processes Fragstats outputs, and merges land-cover and green-space landscape metrics into a cleaned table for further analysis.

## Input Data

Raw raster and shapefile data are not included in this repository.

Please prepare the following files and place them in the `/data` directory.

---

### 1. Neighborhood Boundary Shapefile

**Required filename:**  
`main_buurten selection.shp`

**Source:** BAG data

**Download link:** https://www.3dbag.nl/en/download

**Description:**  Neighborhood boundary shapefile used to clip raster data by neighborhood polygons.

---

### 2. Reclassified LGN Land-cover Raster

**Required filename:**  
`Reclass_LGN22.tif`

**Description:**  Reclassified LGN 2022 land-cover raster.

This raster is used to calculate main land-cover landscape metrics, including class-level and landscape-level metrics of water and built-up areas.

---

### 3. Green-space Classification Raster

**Required filename:**  
`Green_4+1.tif`

**Description:**  Green-space classification raster.

This raster is used to calculate green-space composition and configuration metrics.

## Workflow

The script performs the following steps:

1. Create the project folder structure  
2. Clip the LGN raster by each neighborhood polygon  
3. Clip the green-space raster by each neighborhood polygon  
4. Generate Fragstats batch files (`.fbt`) for clipped rasters  
5. Run Fragstats manually using the generated batch files  
6. Convert Fragstats output files to CSV  
7. Process class-level and landscape-level metrics  
8. Merge LGN and green-space landscape metrics  
9. Clean and rename final variables  

## Manual Fragstats Step

This script generates Fragstats batch files, but Fragstats must be run manually.

After running the raster clipping and `.fbt` generation steps:

1. Open Fragstats  
2. Import the generated `.fbt` files  
3. Run class-level metrics for the LGN raster:
   - `PLAND`
   - `AI`
   - `IJI`
4. Run landscape-level metrics for the LGN raster:
   - `PD`
   - `SHDI`
5. Run class-level metrics for the green-space raster:
   - `PLAND`
   - `AI`
6. Export the Fragstats results into the corresponding output folders  
7. Run the result processing section of the script  

## Main Outputs

Main output folders and files are saved in the `/output` directory.

### Clipped Raster Outputs

```text
output/clipped_rasters/
output/clipped_rasters2/
```

### Fragstats Working Folders

```text
output/4frg_5types/
output/4frg_4green/
```

### Generated Fragstats Batch Files

```text
output/4frg_5types/import_LGN.fbt
output/4frg_4green/import_GREEN.fbt
```

### Processed Fragstats Outputs

```text
output/4frg_5types/LGN_merged_class_land.csv
output/4frg_4green/GREEN_merged_class.csv
```

### Final Outputs

```text
output/form_merged_green_again_all.csv
output/form_merged_green_cleaned.csv
```

The main final output used for further analysis is:

```text
output/form_merged_green_cleaned.csv
```

## Key Metrics

The script processes the following landscape metrics:

| Metric | Description |
|---|---|
| `PLAND` | Percentage of landscape occupied by each class |
| `AI` | Aggregation index |
| `IJI` | Interspersion and juxtaposition index |
| `PD` | Patch density |
| `SHDI` | Shannon diversity index |

## Main Variables

The final cleaned table includes renamed land-cover and green-space variables, such as:

```text
AI_water
AI_built
AI_Dense_Short
AI_Dense_Tall
AI_Sparse_Short
AI_Sparse_Tall
IJI_built
IJI_water
PLAND_water
PLAND_built
PLAND_agriculture
PLAND_Dense_Short
PLAND_Dense_Tall
PLAND_Sparse_Short
PLAND_Sparse_Tall
```

## Requirements

Required Python packages:

```bash
pip install pandas geopandas rasterio
```

Fragstats is also required for calculating landscape metrics.

## Run

```bash
python main_landscape_updated.ipynb
```

## Notes

The script uses the current working directory as the project root:

```python
BASE_DIR = os.getcwd()
```

