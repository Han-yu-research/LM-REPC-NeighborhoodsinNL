## Step 2. Landscape metrics (`landscape.py`)

This step prepares neighborhood-level raster inputs and processes landscape metrics calculated by FRAGSTATS.

### 2.1 Raster clipping and FRAGSTATS input preparation

The script first clips the raster maps by neighborhood boundaries using `buurtcode` as the unique identifier.

Two raster datasets are processed:

- `Reclass_LGN22.tif`: land-use classification raster
- `Green_4+1.tif`: green-space classification raster

The clipped raster files are saved to:

- `output/clipped_rasters/`
- `output/clipped_rasters2/`

The script also automatically generates `.fbt` batch files for FRAGSTATS:

- `output/4frg_5types/import_LGN.fbt`
- `output/4frg_4green/import_GREEN.fbt`

### 2.2 FRAGSTATS processing

The `.fbt` files should be manually imported into FRAGSTATS.

In FRAGSTATS, the following metrics are calculated:

For the LGN raster:

- Class-level metrics: `PLAND`, `AI`, `IJI`
- Landscape-level metrics: `PD`, `SHDI`

For the green-space raster:

- Class-level metrics: `PLAND`, `AI`

After running FRAGSTATS, export the results and place them in the corresponding output folders:

- `output/4frg_5types/`
- `output/4frg_4green/`

### 2.3 Post-processing FRAGSTATS results

The second part of `landscape.py` processes the FRAGSTATS outputs.

It converts `.class` and `.land` files into `.csv` files, extracts `buurtcode`, reshapes class-level metrics into wide format, and merges class-level and landscape-level results.

The generated files include:

- `LGN_merged_class_land.csv`
- `GREEN_merged_class.csv`
- `form_merged_green_again_all.csv`
- `form_merged_green_cleaned.csv`

The final cleaned landscape metrics file is:

```text
output/form_merged_green_cleaned.csv
