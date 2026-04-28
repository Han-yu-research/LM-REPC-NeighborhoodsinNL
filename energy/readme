# LM-REPC-NeighborhoodsinNL
Here is the data for the paper: Measuring the Impact of Urban Morphology on Residential Energy Consumption at the Neighborhood-level in the Netherlands 
# Urban Morphology and Residential Energy Consumption Analysis
This repository contains the full data processing and analysis workflow for examining the relationship between urban morphology indicators and residential energy consumption at the neighborhood level in the Netherlands.

The workflow integrates building-level data, energy label data, neighborhood socio-economic data, and spatial boundary data to construct neighborhood-level analytical datasets for spatial econometric modeling.

---

## Project Structure

```text
project/
├── data/                  # Raw input data (download separately)
├── output/                # Processed datasets
├── weights/               # Spatial weights matrices
├── src/                   # Source code
├── shp/                   # Processed shapefiles
├── main.py                # Main workflow script
├── requirements.txt       # Python dependencies
└── README.md
```

---

## Data Sources

All raw datasets used in this study are publicly available from open-data platforms.
Raw data are not included in this repository.
Please download the datasets below and place them in the `/data` directory.

### 1. Building Identity Dataset

Source: BAG data
Dataset description: Building-level identity and neighborhood matching table
Download link: https://data.3dbag.nl/v20250903/tile_index.fgb

Required filename:"0106_building_indentity_neighborhoodsname.csv"

Access date:/

```text
The version for this research is using the data from the year of 2022, however, from the link provided, you could only download the latest data.
By importing the building data into GIS pro, using the function "identity" to match the buildings and neighborhood boundaries, and following the seclection process to remove the buildings on the edge split by the boundaries and export the table and save it as the above required filename.

```

---

### 2. Energy Label Dataset

Source: Netherlands Enterprise Agency (RVO)
Dataset description: Residential energy label records
Download link: https://nationaalgeoregister.nl/geonetwork/srv/dut/catalog.search#/metadata/8dff9ab0-dc82-4143-866f-2c08450abf61

Required filename:"1205_energy_label_full.csv"

Access date:/

```text
The data downloaded is in the format of WMS, which could be open in GIS pro, and we have the final data after exporting from the GIS pro.
```

---

### 3. Neighborhood Socio-economic Dataset

Source: Statistics Netherlands (CBS)
Dataset description: Neighborhood socio-economic indicators (2022)
Download link:https://www.cbs.nl/nl-nl/maatwerk/2025/13/kerncijfers-wijken-en-buurten-2022

Required filename:

```text
kwb-2022.xlsx
```

Access date:

```text
/
```

---

## Workflow Overview

The workflow consists of the following steps:

### Step 1: Building Dataset Cleaning

* Remove duplicate building records
* Keep the record with the largest area
* Remove records without neighborhood code

Output:

```text
0106_building_indentity_neighborhoodsname_cleaned.csv
0106_building_indentity_neighborhoodsname_cleaned_dropnancode.csv
```

---

### Step 2: Spatial Filtering

* Filter neighborhood boundary shapefile based on valid neighborhood codes

Output:

```text
*_filtered_neighborhood_boundaries.shp
```

---

### Step 3: Energy Label Processing

* Keep essential variables only
* Remove duplicates
* Standardize building identifiers

Output:

```text
1205_energylabel_2022_slim.csv
```

---

### Step 4: Building-Energy Merge

* Merge building records with energy label data

Output:

```text
energylabel_building_identity_neighborhoods.csv
```

---

### Step 5: Building Characteristics Aggregation

Compute neighborhood-level:

* Building count
* Ground area
* Total floor area
* Residential unit count
* Share of high-efficiency buildings (A-label and above)

Output:

```text
neighborhoods_summary_building_characteristics.csv
```

---

### Step 6: Socio-economic Data Processing

* Filter neighborhood-level records
* Select variables
* Clean missing values

Output:

```text
buurt_filtered_table_columns_cleaned.csv
buurt_filtered_table_cleaned_selected.csv
```

---

### Step 7: Energy Consumption Estimation

Compute:

* Electricity consumption
* Gas consumption
* District heating consumption
* Total energy consumption
* Per capita energy consumption intensity

Output:

```text
neighborhoods_energy_table.csv
```

---

### Step 8: Final Dataset Construction

Merge:

* Building characteristics
* Socio-economic indicators
* Energy consumption indicators

Output:

```text
neighborhoods_energy_building_characteristics_merged.csv
```

---

### Step 9: Urban Morphology Indicators (Building Characteristics)

Compute:

* Building Density
* Floor Area Ratio (FAR)
* Open Space Ratio (OSR)
* Housing Density
* Average Floor Number
* Building Intensity

Output:

```text
neighborhoods_energy_building_characteristics_merged_more.csv
```

---

### Step 10: Final Spatial Boundary Export

Export final spatial units for modeling.

Output:

```text
*_filtered_neighborhood_boundaries_final.shp
```

---

## Installation

Create a Python environment and install dependencies:

```bash
pip install -r requirements.txt
```

---

## Dependencies

Main Python packages:

* pandas
* numpy
* geopandas

Install manually if needed:

```bash
pip install pandas numpy geopandas
```

---

## Running the Workflow

Run the main script:

```bash
python main.py
```

The script will automatically:

* create project directories
* clean raw datasets
* merge datasets
* compute indicators
* export final datasets

---

## Output Files

Main analytical outputs:

```text
output/neighborhoods_summary_building_characteristics.csv
output/neighborhoods_energy_table.csv
output/neighborhoods_energy_building_characteristics_merged.csv
output/neighborhoods_energy_building_characteristics_merged_more.csv
shp/*_filtered_neighborhood_boundaries_final.shp
```

---

## Reproducibility

To reproduce the analysis:

1. Download all raw datasets
2. Place them in `/data`
3. Install dependencies
4. Run `main.py`

All intermediate and final datasets will be automatically generated.

---

## Citation

If you use this code, please cite:

[TO BE FILLED AFTER PUBLICATION]

---
