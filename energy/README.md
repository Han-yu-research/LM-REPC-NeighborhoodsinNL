# Energy Module

This module preprocesses residential energy, building, and neighborhood socio-economic data for the study.

The workflow integrates building identity data, energy label records, and neighborhood statistics to generate neighborhood-level energy consumption indicators and part of the urban morphology variables (building characteristics and demographic indicators).

---

## Input Data

All raw datasets used in this study are publicly available from open-data platforms.

Raw datasets are not included in this repository.  
Please download the datasets below and place them in the `/data` directory.

---

### 1. Building Identity Dataset

**Source:** BAG data  
**Dataset description:** Building-level identity and neighborhood matching table  
**Download link:** https://www.3dbag.nl/en/download  

**Required filename:**  
`0106_building_indentity_neighborhoodsname.csv`

**Notes:**  
This study uses the 2022 version of the building dataset.  
The public download link only provides the latest available version.

To reproduce the input file:

1. Download the BAG building dataset  
2. Import the building footprint layer and neighborhood boundary layer into ArcGIS Pro  
3. Perform the Identity tool to spatially match buildings with neighborhood boundaries  
4. Remove edge-split buildings caused by boundary intersections  
5. Export the resulting attribute table as:

`0106_building_indentity_neighborhoodsname.csv`

---

### 2. Energy Label Dataset

**Source:** Netherlands Enterprise Agency (RVO)  
**Dataset description:** Residential energy label records  
**Download link:** https://nationaalgeoregister.nl/geonetwork/srv/dut/catalog.search#/metadata/8dff9ab0-dc82-4143-866f-2c08450abf61  

**Required filename:**  
`1205_energy_label_full.csv`

**Notes:**  
The dataset is provided as a WMS service.

To reproduce the input file:

1. Load the WMS layer into ArcGIS Pro  
2. Export the attribute table as CSV  
3. Save it as:

`1205_energy_label_full.csv`

---

### 3. Neighborhood Socio-economic Dataset

**Source:** Statistics Netherlands (CBS)  
**Dataset description:** Neighborhood socio-economic indicators and energy statistics (2022)  
**Download link:** https://www.cbs.nl/nl-nl/maatwerk/2025/13/kerncijfers-wijken-en-buurten-2022  

**Required filename:**  
`kwb-2022.xlsx`

---

## Workflow Overview

The workflow consists of the following steps:

1. Clean building identity data  
2. Remove duplicated buildings  
3. Remove missing neighborhood codes  
4. Filter neighborhood shapefile  
5. Process energy label data  
6. Merge building and energy label data  
7. Summarize building characteristics at neighborhood level  
8. Clean socio-economic data  
9. Calculate residential energy consumption  
10. Merge energy and building data  
11. Calculate urban morphology indicators (building characteristics and demographic indicators)
12. Export final shapefile  

---

## Main Outputs


All outputs are saved in the `output/` and `shp/` folders.

Main output files:

- `neighborhoods_summary_building_characteristics.csv`
- `neighborhoods_energy_table.csv`
- `neighborhoods_energy_building_characteristics_merged.csv`
- `neighborhoods_energy_building_characteristics_merged_more.csv`

Final shapefile:

- `*_filtered_neighborhood_boundaries_final.shp`

---

## Requirements

Required packages:

```bash
pip install pandas numpy geopandas openpyxl
```

## Run

```bash
python main_energy_updated.ipynb
```

## Final Output

The final dataset used for further analysis:

```text
output/neighborhoods_energy_building_characteristics_merged_more.csv
```
