# Freshwater-beta-diversity
---
## Description
This repository contains the data and R scripts used in the manuscript:
**"Widespread functional loss under freshwater eutrophication"**.

**Review status:** 
This repository is currently shared exclusively for peer review.

## File structure

### `all_scripts/`
This directory contains all R scripts used for data processing, statistical analyses, and figure generation in the manuscript.
- `0.0_calculation_beta`    Functions of calculation of taxonomic and functional β-diversity metrics.

- `1.0_data_compilation_td`    Compile pairwise taxonomic diversity data alongside all predictor variables

- `1.1_data_compilation_fd`    Compile pairwise functional diversity data alongside all predictor variables

- `2.0_brm_model`    Bayesian modeling using **brms**.

- `3.0_plot_figure_2`  
- `3.1_plot_figure_3`  
- `3.2_plot_figure_4`  
- `3.3_plot_figure_S1-S3`  
  Scripts used to generate main-text and supplementary figures.

### `input_file/`
- `dataset_var.csv`  
  Dataset-level explanatory variables used in the models.

- `world_shp/`  
  World map base layers used to generate Figure S1 (global distribution of study sites).

- **Land-use data** were derived from the [*Global Dynamic Land Cover* dataset.](https://land.copernicus.eu/en/products/global-dynamic-land-cover)  
- **Climate data** were obtained from the [*WorldClim v2*](https://www.worldclim.org/data/worldclim21.html).

### `output_figure/`
- This directory contains all figures generated from the analyses, including both main-text and supplementary figures presented in the manuscript.  
  Figures are organized and named to correspond to figure numbers in the manuscript.
---

## Note  
The `raw_data` and the `output_file` directory will be made publicly available upon acceptance of the manuscript.  
At this stage, these data is shared solely for peer review.
