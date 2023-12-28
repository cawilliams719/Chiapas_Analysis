# Chiapas Analysis

### Description
This repository features R code to analyze downscaled cmip6 precipitation data in Chiapas, Mexico.

#### Files:
<b>`file_paths.R`</b>: This internal script is required to run `precip_analysis.R`. This script includes the `working_dir` and `pattern` variables. To replicate users should: 

  1. Create a separate R script in the same working directory as `precip_analysis.R` called `file_paths.R` 
  2. In the file path R script, set variables `working_dir` and `pattern` (ex: working_dir <- "C:/users/file_path/working_dir/ & pattern <- "Chiapas")
  
<b>`precip_analysis.R`</b>: This script calculates total annual precipitation, difference, percent change, and statistical extremes. Inputs include cmip6 precipitation data, the `file_paths.R` script. Outputs include the precipitation analysis, resulting raster images, statistics, and visualizations.

#### R Packages Used:
`stars`, `tidyverse`, `stringr`, `forstringr`, `ggthemes`, `viridis`, `viridisLite`, `cowplot`, `furrr`, `evd`
