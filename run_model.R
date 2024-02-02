# Run Model
# Carlos Dobler & Caroline Williams
# 1/30/2024
# Read in files internally and run precipitation & drought models

# Import Libraries
library(stars)
library(tidyverse)
library(furrr)


# Define working directories
## Set source folder and destination folder paths
from_dir <- "/home/cwilliams/cmip6_directory/kgassert_TEMP/Chiapas/pr"

# Set working directory for data downloads and reading in files
working_dir <- "/home/cwilliams/persistent_disk/chiapas/precip/data/"

# Set project directory
project_dir <- "/home/cwilliams/persistent_disk/chiapas/chiapas_analysis"

# Set results directory
results_dir <- "/home/cwilliams/persistent_disk/chiapas/precip/results/"

## Set file pattern to detect in from_dir
model <- "GFDL-CM4"

# Copy data files from source to working directory
# ########### Uncomment below code chunk to copy files for first time ###########
# ## List and copy files from cmip6_data bucket to persistent drive
# from_dir %>%
#   list.files(full.names = TRUE, pattern = model) %>%
#   map(function(f) {
#     f_gs <- str_replace(f, "/home/cwilliams/cmip6_directory", "gs://cmip6_data")
# 
#     str_glue("gsutil cp {f_gs} {working_dir}") %>%
#       system()
#   })
# 
# ###############################################################################


# Run precipitation and drought analyses
## 1. Run change in precipitation
source("precip_change.R")

## 2. Run precipitation seasonality
source("precip_seasonality.R")

## 3. Run drought
source("drought.R")

# Export Results
## 1. "precip_change.R" exports
#mean_total_annual, change1, chaneg2, [extremes quantile 95 or dist_analysis (GEV_L)?]

## 2. "precip_seasonality.R" exports

## 3. "drought.R" exports