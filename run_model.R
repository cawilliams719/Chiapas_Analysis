# Run Model
# Carlos Dobler & Caroline Williams
# 1/30/2024
# Read in files internally and run precipitation & drought models

# Import Libraries
library(stars)
library(tidyverse)
library(furrr)
library(units)


# Define working directories
source("file_paths.R")
# ## Set source folder and destination folder paths
# from_dir <- "<INSERT_INITIAL_DATA_LOCATION>"
# 
# # Set working directory for data downloads and reading in files
# working_dir <- "<INSERT_CURRENT_DATA_LOCATION>"
# 
# # Set project directory
# project_dir <- "<INSERT_chiapas_analysis_REPO_DIRECTORY>"
# 
# # Set results directory
# results_dir <- "<INSERT_RESULTS_DIRECTORY>"

## Set file pattern to detect in from_dir
model <- "GFDL-CM4"
ssp <- "ssp285"

# Run precipitation and drought analyses
## 1. Run change in precipitation
system.time(source("precip_change.R"))
gc()

## 2. Run precipitation seasonality
source("precip_seasonality.R")

## 3. Run drought
source("drought.R")

# Export Results
## 1. "precip_change.R" exports
#mean_total_annual, change1, chaneg2, [extremes quantile 95 or dist_analysis (GEV_L)?]

## 2. "precip_seasonality.R" exports

## 3. "drought.R" exports
