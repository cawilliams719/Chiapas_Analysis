# Run Model
# Carlos Dobler & Caroline Williams
# 1/30/2024
# Read in files internally and run precipitation & drought models

# Import Libraries
library(stars)
library(tidyverse)
library(furrr)
library(units)


# Set file pattern to detect in from_dir
model <- "UKESM1-0-LL"
ssp <- "ssp245"

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

# Run precipitation and drought analyses
## 1. Run change in precipitation
system.time(source("precip_change.R"))
gc()

## 2. Run precipitation seasonality
system.time(source("precip_seasonality_2.R"))
gc()

## 3. Run drought
system.time(source("drought.R"))
gc()

# Export Results
## Call time vector and write nc spatial routines scripts
source("https://raw.githubusercontent.com/carlosdobler/spatial-routines/master/time_vector.R")
source("https://raw.githubusercontent.com/carlosdobler/spatial-routines/master/write_nc.R")
source("write_nc_mod.R") # modified to export as double data type

## 1. "precip_change.R" exports
### Export Mean Total Annual Precipitation
time_export <- time_vector(st_get_dimension_values(mean_total_annual, "time"), daily = F)
write_nc(mean_total_annual, time_export, str_glue("{results_dir}Chiapas_pr_{model}_{ssp}_mean_total_annual_2001-2060.nc"))

### Export Change in Precipitation (Difference & % Change)
time_export <- PCICt::as.PCICt(as.character("2021-01-01", "2041-01-01"), cal = "gregorian")
write_nc(change, time_export, str_glue("{results_dir}Chiapas_pr_{model}_{ssp}_change_2001-2060.nc"))

### Export Precipitation Extremes
time_export <- PCICt::as.PCICt(as.character("2001-01-01"), cal = "gregorian")
write_nc(dist_analysis, time_export, str_glue("{results_dir}Chiapas_pr_{model}_{ssp}_extremes_2001-2060.nc"))

## 2. "precip_seasonality.R" exports
time_export <- time_vector(st_get_dimension_values(mean_total_annual, "time"), daily = F)
write_nc(l_seas, time_export, str_glue("{results_dir}Chiapas_pr_{model}_{ssp}_seasonality_2001_2001-2060.nc"))


## 3. "drought.R" exports
time_export <- time_vector(st_get_dimension_values(drought, "time"), daily = F)
write_nc_mod(drought, time_export, str_glue("{results_dir}Chiapas_pr_{model}_{ssp}_drought_2001-2060.nc"))

time_export <- time_vector(st_get_dimension_values(extreme_drought, "time"), daily = F)
write_nc(extreme_drought, time_export, str_glue("{results_dir}Chiapas_pr_{model}_{ssp}_extreme_drought_2001-2060.nc"))

# s_gamma_gof <- s_gamma_gof %>% st_set_dimensions(names = c("lon", "lat", "month"))
time_export <- PCICt::as.PCICt(c("2001-01-01", "2001-02-01","2001-03-01",
                                 "2001-04-01", "2001-05-01","2001-06-01",
                                 "2001-07-01", "2001-08-01","2001-09-01",
                                 "2001-10-01", "2001-11-01","2001-12-01"), cal = "gregorian")
# write_nc(s_gamma_gof, time_export, str_glue("{results_dir}test_output/test_run/Chiapas_pr_{model}_{ssp}_drought_2001-2060_write_nc.nc"))
write_nc(s_gamma_gof, time_export, str_glue("{results_dir}Chiapas_pr_{model}_{ssp}_gamma_gof_2001-2060.nc"))
