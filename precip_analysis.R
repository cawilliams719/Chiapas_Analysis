# Chiapas Precipitation Analysis
# Caroline Williams
# 12/7/2023
# Analysis of downscaled cmip6 precipitation data


# Import libraries
library(stars)
library(tidyverse)
library(stringr)
library(forstringr)
library(ggthemes)
library(viridis)
library(viridisLite)
library(cowplot)
library(furrr)
library(evd)

################################################################################
# Read in file paths & file pattern
## To replicate:
## 1. Create a separate R script in the same working directory called "file_paths.R"
## 2. Set variables working_dir and pattern (ex: working_dir <- "C:/users/file_path/working_dir/)
## This file path script can also be used to copy files to/from directories
source("file_paths.R")


# List of cmip6 model names
models <- working_dir %>% list.files(full.names = FALSE, pattern = pattern) %>% 
  str_extract_part("_ssp245", before = TRUE)

# # plan parameters
# plan(multisession, workers = 2)

# Read in precipitation data
precip <- working_dir %>% 
  list.files(full.names = TRUE, pattern = pattern) %>% 
  map(function(f) {
    read_ncdf(f, proxy = FALSE)
  }) 

# Set names of models in stars object list
precip <- setNames(precip, models)

# # plan parameters
# plan(sequential)

# Calculate Annual Precipitation (Total)
## Sum from Daily
total_annual <- precip[1:length(precip)] %>% 
  map(function(f) {
    filter(f, time >= "2001-01-01", time < "2061-01-01") %>% 
      aggregate(by = "1 year", FUN = sum) 
  }) 

# Annual Precipitation (Total) Means per Time Period
## 2001-2020
mean_total_annual_1 <- total_annual[1:length(total_annual)] %>%
  map(function(f) {
    filter(f, time < "2021-01-01") %>% 
      aggregate(by = "20 years", FUN = mean) %>%
      adrop()
  })

## 2021-2040
mean_total_annual_2 <- total_annual[1:length(total_annual)] %>% 
  map(function(f) {
    filter(f, time >"2020-01-01", time < "2041-01-01") %>% 
      aggregate(by = "20 years", FUN = mean) %>%
      adrop()
  })

## 2041-2060
mean_total_annual_3 <- total_annual[1:length(total_annual)] %>% 
  map(function(f) {
    filter(f, time > "2040-01-01") %>% 
      aggregate(by = "20 years", FUN = mean) %>% 
      adrop()
  })

# Change in Precipitation
## Change Period 2-1: (2021-2040) - (2001-2020)
subtracted1 <- mapply("-", mean_total_annual_2, mean_total_annual_1, SIMPLIFY = FALSE)

# Test output is correct
x <- subtracted1[[1]]
y <- mean_total_annual_2[[1]]-mean_total_annual_1[[1]]
x$pr-y$pr # all output pixels are 0 or NA

## Change Period 3-1: (2041-2060) - (2001-2020)
subtracted2 <- mapply("-", mean_total_annual_3, mean_total_annual_1, SIMPLIFY = FALSE)

# Percent Change in Precipitation
## Percent Change Function
perchange <- function(i, j) (((j-i)/i)*100) 

## Percent Change Period 2-1: (2021-2040) - (2001-2020)
perchangeim1 <- mapply(perchange, mean_total_annual_1, mean_total_annual_2, SIMPLIFY = FALSE)

# Test output is correct
a <- perchangeim1[[1]]
b <- (((mean_total_annual_2[[1]]-mean_total_annual_1[[1]])/mean_total_annual_1[[1]])*100)
a$pr-b$pr # all output pixels are 0 or NA

## Percent Change Period 3-1: (2041-2060) - (2001-2020)
perchangeim2 <- mapply(perchange, mean_total_annual_1, mean_total_annual_3, SIMPLIFY = FALSE)


# Precipitation Extremes
# Calculate Annual Precipitation (Maximum) per Time Period
## 2001-2020
max_annual_1 <- precip[1:length(precip)] %>%
  map(function(f) {
    filter(f, time >= "2001-01-01", time < "2021-01-01") %>% 
      aggregate(by = "1 year", FUN = max)
  })

## 2021-2040
max_annual_2 <- precip[1:length(precip)] %>% 
  map(function(f) {
    filter(f, time >="2021-01-01", time < "2041-01-01") %>% 
      aggregate(by = "1 year", FUN = max) 
  })

## 2041-2060
max_annual_3 <- precip[1:length(precip)] %>% 
  map(function(f) {
    filter(f, time >= "2041-01-01", time < "2061-01-01") %>% 
      aggregate(by = "1 year", FUN = max) 
  })

## Quantile
## Shows tails of maximum precipitation per pixel aggregated by time
quantile_99_1 <- max_annual_1[1:length(max_annual_1)] %>% 
  map(function(y) {
    st_apply(y, 2:3, function(q) { 
      quantile(q, 0.99, na.rm = TRUE)
    })
  })

quantile_99_2 <- max_annual_2[1:length(max_annual_2)] %>% 
  map(function(y) {
    st_apply(y, 2:3, function(q) { 
      quantile(q, 0.99, na.rm = TRUE)
    })
  })

quantile_99_3 <- max_annual_3[1:length(max_annual_3)] %>% 
  map(function(y) {
    st_apply(y, 2:3, function(q) { 
      quantile(q, 0.99, na.rm = TRUE)
    })
  })

## Empirical Cumulative Distribution Function (ECDF) Probability
# ECDF shows probability that precipitation is equal to or less than a certain number
ecdf_200_1 <- max_annual_1[1:length(max_annual_1)] %>% # outputs ecdf 200 for each year for both models
  map(function(y) {
    apply(y$pr, 1, function(q) { 
      ecdf(q)(200)
    })
  })

ecdf_200_2 <- max_annual_2[1:length(max_annual_2)] %>% # outputs ecdf 200 for each year for both models
  map(function(y) {
    apply(y$pr, 1, function(q) { 
      ecdf(q)(200)
    })
  })

ecdf_200_3 <- max_annual_3[1:length(max_annual_3)] %>% # outputs ecdf 200 for each year for both models
  map(function(y) {
    apply(y$pr, 1, function(q) { 
      ecdf(q)(200)
    })
  })

## Generalized Extreme Value (GEV) Distribution
gev_1 <- max_annual_1[1:length(max_annual_1)] %>% 
  map(function(y) {
    st_apply(y, 2:3, function(q) { 
      dgev(q, loc=0, scale=1, shape=0)
    })
  })

gev_2 <- max_annual_2[1:length(max_annual_2)] %>% 
  map(function(y) {
    st_apply(y, 2:3, function(q) { 
      dgev(q, loc=0, scale=1, shape=0)
    })
  })

gev_3 <- max_annual_3[1:length(max_annual_3)] %>% 
  map(function(y) {
    st_apply(y, 2:3, function(q) { 
      dgev(q, loc=0, scale=1, shape=0)
    })
  })


# Precipitation Seasonality
## Sum from Daily by wet (May-Oct) and dry (Nov-April) season
## TESTING SEASONAL AGGREGATION - In Progress
### dates corresponding to start of wet and dry season for each year
tm <- seq(as.POSIXct("2001-05-01"), by = "6 months", length.out = 120) # not sure if it aggregates between these dates

### Aggregate by wet and dry season
total_seasonal <- precip[1:length(precip)] %>% 
  map(function(f) {
    filter(f, time >= "2001-01-01", time < "2061-01-01") %>% 
      aggregate(by = tm, FUN = sum) 
  }) 

# Calculate Seasonal Precipitation (Total) per Year per Time Period
## 2001-2020
## 2021-2040
## 2041-2060


################################################################################
# Export Data
## Save workspace to load in visualizations
### Set file path and name for data
precip_results <- str_glue("{project_dir}/data/precip_analysis_results.RData")

### Remove directorys
rm("from_dir", "working_dir", "project_dir")

### Save data   
save.image(file = precip_results)

## Export NetCDF file results
### Create list of filenames to export results to
results_filenames <- str_glue("{results_dir}data/{models}.nc")

### Saves but removes refsys WGS 84 from metadata
mapply(write_mdim, mean_total_annual_1, results_filenames, SIMPLIFY = FALSE)
