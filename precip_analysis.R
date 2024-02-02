# Chiapas Precipitation Analysis
# Caroline Williams
# 12/7/2023
# Analysis of downscaled cmip6 precipitation data


# Import libraries
library(stars)
library(tidyverse)
# library(stringr)
# library(forstringr)
# library(ggthemes)
# library(viridis)
# library(viridisLite)
# library(cowplot)
library(furrr)
# library(evd)

plan(multisession)

################################################################################
# Read in file paths & file pattern
## To replicate:
## 1. Create a separate R script in the same working directory called "file_paths.R"
## 2. Set variables working_dir and pattern (ex: working_dir <- "C:/users/file_path/working_dir/)
## This file path script can also be used to copy files to/from directories
source("file_paths.R")


# List of cmip6 model names
models <- 
  working_dir %>% 
  list.files(pattern = pattern) %>%
  str_split("_", simplify = T) %>% 
  .[,3]


## Percent Change Function
perchange <- function(i, j) (((j-i)/i)*100) 



# Loop through models

models %>% 
  walk(function(m) {
    
    f <- 
      working_dir %>% 
      list.files(full.names = T, pattern = pattern) %>% # pattern = m
      str_subset(m)
    
    time_vector <- 
      f %>% 
      read_ncdf(proxy = T) %>% 
      st_get_dimension_values(3) %>% 
      suppressMessages() %>% 
      str_sub(end = 10)
    
    max_feb <- 
      time_vector[str_sub(time_vector, 6,7) == "02"] %>% # filter feb months
      str_sub(9,10) %>% # extract days
      as.numeric() %>% 
      max()
    
    model_cal <- 
      case_when(max_feb == 30 ~ "360_day",
                max_feb == 29 ~ "gregorian",
                max_feb == 28 ~ "noleap")
    
    time_vector <- 
      PCICt::as.PCICt(time_vector, cal = model_cal)
    
    time_first <- 
      which(time_vector == PCICt::as.PCICt("2001-01-01", cal = model_cal))
    
    time_last <- 
      which.min(time_vector < PCICt::as.PCICt("2061-01-01", cal = model_cal))
    
    # Read in precipitation data
    precip <- 
      f %>% 
      read_ncdf(ncsub = cbind(start = c(1,  1, time_first),
                              count = c(NA, NA, time_last-time_first)),
                proxy = F) %>% 
      suppressMessages()
    
    # Calculate Annual Precipitation (Total)
    ## Sum from Daily
    total_annual <- 
      precip %>% 
      aggregate(by = "1 year", FUN = sum) %>% 
      aperm(c(2,3,1))
    
    # Annual Precipitation (Total) Means per Time Period
    # mean_total_annual
    mean_total_annual <- 
      total_annual %>%
      aggregate(by = "20 years", FUN = mean) %>% 
      aperm(c(2,3,1))
    
    # Change in Precipitation
    ## Change Period 2-1: (2021-2040) - (2001-2020)
    subtracted1 <- slice(mean_total_annual, time, 2) - slice(mean_total_annual, time, 1)
    ## Change Period 3-1: (2041-2060) - (2001-2020)
    subtracted2 <- slice(mean_total_annual, time, 3) - slice(mean_total_annual, time, 1)

    # Percent Change in Precipitation
    ## Percent Change Period 2-1: (2021-2040) - (2001-2020)
    perchangeim1 <- perchange(slice(mean_total_annual, time, 1), slice(mean_total_annual, time, 2))
    ## Percent Change Period 3-1: (2041-2060) - (2001-2020)
    perchangeim2 <- perchange(slice(mean_total_annual, time, 1), slice(mean_total_annual, time, 3))
    
    
    # Precipitation Extremes
    # Calculate Annual Precipitation (Maximum) per Time Period
    max_annual <- 
      precip %>% 
      aggregate(by = "1 year", FUN = max) %>% 
      aperm(c(2,3,1))
    
    max_annual_periods <- 
      list(1:20, 21:40, 41:60) %>% 
      map(~slice(max_annual, time, .x))
    
    
    ## Quantile
    ## Shows tails of maximum precipitation per pixel aggregated by time
    quantile_99 <- 
      max_annual_periods %>% 
      # map(st_apply, c(1,2), quantile, 0.99, na.rm = T)
      map(function(s_1_period) {
        
        s_1_period %>% 
          # st_apply(c(1,2), quantile, 0.99, na.rm = T)
          st_apply(c(1,2), function(q) {
            quantile(q, 0.99, na.rm = T)
          })
        
      })
    
    
    
    ## ALTERNATIVE
    dist_analysis <- 
      max_annual %>% 
      st_apply(c(1,2), function(x) {
        
        if (any(is.na(x))) {
          
          quantiles_95 <- c(NA,NA,NA)
          probs <- c(NA,NA,NA)
          
        } else {
          
          periods <- list(x[1:20], x[21:40], x[41:60])
          
          
          # ********************
          
          # # EMPIRICAL
          # 
          # # levels
          # quantiles_95 <-
          #   periods %>%
          #   map_dbl(quantile, 0.95)
          # 
          # # prob
          # probs <-
          #   periods %>%
          #   map_dbl(function(p) {
          #     1-ecdf(p)(quantiles_95[1])
          #   })
          
          
          # ********************
          
          # # GEV
          # 
          # # parameters
          # params <-
          #   periods %>%
          #   map(~evd::fgev(.x)$estimate %>% unname())
          # 
          # # levels
          # quantiles_95 <-
          #   params %>%
          #   map_dbl(~evd::qgev(0.95, loc = .x[1], scale = .x[2], shape = .x[3]))
          # 
          # # prob
          # probs <-
          #   params %>%
          #   map_dbl(~ 1 - evd::pgev(quantiles_95[1], loc = .x[1], scale = .x[2], shape = .x[3]))
          
          
          # ********************
          
          # GEV WITH L-MOMENTS
          
          # l-moments
          lmoms <- 
            periods %>% 
            map(~lmom::samlmu(.x))
          
          # parameters
          params <- 
            lmoms %>% 
            map(~lmom::pelgev(.x) %>% unname())
          
          # levels
          quantiles_95 <- 
            params %>% 
            map_dbl(~lmom::quagev(0.95, para = c(.x[1], .x[2], .x[3])))
          
          # prob
          probs <- 
            params %>% 
            map_dbl(~ 1 - lmom::cdfgev(quantiles_95[1], para = c(.x[1], .x[2], .x[3])))
          
        }
        
        
        # ***************
        
        results <- 
          c(quantiles_95, probs) %>% 
          set_names(c("q_1", "q_2", "q_3", "p_1", "p_2", "p_3"))
        
        return(results)
        
      },
      FUTURE = T,
      .fname = "stat")
    
    dist_analysis %>% 
      split("stat") %>%
      
      select(6) %>% mapview::mapview()
      
      
      
    
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
