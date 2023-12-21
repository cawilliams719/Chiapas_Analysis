# Chiapas Precipitation Analysis
# Caroline Williams
# 12/7/2023
# Analysis of downscaled cmip6 precipitation data


# Import libraries
library(stars)
library(tidyverse)
library(forstringr)
library(ggthemes)
library(viridis)
library(viridisLite)
library(cowplot)
library(furrr)
library(evd)


################################################################################
# Copy files from cmip6_data to persistent drive
## Set source folder and destination folder paths
from_dir <- "/home/cwilliams/cmip6_directory/kgassert_TEMP/Chiapas/pr"
# Set working directory for data downloads and reading in files
working_dir <- "/home/cwilliams/persistent_disk/chiapas/precip/data/"

## Set file pattern to detect in from_dir
pattern <- "Chiapas_pr_ACCESS"

# ###############################################################################
# ## List and copy files from cmip6_data bucket to persistent drive
# from_dir %>%
#   list.files(full.names = TRUE, pattern = pattern) %>%
#   map(function(f) {
#     f_gs <- str_replace(f, "/home/cwilliams/cmip6_directory", "gs://cmip6_data")
# 
#     str_glue("gsutil cp {f_gs} {working_dir}") %>%
#       system()
#   })
# ###############################################################################

# Read and Pre-process Data
## Set Chiapas extent
# bb = st_bbox(c(xmin = -95.1,
#                ymin = 13.5,
#                xmax = -89.4,
#                ymax = 19.0))

# List of cmip6 model names
models <- working_dir %>% list.files(full.names = FALSE, pattern = pattern) %>% str_extract_part("_ssp245", before = TRUE)

# Read in precipitation data
precip <- working_dir %>% 
  list.files(full.names = TRUE, pattern = pattern) %>% 
  map(function(f) {
    read_ncdf(f, proxy = FALSE)
  }) 

# Set names of models in stars object list
precip <- setNames(precip, models)

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
perchangeim2 <- mapply(perchange, mean_total_annual_1, mean_total_annual_2, SIMPLIFY = FALSE)


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

################################################################################
## TESTING Quantile
## Shows tails
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


## TESTING Probability
# ECDF plot shows probability that precipitation is equal to or less than a certain number
fn <- ecdf(max_annual_1[[1]]$pr) # outputs sum probability for model 1  20 year period
plot(fn)
fn(200)
summary(fn)

ecdf_test <- max_annual_1[1:length(max_annual_1)] %>% # outputs ecdf 200 for each year for both models
  map(function(y) {
    apply(y$pr, 1, function(q) { 
      ecdf(q)(200)
    })
  })

# t3 <- max_annual_1[[1]] %>% st_apply(., 1:3, function(q) {
#   quantile(q, 0.99, na.rm = TRUE) = "200mm"
#   ecdf(q)
#   gev()
# })
# 
# ## TESTING EV Distribution
# pgev(max_annual_1[[1]]$pr)

################################################################################

# Visualize outputs - Testing and updating maps

# TESTING - 2 models x 2 difference maps
s1 <- subtracted1 %>% 
  map(function(f) {
    ggplot() +  
      geom_stars(data = f) + 
      scale_fill_viridis_c(option = "D") +
      coord_equal() +
      theme_map() +
      labs(subtitle = "(2021-2040) - (2001-2020)") 
  })
s2 <- subtracted2 %>% 
  map(function(f) {
    ggplot() +  
      geom_stars(data = f) + 
      scale_fill_viridis_c(option = "D") +
      coord_equal() +
      theme_map() +
      labs(subtitle = "(2041-2060) - (2001-2020)") 
  })

plot_grid(s1[[1]], s1[[2]],
          s2[[1]], s2[[2]],
          nrow = 2, ncol = 2, 
          labels = "AUTO") 


# ## Average Total Precipitation: 2001-2020, 2021-2040, 2041-2060
# ### Change time labels for visuals
# time_labels <- c("2001-2020", "2021-2040", "2041-2060")
# names(time_labels) <- c("2001-01-01", "2021-01-01", "2041-01-01")
# 
# ### Map Annual Total Precipitation
# ggplot() +  
#   geom_stars(data = c(mean_total_annual_1, mean_total_annual_2, mean_total_annual_3)) + 
#   facet_wrap("time", labeller = labeller("time" = time_labels)) +
#   scale_fill_viridis_c(option = "D") +
#   coord_equal() +
#   theme_map() +
#   ncol(3) +
#   nrow(1) +
#   # borders("world", xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
#   labs(title = "Average Total Annual Precipitation") +
#   theme(legend.position = "bottom") +
#   theme(legend.key.width = unit(1.5, "cm"))
# 
# ## Change in Precipitation: 2-1 & 3-1
# ### Change time labels for visuals
# time_labels_change <- c("(2021-2040) - (2001-2020)", "(2041-2060) - (2001-2020)")
# names(time_labels_change) <- c("2020-01-01", "2040-01-01")
# 
# ### Map Change in Precipitation
# ggplot() +  
#   geom_stars(data = c(subtracted1, subtracted2)) + 
#   facet_wrap("time", labeller = labeller("time" = time_labels_change)) +
#   scale_fill_viridis_c(option = "D") +
#   coord_equal() +
#   theme_map() +
#   ncol(2) +
#   nrow(1) +
#   # borders("world", xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
#   labs(title = "Change (mm) in Precipitation") +
#   theme(legend.position = "bottom") +
#   theme(legend.key.width = unit(1.5, "cm"))
# 
# ## Percent Change in Precipitation: 2-1 & 3-1
# ### Map Percent Change in Precipitation
# ggplot() +  
#   geom_stars(data = c(perchangeim1, perchangeim2)) + 
#   facet_wrap("time", labeller = labeller("time" = time_labels_change)) +
#   scale_fill_viridis_c(option = "D") +
#   coord_equal() +
#   theme_map() +
#   ncol(2) +
#   nrow(1) +
#   # borders("world", xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
#   labs(title = "Change (%) in Precipitation") +
#   theme(legend.position = "bottom") +
#   theme(legend.key.width = unit(1.5, "cm"))
