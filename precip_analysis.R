# Chiapas Precipitation Analysis
# Caroline Williams
# 12/7/2023
# Analysis of downscaled cmip6 precipitation data


# Import libraries
library(stars)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(viridis)
library(viridisLite)

################################################################################
# # Copy files from cmip6_data to persistent drive
# ## Set source folder and destination folder paths
# from_dir <- "/home/cwilliams/cmip6_directory/ISIMIP_BASD_data/daily/precipitation/ssp245/"
# to_dir <- "/home/cwilliams/persistent_disk/chiapas/precip/data/"
# 
# ## Set file pattern to detect in from_dir
# pattern <- "pr_ACCESS-CM2_ssp245_basd_0.5deg_20[0-5][0-9].nc"
# 
# ## Save list of file names to copy
# the_files <- list.files(path = from_dir, pattern = pattern)
# 
# ## Copy files to to_dir
# file.copy(from = file.path(from_dir, the_files),
#                     to = file.path(to_dir, the_files))
################################################################################

# Read and Pre-process Data
## Set Chiapas extent
bb = st_bbox(c(xmin = -95.1,
               ymin = 13.5,
               xmax = -89.4,
               ymax = 19.0))
  
## Set working directory to read in files
working_dir <- "/home/cwilliams/persistent_disk/chiapas/precip/data/"

## Set pattern for list.files
model_name_00s <- "pr_ACCESS-CM2_ssp245_basd_0.5deg_200[0-9].nc"
model_name_10s <- "pr_ACCESS-CM2_ssp245_basd_0.5deg_201[0-9].nc"
model_name_20s <- "pr_ACCESS-CM2_ssp245_basd_0.5deg_202[0-9].nc"
model_name_30s <- "pr_ACCESS-CM2_ssp245_basd_0.5deg_203[0-9].nc"
model_name_40s <- "pr_ACCESS-CM2_ssp245_basd_0.5deg_204[0-9].nc"
model_name_50s <- "pr_ACCESS-CM2_ssp245_basd_0.5deg_205[0-9].nc"
model_name_60 <- "pr_ACCESS-CM2_ssp245_basd_0.5deg_2060.nc"

# List files based on directory and file nomenclature pattern
precip_list_00s = list.files(path = working_dir, pattern = model_name_00s, full.names = TRUE)
precip_list_10s = list.files(path = working_dir, pattern = model_name_10s, full.names = TRUE)
precip_list_20s = list.files(path = working_dir, pattern = model_name_20s, full.names = TRUE)
precip_list_30s = list.files(path = working_dir, pattern = model_name_30s, full.names = TRUE)
precip_list_40s = list.files(path = working_dir, pattern = model_name_40s, full.names = TRUE)
precip_list_50s = list.files(path = working_dir, pattern = model_name_50s, full.names = TRUE)
precip_list_60 = list.files(path = working_dir, pattern = model_name_60, full.names = TRUE)

## Read in NetCDF cmip6 data & crop to Chiapas extent
chiapas_2000s <- precip_list_00s %>% read_stars(.) %>% st_crop(.,  bb, crop = TRUE)
chiapas_2010s <- precip_list_10s %>% read_stars(.) %>% st_crop(.,  bb, crop = TRUE)
chiapas_2020s <- precip_list_20s %>% read_stars(.) %>% st_crop(.,  bb, crop = TRUE)
chiapas_2030s <- precip_list_30s %>% read_stars(.) %>% st_crop(.,  bb, crop = TRUE)
chiapas_2040s <- precip_list_40s %>% read_stars(.) %>% st_crop(.,  bb, crop = TRUE)
chiapas_2050s <- precip_list_50s %>% read_stars(.) %>% st_crop(.,  bb, crop = TRUE)
chiapas_2060 <- precip_list_60 %>% read_stars(.) %>% st_crop(.,  bb, crop = TRUE)

## Combine into one stars object
chiapas <- c(chiapas_2000s, chiapas_2010s, chiapas_2020s, chiapas_2030s,
             chiapas_2040s, chiapas_2050s, chiapas_2060)

## Remove initial files
rm(chiapas_2000s, chiapas_2010s, chiapas_2020s, chiapas_2030s, chiapas_2040s, 
   chiapas_2050s, chiapas_2060)

# Calculate Annual Precipitation (Total)
# ## Monthly Mean
# monthly_mean <- chiapas %>% aggregate(., by = "1 month", FUN = mean) %>% 
#   aperm(., c("x", "y", "time")) 
# 
# ## Annual 
# ## Sum of Monthly Mean
# total_annual_mean <- monthly_mean %>% aggregate(., by = "1 year", FUN = sum) %>% 
#   aperm(., c("x", "y", "time")) 

## Sum from Daily
total_annual <- chiapas %>% aggregate(., by = "1 year", FUN = sum) %>% 
  aperm(., c("x", "y", "time")) 

# Annual Means per Time Period
## 2001-2020
mean_total_annual_1 <- total_annual %>% 
  filter(time < "2021-01-01") %>% 
  aggregate(., by = "20 years", FUN = mean) %>% 
  aperm(., c("x", "y", "time")) 

## 2021-2040
mean_total_annual_2 <- total_annual %>% 
  filter(time >"2020-01-01", time < "2041-01-01") %>% 
  aggregate(., by = "20 years", FUN = mean) %>% 
  aperm(., c("x", "y", "time")) 

## 2041-2060
mean_total_annual_3 <- total_annual %>% 
  filter(time > "2040-01-01") %>% 
  aggregate(., by = "20 years", FUN = mean) %>% 
  aperm(., c("x", "y", "time")) 

# Change in Precipitation
## Subtract Function
subtract <- function(x) (x[[2]]-x[[1]]) 

## Change Period 2-1: (2021-2040) - (2001-2020)
mean_total_annual_1_tchange <- st_set_dimensions(mean_total_annual_1, "time", as.POSIXct("2020-01-01")) # Reset time for the first period to 2020 so when aggregating with time period 2 it is only 1 year difference
change1 <- c(mean_total_annual_1_tchange, mean_total_annual_2) # combine two time periods to one star object
st_get_dimension_values(change1, "time") # check dates


subtracted1 <- change1 %>% aggregate(., by = "2 years", FUN = subtract) %>%
  aperm(., c("x", "y", "time"))

## Change Period 3-1: (2041-2060) - (2001-2020)
mean_total_annual_1_tchange2 <- st_set_dimensions(mean_total_annual_1, "time", as.POSIXct("2040-01-01")) # Reset time again
change2 <- c(mean_total_annual_1_tchange2, mean_total_annual_3) # combine two time periods to one star object
st_get_dimension_values(change2, "time") # check dates


subtracted2 <- change2 %>% aggregate(., by = "2 years", FUN = subtract) %>%
  aperm(., c("x", "y", "time"))

# Percent Change in Precipitation
## Percent Change Function
perchange <- function(x) (((x[[2]]-x[[1]])/x[[1]])*100) # double check

## Percent Change Period 2-1: (2021-2040) - (2001-2020)
perchangeim1 <- change1 %>% aggregate(., by = "2 years", FUN = perchange) %>%
  aperm(., c("x", "y", "time"))

## Percent Change Period 3-1: (2041-2060) - (2001-2020)
perchangeim2 <- change2 %>% aggregate(., by = "2 years", FUN = perchange) %>%
  aperm(., c("x", "y", "time"))


# Visualize outputs

## Average Total Precipitation: 2001-2020, 2021-2040, 2041-2060
### Change time labels for visuals
time_labels <- c("2001-2020", "2021-2040", "2041-2060")
names(time_labels) <- c("2001-01-01", "2021-01-01", "2041-01-01")

### Map Annual Total Precipitation
ggplot() +  
  geom_stars(data = c(mean_total_annual_1, mean_total_annual_2, mean_total_annual_3)) + 
  facet_wrap("time", labeller = labeller("time" = time_labels)) +
  scale_fill_viridis_c(option = "D") +
  coord_equal() +
  theme_map() +
  ncol(3) +
  nrow(1) +
  # borders("world", xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
  labs(title = "Average Annual Precipitation") +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(1.5, "cm"))

### Change in Precipitation: 2-1 & 3-1
### Change time labels for visuals
time_labels_change <- c("(2021-2040) - (2001-2020)", "(2041-2060) - (2001-2020)")
names(time_labels_change) <- c("2020-01-01", "2040-01-01")

### Map Change in Precipitation
ggplot() +  
  geom_stars(data = c(subtracted1, subtracted2)) + 
  facet_wrap("time", labeller = labeller("time" = time_labels_change)) +
  scale_fill_viridis_c(option = "D") +
  coord_equal() +
  theme_map() +
  ncol(2) +
  nrow(1) +
  # borders("world", xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
  labs(title = "Change (mm) in Precipitation") +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(1.5, "cm"))

## Percent Change in Precipitation: 2-1 & 3-1
### Map Change in Precipitation
ggplot() +  
  geom_stars(data = c(perchangeim1, perchangeim2)) + 
  facet_wrap("time", labeller = labeller("time" = time_labels_change)) +
  scale_fill_viridis_c(option = "D") +
  coord_equal() +
  theme_map() +
  ncol(2) +
  nrow(1) +
  # borders("world", xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
  labs(title = "Change (%) in Precipitation") +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(1.5, "cm"))
