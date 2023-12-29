# Chiapas Precipitation Analysis
# Caroline Williams
# 12/29/2023
# Visualization of downscaled cmip6 precipitation data


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
# Read in precipitation analysis results
source("file_paths.R")
load("data/precip_analysis_results.RData")

# Visualize outputs
## View and select model for all map precipitation outputs
names(precip) # Examine list of models in console
plot_model <- names(precip[1]) # Select index number corresponding to model of interest or replace with model string


## Average Total Precipitation: 2001-2020, 2021-2040, 2041-2060
### Subtitle time  labels for visuals
time_labels <- c("2001-2020", "2021-2040", "2041-2060")

### Get min/max values for scale limits to get a uniform legend/scale across comparison maps
#### get min
get_min1 <- map(c(mean_total_annual_1[plot_model], 
                  mean_total_annual_2[plot_model], 
                  mean_total_annual_3[plot_model]), 
                function(p) {
                  min(p$pr, na.rm = TRUE)
                }) %>% as.data.frame() %>% min()

#### get_max
get_max1 <- map(c(mean_total_annual_1[plot_model], 
                  mean_total_annual_2[plot_model], 
                  mean_total_annual_3[plot_model]), 
                function(p) {
                  max(p$pr, na.rm = TRUE)
                }) %>% as.data.frame() %>% max()

### Map Annual Total Precipitation
#### Map with legend first to call later separately
mta_legend <- map(c(mean_total_annual_1[plot_model], 
                    mean_total_annual_2[plot_model], 
                    mean_total_annual_3[plot_model]), 
                  function(p) {
                    ggplot() + 
                      geom_stars(data = p) +
                      scale_fill_viridis_c(option = "D", na.value = "transparent", limits = c(get_min1, get_max1)) +
                      labs(fill = "Precipitation") +
                      theme(legend.position = "bottom") +
                      theme(legend.key.width = unit(1.5, 'cm'))
                  })

#### Map without legend that will later have single legend added
mta_nolegend <- map(c(mean_total_annual_1[plot_model], 
                      mean_total_annual_2[plot_model], 
                      mean_total_annual_3[plot_model]), 
                    function(p) {
                      ggplot() + 
                        borders("world", fill = "gray", colour = "black") +
                        geom_stars(data = p) +
                        scale_fill_viridis_c(option = "D", na.value = "transparent", limits = c(get_min1, get_max1)) +
                        coord_equal(xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
                        theme_light() +
                        theme(legend.position = "none")
                    })

#### Add time ranges as labels to no legend map
p_mta_nolegend <- plot_grid(plotlist = mta_nolegend, 
                            nrow = 1, 
                            ncol = 3, 
                            labels = time_labels,
                            label_size = 12) 

#### Add legend and modify position
p_mta_legend <- plot_grid(p_mta_nolegend + theme(legend.position = "none"), 
                          get_legend(mta_legend[[1]] + theme(legend.box.margin = margin(20, 100, 80, 60))),
                          rel_heights = c(0.8, 0.2), 
                          nrow = 2, 
                          ncol = 1)

#### Create new ggplot with overarching title and subtitle
title_gg1 <- ggplot() + 
  labs(title = "Average Total Annual Precipitation", subtitle = str_glue("Model: {plot_model}")) 

#### Final plot for Average Total Annual precipitation
mta_map <- plot_grid(title_gg1, p_mta_legend, ncol = 1, rel_heights = c(0.15, 1))
mta_map


## Change in Precipitation: 2-1 & 3-1
### Subtitle time  labels for visuals
time_labels_change <- c("(2021-2040) - (2001-2020)", "(2041-2060) - (2001-2020)")

### Get min/max values for scale limits to get a uniform legend/scale across comparison maps
#### get min
get_min2 <- map(c(subtracted1[plot_model], 
                  subtracted2[plot_model]), 
                function(p) {
                  min(p$pr, na.rm = TRUE)
                }) %>% as.data.frame() %>% min()

#### get_max
get_max2 <- map(c(subtracted1[plot_model], 
                  subtracted2[plot_model]), 
                function(p) {
                  max(p$pr, na.rm = TRUE)
                }) %>% as.data.frame() %>% max()

### Map Change in Precipitation
#### Map with legend first to call later separately
sub_legend <- map(c(subtracted1[plot_model], 
                    subtracted2[plot_model]), 
                  function(p) {
                    ggplot() + 
                      geom_stars(data = p) +
                      scale_fill_viridis_c(option = "D", na.value = "transparent", limits = c(get_min2, get_max2)) +
                      coord_equal() +
                      labs(fill = "Precipitation") +
                      theme(legend.position = "bottom") +
                      theme(legend.key.width = unit(1.5, 'cm'))
                  })

#### Map without legend that will later have single legend added
sub_nolegend <- map(c(subtracted1[plot_model], 
                      subtracted2[plot_model]), 
                    function(p) {
                      ggplot() + 
                        borders("world", fill = "gray", colour = "black") +
                        geom_stars(data = p) +
                        scale_fill_viridis_c(option = "D", na.value = "transparent", limits = c(get_min2, get_max2)) +
                        coord_equal(xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
                        theme_light() +
                        theme(legend.position = "none")
                    })

#### Add time ranges as labels to no legend map
p_sub_nolegend <- plot_grid(plotlist = sub_nolegend, 
                            nrow = 1, 
                            ncol = 2, 
                            labels = time_labels_change,
                            label_size = 12) 

#### Add legend and modify position
p_sub_legend <- plot_grid(p_sub_nolegend + theme(legend.position = "none"), 
                          get_legend(sub_legend[[1]] + theme(legend.box.margin = margin(0, 0, 0, 0))),
                          rel_heights = c(0.9, 0.1), 
                          nrow = 2, 
                          ncol = 1)

#### Create new ggplot with overarching title and subtitle
title_gg2 <- ggplot() + 
  labs(title = "Change (mm) in Precipitation", subtitle = str_glue("Model: {plot_model}")) 

#### Final plot for Change in Precipitation
sub_map <- plot_grid(title_gg2, p_sub_legend, ncol = 1, rel_heights = c(0.1, 0.75))
sub_map


## Percent Change in Precipitation: 2-1 & 3-1
### Get min/max values for scale limits to get a uniform legend/scale across comparison maps
#### get min
get_min3 <- map(c(perchangeim1[plot_model], 
                  perchangeim2[plot_model]), 
                function(p) {
                  min(p$pr, na.rm = TRUE)
                }) %>% as.data.frame() %>% min()

#### get_max
get_max3 <- map(c(perchangeim1[plot_model], 
                  perchangeim2[plot_model]), 
                function(p) {
                  max(p$pr, na.rm = TRUE)
                }) %>% as.data.frame() %>% max()

### Map Percent Change in Precipitation
#### Map with legend first to call later separately
per_legend <- map(c(perchangeim1[plot_model], 
                    perchangeim2[plot_model]), 
                  function(p) {
                    ggplot() + 
                      geom_stars(data = p) +
                      scale_fill_viridis_c(option = "D", na.value = "transparent", limits = c(get_min3, get_max3)) +
                      coord_equal() +
                      labs(fill = "Precipitation") +
                      theme(legend.position = "bottom") +
                      theme(legend.key.width = unit(1.5, 'cm'))
                  })

#### Map without legend that will later have single legend added
per_nolegend <- map(c(perchangeim1[plot_model], 
                      perchangeim2[plot_model]), 
                    function(p) {
                      ggplot() + 
                        borders("world", fill = "gray", colour = "black") +
                        geom_stars(data = p) +
                        scale_fill_viridis_c(option = "D", na.value = "transparent", limits = c(get_min3, get_max3)) +
                        coord_equal(xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
                        theme_light() +
                        theme(legend.position = "none")
                    })

#### Add time ranges as labels to no legend map
p_per_nolegend <- plot_grid(plotlist = per_nolegend, 
                            nrow = 1, 
                            ncol = 2, 
                            labels = time_labels_change,
                            label_size = 12) 

#### Add legend and modify position
p_per_legend <- plot_grid(p_per_nolegend + theme(legend.position = "none"), 
                          get_legend(per_legend[[1]] + theme(legend.box.margin = margin(0, 0, 0, 0))),
                          rel_heights = c(0.9, 0.1), 
                          nrow = 2, 
                          ncol = 1)

#### Create new ggplot with overarching title and subtitle
title_gg3 <- ggplot() + 
  labs(title = "Change (%) in Precipitation", subtitle = str_glue("Model: {plot_model}")) 

#### Final plot for Percent Change in Precipitation
per_map <- plot_grid(title_gg3, p_per_legend, ncol = 1, rel_heights = c(0.1, 0.75))
per_map


## Quantile Maximum Precipitation: 2001-2020, 2021-2040, 2041-2060
### Get min/max values for scale limits to get a uniform legend/scale across comparison maps
#### get min
get_min4 <- map(c(quantile_99_1[plot_model], 
                  quantile_99_2[plot_model], 
                  quantile_99_3[plot_model]), 
                function(p) {
                  min(p$pr, na.rm = TRUE)
                }) %>% as.data.frame() %>% min()

#### get_max
get_max4 <- map(c(quantile_99_1[plot_model], 
                  quantile_99_2[plot_model], 
                  quantile_99_3[plot_model]), 
                function(p) {
                  max(p$pr, na.rm = TRUE)
                }) %>% as.data.frame() %>% max()

### Map Quantile Maximum Extreme Precipitation
#### Map with legend first to call later separately
qnt_legend <- map(c(quantile_99_1[plot_model], 
                    quantile_99_2[plot_model], 
                    quantile_99_3[plot_model]), 
                  function(p) {
                    ggplot() + 
                      geom_stars(data = p) +
                      scale_fill_viridis_c(option = "D", na.value = "transparent" , limits = c(get_min4, get_max4)) +
                      coord_equal() +
                      labs(fill = "Precipitation") +
                      theme(legend.position = "bottom") +
                      theme(legend.key.width = unit(1.5, 'cm'))
                  })

#### Map without legend that will later have single legend added
qnt_nolegend <- map(c(quantile_99_1[plot_model], 
                      quantile_99_2[plot_model], 
                      quantile_99_3[plot_model]), 
                    function(p) {
                      ggplot() + 
                        borders("world", fill = "gray", colour = "black") +
                        geom_stars(data = p) +
                        scale_fill_viridis_c(option = "D", na.value = "transparent" , limits = c(get_min4, get_max4)) +
                        coord_equal(xlim = c(-95.1, -89.4), ylim = c(13.5, 19.0)) +
                        theme_light() +
                        theme(legend.position = "none")
                    })

#### Add time ranges as labels to no legend map
p_qnt_nolegend <- plot_grid(plotlist = qnt_nolegend, 
                            nrow = 1, 
                            ncol = 3, 
                            labels = time_labels,
                            label_size = 12) 

#### Add legend and modify position
p_qnt_legend <- plot_grid(p_qnt_nolegend + theme(legend.position = "none"), 
                          get_legend(qnt_legend[[1]] + theme(legend.box.margin = margin(20, 100, 80, 60))),
                          rel_heights = c(0.8, 0.2), 
                          nrow = 2, 
                          ncol = 1)

#### Create new ggplot with overarching title and subtitle
title_gg4 <- ggplot() + 
  labs(title = "Quantile Maximum Annual Precipitation", subtitle = str_glue("Model: {plot_model}")) 

#### Final plot for Average Total Annual precipitation
qnt_map <- plot_grid(title_gg4, p_qnt_legend, ncol = 1, rel_heights = c(0.15, 1))
qnt_map
