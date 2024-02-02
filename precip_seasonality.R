# Precipitation Seasonality
# Examine precipitation seasonality and changes in seasonality
# https://github.com/carlosdobler/Project_mixteca_precip_analysis/blob/master/02_precip_analysis.Rmd

# Seasonal Precipitation 
## Testing original code first
chirps_stars <- precip[,60:90, 60:90,1:2]
# Seasonal precipitation for the whole region
seasonal_map <- map_df(seq_along(st_get_dimension_values(chirps_stars, "lon")), function(c){
  map_df(seq_along(st_get_dimension_values(chirps_stars, "lat")), function(r){
    
    m <- chirps_stars %>% 
      slice("lon", c) %>% 
      slice("lat", r) %>% 
      as_tibble() %>% 
      
      mutate(year = year(time)) %>% 
      group_by(year) %>% 
      mutate(doy = row_number()) %>% 
      ungroup()
    
    n <- m %>%
      group_by(doy) %>%
      summarize(daily_mean = mean(precip)) %>%
      mutate(anomaly = daily_mean - mean(m$precip),
             cum_anomaly = cumsum(anomaly)) %>%
      summarize(onset = which.min(cum_anomaly),
                cessation = which.max(cum_anomaly))

    # m %>% 
    #   mutate(seas = case_when(doy >= n$onset & doy < n$cessation ~ "wet",
    #                           TRUE ~ "dry"),
    #          year = ifelse(doy >= n$cessation, year+1, year)) %>% 
    #   group_by(seas, year) %>% 
    #   summarize(precip = sum(precip)) %>% 
    #   mutate(lon = st_get_dimension_values(chirps_stars, "lon")[c],
    #          lat = st_get_dimension_values(chirps_stars, "lat")[r],
    #          c = c,
    #          r = r,
    #          onset = n$onset,
    #          cessation = n$cessation)
    
  })
})

perc_seasonal <- seasonal_map %>%
  group_by(seas) %>%
  summarize(precip = mean(precip)) %>% 
  mutate(perc = precip/sum(precip)*100)

# Mean onset and cessation
seasonal_map %>% 
  ungroup() %>% 
  summarize(onset = as_date(mean(onset)),
            cessation = as_date(mean(cessation))) 