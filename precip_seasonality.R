# Precipitation Seasonality
# Examine precipitation seasonality and changes in seasonality
# https://github.com/carlosdobler/Project_mixteca_precip_analysis/blob/master/02_precip_analysis.Rmd

# Seasonal Precipitation 
## Testing original code first
chirps_stars <- precip[,60:70, 60:80,]
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
      summarize(daily_mean = mean(pr)) %>%
      mutate(anomaly = daily_mean - mean(m$pr),
             cum_anomaly = cumsum(anomaly)) %>%
      summarize(onset = which.min(cum_anomaly),
                cessation = which.max(cum_anomaly))

    m %>%
      mutate(seas = case_when(doy >= n$onset & doy < n$cessation ~ "wet",
                              TRUE ~ "dry"),
             year = ifelse(doy >= n$cessation, year+1, year)) %>%
      group_by(seas, year) %>%
      summarize(pr = sum(pr)) %>%
      mutate(lon = st_get_dimension_values(chirps_stars, "lon")[c],
             lat = st_get_dimension_values(chirps_stars, "lat")[r],
             c = c,
             r = r,
             onset = n$onset,
             cessation = n$cessation)

  })
}) %>% suppressMessages()


perc_seasonal <- seasonal_map %>%
  group_by(seas) %>%
  summarize(pr = mean(pr)) %>% 
  mutate(pr = pr/sum(pr)*100)

# Mean onset and cessation --> produces weird output
seasonal_map %>% 
  ungroup() %>% 
  summarize(onset = as_date(mean(onset)),
            cessation = as_date(mean(cessation))) 

####
m <- seasonal_map %>% 
  group_by(lon, lat, c, r, onset, cessation) %>% 
  dplyr::summarize() %>%
  dplyr::ungroup() %>% 
  
  pmap_df(function(onset, cessation, c, r, lon, lat){
    
    chirps_stars %>% 
      slice("lon", c) %>% 
      slice("lat", r) %>%
      
      as_tibble() %>%
      mutate(daily_mean = mean(pr),
             anomaly = pr - daily_mean,
             year = year(time)) %>%
      group_by(year) %>%
      mutate(cum_anomaly = cumsum(anomaly),
             doy = row_number()) %>%

      mutate(cum_anomaly_onset = ifelse(onset - 60 <= doy & onset + 60 >= doy, cum_anomaly, NA),
             cum_anomaly_cessation = ifelse(cessation - 60 <= doy & cessation + 60 >= doy, cum_anomaly, NA)) %>%

      summarise(onset = which.min(cum_anomaly_onset),
                cessation = which.max(cum_anomaly_cessation)) %>%

      dplyr::summarise(mk_onset = trend::mk.test(onset) %>% .$p.value,
                sen_onset = trend::sens.slope(onset) %>% .$estimates %>% as.numeric(),
                mk_cessation = trend::mk.test(cessation) %>% .$p.value,
                sen_cessation = trend::sens.slope(cessation) %>% .$estimates %>% as.numeric()) %>%

      mutate(lon = lon,
             lat = lat)
    
  })

mk <- m %>% 
  select(mk_onset, mk_cessation, lon, lat) %>% 
  gather(1:2, key = seas, val = mk) %>% 
  mutate(mk = ifelse(mk < 0.1, 1, NA),
         seas = str_sub(seas, start = 4),
         seas = factor(seas, levels = c("onset", "cessation"))) %>% 
  filter(!is.na(mk))

m %>% 
  select(sen_onset, sen_cessation, lon, lat) %>% 
  gather(1:2, key = seas, val = sen) %>% 
  mutate(seas = str_sub(seas, start = 5),
         seas = factor(seas, levels = c("onset", "cessation"))) %>% 
  
  ggplot() +
  geom_raster(aes(x = lon, y = lat, fill = sen)) +
  geom_point(data = mk, aes(x = lon, y = lat), shape = 4, alpha = 0.5) +
  scale_fill_distiller(palette = "PuOr",
                       limits = c(-1.1, 1.1),
                       direction = -1,
                       breaks = c(-1, 0, 1),
                       labels = c("-1 (later)", "0 (no change)", "1 (sooner)"),
                       name = "days/year") +
  # geom_point(data = locs, aes(x = lon, y = lat)) +
  # geom_text(data = locs, aes(x = lon, y = lat, label = loc), vjust = 1, nudge_y = -0.02, nudge_x = 0.05, size = 3) +
  facet_grid(~seas) +
  coord_quickmap()


# Determine Annual vs Biannual Season
harmonic <- geoTS::haRmonics(chirps_stars$pr, numFreq = 1, delta = 0.1)
## If ratio is > 1 then biannual if < 1 then annual regime
ratio <- harmonic$amplitude[2]/harmonic$amplitude[1]


# Calculate Daily Climatological Mean (DOY Averages) per Pixel
## Get mean pixel values aggregated by day of year (DOY) 1-366
precip_doy_mean <-  chirps_stars %>% 
  aggregate(function(x) {
    days = as.factor(yday(x))
    }, mean, na.rm = T) %>% 
  aperm(c(2,3,1))
dimnames(precip_doy_mean)[[3]] <- "doy" # rename dimension

## Get mean pixel value across time series aggregated by DOY
precip_mean <- precip_doy_mean %>% 
  st_apply(c(1:2), function(m) {
    mean(m, na.rm = T)
}) %>% 
  do.call("cbind",.) %>% 
  replicate(366,.) %>% # duplicate by length of precip_doy_mean
  as.vector()

## Add mean precipitation pixels as attribute
precip_anomaly <- precip_doy_mean %>% 
  mutate(pr_mean = precip_mean, 
         anomaly = pr - pr_mean) # subtract pixels for each time slice by mean precip

## Calculate cumulative sum
precip_cum <- precip_anomaly[3,,,] %>% 
  st_apply(c(1:2), function(t) {
    cumsum(t)
  }) %>% 
  aperm(c(2,3,1)) %>% 
  mutate(cum_anomaly = anomaly) %>% 
  select(cum_anomaly)
dimnames(precip_cum)[[3]] <- "doy" # rename dimension
precip_cum <- st_set_dimensions(precip_cum, "doy", # reset doy values
                                values = st_get_dimension_values(
                                  precip_anomaly, 'doy', where = "start"))

## Combine anomaly and cumulative anomaly stars objects
precip_comb <- list(precip_anomaly, precip_cum) %>% 
  do.call("c", .)

# Calculate Onset and Cessation
onset <- precip_comb %>% 
  select(cum_anomaly) %>% 
  st_apply(1:2, function (i) {
    onset <- which.min(i)
  }, 
  .fname = "onset") %>% 
  do.call("cbind",.) %>% 
  replicate(length(st_get_dimension_values(precip_comb, 'doy')),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

cessation <- precip_comb %>% 
  select(cum_anomaly) %>% 
  st_apply(1:2, function (i) {
    cessation <- which.max(i)
  }, 
  .fname = "cessation" ) %>% 
  do.call("cbind",.) %>% 
  replicate(length(st_get_dimension_values(precip_comb, 'doy')),.) %>% # duplicate by length of precip_doy_mean
  as.vector()


precip_seas <- precip_comb %>% 
  mutate(onset = onset,
         cessation = cessation)

## Visuals
df <- precip_seas %>% 
  as.tbl_cube.stars() %>% 
  group_by(doy) %>% 
  summarise(mean_pr = mean(pr),
            mean_anom = mean(anomaly),
            mean_cum_anom = mean(cum_anomaly, na.rm = T),
            mean_onset = mean(onset),
            mean_cessation = mean(cessation)) %>% 
  as.data.frame()

precip_map <- ggplot(df, aes(x = doy, group =2)) +
  geom_line(aes(y = mean_pr), color = "red") +
  geom_line(aes(y = mean_anom), color = "blue") +
  labs(x = "DOY", y = "Summarized Variable of Interest")
precip_map

seasonality <- ggplot(df, aes(x = doy, group =2)) +
  geom_line(aes(y = mean_cum_anom), color = "orange") +
  geom_point(aes(x = mean_onset, y = -460), color = "green") +
  geom_point(aes(x = mean_cessation, y = 97), color = "purple") +
  labs(x = "DOY", y = "Summarized Variable of Interest")
seasonality

## TESTING SECOND PART involving year
### Function to calculate attributes DOY, Year, etc.
attr_calc <- function(s, r) {
  s %>% 
    map(function(i) {
      matrix(i, 
             nrow = nrow(st_get_dimension_values(r, "lat")), 
             ncol = nrow(st_get_dimension_values(r, "lon")))
    }) %>% unlist()
}

## Calculate DOY and year attributes
doy_calc <- attr_calc(yday(st_get_dimension_values(chirps_stars, "time")), chirps_stars)
year_calc <- attr_calc(year(st_get_dimension_values(chirps_stars, "time")), chirps_stars)

## Fill onset and cessation matrix to length of time series
onset_att <- precip_seas[,,,1]$onset %>%   
  replicate(length(st_get_dimension_values(chirps_stars, 'time')),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

cessation_att <- precip_seas[,,,1]$cessation %>%   
  replicate(length(st_get_dimension_values(chirps_stars, 'time')),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

## Add calculated attribtues
precip_ts <- chirps_stars %>% 
  mutate(doy = doy_calc, 
         year = year_calc,
         onset = onset_att,
         cessation = cessation_att,
         seas = case_when(doy >= onset & doy < cessation ~ "wet",
                         TRUE ~ "dry"),
         year = ifelse(doy >= cessation, year+1, year))



