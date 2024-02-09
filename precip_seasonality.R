# Precipitation Seasonality
# Examine precipitation seasonality and changes in seasonality

# Seasonal Precipitation 
## Testing original code first
chirps_stars <- precip[,60:70, 60:80,]


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
  }, keep = T) %>% 
  aperm(c(2,3,1)) %>%
  mutate(cum_anomaly = anomaly) %>% 
  select(cum_anomaly)


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
  .fname = "cessation") %>% 
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
                         TRUE ~ "dry"))#,
         # year = ifelse(doy >= cessation, year+1, year)) # Do we need this year part? It was in the original code

## Get mean pixel value across time series aggregated by DOY
precip_ts_mean <- chirps_stars %>% 
  st_apply(c(1:2), function(m) {
    mean(m, na.rm = T)
  }) %>% 
  do.call("cbind",.) %>% 
  replicate(length(st_get_dimension_values(chirps_stars, 'time')),.) %>% # duplicate by length of main precip time series
  as.vector()

## Add mean precipitation pixels as attribute
precip_ts_anomaly <- precip_ts %>% 
  mutate(pr_mean = precip_ts_mean, 
         anomaly = units::drop_units(pr) - pr_mean) # subtract pixels for each time slice by mean precip



## Calculate cumulative sum
year_positions <- pmatch(unique(year(st_get_dimension_values(chirps_stars, "time"))), 
                         year(st_get_dimension_values(chirps_stars, "time")))
endyear_positions <- year_positions-1
endyear_positions <- c(endyear_positions[-1], dim(chirps_stars)[3] %>% as.numeric())

year_range <- map2(year_positions, endyear_positions, 
     function(y, e){
       seq(y, e)
})

years <- 
  year_range %>% 
  map(~slice(precip_ts_anomaly[8,,,], time, .x))

precip_ts_cum <- years %>% 
  map(function(t){
    st_apply(t, c(1:2), function(t) {
      cumsum(t)
    }, keep = T) %>% 
      aperm(c(2,3,1)) 
  }) %>% do.call("c",.) %>% 
  mutate(cum_anomaly = anomaly) %>% 
  select(cum_anomaly)




## Combine anomaly and cumulative anomaly stars objects
precip_ts_comb <- list(precip_ts_anomaly, precip_ts_cum) %>% 
  do.call("c", .) %>% 
  mutate(cum_anomaly_onset = ifelse(onset - 60 <= doy & onset + 60 >= doy, cum_anomaly, NA),
         cum_anomaly_cessation = ifelse(cessation - 60 <= doy & cessation + 60 >= doy, cum_anomaly, NA))

# Calculate Onset and Cessation
onset_ts <- precip_ts_comb %>% 
  select(cum_anomaly_onset) %>% 
  st_apply(1:2, function (i) {
    onset <- which.min(i)
  }, 
  .fname = "onset_ts", keep = T) %>% 
  do.call("cbind",.) %>% 
  replicate(length(st_get_dimension_values(precip_ts_comb, 'time')),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

cessation_ts <- precip_ts_comb %>% 
  select(cum_anomaly_cessation) %>% 
  st_apply(1:2, function (i) {
    cessation <- which.max(i)
  }, 
  .fname = "cessation_ts", keep = T) %>% 
  do.call("cbind",.) %>% 
  replicate(length(st_get_dimension_values(precip_ts_comb, 'time')),.) %>% # duplicate by length of precip_doy_mean
  as.vector()


precip_ts_seas <- precip_ts_comb %>% 
  mutate(onset_ts = onset_ts,
         cessation_ts = cessation_ts)

onset_stats <- precip_ts_seas[12,,,] %>% 
  st_apply(1:2, function(a){
    mk_onset  <-  trend::mk.test(a) %>% .$p.value
    # sen_onset  <-  trend::sens.slope(onset) %>% .$estimates %>% as.numeric()
  })

cessation_stats <- precip_ts_seas[13,,,] %>% 
  st_apply(1:2, function(a){
     mk_cessation = trend::mk.test(cessation) %>% .$p.value
     # sen_cessation = trend::sens.slope(cessation) %>% .$estimates %>% as.numeric()
  })

