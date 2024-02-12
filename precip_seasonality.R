# Precipitation Seasonality
# Examine precipitation seasonality and changes in seasonality


# Determine Annual vs Biannual Season
# ## Harmonic function to get amplitude
# harmonic <- geoTS::haRmonics(precip$pr, numFreq = 1, delta = 0.1)
# ## If ratio is > 1 then biannual if < 1 then annual regime
# ratio <- harmonic$amplitude[2]/harmonic$amplitude[1]


# Calculate Daily Climatological Mean (DOY Averages) per Pixel
## Get year start and end positions for aggregating throughout script
year_positions <- pmatch(unique(year(st_get_dimension_values(precip, "time"))), 
                         year(st_get_dimension_values(precip, "time")))
endyear_positions <- year_positions-1
endyear_positions <- c(endyear_positions[-1], dim(precip)[3] %>% as.numeric())

year_range <- map2(year_positions, endyear_positions, 
                   function(y, e){
                     seq(y, e)
                   })

precip_periods <-
  list(year_range[1:20] %>% unlist(), 
       year_range[21:40] %>% unlist(),
       year_range[41:60] %>% unlist()) %>%
  map(~slice(precip, time, .x))


## Get mean pixel values aggregated by day of year (DOY) 1-366
precip_doy_mean <-  precip_periods %>% 
  map(function(p) {
    p %>% aggregate(function(x) {
      as.factor(yday(x))
    }, mean) %>% 
      aperm(c(2,3,1))
}) 

### rename dimension
dimnames(precip_doy_mean[[1]])[3] <- "doy" 
dimnames(precip_doy_mean[[2]])[3] <- "doy"  
dimnames(precip_doy_mean[[3]])[3] <- "doy"


## Get mean pixel value across time series aggregated by DOY
precip_mean <- precip_doy_mean %>% 
  map(function(p){
    p %>% st_apply(c(1:2), function(m) {
      mean(m)
    }) %>% 
    do.call("cbind",.)
}) %>% map(function(d) {
  d %>% replicate(366,.) %>%  # duplicate by length of precip_doy_mean
  as.vector()
})

## Add mean precipitation pixels as attribute
precip_anomaly <- map2(precip_doy_mean, precip_mean, function(p, m){
  p %>% mutate(pr_mean = m,
               anomaly = pr - pr_mean) # subtract pixels for each time slice by mean precip
  })

## Calculate cumulative sum
precip_cum <- precip_anomaly %>% 
  map(function(p) {
    p %>% select(anomaly) %>% 
      st_apply(c(1:2), function(t) {
        cumsum(t)
      }, keep = T) %>% 
      aperm(c(2,3,1)) %>%
      mutate(cum_anomaly = anomaly) %>% 
      select(cum_anomaly)
})

## Combine anomaly and cumulative anomaly stars objects
precip_comb <- map2(precip_anomaly, precip_cum, function(a, c) {
  list(a, c) %>% 
    do.call("c", .)
})


# Calculate Onset and Cessation
## Which min and max functions to handle NAs
whichMinNAs <- function(x){
  if(FALSE %in% is.na(x)){
    return(which.min(x))
  } else {
    return(NA)
  }
}

whichMaxNAs <- function(x){
  if(FALSE %in% is.na(x)){
    return(which.max(x))
  } else {
    return(NA)
  }
}


onset <- precip_comb %>% 
  map(function(p) {
    p %>% select(cum_anomaly) %>% 
      st_apply(c(1:2), function(i) {
        whichMinNAs(i)
      }, 
      .fname = "onset") %>% 
      do.call("cbind",.)
  }) %>% map(function(d) {
    d %>% replicate(366,.) %>% # duplicate by length of precip_doy_mean
      as.vector()
  })

cessation <- precip_comb %>% 
  map(function(p) {
    p %>% select(cum_anomaly) %>% 
      st_apply(c(1:2), function(i) {
        whichMaxNAs(i)
      }, 
      .fname = "cessation") %>% 
      do.call("cbind",.)
  }) %>% map(function(d) {
    d %>% replicate(366,.) %>% # duplicate by length of precip_doy_mean
      as.vector()
  })

precip_seas1 <- precip_comb[[1]] %>%
  mutate(onset = onset[[1]],
         cessation = cessation[[1]])
precip_seas2 <- precip_comb[[2]] %>%
  mutate(onset = onset[[2]],
         cessation = cessation[[2]])
precip_seas3 <- precip_comb[[3]] %>%
  mutate(onset = onset[[3]],
         cessation = cessation[[3]])

precip_seas <- list(precip_seas1, precip_seas2, precip_seas3)

# COME BACK TO: pmap not working so using above method for now
# precip_seas <- 
#   pmap(precip_comb, onset, cessation, 
#        function(m, o, c) {
#          m %>% mutate(onset = as.vector(o),
#                       cessation = as.vector(c))
#   })  


## Visuals
dfmap <- function(stars) {
  stars %>% 
    as.tbl_cube.stars() %>% 
    group_by(doy) %>% 
    summarise(mean_pr = mean(pr, na.rm = T),
              mean_anom = mean(anomaly, na.rm = T),
              mean_cum_anom = mean(cum_anomaly, na.rm = T),
              mean_onset = mean(onset, na.rm = T),
              mean_cessation = mean(cessation, na.rm = T)) %>% 
    as.data.frame()
}
  
df1 <- dfmap(precip_seas[[1]])
df2 <- dfmap(precip_seas[[2]])
df3 <- dfmap(precip_seas[[3]])

precip_map <- function(df) {
  ggplot(df, aes(x = doy, group =2)) +
    geom_line(aes(y = mean_pr), color = "red") +
    geom_line(aes(y = mean_anom), color = "blue") +
    labs(x = "DOY", y = "Precipitation")
}

p1 <- precip_map(df1)
p2 <- precip_map(df2)
p3 <- precip_map(df3)

cowplot::plot_grid(p1, p2, p3, ncol = 1)

seasonality <- function(df, o, c) {
  ggplot(df, aes(x = doy, group =2)) +
    geom_line(aes(y = mean_cum_anom), color = "orange") +
    geom_point(aes(x = mean_onset, y = o), color = "green") +
    geom_point(aes(x = mean_cessation, y = c), color = "purple") +
    labs(x = "DOY", y = "Summarized Variable of Interest")
}
s1 <- seasonality(df1, min(df1$mean_cum_anom), max(df1$mean_cum_anom))
s2 <- seasonality(df2, min(df2$mean_cum_anom), max(df2$mean_cum_anom))
s3 <- seasonality(df3, min(df3$mean_cum_anom), max(df3$mean_cum_anom))

cowplot::plot_grid(s1, s2, s3, ncol = 1)

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
doy_calc <- attr_calc(yday(st_get_dimension_values(precip, "time")), precip)
year_calc <- attr_calc(year(st_get_dimension_values(precip, "time")), precip)

## Fill onset and cessation matrix to length of time series
onset_att1 <- precip_seas1[,,,1]$onset %>%   
  replicate(length(st_get_dimension_values(precip_periods[[1]], "time")),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

onset_att2 <- precip_seas2[,,,1]$onset %>%   
  replicate(length(st_get_dimension_values(precip_periods[[2]], "time")),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

onset_att3 <- precip_seas3[,,,1]$onset %>%   
  replicate(length(st_get_dimension_values(precip_periods[[3]], "time")),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

cessation_att1 <- precip_seas1[,,,1]$cessation %>%   
  replicate(length(st_get_dimension_values(precip_periods[[1]], "time")),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

cessation_att2 <- precip_seas2[,,,1]$cessation %>%   
  replicate(length(st_get_dimension_values(precip_periods[[2]], "time")),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

cessation_att3 <- precip_seas3[,,,1]$cessation %>%   
  replicate(length(st_get_dimension_values(precip_periods[[3]], "time")),.) %>% # duplicate by length of precip_doy_mean
  as.vector()

onset_att <- list(onset_att1, onset_att2, onset_att3) %>% unlist()
cessation_att <- list(cessation_att1, cessation_att2, cessation_att3) %>% unlist()

## Add calculated attributes
precip_ts <- precip %>%
  mutate(doy = doy_calc,
         year = year_calc,
         onset = onset_att,
         cessation = cessation_att,
         seas = case_when(doy >= onset & doy < cessation ~ "wet",
                         TRUE ~ "dry"))#,
         # year = ifelse(doy >= cessation, year+1, year)) # Do we need this year part? It was in the original code

## Get mean pixel value across time series aggregated by DOY
precip_ts_mean <- precip_periods %>% map(m) %>%
  st_apply(c(1:2), function(m) {
    mean(m, na.rm = T)
  }) %>%
  do.call("cbind",.) %>%
  replicate(length(st_get_dimension_values(precip_periods[[1]], "time")),.) %>% # duplicate by length of main precip time series
  as.vector()


precip_ts_mean <- precip_periods %>% 
  map(function(p){
    p %>% st_apply(c(1:2), function(m) {
      mean(m)
    }) %>% 
      do.call("cbind",.)
  }) %>% map(function(d) {
    d %>% replicate(length(st_get_dimension_values(precip_periods[[3]], "time")),.) %>%  # duplicate by length of precip_doy_mean
      as.vector()
  })

## From here below, need to continue updating for each period
## Add mean precipitation pixels as attribute
precip_ts_anomaly <- precip_ts %>%
  mutate(pr_mean = precip_ts_mean,
         anomaly = units::drop_units(pr) - pr_mean) # subtract pixels for each time slice by mean precip


## Calculate cumulative sum
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

