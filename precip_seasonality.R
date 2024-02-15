# Precipitation Seasonality
# Examine precipitation seasonality and changes in seasonality


# Determine Annual vs Biannual Season
# ## Harmonic function to get amplitude
# harmonic <- geoTS::haRmonics(precip$pr, numFreq = 1, delta = 0.1)
# ## If ratio is > 1 then biannual if < 1 then annual regime
# ratio <- harmonic$amplitude[2]/harmonic$amplitude[1]


# Calculate Daily Climatological Mean (DOY Averages) per Pixel
subset <- precip[,80:81,70:71,] # outputs NAs with larger file but works with smaller subset


# Calculate Daily Climatological Mean (DOY Averages) per Pixel
## Get year start and end positions for aggregating throughout script
year_positions <- pmatch(unique(year(st_get_dimension_values(subset, "time"))), 
                         year(st_get_dimension_values(subset, "time")))
endyear_positions <- year_positions-1
endyear_positions <- c(endyear_positions[-1], dim(subset)[3] %>% as.numeric())

year_range <- map2(year_positions, endyear_positions, 
                   function(y, e){
                     seq(y, e)
                   })

precip_periods <-
  list(year_range[1:20] %>% unlist(), 
       year_range[21:40] %>% unlist(),
       year_range[41:60] %>% unlist()) %>%
  map(~slice(subset, time, .x))


system.time(testrf_periods <- precip_periods %>% map(function(p) {
  p %>% st_apply(c(1:2), function(x) {
      

    pr <- as.vector(x)
    t <- st_get_dimension_values(p, "time")
    
    df <- cbind.data.frame(pr,t) %>% 
      mutate(year = year(t), doy = yday(t)) 
    
    s <- df %>% group_by(doy) %>%
      reframe(pr_mean = mean(pr), by = "doy") %>%
      mutate(anomaly = pr_mean-mean(df$pr),
             cum_anomaly = cumsum(anomaly)) %>%
      reframe(onset = which.min(cum_anomaly),
              cessation = which.max(cum_anomaly), by = "doy")
    
    timing <- df %>% 
      mutate(o= s$onset %>% replicate(length(st_get_dimension_values(p, "time")), .), 
             c= s$cessation %>% replicate(length(st_get_dimension_values(p, "time")), .),
             seas = case_when(doy >= o & doy < c ~ "wet",
                              TRUE ~ "dry"),
             year = ifelse(doy >= c, year+1, year)) %>%
      group_by(seas, year) %>%
      summarize(pr = sum(pr)) %>% 
      mutate(onset = s$onset %>% replicate(length(seas), .),
             cessation = s$cessation %>% replicate(length(seas), .))
    
    # perc <- timing %>% 
    #   mutate(seas_perc = seas) %>% 
    #   group_by(seas_perc) %>%
    #   summarize(pr = mean(pr)) %>%
    #   mutate(pr = pr/sum(pr)*100)
    # 
    # dates <- timing %>%
    #   ungroup() %>%
    #   summarize(onset = as_date(mean(onset)),
    #             cessation = as_date(mean(cessation)))
    # 
    # return(c(timing, perc, dates))
    
    s <- df %>%
      mutate(pr_mean = mean(pr),
             anomaly = pr - pr_mean) %>% 
      group_by(year) %>%
      mutate(cum_anomaly = cumsum(anomaly),
             doy = row_number()) %>%
      mutate(cum_anomaly_onset = ifelse(as.numeric(timing$onset) - 60 <= doy & as.numeric(timing$onset) + 60 >= doy, cum_anomaly, NA),
             cum_anomaly_cessation = ifelse(as.numeric(timing$cessation) - 60 <= doy & as.numeric(timing$cessation) + 60 >= doy, cum_anomaly, NA)) %>%
      summarise(onset = which.min(cum_anomaly_onset),
                cessation = which.max(cum_anomaly_cessation))
    
  }) %>% suppressMessages()
}))

testrf_periods[[1]]$pr[1,2]

## Output format - stars object with lat/lon dimensions. Time is show as each pixel has a data frame where rows = time
## Not sure how to transform the result into something that is easy to work with.