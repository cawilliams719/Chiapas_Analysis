# Precipitation Seasonality
# Examine precipitation seasonality and changes in seasonality


# Determine Annual vs Biannual Season
# ## Harmonic function to get amplitude
# harmonic <- geoTS::haRmonics(precip$pr, numFreq = 1, delta = 0.1)
# ## If ratio is > 1 then biannual if < 1 then annual regime
# ratio <- harmonic$amplitude[2]/harmonic$amplitude[1]


# Calculate Daily Climatological Mean (DOY Averages) per Pixel
subset <- precip[,,,1:11] # ouptus NAs with larger file but works with smaller subset

test <- subset %>% 
  st_apply(c(1:2), function(x) {
    
    pr <- as.vector(x)
    t <- st_get_dimension_values(subset, "time")
    
    df <- cbind.data.frame(pr,t) %>% 
      mutate(year = year(t), doy = yday(t)) 
    
    s <- df %>% group_by(doy) %>%
      summarize(pr_mean = mean(pr)) #%>%
      # mutate(anomaly = pr_mean-mean(df$pr),
      #        cum_anomaly = cumsum(anomaly)) %>%
      # summarize(onset = which.min(cum_anomaly),
      #           cessation = which.max(cum_anomaly))

    # df %>%
    #   mutate(seas = case_when(doy >= s$onset & doy < s$cessation ~ "wet",
    #                           TRUE ~ "dry"),
    #          year = ifelse(doy >= s$cessation, year+1, year)) %>%
    #   group_by(seas, year) %>%
    #   summarize(pr = sum(pr)) %>%
    #   mutate(onset = s$onset,
    #          cessation = s$cessation)
    
})%>% suppressMessages()

test$pr[3,2]

