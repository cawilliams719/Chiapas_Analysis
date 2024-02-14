# Precipitation Seasonality
# Examine precipitation seasonality and changes in seasonality


# Determine Annual vs Biannual Season
# ## Harmonic function to get amplitude
# harmonic <- geoTS::haRmonics(precip$pr, numFreq = 1, delta = 0.1)
# ## If ratio is > 1 then biannual if < 1 then annual regime
# ratio <- harmonic$amplitude[2]/harmonic$amplitude[1]


# Calculate Daily Climatological Mean (DOY Averages) per Pixel
subset <- precip#[,60:70,60:70,1:720] # outputs NAs with larger file but works with smaller subset
# 
# test <- subset %>% 
#   st_apply(c(1:2), function(x) {
#     
#     pr <- as.vector(x)
#     t <- st_get_dimension_values(subset, "time")
#     
#     df <- cbind.data.frame(pr,t) %>% 
#       mutate(year = year(t), doy = yday(t)) 
#     
#     s <- df %>% group_by(doy) %>%
#       summarize(pr_mean = mean(pr)) %>%
#       mutate(anomaly = pr_mean-mean(df$pr),
#              cum_anomaly = cumsum(anomaly)) %>%
#       summarize(onset = which.min(cum_anomaly),
#                 cessation = which.max(cum_anomaly))
# 
#     df %>%
#       mutate(seas = case_when(doy >= s$onset & doy < s$cessation ~ "wet",
#                               TRUE ~ "dry"),
#              year = ifelse(doy >= s$cessation, year+1, year)) %>%
#       group_by(seas, year) %>%
#       summarize(pr = sum(pr)) %>%
#       mutate(onset = s$onset,
#              cessation = s$cessation)
#     
# }) %>% suppressMessages()
# 
# test$pr[2,3]

system.time(testrf <- subset %>% 
  st_apply(c(1:2), function(x) {
    
    pr <- as.vector(x)
    t <- st_get_dimension_values(subset, "time")
    
    df <- cbind.data.frame(pr,t) %>% 
      mutate(year = year(t), doy = yday(t)) 
    
    s <- df %>% group_by(doy) %>%
      reframe(pr_mean = mean(pr), by = "doy") %>%
      mutate(anomaly = pr_mean-mean(df$pr),
             cum_anomaly = cumsum(anomaly)) %>%
      reframe(onset = which.min(cum_anomaly),
                cessation = which.max(cum_anomaly), by = "doy")
    
    timing <- df %>% 
      mutate(o= s$onset %>% replicate(length(st_get_dimension_values(subset, "time")), .), 
             c= s$cessation %>% replicate(length(st_get_dimension_values(subset, "time")), .),
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

  }) %>% suppressMessages())

testrf$pr[2,3]
# testrf$pr[60,70]

# testrf %>% 
#   st_apply(c(1:2), function(x){
#     head(x$pr) 
#     # %>% group_by(seas) %>%
#     #     summarize(pr = mean(pr)) %>%
#     #     mutate(pr = pr/sum(pr)*100)
#   })
# 
# testrf %>% split("pr")
