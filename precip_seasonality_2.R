
# New seasonality script

library(stars)
library(tidyverse)
library(furrr)
library(units)

plan(multisession)


f <- 
  str_glue("{working_dir}") %>% 
  list.files(full.names = T, pattern = "pr") 
#   working_dir %>% 
#   fs::dir_ls() %>% 
#   first()

# time vector
time_vector <- 
  f %>% 
  read_ncdf(proxy = T) %>% 
  st_get_dimension_values(3) %>% 
  suppressMessages() %>% 
  str_sub(end = 10)

# what calendar?
max_feb <- 
  time_vector[str_sub(time_vector, 6,7) == "02"] %>% # filter feb months
  str_sub(9,10) %>% # extract days
  as.numeric() %>% 
  max()

model_cal <- 
  case_when(max_feb == 30 ~ "360_day",
            max_feb == 29 ~ "gregorian",
            max_feb == 28 ~ "noleap")

# reformat time_vector as PCICt object
time_vector <- 
  PCICt::as.PCICt(time_vector, cal = model_cal)


# Import data into a list
# One time-period per element

l_s <- 
  c(2001, 2021, 2041) %>% 
  future_map(function(year_i) { # import in parallel
    
    time_first <- 
      which(time_vector == PCICt::as.PCICt(str_glue("{year_i}-01-01"), cal = model_cal))
    
    time_last <- 
      which.min(time_vector < PCICt::as.PCICt(str_glue("{year_i + 20}-01-01"), cal = model_cal))
    
    time_vector_sub <- 
      time_vector[time_first:(time_last-1)]
    
    
    # load data
    s <- 
      f %>% 
      read_ncdf(ncsub = cbind(start = c(1,  1, time_first),
                              count = c(NA, NA, time_last-time_first)), # import only necessary time steps
                proxy = F) %>% 
      suppressMessages() %>%
      st_set_dimensions(3, time_vector_sub)
    
    return(s)
  })


# Obtain seasonality
l_seas <- 
  l_s %>% 
  map(function(s) {
    
    tictoc::tic()
    
    doy <- 
      st_get_dimension_values(s, 3) %>% 
      yday() %>% 
      suppressWarnings()
    
    ss <- 
      st_apply(s, c(1,2), function(x) {
        
        if (any(is.na(x))) {
          
          res <- 
            c(onset = NA,
              cessation = NA)
          
        } else {
          
          m <-
            tibble(precip = x,
                   doy = doy)
          
          n <-
            m %>%
            group_by(doy) %>%
            summarize(daily_mean = mean(precip)) %>%
            mutate(anomaly = daily_mean - mean(m$precip),
                   cum_anomaly = cumsum(anomaly)) %>%
            summarize(onset = which.min(cum_anomaly),
                      cessation = which.max(tail(cum_anomaly, -onset)) + onset) # cessation always after onset
          
          res <-
            c(onset = n$onset,
              cessation = n$cessation)
          
        }
        
        return(res)
        
      },
      .fname = "seas",
      FUTURE = F) %>% 
      aperm(c(2,3,1))
    
    tictoc::toc()
    
    return(ss)
    
  })


# Convert to single stars obj
l_seas <- do.call(c, c(l_seas, along = "time_period")) %>% 
  split("seas")





