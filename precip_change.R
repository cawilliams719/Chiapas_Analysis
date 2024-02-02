# Change in Precipitation & Extremes
# Examine change in precipitation between future scenarios and baseline and extremes


# Read and Process Precipitation Data
## Select file based on model name
precip_file <- 
  working_dir %>% 
  list.files(full.names = T, pattern = model) # model defined in run_model

## Extract time as vector
time_vector <- 
  precip_file %>% 
  read_ncdf(proxy = T) %>% 
  st_get_dimension_values(3) %>% 
  suppressMessages() %>% 
  str_sub(end = 10)

## Get the maximum date used in February based on time vector
max_feb <- 
  time_vector[str_sub(time_vector, 6,7) == "02"] %>% # filter feb months
  str_sub(9,10) %>% # extract days
  as.numeric() %>% 
  max()

## Define type of calendar used based on the maximum February date
model_cal <- 
  case_when(max_feb == 30 ~ "360_day",
            max_feb == 29 ~ "gregorian",
            max_feb == 28 ~ "noleap")

## Format time vector as PCICt
time_vector <- 
  PCICt::as.PCICt(time_vector, cal = model_cal)

## Define first time step
time_first <- 
  which(time_vector == PCICt::as.PCICt("2001-01-01", cal = model_cal))

## Define final time step 
time_last <- 
  which.min(time_vector < PCICt::as.PCICt("2061-01-01", cal = model_cal))

## Read in precipitation subletting by time
precip <- 
  precip_file %>% 
  read_ncdf(ncsub = cbind(start = c(1,  1, time_first),
                          count = c(NA, NA, time_last-time_first)),
            proxy = F) %>% 
  suppressMessages()

# Calculate Annual Precipitation (Total)
## Sum from Daily
total_annual <-
  precip %>%
  aggregate(by = "1 year", FUN = sum) %>%
  aperm(c(2,3,1))

# Annual Precipitation (Total) Means per Time Period
# mean_total_annual
mean_total_annual <-
  total_annual %>%
  aggregate(by = "20 years", FUN = mean) %>%
  aperm(c(2,3,1))

# Change in Precipitation
## Define percent change function
perchange <- 
  function(i, j) {(((j-i)/i)*100)}

## Change Period 2-1: (2021-2040) - (2001-2020)
subtracted1 <- slice(mean_total_annual, time, 2) - slice(mean_total_annual, time, 1)
## Change Period 3-1: (2041-2060) - (2001-2020)
subtracted2 <- slice(mean_total_annual, time, 3) - slice(mean_total_annual, time, 1)

# Percent Change in Precipitation
## Percent Change Period 2-1: (2021-2040) - (2001-2020)
perchangeim1 <- perchange(slice(mean_total_annual, time, 1), slice(mean_total_annual, time, 2))
## Percent Change Period 3-1: (2041-2060) - (2001-2020)
perchangeim2 <- perchange(slice(mean_total_annual, time, 1), slice(mean_total_annual, time, 3))

## Define combine subtracted and percent result function
change_merge <-
  function(s, p) {
  sub <- s %>% rename(sub = pr)
  perc <- p %>% rename(perc= pr)
  l <- list(sub, perc) %>% 
    do.call("c", .)
  return(l)
}

## Combine subtracted and percent change results by period
change1 <- change_merge(subtracted1, perchangeim1)
change2 <- change_merge(subtracted2, perchangeim2)

# Precipitation Extremes
# Calculate Annual Precipitation (Maximum) per Time Period
max_annual <-
  precip %>%
  aggregate(by = "1 year", FUN = max) %>%
  aperm(c(2,3,1))

max_annual_periods <-
  list(1:20, 21:40, 41:60) %>%
  map(~slice(max_annual, time, .x))

## Quantile
## Shows tails of maximum precipitation per pixel aggregated by time
quantile_95 <-
  max_annual_periods %>%
  # map(st_apply, c(1,2), quantile, 0.99, na.rm = T)
  map(function(s_1_period) {
    s_1_period %>%
      # st_apply(c(1,2), quantile, 0.99, na.rm = T)
      st_apply(c(1,2), function(q) {
        quantile(q, 0.95, na.rm = T)
      })

  })

## ALTERNATIVE
dist_analysis <-
  max_annual %>%
  st_apply(c(1,2), function(x) {
    if (any(is.na(x))) {

      quantiles_95 <- c(NA,NA,NA)
      probs <- c(NA,NA,NA)

    } else {

      periods <- list(x[1:20], x[21:40], x[41:60])

      # ********************

      # # EMPIRICAL
      # 
      # # levels
      # quantiles_95 <-
      #   periods %>%
      #   map_dbl(quantile, 0.95)
      # 
      # # prob
      # probs <-
      #   periods %>%
      #   map_dbl(function(p) {
      #     1-ecdf(p)(quantiles_95[1])
      #   })


      # ********************

      # # GEV
      # 
      # # parameters
      # params <-
      #   periods %>%
      #   map(~evd::fgev(.x)$estimate %>% unname())
      # 
      # # levels
      # quantiles_95 <-
      #   params %>%
      #   map_dbl(~evd::qgev(0.95, loc = .x[1], scale = .x[2], shape = .x[3]))
      # 
      # # prob
      # probs <-
      #   params %>%
      #   map_dbl(~ 1 - evd::pgev(quantiles_95[1], loc = .x[1], scale = .x[2], shape = .x[3]))


      # ********************

      # GEV WITH L-MOMENTS

      # l-moments
      lmoms <-
        periods %>%
        map(~lmom::samlmu(.x))

      # parameters
      params <-
        lmoms %>%
        map(~lmom::pelgev(.x) %>% unname())

      # levels
      quantiles_95 <-
        params %>%
        map_dbl(~lmom::quagev(0.95, para = c(.x[1], .x[2], .x[3])))

      # prob
      probs <-
        params %>%
        map_dbl(~ 1 - lmom::cdfgev(quantiles_95[1], para = c(.x[1], .x[2], .x[3])))

    }

    # ***************
    
    results <-
      c(quantiles_95, probs) %>%
      set_names(c("q_1", "q_2", "q_3", "p_1", "p_2", "p_3"))

    return(results)
  },
  FUTURE = T,
  .fname = "stat")

dist_analysis <- dist_analysis %>%
  split("stat")

