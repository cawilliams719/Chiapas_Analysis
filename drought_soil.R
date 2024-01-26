
library(tidyverse)
library(stars)
library(furrr)
library(units)

plan(multisession)


source("file_paths.R")


# iterate over models
model <- "GFDL-CM4"


# download vars to working_dir
# download.R



# time vector
# f <- 
#   working_dir %>% 
#   fs::dir_ls() %>% 
#   first()
# 
# time_vector <- 
#   f %>% 
#   read_ncdf(proxy = T) %>% 
#   st_get_dimension_values(3) %>% 
#   suppressMessages() %>% 
#   str_sub(end = 10)
# 
# max_feb <- 
#   time_vector[str_sub(time_vector, 6,7) == "02"] %>% # filter feb months
#   str_sub(9,10) %>% # extract days
#   as.numeric() %>% 
#   max()
# 
# model_cal <- 
#   case_when(max_feb == 30 ~ "360_day",
#             max_feb == 29 ~ "gregorian",
#             max_feb == 28 ~ "noleap")
# 
# time_vector <- 
#   PCICt::as.PCICt(time_vector, cal = model_cal)
# 
# time_first <- 
#   which(time_vector == PCICt::as.PCICt("1971-01-01", cal = model_cal))
# 
# time_last <- 
#   which.min(time_vector < PCICt::as.PCICt("2061-01-01", cal = model_cal))
# 


# load data
ff <- 
  working_dir %>% 
  fs::dir_ls()

l_s <- 
  ff %>% 
  map(read_ncdf, proxy = F)

time_vector_m <- 
  seq(as_date("1991-01-01"), as_date("2060-12-01"), by = "1 month")

ss <- 
  do.call(c, c(l_s, along = "time")) %>% 
  st_set_dimensions(3, values = time_vector_m)




# calculate drought

s_drought_perc <- 
  ss %>% 
  st_apply(c(1,2), function(x) {
    
    if (any(is.na(x))) {
      
      rep(NA, length(x))
      
    } else {
      
      x %>% 
        zoo::rollsum(k = 3, na.pad = T, align = "right") %>%
        matrix(ncol = 12, byrow = T) %>% 
        apply(2, FUN = function(u) {
          
          f_ecdf <- ecdf(u[1:30]) # baseline
          
          f_ecdf(u)
          
        }) %>% 
        t() %>% 
        as.vector()
      
    }
  },
  FUTURE = T,
  .fname = "time") %>% 
  aperm(c(2,3,1)) %>% 
  st_set_dimensions(3, values = time_vector_m)


prob_extr_drought <- 
  s_drought_perc %>% 
  filter(year(time) >= 2001) %>% 
  aggregate(by = "20 years", FUN = function(x) mean(x <= 0.1)) %>% 
  aperm(c(2,3,1))












