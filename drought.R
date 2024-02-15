# Drought
# Examine drought


# plan(multisession)

# time vector
f <- 
  working_dir %>% 
  fs::dir_ls() %>% 
  first()

time_vector <- 
  f %>% 
  read_ncdf(proxy = T) %>% 
  st_get_dimension_values(3) %>% 
  suppressMessages() %>% 
  str_sub(end = 10)

max_feb <- 
  time_vector[str_sub(time_vector, 6,7) == "02"] %>% # filter feb months
  str_sub(9,10) %>% # extract days
  as.numeric() %>% 
  max()

model_cal <- 
  case_when(max_feb == 30 ~ "360_day",
            max_feb == 29 ~ "gregorian",
            max_feb == 28 ~ "noleap")

time_vector <- 
  PCICt::as.PCICt(time_vector, cal = model_cal)

time_first <- 
  which(time_vector == PCICt::as.PCICt("1971-01-01", cal = model_cal))

time_last <- 
  which.min(time_vector < PCICt::as.PCICt("2061-01-01", cal = model_cal))



# load data
ff <- 
  working_dir %>% 
  fs::dir_ls()

ff <- 
  ff %>%
  str_subset("_tas_", negate = T)

vars <- 
  ff %>% 
  str_split("_") %>% 
  map_chr(~.x[[5]])

l_s <- 
  ff %>% 
  set_names(vars) %>%
  imap(function(f, i) {
    
    print(str_glue("Importing {i}"))
    
    s <- 
      read_ncdf(f,
                ncsub = cbind(start = c(1,  1, time_first),
                              count = c(NA, NA, time_last-time_first)),
                proxy = F) %>% 
      suppressMessages() %>%
      setNames("v") %>% 
      st_set_dimensions(3, time_vector[time_first:(time_last-1)])
    
    
    if (i == "rsds") {
      s <- 
        s %>% 
        mutate(v = v %>% units::set_units(MJ/d/m^2))
    } else if (str_detect(i, "tas")) {
      s <- 
        s %>% 
        mutate(v = v %>% units::set_units(degC))
    }
    
    
    fun_agg <- case_when(i == "pr" ~ "sum",
                         i %in% c("rsds", "tasmax", "tasmin") ~ "mean")
    
    a <-
      aggregate(s, by = "1 month", FUN = eval(parse(text = fun_agg))) %>%
      aperm(c(2,3,1))
    
    a <- a %>% setNames(i)
    
    return(a)
    
  })


write_rds(l_s, str_glue("{working_dir}/l_s.rds"))
read_rds(str_glue("{working_dir}/l_s.rds")) -> l_s


# calculate pet

ss <- 
  do.call(c, c(l_s[which(names(l_s) != "pr")], 
               along = "vars"))


s_pet <- 
  ss %>% 
  st_apply(c(1,2), function(d) {
    
    if (any(is.na(d[,1]))) {
      rep(NA, nrow(d))
    } else {
      
      SPEI::hargreaves(Tmin = d[,3],
                       Tmax = d[,2],
                       Ra = d[,1],
                       verbose = F)
    }
  }, 
  .fname = "time") %>% 
  aperm(c(2,3,1))

st_dimensions(s_pet) <- st_dimensions(l_s$pr)


# calculate drought
s_wb <- l_s$pr - s_pet

time_vector_m <- st_get_dimension_values(s_wb, 3) %>% as_date()

int_window <- 3

s_drought_perc <- 
  s_wb %>% 
  st_apply(c(1,2), function(x) {
    
    if (any(is.na(x))) {
      
      rep(NA, length(x))
      
    } else {
      
      x %>% 
        zoo::rollsum(k = int_window, na.pad = T, align = "right") %>%
        matrix(ncol = 12, byrow = T) %>% 
        apply(2, FUN = function(u) {
          
          f_ecdf <- ecdf(u[1:50]) # baseline
          
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

prob_extr_drought_perc <- 
  s_drought_perc %>% 
  filter(year(time) >= 2001) %>% 
  aggregate(by = "20 years", FUN = function(x) mean(x <= 0.1)) %>% 
  aperm(c(2,3,1))


s_drought_perc_nr <- 
  s_wb %>% 
  st_apply(c(1,2), function(x) {
    
    if (any(is.na(x))) {
      
      rep(NA, length(x))
      
    } else {
      
      x %>% 
        # zoo::rollsum(k = int_window, na.pad = T, align = "right") %>%
        matrix(ncol = 12, byrow = T) %>% 
        apply(2, FUN = function(u) {
          
          f_ecdf <- ecdf(u[1:50]) # baseline
          
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

prob_extr_drought_perc_nr <- 
  s_drought_perc_nr %>% 
  filter(year(time) >= 2001) %>% 
  aggregate(by = "20 years", FUN = function(x) mean(x <= 0.1)) %>% 
  aperm(c(2,3,1))



# ignore

{
  s_drought_gamma <- 
    s_wb %>% 
    st_apply(c(1,2), function(x) {
      
      if (any(is.na(x))) {
        
        rep(NA, length(x))
        
      } else {
        
        x %>% 
          zoo::rollsum(k = int_window, na.pad = T, align = "right") %>% 
          matrix(ncol = 12, byrow = T) %>% 
          apply(2, FUN = function(u) {
            
            u_bl <- u[1:30]
            
            u_min <- min(u_bl, na.rm = T) - 0.00001
            
            u1 <- u_bl - u_min
            
            lmoms <-
              lmom::samlmu(u1)
            
            params <-
              lmom::pelgam(lmoms)
            
            lmom::cdfgam(u - u_min, params)
            
          }) %>% 
          t() %>% 
          as.vector()
      }
    },
    .fname = "time") %>% 
    aperm(c(2,3,1)) %>% 
    st_set_dimensions(3, values = time_vector_m)
  
  
  # how well gamma fits obs data?
  s_gamma_gof <- 
    s_wb %>% 
    st_apply(c(1,2), function(x) {
      
      if (any(is.na(x))) {
        
        rep(NA, 12)
        
      } else {
        
        x %>% 
          zoo::rollsum(k = int_window, na.pad = T, align = "right") %>% 
          matrix(ncol = 12, byrow = T) %>% 
          apply(2, FUN = function(u) {
            
            u_bl <- u[1:30]
            
            u_min <- min(u_bl, na.rm = T) - 0.00001
            
            u1 <- u_bl - u_min
            
            lmoms <-
              lmom::samlmu(u1)
            
            params <-
              lmom::pelgam(lmoms)
            
            na <- which(is.na(u_bl))
            
            if (length(na) < 1) na <- 9999
            
            cor.test(
              # method = "spearman",
              sort(u_bl[-na]),
              lmom::quagam(seq(0.01, 0.99, length.out = length(u_bl)),
                           para = params) %>%
                .[-na]
              #{. + u_min}
            )$estimate
            
          }) %>% 
          as.vector()
      }
    },
    .fname = "month") %>% 
    aperm(c(2,3,1))
  }

int_window <- 12

s_drought_spei <- 
  s_wb %>% 
  st_apply(c(1,2), function(x) {
    
    if (any(is.na(x))) {
      
      rep(NA, length(x))
      
    } else {
      
      SPEI::spei(x, scale = int_window, ref.end = c(50,12), verbose = F)$fitted
      
    }
  },
  FUTURE = T,
  .fname = "time") %>% 
  aperm(c(2,3,1)) %>% 
  st_set_dimensions(3, values = time_vector_m)

prob_extr_drought_spei <- 
  s_drought_spei %>% 
  filter(year(time) >= 2001) %>% 
  aggregate(by = "20 years", FUN = function(x) mean(x <= -1.25)) %>% 
  aperm(c(2,3,1))







# foo <- units::set_units(1, W/m^2) # == J/s/m^2 ( rsds = srad * dayl / 86400 )
# foo %>% units::set_units(MJ/d/m^2)

