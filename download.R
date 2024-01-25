

# "/mnt/bucket_cmip6/kgassert_TEMP/Chiapas/pr/" %>% 
#   fs::dir_ls(regexp = "ssp245") %>% 
#   future_walk(function(f) {
#     
#     f_gs <- str_replace(f, "/mnt/bucket_cmip6", "gs://cmip6_data")
#     
#     "gsutil cp {f_gs} /mnt/pers_disk/chiapas_cmip6_basd_pr" %>%
#       str_glue() %>% 
#       system(ignore.stdout = T, ignore.stderr = T)
#     
#   })


"/mnt/bucket_cmip6/kgassert_TEMP/Chiapas/" %>% 
  fs::dir_ls() %>% 
  map(fs::dir_ls) %>% 
  map(str_subset, "GFDL-CM4") %>% 
  map(str_subset, "ssp245") %>% 
  walk(function(f) {
    
    f_gs <- str_replace(f, "/mnt/bucket_cmip6", "gs://cmip6_data")
    
    str_glue("gsutil cp {f_gs} /mnt/pers_disk/chiapas_cmip6_basd") %>% 
      system()
    
  })
    