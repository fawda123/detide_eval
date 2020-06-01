library(lubridate)
library(tidyverse)
library(WtRegDO)
library(doParallel)
library(patchwork)
library(here)

wingrds <- crossing(
    tibble(flnm = c('APNERR', 'SAPDC'), tz = c('America/New_York', 'America/Jamaica'), lat = c(29.75, 31.39), long = c(-85, -81.28)),
    daywin = c(1, 3, 6, 9, 12),
    hrswin = c(1, 3, 6, 9, 12), 
    tidwin = c(0.2, 0.4, 0.6, 0.8, 1)
  ) 

ncores <- detectCores()  
registerDoParallel(cores = ncores - 1)
strt <- Sys.time()

foreach(i = 1:nrow(wingrds), .packages = c('WtRegDO', 'here', 'dplyr', 'lubridate')) %dopar% {
  
  sink('log.txt')
  cat(i, 'of', nrow(wingrds), '\n')
  print(Sys.time()-strt)
  sink()
  
  flnm <- wingrds[[i, 'flnm']]
  tz <- wingrds[[i, 'tz']]
  lat <- wingrds[[i, 'lat']]
  long <- wingrds[[i, 'long']]
  dy <- wingrds[[i, 'daywin']]
  hr <- wingrds[[i, 'hrswin']]
  td <- wingrds[[i, 'tidwin']]
  
  # data
  dat <- read.csv(here('data/raw/', paste0(flnm, '.csv'))) %>% 
    mutate(DateTimeStamp = mdy_hm(DateTimeStamp, tz = tz)) %>% 
    na.omit %>% 
    unique
    
  # weighted regression, optimal window widths for SAPDC from the paper
  wtreg_res <- wtreg(dat, parallel = F, wins = list(dy, hr, td), progress = F, 
                     tz = tz, lat = lat, long = long)
  
  # estimate ecosystem metabolism using observed DO time series
  metab_obs <- ecometab(wtreg_res, DO_var = 'DO_obs', tz = tz, 
                        lat = lat, long = long)
  
  # estimate ecosystem metabolism using detided DO time series
  metab_dtd <- ecometab(wtreg_res, DO_var = 'DO_nrm', tz = tz, 
                        lat = lat, long = long)
  
  # file base name
  flbs <- paste(flnm, dy, hr, td, sep = '_')
  
  # save do
  fldo <- paste0('DO_', flbs) 
  assign(fldo, wtreg_res)
  save(list = fldo, file = paste0('data/', fldo, '.RData'), compress = 'xz')
  
  # save metab, obs do
  flmetobs <- paste0('metobs_', flbs) 
  assign(flmetobs, metab_obs)
  save(list = flmetobs, file = paste0('data/', flmetobs, '.RData'), compress = 'xz')
  
  # save metab, dtd do
  flmetdtd <- paste0('metdtd_', flbs) 
  assign(flmetdtd, metab_dtd)
  save(list = flmetdtd, file = paste0('data/', flmetdtd, '.RData'), compress = 'xz')
  
}
