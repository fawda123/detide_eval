library(lubridate)
library(tidyverse)
library(WtRegDO)
library(doParallel)
library(patchwork)
library(here)
library(janitor)
library(readxl)

# piermont data prep ------------------------------------------------------

met <- read.table('data/raw/piermont_met_all_2017.txt', sep = '\t', header = F) %>% 
  rename(
    DateTimeStamp = V1, 
    PAR = V2, 
    rain = V3, 
    ATemp = V4, 
    RH = V5,
    WD = V6,
    WSpd = V7,
    Gust = V8,
    blankT = V9, 
    blankRH = V10, 
    BP = V11
  ) %>% 
  select(DateTimeStamp, ATemp, BP, WSpd) %>% 
  mutate(
    DateTimeStamp = mdy_hms(DateTimeStamp, tz = 'America/Jamaica')
  ) %>% 
  arrange(DateTimeStamp)

wq <- read.table('data/raw/piermont_ysi_2017correct.txt', sep = '\t') %>% 
  rename(
    DateTimeStamp = V1,
    BP = V2,
    Temp = V3,
    Cond = V4,
    Sal = V5,
    depth = V6,
    pH = V7,
    pHvolt= V8,
    NTU = V9,
    blank = V10,
    O2sat = V11,
    time = V12,
    DO_obs = V13,
    battery = V14,
    depthcorr = V15
  ) %>% 
  mutate(
    DateTimeStamp = mdy_hm(DateTimeStamp, tz = 'America/Jamaica')
    ) %>% 
  select(DateTimeStamp, Temp, Sal, DO_obs, Tide = depthcorr) %>% 
  na.omit() %>% 
  arrange(DateTimeStamp) %>% 
  group_by(DateTimeStamp) %>% 
  summarise_all(mean, na.rm = T)

PIERMO2017 <- inner_join(wq, met, by = 'DateTimeStamp') %>% 
  arrange(DateTimeStamp) %>% 
  mutate(
    DateTimeStamp = format(DateTimeStamp, '%m/%d/%Y %H:%M')
  )

write.csv(PIERMO2017, 'data/raw/PIERMO2017.csv', row.names = F)

# process wtregdo ---------------------------------------------------------

wingrds <- crossing(
    tibble(flnm = c('APNERR2012', 'APNERR', 'APNERR2018', 'APNERR2020', 'HUDNERR', 'SAPDC', 'PIERMO', 'PIERMO2017'), 
           tz = c('America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica'), 
           lat = c(29.75, 29.75, 29.75, 29.75, 42.017, 31.39, 41.04, 41.04), 
           long = c(-85, -85, -85, -85, -73.915, -81.28, -73.90, -73.90)
           ),
    daywin = c(1, 3, 6, 9, 12),
    hrswin = c(1, 3, 6, 9, 12), 
    tidwin = c(0.2, 0.4, 0.6, 0.8, 1)
  ) 

# use this to filter out new files from the grid
wingrds <- wingrds %>%
  filter(flnm %in% 'APNERR2012')

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

# optimal window widths from results above --------------------------------

# devtools::load_all('../WtRegDO/')

mins <- list.files('data/', pattern = '^metdtd_') %>% 
  tibble(fl = .) %>% 
  mutate(
    fl = gsub('\\.RData$', '', fl),
    nm = gsub('^metdtd_(.*)_.*_.*_.*$', '\\1', fl)
  ) %>% 
  group_by(nm) %>% 
  nest() %>% 
  mutate(
    metobs = purrr::map(nm, function(nm){
      
      flobs <- list.files('data', paste0('metobs_', nm))[1]
      flobs <- gsub('\\.RData$', '', flobs)
      load(file = paste0('data/', flobs, '.RData'))
      get(flobs)
      
    }), 
    obseval = purrr::map(metobs, function(x) enframe(meteval(x, all = F)) %>% unnest(value))
  ) %>% 
  unnest(data) %>% 
  ungroup %>% 
  mutate(
    metdtd = purrr::pmap(list(metobs, fl), function(metobs, fl){
      
      cat(fl, '\n')
      load(file = paste0('data/', fl, '.RData'))
      dat <- get(fl)
      
      return(dat)
      
    }), 
    dtdeval = purrr::map(metdtd, function(x) enframe(meteval(x, all = F)) %>% unnest(value))#,
  ) %>% 
  mutate(
    minano = purrr::pmap(list(metobs, metdtd), objfun, vls = c('anomPg', 'anomRt')),
    minall = purrr::pmap(list(metobs, metdtd), objfun),
    fl = gsub('^metdtd_', '', fl)
  ) %>% 
  select(-metobs, -metdtd, -nm) %>%
  separate(fl, c('nm', 'dy', 'hr', 'td'), sep = '_') %>% 
  unnest(c('minano', 'minall')) %>%
  gather('obj', 'est', minano, minall) %>% 
  gather('evaltyp', 'eval', obseval, dtdeval) %>% 
  unnest('eval') %>% 
  spread(name, value) %>% 
  group_by(nm, obj, evaltyp)

soln <- mins %>% 
  filter(est == min(est)) %>% 
  mutate(
    evaltyp = factor(evaltyp, levels = c('obseval', 'dtdeval'), labels = c('Observed', 'Detided')),
    obj = factor(obj, levels = c('minall', 'minano'), labels = c('All', 'Anomalous'))
  ) %>% 
  arrange(nm, obj, evaltyp) %>% 
  filter(dy == dy[1] & hr == hr[1] & td == td[1]) %>% 
  ungroup() %>% 
  select(nm, dy, hr, td, Opt = obj, Val = est, Input = evaltyp, `Pg Mean` = meanPg, `Pg SD` = sdPg, `Pg Anom` = anomPg, `Rt Mean` = meanRt, `Rt SD` = sdRt, `Rt Anom` = anomRt)

save(soln, file = 'data/soln.RData', compress = 'xz')
