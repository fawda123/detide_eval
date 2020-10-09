library(tidyverse)
library(lubridate)
library(WtRegDO)
library(foreach)
library(doParallel)

# flnm <- 'PIERMO2017'
# tz <- 'America/Jamaica'
# lat <- 41.04
# long <- -73.90
# 
# # data
# dat <- read.csv(paste0('data/raw/', flnm, '.csv')) %>%
#   mutate(DateTimeStamp = mdy_hm(DateTimeStamp, tz = tz)) %>%
#   na.omit %>%
#   unique
# 
# ncores <- detectCores()
# cl <- makeCluster(ncores - 1)
# registerDoParallel(cl)
# 
# out <- winopt(dat, tz = tz, lat = lat, long = long, wins = list(6, 6, 0.5), 
#               control = list(factr = 1e+08, parscale = c(100, 200, 100), maxit = 100), parallel = T, vls = c('anomPg', 'anomRt'))
# 

optgrds <- tibble(
  flnm = c('APNERR', 'APNERR2018', 'APNERR2020', 'HUDNERR', 'SAPDC', 'PIERMO', 'PIERMO2017'), 
  tz = c('America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica'), 
  lat = c(29.75, 29.75, 29.75, 42.017, 31.39, 41.04, 41.04), 
  long = c(-85, -85, -85, -73.915, -81.28, -73.90, -73.90)
  )

ncores <- detectCores()  
registerDoParallel(cores = ncores - 1)
strt <- Sys.time()

opts <- vector('list', length = nrow(optgrds))
names(opts) <- optgrds$flnm

for(i in 1:nrow(optgrds)){
  
  flnm <- optgrds[[i, 'flnm']]
  tz <- optgrds[[i, 'tz']]
  lat <- optgrds[[i, 'lat']]
  long <- optgrds[[i, 'long']]

  cat(flnm, '\n\n')
  
  # data
  dat <- read.csv(paste0('data/raw/', flnm, '.csv')) %>% 
    mutate(DateTimeStamp = mdy_hm(DateTimeStamp, tz = tz)) %>% 
    na.omit %>% 
    unique
  
  out <- winopt(dat, tz = tz, lat = lat, long = long, wins = list(6, 6, 0.5), 
                control = list(factr = 1e+07, parscale = c(100, 200, 100), maxit = 100), parallel = T, vls = c('anomPg', 'anomRt'))
  
  opts[[flnm]] <- out
  
}

# test when optim is done -------------------------------------------------

optgrds <- tibble(
  flnm = c('APNERR', 'APNERR2018', 'APNERR2020', 'HUDNERR', 'SAPDC', 'PIERMO', 'PIERMO2017'), 
  tz = c('America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica', 'America/Jamaica'), 
  lat = c(29.75, 29.75, 29.75, 42.017, 31.39, 41.04, 41.04), 
  long = c(-85, -85, -85, -73.915, -81.28, -73.90, -73.90),
  wins = list(
    APNERR = list(4.45114875821087, 8.90461094963134, 0.876024796925915), 
    APNERR2018 = list(6, 6, 0.5), 
    APNERR2020 = list(6, 6, 0.5), 
    HUDNERR = list(12, 0.1, 1), 
    SAPDC = list(6, 6, 1), 
    PIERMO = list(12, 12, 1), 
    PIERMO2017 = list(12, 7.02594646999893, 1)
    )
)

cmps <- vector('list', length = nrow(optgrds))
names(cmps) <- optgrds$flnm

for(i in 1:nrow(optgrds)){
  
  flnm <- optgrds[[i, 'flnm']]
  tz <- optgrds[[i, 'tz']]
  lat <- optgrds[[i, 'lat']]
  long <- optgrds[[i, 'long']]
  wins <- optgrds[[i, 'wins']][[1]]
  
  cat(flnm, '\n\n')
  
  # data
  dat <- read.csv(paste0('data/raw/', flnm, '.csv')) %>% 
    mutate(DateTimeStamp = mdy_hm(DateTimeStamp, tz = tz)) %>% 
    na.omit %>% 
    unique

  wtreg_res <- wtreg(dat, wins = wins, progress = T, parallel = F,
                     tz = tz, lat = lat, long = long)


  metab_obs <- ecometab(wtreg_res, DO_var = 'DO_obs', tz = tz,
                        lat = lat, long = long)
  metab_dtd <- ecometab(wtreg_res, DO_var = 'DO_nrm', tz = tz,
                        lat = lat, long = long)

  obj <- objfun(metab_obs, metab_dtd, vl = c('anomPg', 'anomRt'))
  
  evalobs <- meteval(metab_obs, all = F)
  evaldtd <- meteval(metab_dtd, all = F)
 
  out <- data.frame(rbind(evalobs, evaldtd))
  out$obj <- obj
  
  cmps[[flnm]] <- out
  
}
