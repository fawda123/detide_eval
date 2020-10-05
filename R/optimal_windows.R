library(tidyverse)
library(lubridate)
library(WtRegDO)
library(foreach)
library(doParallel)
# devtools::load_all('../WtRegDO/')

ncores <- detectCores()  
cl <- makeCluster(ncores)
registerDoParallel(cl)

flnm <- 'APNERR2020'
tz <- 'America/Jamaica'
lat <- 29.75
long <- -85
dy <- 12
hr <- 6
td <- 0.5

# data
dat <- read.csv(paste0('data/raw/', flnm, '.csv')) %>%
  mutate(DateTimeStamp = mdy_hm(DateTimeStamp, tz = tz)) %>%
  na.omit %>%
  unique

# weighted regression, optimal window widths for SAPDC from the paper
wtreg_res <- wtreg(dat, wins = list(12, 6, 0.5), progress = T, parallel = T,
                   tz = tz, lat = lat, long = long)

# estimate ecosystem metabolism using observed DO time series
metab_obs <- ecometab(wtreg_res, DO_var = 'DO_obs', tz = tz,
                      lat = lat, long = long)

wins_in <- list(12, 6, 0.5)

control <- list(factr = 1e7, parscale = c(50, 100, 50))
lower <- c(0.1, 0.1, 0.1)
upper <- c(12, 12, 1)
 
strt <- Sys.time()

fun_in <- function(wins_in){

  txt <- unlist(wins_in)
  print(Sys.time() - strt)
  cat(txt, '\n')
  
  wtreg_res <- wtreg(dat, parallel = T, wins = wins_in, progress = T,
                     tz = tz, lat = lat, long = long)
  
  metab_dtd <- ecometab(wtreg_res, DO_var = 'DO_nrm', tz = tz,
                        lat = lat, long = long)
  
  out <- objfun(metab_obs, metab_dtd, vls = c('anomPg', 'anomRt'))
  
  cat(out, '\n\n')
  
  return(out)
  
}

out <- optim(
  wins_in, 
  fun_in, 
  method = 'L-BFGS-B', 
  lower = lower, 
  upper = upper, 
  control = control
)

stopCluster(cl)
  

ncores <- detectCores()  
cl <- makeCluster(ncores)
registerDoParallel(cl)

flnm <- 'APNERR2020'
tz <- 'America/Jamaica'
lat <- 29.75
long <- -85
dy <- 12
hr <- 6
td <- 0.5

# data
dat <- read.csv(paste0('data/raw/', flnm, '.csv')) %>%
  mutate(DateTimeStamp = mdy_hm(DateTimeStamp, tz = tz)) %>%
  na.omit %>%
  unique

wins <- list(12, 21.57143, 1.747619)
wins <- list(13.86532, 20.6448, 1.6704)
# weighted regression, optimal window widths for SAPDC from the paper
wtreg_res <- wtreg(dat, wins = list(12, 12, 0.95), progress = T, parallel = T,
                   tz = tz, lat = lat, long = long)

metab_dtd <- ecometab(wtreg_res, DO_var = 'DO_nrm', tz = tz,
                      lat = lat, long = long)

meteval(metab_dtd, all = F)
