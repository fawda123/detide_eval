library(tidyverse)
library(lubridate)
library(WtRegDO)

dat <- read.csv('data/raw/APNERR2020.csv') %>% 
  mutate(
    DateTimeStamp = mdy_hm(DateTimeStamp, tz = 'America/Jamaica'), 
    mo = month(DateTimeStamp), 
    hr = hour(DateTimeStamp)
  )

flnm <- 'APNERR2020.csv'
tz <- 'America/Jamaica'
lat <- 29.75
long <- -85
dy <- 12
hr <- 6
td <- 0.5

# data
dat <- read.csv('data/raw/', paste0(flnm, '.csv')) %>% 
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

# the value returned by this function should be maximized
objfun <- function(metab_obs, metab_dtd, vls = c('meanPg', 'sdPg', 'anomPg', 'meanRt', 'sdRt', 'anomRt')){

  eval <- list(
      obseval = enframe(meteval(metab_obs)),
      dtdeval = enframe(meteval(metab_dtd))
    ) %>% 
    enframe('metdat', 'fitdat') %>% 
    unnest('fitdat') %>% 
    filter(name %in% vls) %>% 
    unnest('value') %>% 
    spread(metdat, value) %>% 
    mutate(
      perdif = (abs(obseval - dtdeval)) / ((obseval + dtdeval) / 2), 
      perdif = case_when(
        name %in% c('meanPg', 'meanRt') ~ 1 / perdif,
        T ~ perdif
      )
    )
  
  est <- sum(eval$perdif)
  
  return(est)

}

obj_fun(metab_obs, metab_dtd)

# apply objective function to already processed results
metab_obs <- metobs_SAPDC_12_12_0.2

mins <- list.files('data/', pattern = 'metdtd_SAPDC') %>% 
  tibble(fl = .) %>% 
  mutate(fl = gsub('\\.RData$', '', fl)) %>% 
  group_by(fl) %>% 
  nest() %>% 
  mutate(
    data = purrr::pmap(list(fl), function(fl){
      
      cat(fl, '\n')
      load(file = paste0('data/', fl, '.RData'))
      dat <- get(fl)
      objfun(metab_obs, dat)
      
    })
  ) %>% 
  unnest(data)

toplo <- mins %>% 
  mutate(fl = gsub('^metdtd\\_SAPDC\\_', '', fl)) %>% 
  separate(fl, c('dy', 'hr', 'tid'), sep = '_') %>% 
  mutate_if(is.character, as.numeric) %>% 
  mutate(data = pmax(-0.1, data))

ggplot(toplo, aes(x = factor(dy), y = factor(hr), fill = data)) + 
  geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  facet_grid(~tid)


  