library(tidyverse)
library(lubridate)
library(WtRegDO)

# dat <- read.csv('data/raw/APNERR2020.csv') %>% 
#   mutate(
#     DateTimeStamp = mdy_hm(DateTimeStamp, tz = 'America/Jamaica'), 
#     mo = month(DateTimeStamp), 
#     hr = hour(DateTimeStamp)
#   )
# 
# flnm <- 'APNERR2020.csv'
# tz <- 'America/Jamaica'
# lat <- 29.75
# long <- -85
# dy <- 12
# hr <- 6
# td <- 0.5
# 
# # data
# dat <- read.csv('data/raw/', paste0(flnm, '.csv')) %>% 
#   mutate(DateTimeStamp = mdy_hm(DateTimeStamp, tz = tz)) %>% 
#   na.omit %>% 
#   unique
# 
# # weighted regression, optimal window widths for SAPDC from the paper
# wtreg_res <- wtreg(dat, parallel = F, wins = list(dy, hr, td), progress = F, 
#                    tz = tz, lat = lat, long = long)
# 
# # estimate ecosystem metabolism using observed DO time series
# metab_obs <- ecometab(wtreg_res, DO_var = 'DO_obs', tz = tz, 
#                       lat = lat, long = long)
# 
# # estimate ecosystem metabolism using detided DO time series
# metab_dtd <- ecometab(wtreg_res, DO_var = 'DO_nrm', tz = tz, 
#                       lat = lat, long = long)


toplo <- mins %>% 
  filter(obj == 'minano') %>% 
  filter(nm == 'SAPDC') %>% 
  mutate(est = pmax(350, est))

ggplot(toplo, aes(x = factor(dy), y = factor(hr), fill = est)) + 
  geom_tile() + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  facet_grid(nm~td)


  