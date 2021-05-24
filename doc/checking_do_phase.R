# library(WtRegDO)
devtools::load_all('../WtRegDO')
library(tidyverse)
library(here)
library(lubridate)
library(patchwork)

# setup parallel backend
library(doParallel)
registerDoParallel(cores = detectCores() - 1)

flnm <- 'APNERR2012'
tz <- 'America/Jamaica'
lat <- 29.75
long <- -85

# data
dat <- paste0('data/raw/', flnm, '.csv') %>% 
  read.csv() %>%
  mutate(DateTimeStamp = mdy_hm(DateTimeStamp, tz = tz)) %>%
  na.omit %>%
  unique %>% 
  filter(DateTimeStamp < date('2012-08-01') & DateTimeStamp > date('2012-03-01'))

# weighted regression, optimal window widths for SAPDC from the paper
str <- Sys.time()
wtreg_res <- wtreg(dat, wins = list(2, 6, 0.8), progress = F, parallel = T,
                   tz = tz, lat = lat, long = long)
Sys.time() - str

metest <- ecometab(wtreg_res, tz = tz, DO_var = 'DO_nrm', lat = lat, long = long)
perf <- meteval(metest, all = F) %>% 
  data.frame %>% 
  round(1)

toplo <- wtreg_res %>% 
  select(DateTimeStamp, DO_obs, DO_prd, DO_nrm, totpar, Tide) %>%
  mutate_at(vars(matches('DO')), function(x) (x * 1 / 0.032)) %>% 
  mutate(
    Date = date(DateTimeStamp),
    utc_time = as.POSIXlt(DateTimeStamp, tz = 'UTC'),
    sun_angle = oce::sunAngle(utc_time, long, lat)$altitude
  ) %>% 
  group_by(Date) %>% 
  mutate(
    donrmmx = DateTimeStamp[which.max(DO_nrm)],
    sunmx = DateTimeStamp[which.max(sun_angle)]
  ) %>% 
  ungroup() %>%
  select(-utc_time, -sun_angle) %>%
  gather('var', 'val', -Date, -DateTimeStamp, -donrmmx, -sunmx) %>% 
  mutate(
    var = case_when(
      var == 'DO_nrm' ~ 'DO dtd (mmol/m3)', 
      var == 'DO_obs' ~ 'DO obs (mmol/m3)', 
      var == 'DO_prd' ~ 'DO prd (mmol/m3)', 
      var == 'totpar' ~ 'PAR (W/m2)' ,
      var == 'Tide' ~ 'Tide (m)'
    )
  )

pthm <- theme_bw() +
  theme(
    axis.title.y = element_blank(), 
    strip.placement = 'outside', 
    legend.position = 'none', 
    strip.background = element_blank(), 
    axis.title.x = element_blank(), 
    panel.grid = element_blank(), 
    strip.text.y = element_text(size = 8)
  )

toplo1 <- toplo %>% 
  filter(Date >= date('2012-03-18') & Date <= date('2012-03-20'))

p1 <- ggplot(toplo1, aes(x = DateTimeStamp, y = val, color = var)) + 
  geom_line() + 
  geom_vline(aes(xintercept = sunmx, colour = 'PAR (W/m2)')) +
  geom_vline(aes(xintercept = donrmmx, colour = 'DO dtd (mmol/m3)')) +
  facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y') + 
  pthm + 
  labs(
    title = "Model with sine save, 2, 6, 0.8",
    subtitle = 'DO/PAR in phase'
  )

toplo2 <- toplo %>% 
  filter(Date >= date('2012-06-01') & Date <= date('2012-06-03'))

p2 <- ggplot(toplo2, aes(x = DateTimeStamp, y = val, color = var)) + 
  geom_line() + 
  geom_vline(aes(xintercept = sunmx, colour = 'PAR (W/m2)')) + 
  geom_vline(aes(xintercept = donrmmx, colour = 'DO dtd (mmol/m3)')) +
  facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y') + 
  pthm +  
  labs(
    subtitle = 'DO lagging PAR'
  )

p <- p1 + p2 + plot_layout(ncol = 2)

png('~/Desktop/orignalmod.png', res = 300, units = 'in', height = 6, width = 7)
print(p)
dev.off()
