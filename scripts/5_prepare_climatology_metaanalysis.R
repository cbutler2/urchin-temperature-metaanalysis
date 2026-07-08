library(pacman)
p_load(tidyverse, ncdf4, lubridate)

# Load sst file from jiaxin --------------------------------------------

sst2.nc <- nc_open('data/metaanalysis_sst_1984_2023.nc')

print(sst2.nc)

time <- ncvar_get(sst2.nc, 'time')
origin_lat <- ncvar_get(sst2.nc, 'origin_lat')
origin_lon <- ncvar_get(sst2.nc, 'origin_lon')
lat <- ncvar_get(sst2.nc, 'lat')
lon <- ncvar_get(sst2.nc, 'lon')
sst2 <- ncvar_get(sst2.nc, 'sst')


sst2.df <- as.data.frame(sst2)

sst2.df <- t(sst2.df) # swap columns to rows and rows to columns

sst2.df <- as.data.frame(sst2.df)

sst2.df <- cbind(origin_lat, origin_lon, lat, lon, sst2.df)

names(sst2.df) <- c('Latitude', 'Longitude', 'sst_lat', 'sst_lon', time)

sst2 <- sst2.df %>% 
  pivot_longer(cols = 5:14613, names_to = "day", values_to = "sst") %>% 
  mutate(day = as.numeric(day),
         date = as_date(day, origin = ymd('1978-01-01')))

# repeat this for additional sites in second extraction

sst1.nc <- nc_open('data/raw/followup_sst_1984_2023.nc')

print(sst1.nc)

time <- ncvar_get(sst1.nc, 'time')
origin_lat <- ncvar_get(sst1.nc, 'origin_lat')
origin_lon <- ncvar_get(sst1.nc, 'origin_lon')
lat <- ncvar_get(sst1.nc, 'lat')
lon <- ncvar_get(sst1.nc, 'lon')
sst1 <- ncvar_get(sst1.nc, 'sst')


sst1.df <- as.data.frame(sst1)

sst1.df <- t(sst1.df) # swap columns to rows and rows to columns

sst1.df <- as.data.frame(sst1.df)

sst1.df <- cbind(origin_lat, origin_lon, lat, lon, sst1.df)

names(sst1.df) <- c('Latitude', 'Longitude', 'sst_lat', 'sst_lon', time)

# I only want one of these sites - the 3rd row (Latitude == 59.629)
sst1.df <- sst1.df[3, ]

sst1 <- sst1.df %>% 
  pivot_longer(cols = 5:14613, names_to = "day", values_to = "sst") %>% 
  mutate(day = as.numeric(day),
         date = as_date(day, origin = ymd('1978-01-01')))

# join to other dataset
sst2 <- rbind(sst2, sst1)

# convert lon values to standard -180 to +180 (not 0-360)
sst2 <- sst2 %>% 
  mutate(sst_lon = case_when(sst_lon >= 180 ~ sst_lon-360,
                             sst_lon < 180 ~ sst_lon,
                             TRUE ~ 9999),
         Longitude = case_when(Longitude >= 180 ~ Longitude-360,
                               Longitude < 180 ~ Longitude,
                               TRUE ~ 9999))

# when r handles the matching between data sets, it seems to round the lat/lon 
# coords to slightly different values and consequently there are many missing 
# values in the match. By rounding the 3 decimal places before the match, most of 
# these are avoided and just have to manually match one site

sst2 <- sst2 %>% 
  mutate(Latitude = round(Latitude, 3),
         Longitude = round(Longitude, 3))


### calculate climatology at each site
sst2 <- sst2 %>%
  mutate(dayMonth = format(date, "%d-%m"),
         year = format(date, "%Y"))


sst2_sum1  <-  sst2 %>%
  group_by(Latitude, Longitude, dayMonth) %>% 
  dplyr::summarise(N_daily = length(sst),
                   mean_daily = mean(sst),
                   median_daily = median(sst),
                   sd_daily = sd(sst),
                   se_daily = sd_daily / sqrt(N_daily), 
                   percentile_1 = quantile(sst, 0.01), 
                   percentile_99 = quantile(sst, 0.99))

saveRDS(sst2_sum1, 'data/climatology_meta-analysis.rds')
