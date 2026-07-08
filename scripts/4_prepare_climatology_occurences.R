library(pacman)
p_load(tidyverse, ncdf4, lubridate)


# Load required data sets -----------------------------------------------

## Read in occurrence data
urch_occ <- readRDS('data/urchin_occurrence_records.rds')

summary(urch_occ)


## read in unique locations key
unique_locations <- readRDS('data/urchin_occurrence_records_unique_grouped.rds')

summary(unique_locations)


## Read in SST data for occurrence sites
sst1.nc <- nc_open('data/urchin_sst_1984_2023.nc')

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

# convert lon values to standard -180 to +180 (not 0-360)
sst1.df <- sst1.df %>% 
  mutate(lon = case_when(lon >= 180 ~ lon-360,
                             lon < 180 ~ lon,
                             TRUE ~ 9999),
         origin_lon = case_when(origin_lon >= 180 ~ origin_lon-360,
                                origin_lon < 180 ~ origin_lon,
                                    TRUE ~ 9999))

names(sst1.df) <- c('group_lat', 'group_lon', 'sst_lat', 'sst_lon', time) 

# repeat for additional sites from second extraction
sst2.nc <- nc_open('data/followup_sst_1984_2023.nc')

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

# select only the required two sites
sst2.df <- sst2.df[c(1,2), ]

# convert lon values to standard -180 to +180 (not 0-360)
sst2.df <- sst2.df %>% 
  mutate(lon = case_when(lon >= 180 ~ lon-360,
                         lon < 180 ~ lon,
                         TRUE ~ 9999),
         origin_lon = case_when(origin_lon >= 180 ~ origin_lon-360,
                                origin_lon < 180 ~ origin_lon,
                                TRUE ~ 9999))

names(sst2.df) <- c('group_lat', 'group_lon', 'sst_lat', 'sst_lon', time) 

# join them together
sst1.df <- rbind(sst1.df, sst2.df)

# when r handles the matching between data sets, it seems to round the lat/lon 
# coords to slightly different values and consequently there are many missing 
# values in the match. By rounding the 3 decimal places before the match, most of 
# these are avoided and just have to manually match one site

sst1.df <- sst1.df %>% 
  mutate(group_lat = round(group_lat, 3),
         group_lon = round(group_lon, 3))

unique_locations <- unique_locations %>% 
  mutate(group_lat = round(group_lat, 3),
         group_lon = round(group_lon, 3))

# Join SST stats to occurrence data set ------------

sst1 <- sst1.df %>% 
  pivot_longer(cols = 5:14613, names_to = "day", values_to = "sst") %>% 
  mutate(day = as.numeric(day),
         date = as_date(day, origin = ymd('1978-01-01')))

sst1 <- sst1 %>%
  mutate(dayMonth = format(date, "%d-%m"),
         year = format(date, "%Y"))

# summarise data into a climatology for each location
sst1_sum1 <- sst1 %>% 
  group_by(group_lat, group_lon, dayMonth) %>% 
  summarise(N_daily = length(sst),
            mean_daily = mean(sst),
            median_daily = median(sst),
            sd_daily = sd(sst),
            se_daily = sd_daily / sqrt(N_daily), 
            percentile_1_daily = quantile(sst, 0.01), 
            percentile_99_daily = quantile(sst, 0.99))

# summarise to get the mean, median 1st, and 99th percentiles for each location
sst1_sum2 <- sst1_sum1 %>% 
  group_by(group_lat, group_lon) %>% 
  summarise(N = length(mean_daily),
            mean = mean(mean_daily),
            median = median(median_daily),
            sd = sd(mean_daily),
            se = sd / sqrt(N),
            percentile_1_med = quantile(median_daily, 0.01), 
            percentile_99_med = quantile(median_daily, 0.99),
            percentile_1_mean = quantile(mean_daily, 0.01), 
            percentile_99_mean = quantile(mean_daily, 0.99))


# join sst data to unique location key using group_lat and group_lon

occ.sst <- unique_locations %>% 
  left_join(sst1_sum2, by = join_by(group_lat, group_lon), relationship = 'many-to-one')

# there are several (424) non-unique sites in this list still, so will remove these locations
occ.sst <- occ.sst %>% 
  distinct(latitude, longitude, group_lat, group_lon, N, mean, sd, se, median, 
           percentile_1_med, percentile_99_med, percentile_1_mean, percentile_99_mean)

occ.sst <- occ.sst %>% 
  mutate(latitude_rounded = round(latitude, 6),
         longitude_rounded = round(longitude, 6))

unique1 <- occ.sst %>% 
  distinct(latitude_rounded, longitude_rounded)

urch_occ <- urch_occ %>% 
  mutate(latitude_rounded = round(latitude, 6),
         longitude_rounded = round(longitude, 6))

unique2 <- urch_occ %>% 
  distinct(latitude_rounded, longitude_rounded)

### join to urch_occ
urch_sst <- urch_occ %>% 
  left_join(occ.sst, by = join_by(latitude_rounded, longitude_rounded), relationship = 'many-to-one')


# manually match two sites - using survey ID instead of lat/long as it's being a pain
urch_sst[c(which(urch_sst$survey_id == 923400690)), 49:61] <- 
  occ.sst[c(which(occ.sst$latitude_rounded < -28.162 & occ.sst$latitude_rounded > -28.163)), 1:13]

urch_sst[c(which(urch_sst$survey_id == 923400004)), 49:61] <- 
  occ.sst[c(which(occ.sst$latitude_rounded < -30.478 & occ.sst$latitude_rounded > -30.481)), 1:13]

urch_sst[c(which(urch_sst$survey_id == 923400005)), 49:61] <- 
  occ.sst[c(which(occ.sst$latitude_rounded < -30.478 & occ.sst$latitude_rounded > -30.481)), 1:13]


summary(urch_sst) 


# Save occurrence data matched to SST -------------------------------------


saveRDS(urch_sst, 'data/urchin_occurrence_sst.rds')



