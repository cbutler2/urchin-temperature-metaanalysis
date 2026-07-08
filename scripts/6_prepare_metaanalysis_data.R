library(pacman)
p_load(tidyverse, nplyr)


# Load required datasets --------------------------------------------------
d1 <- read_rds('data/UrchinReview_CompiledData_clean.rds')


# Add climatology data ----------------------------------------------------

sst <- readRDS('data/prepared/climatology_meta-analysis.rds')

sst_summary  <-  sst %>%
  group_by(Latitude, Longitude) %>% 
  dplyr::summarise(N = length(mean_daily),
                   mean = mean(mean_daily),
                   median = median(median_daily),
                   sd = sd(mean_daily),
                   se = sd / sqrt(N),
                   percentile_1_mean = quantile(mean_daily, 0.01), 
                   percentile_34_mean = quantile(mean_daily, 0.34), 
                   percentile_66_mean = quantile(mean_daily, 0.66), 
                   percentile_99_mean = quantile(mean_daily, 0.99))

d1 <- d1 %>% 
  mutate(Latitude = round(Latitude, 3),
         Longitude = round(Longitude, 3))

d1 <- d1 %>% 
  left_join(sst_summary, by = join_by(Latitude, Longitude), relationship = 'many-to-one')

d1 <- d1 %>% 
  dplyr::rename(sd_sst = sd.y, 
                se_sst = se.y,
                sd = sd.x,
                se = se.x)

summary(d1)

# Calculate range index ---------------------------------------------------

# load occurrence data
urch_sst <- readRDS('data/prepared/urchin_occurrence_sst.rds')

### calculate thermal distribution metrics
thermal_dist <- urch_sst %>% 
  group_by(species) %>% 
  dplyr::summarise(species_N = length(mean),
                   species_mean = mean(mean),
                   species_median = median(median),
                   median_max = max(median),
                   median_min = min(median)) %>% 
  mutate(species_midpoint_median = (median_min + median_max)/2, #get standard deviation
         species_median_range = median_max-median_min) %>% 
  select(species, species_N, species_mean, species_median,
         species_midpoint_median, species_median_range)

# join to meta-analysis data
d1 <- d1 %>% 
  left_join(thermal_dist, by = join_by(Species == species), relationship = 'many-to-one')

d1 <- d1 %>% 
  mutate(Range_Index = 2*(median - species_midpoint_median)/species_median_range)

d1 <- d1 %>% 
  mutate(Range_Position = factor(case_when(Range_Index < -0.34 ~ 'Cool edge',
                                    Range_Index >= -0.34 & Range_Index <= 0.34 ~ 'Central',
                                    Range_Index > 0.34 ~ 'Warm edge')))


summary(d1)


# Calculate inverse of mortality -------------------------------------------

d1 <- d1 %>% 
  mutate(Response_Value = case_when(Response_Type2 == 'Abnormal/Dead' & Response_Value_Unit == 'percent success' ~ 100 - Response_Value,
                                    Response_Type2 == 'Arrested' & Response_Value_Unit == '%' ~ 100 - Response_Value,
                                    Response_Type2 == 'Dead/Arrested' & Response_Value_Unit == 'percent success' ~ 100 - Response_Value,
                                    Response_Type2 == 'Mortality' & Response_Value_Unit == '%' ~ 100 - Response_Value,
                                    Response_Type2 == 'Mortality' & Response_Value_Unit == 'percent' ~ 100 - Response_Value,
                                    Response_Type2 == 'dead' & Response_Value_Unit == 'No. individuals (out of 15)' ~ 15 - Response_Value,
                                    TRUE ~ Response_Value))

d1 <- d1 %>% 
  mutate(Response_Type2 = case_when(Response_Type2 == 'Abnormal/Dead' & Response_Value_Unit == 'percent success' ~ 'Survival',
                                    Response_Type2 == 'Arrested' & Response_Value_Unit == '%' ~ 'Survival',
                                    Response_Type2 == 'Dead/Arrested' & Response_Value_Unit == 'percent success' ~ 'Survival',
                                    Response_Type2 == 'Mortality' & Response_Value_Unit == '%' ~ 'Survival',
                                    Response_Type2 == 'Mortality' & Response_Value_Unit == 'percent' ~ 'Survival',
                                    Response_Type2 == 'dead' & Response_Value_Unit == 'No. individuals (out of 15)' ~ 'Survival',
                                    TRUE ~ Response_Type2))



# Calculate inverse of abnormal -------------------------------------------

d1 <- d1 %>% 
  mutate(Response_Value = case_when(Response_Type2 == 'Abnormal' & Response_Value_Unit == 'No. individuals (out of 15)' ~ 15 - Response_Value,
                                    Response_Type2 == 'Abnormal plutei' & Response_Value_Unit == '%' ~ 100 - Response_Value,
                                    Response_Type2 == 'Abnormal prism' & Response_Value_Unit == '%' ~ 100 - Response_Value,
                                    TRUE ~ Response_Value))

d1 <- d1 %>% 
  mutate(Response_Type2 = case_when(Response_Type2 == 'Abnormal' & Response_Value_Unit == 'No. individuals (out of 15)' ~ 'Normal development',
                                    Response_Type2 == 'Abnormal plutei' & Response_Value_Unit == '%' ~ 'Normal development',
                                    Response_Type2 == 'Abnormal prism' & Response_Value_Unit == '%'  ~ 'Normal development',
                                    TRUE ~ Response_Type2))


d1 %>% 
  distinct(Rayyan_System_ID, Study_ID)



# Add columns for contrasts -----------------------------------------------

contrasts <- d1 %>% 
  filter(!is.na(Contrast)) # 175 obs

levels(contrasts$Contrast)

d1 <- d1 %>% 
  mutate(Feed_Regime = factor(case_when(Contrast == 'Fed' ~ 'High',
                                        Contrast == 'Starved' ~ 'Low',
                                        Contrast == 'High food*Winter' ~ 'High',
                                        Contrast == 'Low food*Winter' ~ 'Low',
                                        Contrast == 'High food*Summer' ~ 'High',
                                        Contrast == 'Low food*Summer' ~ 'Low',
                                        is.na(Contrast) ~ 'Unknown',
                                        TRUE ~ 'Unknown')))

d1 <- d1 %>% 
  mutate(Urchin_Sex = factor(case_when(Contrast == 'Male' ~ 'Male',
                                       Contrast == 'Female' ~ 'Female',
                                       is.na(Contrast) ~ 'Unknown',
                                       TRUE ~ 'Unknown')))


# Convert experimental durations to same scale ---------------

# for all instances where duration is stated as a range, select the shortest duration
d1 <- d1 %>% 
  mutate(Experimental_Duration = case_when(Experimental_Duration == '4-Jun' ~ '4',
                                           TRUE ~ Experimental_Duration)) %>% 
  mutate(Experimental_Duration = as.numeric(Experimental_Duration))

# convert to same scale

d1 <- d1 %>% 
  mutate(Experimental_Duration = case_when(Experimental_Duration_Unit == 'min' ~ Experimental_Duration/1440,
                                           Experimental_Duration_Unit == 'hour' ~ Experimental_Duration/24,
                                           Experimental_Duration_Unit == 'day' ~ Experimental_Duration,
                                           Experimental_Duration_Unit == 'week' ~ Experimental_Duration*7,
                                           Experimental_Duration_Unit == 'month' ~ Experimental_Duration*30.4)) 


# Add time of year --------------------------------------------------------

d1 <- d1 %>% 
  mutate(Season_Estimate = factor(case_when(Control_Temp < percentile_34_mean ~ 'Winter',
                                            Control_Temp >= percentile_34_mean & Control_Temp <= percentile_66_mean ~ 'Autumn/Spring',
                                            Control_Temp > percentile_66_mean ~ 'Summer')))



# Calculate sd and effect size  ----------------------------------------------------

# for all instances where N is given as a range, select the lowest value
# select only mean response types (not median etc)
d1 <- d1 %>%
  filter(Response_Value_Type == 'mean' | is.na(Response_Value_Type)) %>%
  filter(!is.na(Replication_Treatment)) %>%
  mutate(Replication_Treatment = case_when(Replication_Treatment == '14-17' ~ '14',
                                           Replication_Treatment == '3-Jun' ~'3',
                                           Replication_Treatment == '5-Dec' ~'5',
                                           Replication_Treatment == '5-Nov' ~'5',
                                           Replication_Treatment == '6-Jan' ~ '6',
                                           Replication_Treatment == '9-Nov' ~ '9',
                                           Replication_Treatment == 'May-15' ~'5',
                                           Replication_Treatment == 'variable (7-17 aprox)' ~ '7',
                                           TRUE ~ Replication_Treatment)) %>%
  mutate(Replication_Treatment = as.numeric(as.character(Replication_Treatment)))


# calculate sd
d1 <- d1 %>% 
  mutate(sd = case_when(!is.na(se) ~ se*sqrt(Replication_Treatment),
                        is.na(se) & is.na(sd) ~ NA,
                        is.na(se) & !is.na(sd) ~ sd))

d1.nest <- d1 %>% 
  nest(Study_Nest = -c(Rayyan_System_ID, Study_ID))

# add columns for control mean, sd and n within each study 
for (i in 1:nrow(d1.nest)){
  d1.nest[[3]][[i]]$control_mean <-  d1.nest[[3]][[i]][which(d1.nest[[3]][[i]]$Control_Treatment == 'Control'), 35] # 40 = Response_Value
  d1.nest[[3]][[i]]$control_sd <-  d1.nest[[3]][[i]][which(d1.nest[[3]][[i]]$Control_Treatment == 'Control'), 36] # 41 = sd
  d1.nest[[3]][[i]]$control_n <-  d1.nest[[3]][[i]][which(d1.nest[[3]][[i]]$Control_Treatment == 'Control'), 18] # 17 = Replication_Treatment
  
}


d1 <- d1.nest %>% 
  unnest(Study_Nest)


# for some reason the above indexing is creating columns of 1 x 1 tibbles. Need to convert to numeric

d1 <- d1 %>% 
  mutate(control_mean = as.numeric(as.character(unlist(pull(d1[72])))),
         control_sd = unlist(pull(d1[73])), 
         control_n = unlist(pull(d1[74])))


# d1 <- d1 %>% 
#   mutate(Response_Value = as.numeric(as.character(Response_Value)))


# calculate effect size
d1 <- d1 %>% 
  mutate(sd_pooled = sqrt(((Replication_Treatment-1)*(sd^2) + 
                             (control_n-1)*(control_sd^2))/
                            (Replication_Treatment + control_n -2)),
         j = 1-(3/(4*(Replication_Treatment + control_n - 2)-1)),
         effect_size = ((Response_Value - control_mean)/sd_pooled)*j,
         es_variance = ((Replication_Treatment + control_n)/(Replication_Treatment*control_n)) + 
           (effect_size^2/(2*(Replication_Treatment + control_n))))


summary(d1)

d1 <- d1 %>% 
  filter(!is.na(effect_size)) %>% 
  droplevels() 


# Add thermal anomaly ------------------------------------------------------

d1 <- d1 %>% 
  nest(Study_Nest = -c(Rayyan_System_ID, Study_ID)) %>% 
  nest_mutate(Study_Nest, Temp_Anomaly = Treatment_Temp - Control_Temp) %>% 
  unnest(Study_Nest)

summary(d1)



# Organise dataframe columns ------------------------------------------------------

d1 <- d1 %>% 
  dplyr::rename(sst_N = N,
                mean_sst = mean,
                median_sst = median) %>% 
  nest(Thermal_Treat_Details = c(Ambient_Temp:Adjustment_Period_Unit),
       Variance_Measurements = c(se, Temp_se),
       Effect_Size_Parameters = c(control_mean:j, es_variance),
       Depth_Size = c(Depth_Specific_Collection, Urchin_TD_mm, Urchin_biomass_g),
       Other_Details = c(DOI),
       Site_SST_Metrics = c(sst_N, mean_sst, median_sst, sd_sst, se_sst, percentile_1_mean, percentile_34_mean, percentile_66_mean, percentile_99_mean),
       Thermal_Dist_Metrics = c(species_N, species_mean, species_median, species_midpoint_median, species_median_range)) %>% 
  select(Rayyan_System_ID:Depth_Range_Collection, Depth_Size, Species:Development_Stage, Replication_Treatment, Replication_Experimental, Thermal_Treatment, Thermal_Treat_Details, 
         Experimental_Duration, Response_ID, Response_Type2, Response_Type_Category, Response_Type_Level, Control_Treatment, Control_Temp:Response_Value_Type, Variance_Measurements,
         effect_size, Effect_Size_Parameters, Temp_Anomaly, Number_Treatments, Feed_Regime, Urchin_Sex, Contrast, Time_Of_Year, Season_Estimate, Comments, Other_Details, Range_Index, Range_Position, Site_SST_Metrics, Thermal_Dist_Metrics)


d1 <- d1 %>% 
  mutate(Study_ID = factor(Study_ID),
         Response_ID = factor(Response_ID),
         Response_Type2 = factor(Response_Type2))

glimpse(d1)

# Save as RDS -------------------------------------------------------------

saveRDS(d2, 'data/UrchinReview_CompiledData_spatial.rds')

