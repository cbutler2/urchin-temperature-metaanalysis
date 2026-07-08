library(pacman)
p_load(plyr, tidyverse, rgbif, sf, CoordinateCleaner, rnaturalearthdata, countrycode, ncdf4, FNN)


# Preliminaries -----------------------------------------------------------

##  Basic plot formatting
theme_format <- 
  theme_bw()+
  theme(axis.text.x  = element_text(vjust=0.5,size=12, colour = "black")) +
  theme(axis.text.y  = element_text(size=12, colour = "black")) +
  theme(axis.title.x = element_text(size=18, colour = "black")) +
  theme(axis.title.y = element_text(size=18, colour = "black")) +
  theme(axis.ticks = element_line(colour="black")) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major=element_blank())


## World map
worldmap<-map_data("world")

## load de-bugged gbif functions
source("scripts/gbif_functions.R")



# Load GBIF data --------------------------------------------------

gbif_urch <- read_rds('data/raw/gbif_urchins_2023.rds')

gbif_urch <- gbif_urch %>% 
  select(gbifID, species, countryCode, decimalLatitude, decimalLongitude, eventDate, year, 
         depth, basisOfRecord, individualCount, organismQuantityType, taxonRank, 
         family, institutionCode, datasetName, occurrenceStatus, georeferenceVerificationStatus, 
         coordinateUncertaintyInMeters, coordinatePrecision, depthAccuracy, issue) %>% 
  mutate_if(is.character, as.factor)

glimpse(gbif_urch)
summary(gbif_urch)



# Initial inspections ------------------------------------------------------

# plot each species distribution
gbif_urch_sf <- st_as_sf(gbif_urch, coords = c('decimalLongitude', 'decimalLatitude'))

gbif_urch_sf <- st_set_crs(gbif_urch_sf, 4326)


urchin_species <- unique(gbif_urch$species)
# create empty list to store plots in
species_dist_maps <- list()

for (i in 1:length(urchin_species)) {

# subset species
sub <- gbif_urch_sf %>%
  filter(species == urchin_species[i])

# plot occurence data
plot <- ggplot() +
    geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), fill="#D6DBDF", alpha=1) +
    geom_sf(data=sub, size=3, shape=21, aes(fill=countryCode)) +
    theme_format +
    labs(title = urchin_species[i], x = 'Longitude', y = 'Latitude')

# add plot to list
species_dist_maps[[i]] <- plot


}

names(species_dist_maps) <- c('sdroebachiensis',  'plividus', 'alixula', 'herythrogramma', 'crodgersii')


species_dist_maps


# Filter data -----------------------------------------------------------

# filter out records based on their type
# see https://docs.gbif.org/course-data-use/en/basis-of-record.html for basisOfRecord definitions
gbif_urch <- gbif_urch %>% 
  mutate(countryCode = factor(countrycode(sourcevar = countryCode, origin = 'iso2c', destination = 'iso3c'))) %>% 
  filter(basisOfRecord != 'FOSSIL_SPECIMEN' & 
           basisOfRecord != 'MATERIAL_SAMPLE' &
           basisOfRecord != 'PRESERVED_SPECIMEN' &
           basisOfRecord != 'LIVING_SPECIMEN') %>% 
  droplevels()


# filter out records based on year
gbif_urch <- gbif_urch %>% 
  filter(year >= 1950) %>% 
  droplevels() 

# filter out records with high uncertainty
gbif_urch <- gbif_urch %>%
  filter(coordinateUncertaintyInMeters <= 10000 | 
           is.na(coordinateUncertaintyInMeters)) %>% 
  droplevels() 

# filter out records based on depth
gbif_urch <- gbif_urch %>% 
  filter(depth < 30 | 
           is.na(depth)) %>% 
  droplevels()




# Run CoordinateCleaner tests ---------------------------------------------

# run tests on each subset of species
urchin_species <- unique(gbif_urch$species)

# create empty list to store data in
gbif_urch_tested <- tibble()

# load buffland gazeteer
data("buffland")
buffland <- terra::vect(buffland) # this is a 1 degree buffered vector of global landmass (i.e. will exclude only those points in open ocean)

for (i in 1:length(urchin_species)) {
  
  # subset species  
  sub <- gbif_urch %>% 
    filter(species == urchin_species[i]) 
  
  # run outliers test
  sub <- sub %>%
    mutate(val = cc_val(., value = "flagged"), 
           equ = cc_equ(., value = "flagged"), 
           zero = cc_zero_2(., value = "flagged"), 
           cap = cc_cap(., value = "flagged", buffer = 1000), # reduced buffer from 10000 as several capitals are near the coast
           cen = cc_cen(., value = "flagged"), # default buffer is 1km, which should be reasonable even for small Island country like Malta
           gbf = cc_gbif(., value = "flagged"), 
           inst = cc_inst(., value = "flagged"), # default buffer is 100m
           urb = cc_urb(., value = "flagged"), 
           coun2 = cc_coun_2(., value = "flagged", iso3 = "countryCode", 
                             buffer = 20000),  
           seas = cc_sea(., ref = buffland, value = "flagged"), 
           dupl = cc_dupl(., additions = c("species", "decimalLatitude", "decimalLongitude", 
                                           "eventDate", "depth", "individualCount", "basisOfRecord", 
                                           "institutionCode", "datasetName"), value = "flagged"), 
           outl = cc_outl(., value = "flagged")) # NOTE: method = distance may be more accurate - see section below regarding results from different methods    
  
  gbif_urch_tested <- rbind(gbif_urch_tested, sub)
  
  
}


saveRDS(gbif_urch_tested, 'data/raw/gbif_urchins_2023_CCtests.rds')


# Inspect flagged points --------------------------------------------------

gbif_urch <- read_rds('data/raw/gbif_urchins_2023_CCtests.rds')

flagged_points <- gbif_urch %>% 
  pivot_longer(cols = val:outl, names_to = 'test', values_to = 'flagged') %>% 
  mutate(test = factor(test))  

flagged_points <- flagged_points %>% 
  filter((test == 'val' & flagged == FALSE) |
           (test == 'equ' & flagged == FALSE) |
           (test == 'zero' & flagged == FALSE) |
           (test == 'cap' & flagged == FALSE) |
           (test == 'cen' & flagged == FALSE) |
           (test == 'gbf' & flagged == FALSE) |
           (test == 'inst' & flagged == FALSE) |
           (test == 'urb' & flagged == FALSE) |
           (test == 'coun2' & flagged == FALSE) |
           (test == 'seas' & flagged == FALSE) |
           (test == 'outl' & flagged == FALSE) |
           (test == 'dupl' & flagged == FALSE)) 

unique(flagged_points[c('test', 'flagged')])

print(flagged_points %>% 
        count(species, test), n = Inf)

# add CRS
flagged_points_sf <- flagged_points %>% 
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) 

cc_tests <- unique(flagged_points_sf$test)

# create empty list to store plots in
flagged_points_maps <- list()

for (i in 1:length(cc_tests)) {
  
  # subset species  
  sub <- flagged_points_sf %>% 
    filter(test == cc_tests[i]) 
  
  # plot occurence data
  plot <- ggplot() +
    geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), fill="#D6DBDF", alpha=1) +
    geom_sf(data=sub, size=3, aes(colour = species, shape = species)) +
    theme_format +
    labs(title = paste('flagged_', cc_tests[i], sep = ""), x = 'Longitude', y = 'Latitude') 
  
  # add plot to list
  flagged_points_maps[[i]] <- plot
  
  
}

names(flagged_points_maps) <- c('coun2', 'seas', 'dupl', 'urb', 'inst', 'outl', 'zero2', 'cen', 'cap')


flagged_points_maps 


# Clean data --------------------------------------------------------------

gbif_clean <- gbif_urch %>% 
  filter(val == TRUE & # no flagged values
           equ == TRUE & # no flagged values
           zero == TRUE &
           cap == TRUE & 
           cen == TRUE &  
           gbf == TRUE & # no flagged values
           inst == TRUE &
           seas == TRUE & 
           dupl == TRUE) 


# Manually remove outliers based on Lawrence (2020) Sea Urchins: Biology and Ecology -------------------------

gbif_clean <- gbif_clean %>% 
  dplyr::mutate(rowid = row_number()) %>% 
  dplyr::select(rowid, everything())


# remove arbacia outliers
gbif_clean <- gbif_clean %>% 
  filter(!(species == 'Arbacia lixula' & countryCode == 'USA') &
           !(species == 'Arbacia lixula' & countryCode == 'AUS') &
           !(species == 'Arbacia lixula' & countryCode == 'NZL') &
           !(species == 'Arbacia lixula' & countryCode == 'ARE') &
           !(species == 'Arbacia lixula' & countryCode == 'KWT') &
           !(species == 'Arbacia lixula' & countryCode == 'EGY') &
           !(species == 'Arbacia lixula' & countryCode == 'FRA' & decimalLatitude > 45))

# 9 records


# remove centrostephanus outliers
gbif_clean <- gbif_clean %>% 
  filter(!(species == 'Centrostephanus rodgersii' & countryCode == 'CHL')) %>% 
  filter(!(rowid == 22949))

# remove paracentotus outliers
gbif_clean <- gbif_clean %>% 
  filter(!(species == 'Paracentrotus lividus' & countryCode == 'REU') &
           !(species == 'Paracentrotus lividus' & countryCode == 'CAN') &
           !(species == 'Paracentrotus lividus' & countryCode == 'USA') &
           !(species == 'Paracentrotus lividus' & countryCode == 'NCL') &
           !(species == 'Paracentrotus lividus' & countryCode == 'GMB') &
           !(species == 'Paracentrotus lividus' & countryCode == 'SLE') &
           !(species == 'Paracentrotus lividus' & countryCode == 'IDN') &
           !(species == 'Paracentrotus lividus' & countryCode == 'THA')) 

# remove strongylocentotus outliers
gbif_clean <- gbif_clean %>% 
  filter(!(species == 'Strongylocentrotus droebachiensis' & countryCode == 'MLT') &
           !(species == 'Strongylocentrotus droebachiensis' & countryCode == 'GRC') &
           !(species == 'Strongylocentrotus droebachiensis' & countryCode == 'FRA') &
           !(species == 'Strongylocentrotus droebachiensis' & countryCode == 'SGP')) %>% 
  filter(!(rowid == 5072) & 
           !(rowid == 1349))

gbif_clean_sf <- gbif_clean %>%
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326)

urchin_species <- unique(gbif_clean$species)

# create empty list to store plots in
species_dist_maps_cleaned <- list()

for (i in 1:length(urchin_species)) {

  # subset species
  sub <- gbif_clean_sf %>%
    filter(species == urchin_species[i])

  # plot occurence data
  plot <- ggplot() +
    geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), fill="#D6DBDF", alpha=1) +
    geom_sf(data=sub, size=3, shape=3, aes(colour=countryCode)) +
    theme_format +
    labs(title = urchin_species[i], x = 'Longitude', y = 'Latitude')

  # add plot to list
  species_dist_maps_cleaned[[i]] <- plot


}

names(species_dist_maps_cleaned) <- c('sdroebachiensis', 'plividus', 'alixula', 'herythrogramma', 'crodgersii')


species_dist_maps_cleaned



# Load RLS data -----------------------------------------------------------


rls_inv <- read_csv('data/ep_m2_inverts_ALL_202301.csv')

rls_inv <- rls_inv %>% 
  select(!diver) %>% 
  mutate(across(where(is.character), as.factor)) %>% 
  mutate(survey_id = factor(survey_id),
         method = factor(method),
         block = factor(block))

glimpse(rls_inv)
summary(rls_inv)

rls_urch <- rls_inv %>% 
  filter(species_name == 'Centrostephanus rodgersii' |
           species_name == 'Heliocidaris erythrogramma' |
           species_name == 'Paracentrotus lividus' |
           species_name == 'Strongylocentrotus droebachiensis' |
           species_name == 'Arbacia lixula') %>% 
  droplevels()

summary(rls_urch)  

# remove the centro point from WA
rls_urch <- rls_urch[-which(rls_urch$longitude < 130 & 
                              rls_urch$species_name == 'Centrostephanus rodgersii' & 
                              rls_urch$location == 'Jurien'), ]

# Plot maps

rls_urch_sf <- rls_urch %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)

urchin_species <- unique(rls_urch$species_name)

# create empty list to store plots in
rls_occurence_maps <- list()

for (i in 1:length(urchin_species)) {

  # subset species
  sub <- rls_urch_sf %>%
    filter(species_name == urchin_species[i])

  # plot occurence data
  plot <- ggplot() +
    geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), fill="#D6DBDF", alpha=1) +
    geom_sf(data=sub, size=3, shape=3, aes(colour=country)) +
    theme_format +
    labs(title = urchin_species[i], x = 'Longitude', y = 'Latitude')


  # add plot to list
  rls_occurence_maps[[i]] <- plot


}

names(rls_occurence_maps) <- c('herythrogramma', 'crodgersii', 'sdroebachiensis', 'alixula', 'plividus')


rls_occurence_maps



# Join datasets together --------------------------------------------------

gbif_clean <- gbif_clean %>% 
  mutate(program = factor('GBIF')) %>% 
  mutate(country = factor(countrycode(sourcevar = countryCode, origin = 'iso3c', destination = 'country.name'))) %>%
  dplyr::rename(latitude = decimalLatitude,
                longitude = decimalLongitude,
                date = eventDate) %>% 
  select(!val:outl) %>% 
  select(!countryCode)

rls_urch <-  rls_urch %>% 
  dplyr::rename(species = species_name,
                date = survey_date)
  
urch_occ <- join(rls_urch, gbif_clean, type = 'full')


# saved rds of GBIF and RLS occurrence data --------------------------------


saveRDS(urch_occ, 'data/urchin_occurrence_records.rds')
















