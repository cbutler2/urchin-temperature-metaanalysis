library(pacman)
p_load(tidyverse, rgbif, usethis)


# create dataframe of required species
urchin_species <- as.data.frame(c('Centrostephanus rodgersii', 'Heliocidaris erythrogramma', 
                                  'Strongylocentrotus droebachiensis', 'Paracentrotus lividus', 
                                  'Arbacia lixula'))
names(urchin_species) <- 'Taxon_name'

# extract taxon keys for each species from GBIF
urchin_keys <- 
  urchin_species %>% 
  pull("Taxon_name") %>% # get species name
  name_backbone_checklist()  %>% # match to backbone
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) # get the GBIF taxonkeys



# submit download query GBIF database for selected records 
urchin_query <- occ_download(pred_in("taxonKey", urchin_keys),
                                 #We will keep presence data only
                                 pred("occurrenceStatus", "PRESENT"),
                                 #We will keep georeferenced records only
                                 pred("hasCoordinate", TRUE),
                                 format = )


# check the status of our download query
occ_download_wait(urchin_query)

# now download the dataset

#This section downloads a zip file to disk
gbif_urch <- occ_download_get(key = urchin_query, path = "data/raw/",
                               overwrite = T) %>% 
  #This section loads zip file into memory
  occ_download_import(gbif_urch, path = "data/raw/",
                      na.strings = c("", NA))



saveRDS(gbif_urch, 'data/gbif_urchins_2023.rds')