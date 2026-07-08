library(pacman)
p_load(tidyverse, sf, ncdf4, FNN)


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

worldmap_sf <- worldmap %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

# Load data ---------------------------------------------------------------

urch_occ <- read_rds('data/urchin_occurrence_records.rds')

## get unique locations
unique_locations <- urch_occ %>% 
  distinct(program, country, site_name, longitude, latitude) 

unique_locations_sf <- unique_locations %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


## the below code creates a grid to match that of the SST dataset, and then matches each 
## occurrence point to the nearest grid point

# Make grid based off sst grid centre point -------------------------------

# create centre point of grid
centre_point <- data.frame(lon = -0.125, lat = 0.125)

centre_point_sf <- centre_point %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

## create outside corners of grid
box <- data.frame(lon = c(-180, 180, 180, -180), lat = c(90, 90, -90, -90)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

ggplot() +
  geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), fill="#D6DBDF", alpha=1) +
  geom_sf(data = box, colour = 'blue') +
  theme_format

bbox <- box%>% 
  st_bbox() %>% 
  st_as_sfc()

## fill in grid
grid <- st_make_grid(bbox, cellsize = c(0.25, 0.25), what = 'centers')

grid_coords <- as.data.frame(st_coordinates(grid))

## what are coordinates of the centre point?
zeros <- grid_coords[c(which(grid_coords$X == -0.125 & grid_coords$Y < 1.5)), ]
zeros <- zeros[c(which(zeros$Y > -1.5)), ]

## plot grid
ggplot() +
  geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), fill="#D6DBDF", alpha=1) +
  geom_point(data = grid_coords, aes(x = X, y = Y), colour = 'blue') +
  geom_point(data = centre_point, aes(x = lon, y = lat), colour = 'red') +
  coord_sf(xlim = c(-5, 5), ylim = c(-5, 5)) +
  theme_format


# Match occurrence point to SST grid --------------------------------------

## calculate the nearest grid point to each occurrence point
unique_locations_coords <- unique_locations %>% 
  dplyr::select(longitude, latitude) %>% 
  dplyr::rename(X = longitude, Y = latitude)

nn <- get.knnx(grid_coords, unique_locations_coords, k = 1) # uses nearest neighbour algorithm

unique_locations <- unique_locations %>% 
  mutate(group_lat = grid_coords$Y[nn$nn.index],
         group_lon = grid_coords$X[nn$nn.index])

unique_locations_grouped <- unique_locations %>% 
  distinct(group_lat, group_lon)

## plot the original and grouped points

ggplot() +
  geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), fill="#D6DBDF", alpha=1) +
  geom_point(data = unique_locations, aes(x = longitude, y = latitude), colour = 'blue') +
  geom_point(data = unique_locations, aes(x = group_lon, y = group_lat), colour = 'red') +
  coord_sf(xlim = c(151, 152), ylim = c(-35, -33)) +
  theme_format


# Save as rds -------------------------------------------------------------


saveRDS(unique_locations, 'data/urchin_occurrence_records_unique_grouped.rds')
