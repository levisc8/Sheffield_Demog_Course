# Function to recreate Valentin's beautiful world maps
# Now, Valentin's map

# Load the data that he has already created. It contains some labels and objects 
# used to create different features within the map
# NE_box creates the bounding box for the world, NE_graticules creates
# lat/long lines, NE_countries creates the country polygons, and 
# lbl.x and lbl.y have lat/long labels. NE_places is actually unneeded, so I'll
# remove that next.

library(dplyr)
library(stringr)
library(rgdal)
library(ggplot2)

VSWorldMap <- function(){
  load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
  
  rm(NE_places)
  
  # Next, we'll need to reproject that data to get it into the right format
  PROJ <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
  
  # Reproject the countries, lat/long lines, and bounding box
  NE_countries.prj  <- spTransform(NE_countries, CRSobj = PROJ)
  NE_graticules.prj <- spTransform(NE_graticules, CRSobj = PROJ)
  NE_box.prj        <- spTransform(NE_box, CRSobj = PROJ)
  

  # __ project long-lat coordinates columns for data frames 
  # (two extra columns with projected XY are created)
  prj.coord <- project(cbind(lbl.Y$lon, lbl.Y$lat), proj = PROJ)
  lbl.Y.prj <- cbind(prj.coord, lbl.Y)
  names(lbl.Y.prj)[1:2] <- c("X.prj","Y.prj")
  
  prj.coord <- project(cbind(lbl.X$lon, lbl.X$lat), proj = PROJ)
  lbl.X.prj <- cbind(prj.coord, lbl.X)
  names(lbl.X.prj)[1:2] <- c("X.prj","Y.prj")
  
  
  # Now, to the mapping part!
  ggplot() +
    # Add countries
    geom_polygon(data = NE_countries.prj, 
                 aes(long,lat, group = group), 
                 colour = "gray70", fill = "gray90", size = .25) +
    # Create a layer for the bounding box
    geom_polygon(data = NE_box.prj, 
                 aes(x = long, y = lat), 
                 colour = "black", fill = "transparent", size = .25) +
    # Add lat and longitude lines
    geom_path(data = NE_graticules.prj, 
              aes(long, lat, group = group), 
              linetype = "dotted", colour = "grey50", size = .25) +
    
    # add labels - latitude and longitude
    geom_text(data = lbl.Y.prj, # latitude
              aes(x = X.prj, y = Y.prj, label = lbl), 
              colour = "grey50", size = 2) +
    geom_text(data = lbl.X.prj, # longitude
              aes(x = X.prj, y = Y.prj, label = lbl), 
              colour = "grey50", size = 2)
}
