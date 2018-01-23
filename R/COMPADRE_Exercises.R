# Compadre Exercises for Sheffield Demography Course

rm(list = ls())

# Script that contains the working directories we need.
# Load the data
source('R/WD_Paths.R')

setwd(data.wd)

load('23.1.18/Big data/COMPADRE_v.4.0.1.RData')

setwd(old.wd)

# Remove negative matrices. See lines 17-28 of: 
# https://github.com/jonesor/Rcompadre/blob/master/data-raw/createDataSets.R

NegativeMatrices <- logical(length(compadre$mat))
for(i in seq_len(dim(compadre$metadata)[1])) {
  NegativeMatrices[i] <- any(compadre$mat[[i]]$matA < 0 |
                               compadre$mat[[i]]$matU < 0,
                             na.rm = TRUE)
}

# 24 bad ones. Someone should bug Eelke about having one of them
length(which(NegativeMatrices == TRUE))

# Remove negative matrices from Compadre
compadre <- list(metadata = compadre$metadata[!NegativeMatrices, ],
                 matrixClass = compadre$matrixClass[!NegativeMatrices],
                 mat = compadre$mat[!NegativeMatrices],
                 version = compadre$version)

# Create a sample data set. Using old version of Rage (Mage) and
# skipping the currently rather buggy RCompadre package
library(Mage)
library(popbio)
CommonTrees <- subsetDB(compadre, Genus == 'Quercus' | 
                   Genus == 'Fagus')

# Plot a life cycle for a random tree
selector <- base::sample(1:dim(CommonTrees$metadata)[1],
                         1)

MatA <- CommonTrees$mat[[selector]]$matA
plotLifeCycle(MatA, title = CommonTrees$metadata$SpeciesAccepted[selector])

# Obtain lambda, inertia, and reactivity for each species
library(popdemo)
Output <- list(Species = CommonTrees$metadata$SpeciesAccepted,
               Lambda = rep(NA, dim(CommonTrees$metadata)[1]),
               GenTime = rep(NA, dim(CommonTrees$metadata)[1]),
               Lat = CommonTrees$metadata$Lat,
               Lon = CommonTrees$metadata$Lon)

for(i in seq_len(dim(CommonTrees$metadata)[1])) {
  
  # Set up A matrix and population vector
  MatA <- CommonTrees$mat[[i]]$matA
  PopVec <- runif(dim(MatA)[1])
  PopVec <- PopVec / sum(PopVec)
  
  Output$Lambda[i] <- max(Re(eigen(MatA)$values))
  Output$GenTime[i] <- generation.time(MatA)
}

hist(Output$Lambda)
hist(Output$GenTime)

# Create a map of the geographic distribution
source('R/VS_World_Maps.R')

Output <- data.frame(Output[c(1:2, 4:5)])
Coords <- project(cbind(Output$Lon, Output$Lat),
                  proj =  "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" ) %>%
 data.frame() %>%
 cbind(., Output$Lambda, Output$Species) %>%
  setNames(c("Lon", "Lat", "Lambda", "Species"))

VSWorldMap() + 
  theme_void() + 
  geom_point(data = Coords,
             aes(x = Lon, 
                 y = Lat,
                 color = Lambda,
                 shape = Species))

# Age from stage exercise
# First, isolate the matrix of interest

MatU <- CommonTrees$mat[[3]]$matU
MatF <- CommonTrees$mat[[3]]$matF

test <- makeLifeTable(MatU, matF = MatF, nSteps = 10)



  


