#########################################################################
### Exploratory data analysis of survival 
#########################################################################
rm(list=ls(all=TRUE)); graphics.off();
require(mgcv); 

old.wd <- 'C:/Users/sl13sise/Dropbox/MLU/Knight Lab/Sheffield_Demog_Course/'
# 
# ## Provide directory information
# root = "~"; # This should work on all Macs
# # root = "C:/Users/me"; #something like this should work on a PC 
# 
# ## set working direcetory to the right folder 
# setwd(root); 
# setwd("Dropbox/ATSC 2018 participant folder/24.1.18/Density Dependence") 

setwd(paste0(old.wd, 'R/Density_Dependence'))


## Load the function that reads in the data
source("LoadData.R") 

## Pick which species to analyze 
sppList = c("ARTR","HECO","POSE","PSSP"); 
doSpp = "ARTR"; 

################  Read in the data #######################
setwd(paste0(getwd(), '/',doSpp)); # move into the data folder 
data<-fetchRingData(doSpp=doSpp,surv=TRUE)  # if SURV=FALSE, you get growth data 
setwd(".."); # move back up to script folder 

data$year<-as.factor(data$year)
distances=c(seq(2,20,by=2),seq(25,50,by=5), seq(60,150, by = 10)); # outer radius of first 16 rings 
View(data); 

################################################################### 
# Fit a null model without competition 
###################################################################
fitNull=glm(survives~logarea+Group,data=data,family="binomial"); 

############################################################################### 
# Fit a model with competition limited to conspecifics only, first 5 rings 
###############################################################################
# which columns are conspecific area? 
ARTRcols=which(substr(names(data),1,4)=="ARTR"); 

# pull them out 
ARTRrings = data[,ARTRcols]; 
# View(ARTRrings); 

####### compute percent cover in the first 5 rings 
k=5; ARTR.k = ARTRrings[,1:k]; 
if(k>1) {ARTR.k = rowSums(ARTR.k)} 
ARTR.k=100*ARTR.k/(pi*distances[k]^2); # convert to percent cover 
hist(ARTR.k); 

###### use that as a covariate in predicting survival 
fit5a=glm(survives ~ logarea + Group + ARTR.k, data=data,family="binomial"); 
fit5b=glm(survives ~ logarea + Group + sqrt(ARTR.k), data=data,family="binomial"); 
AIC(fit5a,fit5b); # likes the sqrt 

# Is sqrt good? A spline should collapse to a line. And, it does. 
fit5c = gam(survives ~ logarea + Group + s(sqrt(ARTR.k)), data=data,family="binomial"); 
# plot(fit5c); 

##### Is competition important? 
summary(fit5b); #yes 
AIC(fitNull,fit5b); #yes 
plot(fitNull$fitted,fit5b$fitted,xlab="Predicted survival ignoring competition",
ylab="Predicted survival including competition");  
abline(0,1,col="blue",lty=2); 


# Q 1 from DensityDependence Slideshow
# allocate memory
output <- list(Dist = distances[-1],
               AIC = rep(NA, length(ARTRcols) -1))

# Loop across possible combinations

for(i in 2:length(ARTRcols)) { 
  
  
  kEnd <- i
  
  # Select rings based on iterator. convert to % cover
  ARTR.k = ARTRrings[ ,1:kEnd]
  ARTR.k = rowSums(ARTR.k)
  ARTR.k = 100*ARTR.k/(pi*distances[i]^2)
  
  # Fit your model!
  fit5d = glm(survives ~ logarea + Group + sqrt(ARTR.k),
              data=data,
              family="binomial")
  
  # Extract the AIC
  output$AIC[i] <- AIC(fit5d)  
  
  
}

PlotData <- data.frame(Dist = output$Dist[!is.na(output$Dist)],
                       AIC = output$AIC[!is.na(output$AIC)])


par(mfrow = c(2,2))
plot(AIC ~ Dist, data = PlotData,
     main = 'ARTR',
     type = 'l')

# Look at the actual radius
PlotData$Dist[which.min(PlotData$AIC)]

# Now, we'll add in other species one by one. This is Q2A and B
ARTR.IDX <- which.min(PlotData$AIC)

BestARTRCols <- ARTRcols[1:ARTR.IDX]
ARTR.k = ARTRrings[ ,BestARTRCols]
ARTR.k = rowSums(ARTR.k)
ARTR.k = 100*ARTR.k/(pi*distances[ARTR.IDX]^2)

AddSpp <- sppList[-c(1)]

rm(i)
for(Spp in unique(AddSpp)) { 
  
  output <- list(Dist = distances[-1],
                 AIC = rep(NA, 25))
  
  # Find the columns that contain the data for the species we want to add
  # to the model
  SppCols <- which(grepl(Spp, names(data)))
  
  # pull out the data
  SppData <- data[ ,SppCols]
  
  for(i in 2:26) {
    kEnd <- i
    
    # Select rings based on iterator. convert to % cover
    Spp.k = SppData[ ,1:kEnd]
    Spp.k = rowSums(Spp.k)
    Spp.k = 100*Spp.k/(pi*distances[i - 1]^2)
    
    # Fit your model!
    fit5f = glm(survives ~ logarea + Group + sqrt(ARTR.k)  + sqrt(Spp.k),
                data=data,
                family="binomial")
    
    # Extract the AIC
    output$AIC[i-1] <- AIC(fit5f)
    
  }
  # Plot the data and move on to the next species
  PlotData <- as.data.frame(output)
  
  plot(AIC ~ Dist, data = PlotData,
       main = paste0('ARTR + ', Spp),
       type = 'l')
  
}

# Q 3: 

First6Rings <- which(str_detect( names(data), ".0.2|2.4|4.6|6.8|8.10|10.12"))
badCols <- c(19, 31, 45, 57, 71, 83, 97, 109, 123, 135)

First6Rings <- First6Rings[!First6Rings %in% badCols]

RingData <- data[ ,First6Rings]

for(Spp in unique(sppList)) {

  
  assign(paste0(Spp, '.k'), RingData[ ,which(grepl(Spp, names(RingData)))])
  assign(paste0(Spp, '.k'), rowSums(eval(parse(text = paste0(Spp, '.k')))))
  
}

ARTR.k <- 100*ARTR.k/(pi*12^2)
POSE.k <- 100 * ARTR.k / pi * 12 ^ 2
HECO.k <- 100 * ARTR.k / pi * 12 ^ 2
PSSP.k <- 100 * ARTR.k / pi * 12 ^ 2

BigModel <- glm(survives ~ logarea + Group + 
                  sqrt(ARTR.k)  + sqrt(POSE.k) + 
                  sqrt(HECO.k) + sqrt(PSSP.k),
                data=data,
                family="binomial")

AIC(BigModel, fit5b)


# Sequences of radii and alpha values to choose from
Radii <- c(seq(1, 19, 2), seq(22.5, 47.5, 5), seq(55, 145, 10))
alpha <- seq(1e1, 1e2, 1e1)

# Function to weight columns with exponential decay
expDecay <- function(r, a, col, squareR = FALSE) {
  if(squareR) {
    out <- exp(-a * r^2) * col
  } else {
    out <- exp(-a * r) * col
  }
  return(out)
}

# Initialize list to save outputs
OutData <- list()

for(spp in unique(sppList)){
  Data <- data[ ,grepl(spp, names(data))]
  for(i in alpha) {
    for(r in seq_len(length(Radii))){
      OutData[[paste(Spp, i, r, sep = '_')]] <- expDecay(Radii[r], i, Data[ ,r])
      
    }
  }
}

DecayData <- as.data.frame(OutData)

#####################################################################################
# For use in the Exercise: load the data about other-species neighbors of ARTR plants
#####################################################################################
HECOcols=which(substr(names(data),1,4)=="HECO"); 
HECOrings = data[,HECOcols]; 

POSEcols=which(substr(names(data),1,4)=="POSE"); 
POSErings = data[,POSEcols]; 

PSSPcols=which(substr(names(data),1,4)=="PSSP"); 
PSSPrings = data[,PSSPcols]; 



