fetchRingData<-function(doSpp,surv=TRUE){
  
  
  ##############################################################
  #---------------- load data for growth 
  ##############################################################
  if(surv==FALSE){ 
  
  # load growth data
  groD=read.csv(file="growDnoNA.csv")
  D=groD[groD$allEdge==0,]; 
  D$logarea.t0=log(D$area.t0)
  D$logarea.t1=log(D$area.t1)
  D$quad=as.character(D$quad)
  D=D[order(D$X),]
  
  ##########################################################
  # Read in data on neighbors 
  ##########################################################
  ringD <- read.csv(paste(doSpp,"_nbhood_rings.csv",sep=""))
  ringD$year<-as.factor(ringD$year)
  
  # merge D with ringD (D contains fewer rows)
  D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
  D=D[order(D$X),]
  row.names(D) <- NULL  
  
  #-- ANNULUS DATA: all species 
  allCols= which(substr(names(D),1,4)%in%sppList); 
  propCols= which(substr(names(D),1,4)=="prop"); 
  sppData<-data.frame(logarea.t1=D$logarea.t1,quad=as.factor(D$quad),year=D$year,ID=D$trackID,age=D$age,
            logarea.t0=D$logarea.t0,Group=as.factor(D$Group),as.matrix(D[,allCols]),as.matrix(D[,propCols])) 
  
  ##############################################################
  #---------------- or, load data for survival 
  ##############################################################
    }else{
    
    survDfile="survD.csv"
    survD=read.csv(file=survDfile)
    D=survD[survD$allEdge==0,];
    D$year=D$year
    D$logarea=log(D$area)
    D$quad=as.character(D$quad)
    D1=D=D[order(D$X),];
    
  ##########################################################
  # Read in data on neighbors 
  ##########################################################
    ringD <- read.csv(paste(doSpp,"_nbhood_rings.csv",sep=""))
    ringD$year<-ringD$year
    
    # merge D with ringD (D contains fewer rows)
    D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
    D=D[order(D$X),]
    row.names(D) <- NULL  
    
    # -- ANNULUS DATA: all species  
    allCols= which(substr(names(D),1,4)%in%sppList); 
    propCols= which(substr(names(D),1,4)=="prop"); 
    sppData<-data.frame(survives=D$survives,age=D$age,ID=D$trackID, year=D$year, logarea=D$logarea, 
            Group=as.factor(D$Group), seedling=D$seedling, quad=D$quad, as.matrix(D[,allCols]),as.matrix(D[,propCols]));  
    }

return(sppData)
}

