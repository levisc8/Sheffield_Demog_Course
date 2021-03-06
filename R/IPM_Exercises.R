## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM to illustrate the construction of an IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
set.seed(53241986)


source('C:/Users/sl13sise/Dropbox/MLU/Knight Lab/Sheffield_Demog_Course/R/WD_Paths.R')
## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","C:/Users/sl13sise","~"); 
setwd(paste(root,"/Dropbox/ATSC 2018 participant folder/23.1.18/Rees/R code and Data",sep="")); 

source("Monocarp Demog Funs.R");

source("Standard Graphical Pars.R");

source('MatrixImage.R')

# Set simulation parameters
init.pop.size <- 250
n.yrs <-100

source("Monocarp Simulate IBM.R") 
cat(pop.size.t,"\n")

## trim an initial transient off the simulation
sim.data <- sim.data[sim.data$yr > 10,]
sim.data$yr <- sim.data$yr-10


# Extract 1000 observations to use as our data set

sample.index <- sample(1:nrow(sim.data), size=1000, replace=FALSE)
sim.data <- sim.data[sample.index,]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Fit statistical models to simulated data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Let's use sim data to construct an IPM, assuming we know the order of events
## first a bit of house keeping for plotting - chop up the sizes into classes and sort the data by size

sim.data <- transform(sim.data, z.classes=cut(z, 20))
sim.data <- sim.data[order(sim.data$z),]

##
## fit the functions and plot some graphs
##

set_graph_pars("panel4"); 

## 1 - growth

plot(z1 ~ z, data = sim.data,
     pch=16, cex=0.25,
     xlab=expression("Size t, "*italic(z)),
     ylab=expression("Size t+1, "*italic(z)*minute))

mod.Grow <- lm(z1 ~ z, data = sim.data)

abline(mod.Grow, col="red")

add_panel_label("a")

## 2 - flowering

Repr.ps <- summaryBy(z + Repr ~ z.classes, data = sim.data)

plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, 
xlab=expression("Size t, "*italic(z)),ylab="Probability of flowering")

mod.Repr <- glm(Repr ~ z, family = binomial, data = sim.data)

lines(fitted(mod.Repr) ~ z, data = sim.data, col = "red")

add_panel_label("b")

## 3 - survival

sim.data.noRepr <- subset(sim.data, Repr==0)

surv.ps <- summaryBy(z + Surv ~ z.classes, data = sim.data.noRepr, na.rm = TRUE)

plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19,
xlab=expression("Size t, "*italic(z)),ylab="Probability of survival")

mod.Surv <- glm(Surv ~ z, family = binomial, data = sim.data.noRepr)

lines(fitted(mod.Surv) ~ z, data = sim.data.noRepr, col = "red")

add_panel_label("c")

## 4 - seed production

sim.data.Repr <- subset(sim.data, Repr==1)

plot(log(Seeds) ~ z, data = sim.data.Repr, pch = 19,
xlab=expression("Size t, "*italic(z)), ylab="Seeds production")

mod.Seeds <- glm(Seeds ~ z, family=poisson, data=sim.data.Repr)

abline(mod.Seeds,col="red")

add_panel_label("d")

# dev.copy2eps(file="~/Repos/ipm_book/c2/figures/OenotheraDemog.eps");

## 5 - recruit size

sim.data.Rec <- subset(sim.data, age==0)

mod.Rcsz <- lm(z ~ 1 , sim.data.Rec)


## 6 - establishment probability
## WE SHOULD USE sim.data TO DO THIS CALC

p.r.est <- as.numeric(Recr)/sum(Seeds, na.rm=TRUE)

##
## Finally, store the estimated parameters
##

m.par.est <- c(surv    = coef(mod.Surv),
               flow    = coef(mod.Repr),
               grow    = coef(mod.Grow),
               grow.sd = summary(mod.Grow)$sigma,
               rcsz    = coef(mod.Rcsz),
               rcsz.sd = summary(mod.Rcsz)$sigma,
               seed    = coef(mod.Seeds),
               p.r     = p.r.est)

# par.names <- names(m.par.est)
# names(m.par.est) <- sub("(Intercept)", "int", par.names, fixed=TRUE)

names(m.par.est) <- names(m.par.true)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Construct Kernels and projection population size
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nBigMatrix <- 250;
min.size <- with(sim.data, min(z))
min.size
max.size <- with(sim.data, max(z))
max.size

#so let's set the lower and upper limits for the size range at -2.65 and 4.5
#so slightly smaller/bigger than observed

IPM.true <- mk_K(nBigMatrix, m.par.true, -2.65, 4.5)

IPM.est <- mk_K(nBigMatrix, m.par.est, -2.65, 4.5)

LambdaTrue <- Re(eigen(IPM.true$K)$values[1])

LambdaEstimated <- Re(eigen(IPM.est$K)$values[1])

LambdaTrue
LambdaEstimated


# Q1: Vary integration steps. Using a for loop to calculate
# how lambda converges

MeshpSeq <- round(seq(50, 1000, length.out = 25))
Output <- list(nMeshpts = MeshpSeq,
               LambdaEst = rep(NA, 25))



par(mfrow = c(1,1))

for(i in seq_len(25)) {
  # Iterate over sequence of meshpoint #s
  nMeshP <- MeshpSeq[i]
  
  # Create new IPM
  IPM.It.est <- mk_K(nMeshP, m.par.est, -2.65, 4.5)
  
  Lambda.est <- Re(eigen(IPM.It.est$K)$values[1])
  Output$LambdaEst[i] <- Lambda.est
}

Output <- data.frame(Output)

plot(LambdaEst ~ nMeshpts, data = Output,
     xlab = '# of Meshpoints',
     ylab = 'Estimated Lambda')
abline(h = LambdaTrue, col = 'red',
       lty = 3)

# Looks like the estimated lambdas converge around 250 - 300 meshpts.
# It really isn't varying very much at all though. From Mark:
# This will matter when the standard deviation around a given
# vital rate function is very small (i.e. your recruit size
# distribution is super narrow). At that point, you need
# very small bins and probably should use sparse matrices
# to not kill your computer

fit.pop.growth <- lm(log(pop.size.t)~c(1:yr))

exp(coef(fit.pop.growth)[2])

meshpts <- IPM.true$meshpts

w.est <- Re(eigen(IPM.est$K)$vectors[,1])
stable.z.dist.est <- w.est / sum(w.est)
mean.z.est <- sum(stable.z.dist.est*meshpts)
mean.z.est # Mean size of stable size distribution

w.true <- Re(eigen(IPM.true$K)$vectors[,1])
stable.z.dist.true <- w.true / sum(w.true)
mean.z.true <- sum(stable.z.dist.true*meshpts)
mean.z.true 

wb.est <- p_bz(meshpts,m.par.est)*w.est
stable.flowering.dist.est <- wb.est/ sum(wb.est)
mean.flowering.z.est <- sum(stable.flowering.dist.est*meshpts)
mean.flowering.z.est


wb.true <- p_bz(meshpts,m.par.true)*w.true
stable.flowering.dist.true <- wb.true/ sum(wb.true)
mean.flowering.z.true <- sum(stable.flowering.dist.true*meshpts)
mean.flowering.z.true

matrix.image(IPM.true$K ^ 0.1, xlab = 'Size (T)', ylab = 'Size (T+1)',
             do.contour = TRUE, do.legend = TRUE)
matrix.image(IPM.est$K ^ 0.1, xlab = 'Size (T)', ylab = 'Size (T+1)',
             do.contour = TRUE, do.legend = TRUE)

set_graph_pars("panel4"); 

## 1 - plot population density versus time...

plot(1:yr, log(pop.size.t[1:yr]), type="l",xlab="Time",ylab="Population size")
mod.pop.growth <- lm(log(pop.size.t) ~ seq.int(1,yr))
abline(mod.pop.growth, col="blue")

add_panel_label("a")

## ...roughly linear for log(Nt) vs time so exponential growth


## 2 - plot mean size versus time...
plot(1:yr, mean.z.t[1:yr], type="l",xlab="Time",ylab="Mean plant size")
abline(h=mean.z.true, col="red")
add_panel_label("b")

## ...mean size seems to settle down after a while

## 3 - plot mean flowering size versus time...
plot(1:yr, mean.fl.z.t[1:yr], type="l",xlab="Time",ylab="Mean flowering plant size")
abline(h=mean.flowering.z.true, col="red")
add_panel_label("c")

## ...mean flowering size seems to settle down after a while

## 4 - plot of density estimates at time 50 and the end 
plot(density(sim.data$z),main="",ylim=c(0,0.4), xlab="Plant size")
#lines(density(sim.data$z))
lines(IPM.est$meshpts,stable.z.dist.est/diff(IPM.est$meshpts)[1],col="red")
add_panel_label("d")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## End of code
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Q2: Refit survival model using a GAM, modify s_z() and recalculate
# lambda

source('C:/Users/sl13sise/Dropbox/MLU/Knight Lab/Sheffield_Demog_Course/R/IPM_Exercises_Functions.R')

library(mgcv)

NewSurvModel <- gam(Surv ~ s(z), family = binomial(),
                    data = sim.data)

summary(NewSurvModel)


xx <- seq(min(sim.data$z, na.rm = TRUE),
          max(sim.data$z, na.rm = TRUE),
          length.out = length(NewSurvModel$fitted.values))

plot(Surv.mean ~ z.mean, data =surv.ps)
lines(xx,
      predict(NewSurvModel, data.frame(z = xx), type = 'response'))


NewIPMEst <- mk_K(250,
                  m.par.est, 
                  NewSurvModel,
                  -2.65, 4.5)

matrix.image(NewIPMEst$K ^ 0.1,
             main = 'Estimated IPM',
             xlab = 'size (T)',
             ylab = 'size (T+1)')
matrix.image(IPM.true$K ^ 0.1, 
             main = 'True Data',
             xlab = 'size (T)',
             ylab = 'size (T+1)')

# dev.copy2eps(file="~/Repos/ipm_book/c2/figures/OenotheraSim.eps");

# Q3: Calculate lambda and W by iteration and compare results 
# to eigen
meshpts <- IPM.true$meshpts

# calculate bin width
h <- diff(meshpts)[1]

# assume a uniform distribution of size frequencies
w <- rep(1/(nBigMatrix*h),nBigMatrix)



for(i in 1:100){

  # iteratively show that w converges to the dominant
  # eigenvector
	w <- NewIPMEst$K %*% w
	
	# ibid, but with lambda/dominant eigenvalue
	lambda <- sum(w*h)
	
	# normalize and account for integration
	w <- w/sum(w*h)
	
}

cat("lambda by iteration:  ",
    lambda,
    " \ndominant eigenvalue: ",
    Re(eigen(IPM.est$K)$values[1]),
    '\n')

plot(w/sum(w),stable.z.dist.est)
abline(0,1,col="red")

plot(w/sum(w))

# Q4: Calculate variance in flowering sizes
wb.true <- p_bz(meshpts,m.par.true)*w.true
stable.flowering.dist.true <- wb.true / sum(wb.true)

# This is mean of log size
mean.flowering.z.true <- sum(stable.flowering.dist.true*meshpts)
mean.flowering.z.true

# mean of actual size is ...
# mean.flowering.z.true <- sum(stable.flowering.dist.true* exp(meshpts))
# mean.flowering.z.true

# Variance = (E(x^2) - E(x)) ^2
# x =  (meshpts) * stable.flowering.dist.true
# Variance attempt:
var <- sum(stable.flowering.dist.true * meshpts^2) - mean.flowering.z.true ^ 2
sd <- sqrt(var)

plot(stable.flowering.dist.true ~ meshpts, xlim = c(-3, 5))
abline(v = mean.flowering.z.true, col = 'red', lty = 2)
abline(v = mean.flowering.z.true - 2*sd, col = 'blue', lty = 3)
abline(v = mean.flowering.z.true + 2*sd, col = 'blue', lty = 3)

#calculate lambda for a range of integration steps

ns      <- seq(10,500,by=10)
lambdas <- rep(NA,length(ns))

for(i in 1:length(ns)){
	
	IPM.est     <- mk_K(ns[i], m.par.est, -2.65, 4.5)
	lambdas[i]  <- Re(eigen(IPM.est$K,only.values=TRUE)$values[1]) 
	
}

set_graph_pars("panel2");

plot(ns,lambdas,type="l",xlab="Number integration steps",ylab="lambda")

error <- (lambdas-lambdas[length(ns)])/lambdas[length(ns)]*100

plot(ns, error,xlab="Number integration steps",ylab="Error",type="l")

round(lambdas,3)
round(error,4)



nBigMatrix <- 100;
IPM.est   <- mk_K(nBigMatrix, m.par.est,
                  survModel = NewSurvModel,
                  -2.65, 4.5)

P <- IPM.est$P
F <- IPM.est$F
N <- solve((diag(nBigMatrix)-P))
R <- F %*% N


R0 <- Re(eigen(R,only.values=TRUE)$values[1])

R0.exp <- rep(NA,50)

sum <- diag(nBigMatrix)

Pa <- P


for(i in 1:50){
	
	sum <- sum + Pa
	
	R <- F %*% sum
	
	R0.exp[i] <- Re(eigen(R,only.values=TRUE)$values[1])
	
	Pa <- Pa %*% P
	
	cat("percentage error ",100*(R0.exp[i]-R0)/R0,"\n")
	
}
	

#refit models with GAM
library(mgcv)

mod.Surv.gam <- gam(Surv ~ s(z), family = binomial, data = sim.data.noRepr)

summary(mod.Surv.gam)

s_z <- function(z, m.par)
{
    p <- predict(mod.Surv.gam,newdata=data.frame(z=z),type="response")       # logistic transformation to probability
    return(p)
}

IPM.est <- mk_K(nBigMatrix, m.par.est, -2.65, 4.5)
Re(eigen(IPM.est$K)$values[1])







