
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

# Perturbation analysis!!!!!!!!
# Sam Levin 1/25/18

# Exercise 1. Full kernel perturbation
Pert <- 1.0001

IpmKMat <- IPM.true$K
h <- IPM.true$meshpts[2] - IPM.true$meshpts[1]

Lambda <- Re(eigen(IpmKMat)$values[1])
LeftEV <- Re(eigen(t(IpmKMat))$vectors[ ,1])
RightEV <- Re(eigen(IpmKMat)$vectors[ ,1])
 
SensKmat <- outer(LeftEV, RightEV, '*') / sum(LeftEV * RightEV * h)
ElasKmat <- ((IpmKMat / h) / Lambda) * SensKmat
sum(ElasKmat) * h ^ 2

matrix.image(ElasKmat, xlab = 'Size (t)', ylab = 'Size (t+1)',
             do.contour = TRUE)

# Exercise 2: Subkernel perturbation
Pmat <- IPM.true$P
Fmat <- IPM.true$F

ElasPMat <- ((Pmat / h) / Lambda) * SensKmat
ElasFMat <- ((Fmat / h) / Lambda) * SensKmat

sum(sum(ElasPMat), sum(ElasFMat)) * h ^ 2
sum(ElasPMat) * h^2
sum(ElasFMat) * h^2

# Exercise 3: parameter level perturbation
# Create vector to store elasticity values
ElasParams <- m.par.true
ElasParams[] <- NA

for(k in seq_len(length(m.par.true))) {
  mNew <- m.par.true
  mNew[k] <- mNew[k] + abs(mNew[k]) * 0.0001
  IpmNew <- mk_K(nBigMatrix, mNew, -2.65, 4.5)$K
  LambdaNew <- Re(eigen(IpmNew)$values[1])
  ElasParams[k] <- ((LambdaNew - Lambda) / Lambda) / 0.0001
}

par(mar = c(4, 7, 1, 1))
barplot(ElasParams, horiz = TRUE,
        las = 1)

# Steve: These are telling you how the uncertainty in your parameter estimates
# are influencing your estimates of population growth. Higher elasticities 
# may indicate that you need more power to fit an adequate vital rate 
# regression
