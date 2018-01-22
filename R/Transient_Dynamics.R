# Transient dynamics w/ popdemo
# Based on vignette from Iain Stott

rm(list = ls())

# load in the desert tortoise data set
library(popdemo)
data("Tort")


# Create a starting population vector and scale it to unity
TortVec <- runif(8)
TortVec <- TortVec / sum(TortVec)


TortProj <- project(Tort, TortVec, time = 100)

plot(TortProj)

# This time, return the the projected size distributions
# as well
TortProjVecs <- project(Tort, 
                        TortVec, 
                        time = 100,
                        return.vec = TRUE)

plot(TortProjVecs$vec[,1],
     ylim = c(0,1.1),
     type = 'l',
     ylab = 'Population Size/Density',
     xlab = 'Projection Year',
     main = 'Stages vs time')
legtext <- c(attr(TortProjVecs$vec, 'dimnames')[[2]][1])
legcol <- c(1)
legty <- 1

for(i in 2:dim(TortProjVecs$vec)[2]){
  
  lines(TortProjVecs$vec[ ,i],
       col = i, lty = i)
  legtext <- c(legtext, attr(TortProjVecs$vec, 'dimnames')[[2]][i])
  legcol <- c(legcol, i)
  legty <- c(legty, i)
}

legend('topright', legtext, lty = legty, col = legcol)

# On to section 4 of the vignette!
ss <- eigs(Tort, 'ss')

# Standardizing matrices so we can compare across
# populations. These are called transient indices. One standardization
# is to standardize to the initial population vector so it sums to 1.
# We actually did that in the first run above. 

# Another way to standardize transient dynamics to re-scale matA so that
# it's dominant eigenvalue = 1. You can do this by dividing the 
# matA by lambda
TortP1 <- project(Tort, TortVec, 
                  time = 100,
                  standard.A = TRUE,
                  standard.vec = TRUE)

TortSS <- project(Tort, ss,
                  time = 100,
                  standard.A = TRUE,
                  standard.vec = TRUE)
plot(TortP1, log = 'y')
lines(TortSS, lty = 2)

# Inertia and reactivity are two measures of transient dynamics.
# Reactivity measures the population growth rate between T=0 and T=1
# relative to the prediction by the asymptotic population growth rate.
# Inertia does the same thing, but at T = infinity

react <- reac(Tort, TortVec)
inert <- inertia(Tort, TortVec)

points(c(1, 100), c(react, inert), pch = 3, col = 'blue')

# Next, create vectors that amplify and attenuate transient dynamics
TortAmp <- c(1,1,2,3,5,8,13,21)
TortAtt <- rev(TortAmp)

TortVec3 <- cbind(RAND = TortVec,
                  ATT = TortAmp,
                  ATT = TortAtt)

# Do the projections for each starting vector
TortProj3 <- project(Tort, 
                     TortVec3,
                     time = 100,
                     standard.A = TRUE,
                     standard.vec = TRUE)

plot(TortProj3, log = 'y')
lines(TortSS, col = 'blue', lty = 2)

# Calculate reactivity and inertia for our new 
# starting vectors
reacAtt <- reac(Tort, TortAtt)
reacAmp <- reac(Tort, TortAmp)
InertAtt <- inertia(Tort, TortAtt)
InertAmp <- inertia(Tort, TortAmp)

# add all these points to the plot
points(c(rep(1, 3), rep(100, 3)),
       c(react, reacAtt, reacAmp,
         inert, InertAtt, InertAmp),
       pch = 3, col = 'red')

reacN <- reac(Tort, TortVec, return.N = TRUE)

# Calculate Kreiss bounds for the matrix.
# these tell you the maximum and minimum values for inertia and 
# reactivity

rlwr <- reac(Tort, bound = 'lower')
rupr <- reac(Tort, bound = 'upper')
ilwr <- inertia(Tort, bound = 'lower')
iupr <- inertia(Tort, bound = 'upper')

plot(TortProj3, bounds = TRUE)
points(c(rep(1, 5), rep(100, 5)),
       c(react, reacAmp, reacAtt, rlwr, rupr,
         inert, InertAmp, InertAtt, ilwr, iupr),
       pch = 3, col = 'blue')
# Another way of plotting this is based on the likelihood of a given
# initial vector. Apparently, project(vector = 'diri') does this for you.
# Example below

TortDiri <- project(Tort, vector = 'diri', time = 100,
                    standard.A = TRUE)
plot(TortDiri, plottype = 'shady', log = 'y', bounds = TRUE)
lines(1:101, TortSS, type = 1, col = 'blue')
