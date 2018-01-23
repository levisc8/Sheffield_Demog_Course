# Stochastic MPMs Exercises

rm(list = ls())
library(popdemo)

# Load data
data(Pbear)

PbearVec <- c(0.106, 0.068, 0.106, 0.461, 0.151, 0.108)


# project() will do stochastic projections if you pass it a list
# of matrices. Didn't know that before, but good to know now!
PbearProject <- project(Pbear, PbearVec, time = 100)

plot(PbearProject, log = 'y')

# Next, we can implement this as Markov process. We create a 
# matrix of probabilities that represents the probability choosing
# choosing matrix Y at T+1 given that we selected matrix X at T.
# This employs a uniform matrix (i.e. every entry has the same value)
ASeq <- matrix(rep(0.2, 25), 5, 5)
PbearProjectM <- project(Pbear, PbearVec, Aseq = ASeq, time = 100)
plot(PbearProjectM, log = 'y')

# Next, construct a Markov chain matrix with unequal probabilities. 
# this makes it more or less likely to select a given matrix
# describing transitions based on how frequent we know
# good or bad year occur
p <- 0.5
ASeq2 <- matrix(rep(c((1 - p) / 3, (1 - p) / 3, (1 - p) / 3,
                      p / 2, p / 2), 5), 5, 5)

PbearProjMBad <- project(Pbear, PbearVec, Aseq = ASeq2,
                         time = 100)

plot(PbearProjMBad, log = 'y')

# Make them even worse for the polar bear. Climate change is projected
# to make bad years more frequent, so perhaps this is realistic
pVeryBad <- 0.8
ASeq2 <- matrix(rep(c((1 - pVeryBad) / 3, 
                      (1 - pVeryBad) / 3,
                      (1 - pVeryBad) / 3,
                      pVeryBad / 2, 
                      pVeryBad / 2),
                    5),
                5, 5)
PbearProjVBad <- project(Pbear, PbearVec, Aseq = ASeq2,
                         time = 100)
plot(PbearProjVBad)

# Now, we simulate a system in which there are feebacks. A good year for sea ice
# increases the probability that there will be another good year for ice cover
# in the future. On the other hand, a bad year for ice cover increases the 
# probability that another one will follow

pBadGood <- 0.2
pBadBad <- 0.8

ASeq <- matrix(c(rep(c((1 - pBadGood) / 3,
                       (1 - pBadGood) / 3,
                       (1 - pBadGood) / 3,
                       pBadGood / 2,
                       pBadGood / 2),
                     3),
                 rep(c((1 - pBadBad) / 3,
                       (1 - pBadBad) / 3, 
                       (1 - pBadBad) / 3,
                       pBadBad / 2,
                       pBadBad / 2), 
                     2)),
               5, 5)
# Now project!
PbearProjPosFB <- project(Pbear, PbearVec,
                          Aseq = ASeq, time = 100)
plot(PbearProjPosFB, log = 'y')

# Now, increase the probability of going from good to bad
pBadGood <- 0.5
ASeq <- matrix(c(rep(c((1 - pBadGood) / 3,
                       (1 - pBadGood) / 3,
                       (1 - pBadGood) / 3,
                       pBadGood / 2,
                       pBadGood / 2),
                     3),
                 rep(c((1 - pBadBad) / 3,
                       (1 - pBadBad) / 3, 
                       (1 - pBadBad) / 3,
                       pBadBad / 2,
                       pBadBad / 2), 
                     2)),
               5, 5)
# Now project!
PbearProjPosFB <- project(Pbear, PbearVec,
                          Aseq = ASeq, time = 100)
plot(PbearProjPosFB, log = 'y')

# Finally, use the matrices in the fixed sequence they occurred in
ASeq <- rep(1:5, 10)

PbearProjFixed <- project(Pbear, PbearVec, Aseq = ASeq, time = 100)
plot(PbearProjFixed, log = 'y')
