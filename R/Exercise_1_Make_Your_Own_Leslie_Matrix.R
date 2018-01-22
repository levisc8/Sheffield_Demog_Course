# Exercise 1. Create an age-based matrix model for an organism.
# I am using C. edulis with a seedbank and 4 age classes.
# I am going to make some assumptions about vital rates
# and limiting the seedbank to 1 year

rm(list = ls())


library(popbio)

Gi <- 0.1 
Gsb <- 0.2
# Induce some mortality with seedlings. Seems like this happens in nature
S0 <- 0.4
# Decrease mortality with age
S1 <- 0.55
S2 <- 0.75
S3 <- 0.999
# This is tricky to decide on. My guess is they make way more than this
# each year, but we'll see what it looks like
Fec <- 50

# N0 <- c(50, 20, 30, 25, 10)
N0 <- c(15, 5, 4, 3, 2)
LeslieMatA <- matrix(c(0, 0, 0, 0, Fec * (1 - Gi),
                 Gsb, 0, 0, 0, Fec * Gi,
                 0, S0, 0, 0, 0,
                 0, 0, S1, 0, 0,
                 0, 0, 0, S2, S3),
               nrow = 5,
               byrow = TRUE)

LeslieOutput <- eigen.analysis(LeslieMatA)

# Iterative analysis of Leslie matrix models
LeslieAges <- matrix(0, 50, 5)

LeslieAges[1, ] <- N0

for(i in 2:dim(LeslieAges)[1]) {
  
  LeslieAges[i, ] <- LeslieMatA %*% LeslieAges[i - 1, ]
}
  
# Next, create a Lefkovitch stage matrix with no restrictions on
# classes that you can move between
LefkovitchMatA <- matrix(c(0, 0, 0, 0, Fec * (1 - Gi),
                           Gsb, 0, 0, 0, Fec * Gi,
                           0, 0.2, 0, 0, 0, 
                           0, 0.1, 0.4, 0.35, 0.39,
                           0, 0.1, 0.15, 0.4, 0.6),
                         nrow = 5, 
                         byrow = TRUE)
colSums(LeslieMatA)
colSums(LefkovitchMatA)

LefkovitchOutput <- eigen.analysis(LefkovitchMatA)

# Iterative analysis of Lefkovitch matrix models
LefkovitchAges <- matrix(0, 50, 5)

LefkovitchAges[1, ] <- N0

for(i in 2:dim(LefkovitchAges)[1]) {
  
  LefkovitchAges[i, ] <- LefkovitchMatA %*% LefkovitchAges[i - 1, ]
}

par(mfrow = c(1, 2))
matplot(log(LeslieAges[1:15, ]),
        type = 'l',
        main = 'Leslie Matrix')
matplot(log(LefkovitchAges[1:15, ]),
        type = 'l', 
        main = 'Lefkovitch Matrix')
