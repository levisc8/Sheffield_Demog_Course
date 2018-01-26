# Age x Size IPM
# Presentation available at: 
# http://rpubs.com/dzchilds/age-size-ipm-exercise

# Further details stored in :
# C:/Users/sl13sise/Dropbox/ATS 2018 participant folder/25.1.18/Soay_Age_Size/general_ipm.html

# Sam Levin
# 1/25/18

rm(list = ls())
library(car)
library(dplyr)
library(Matrix)

source('R/WD_Paths.R')

# Load data, reset working directory
setwd(data.wd)
data <- read.csv('25.1.18/Soay_Age_Size/soay_demog_data.csv',
                 stringsAsFactors = FALSE) %>%
  as_tibble()

source('23.1.18/Rees/R code and Data/MatrixImage.R')

setwd(old.wd)

# Variable explanations
# z = size(t) (log)
# a = age(?)
# Surv = survival
# z1 = size(t+1) (log)
# Repr = reproductive
# sex: 1 = female
# Recr = successful recruitment?
# Rcsz = recruit size if successful

# Growth!
GrowData <- filter(data, !is.na(z1))
GrowModA <- lm(z1 ~ z + a, 
              data = GrowData)
GrowModZ <- lm(z1 ~ z,
               data = GrowData)

# Survival
if(any(is.na(data$Surv))) {
  SurvData <- filter(data, !is.na(Surv))
} else {
  SurvData <- data
}

SurvModA <- glm(Surv ~ z + a, 
                family = binomial(),
                data = SurvData)

SurvModZ <- glm(Surv ~ z, 
                family = binomial(),
                data = SurvData)

# Reproduction probability
ReproData <- filter(data, Surv == 1 & a > 0)
ReproModA <- glm(Repr ~ z + a,
                 family = binomial(),
                 data = ReproData)

ReproModZ <- glm(Repr ~ z,
                 family = binomial(),
                 data = ReproData)

# Recruit survival probability
RecrData <- filter(data, !is.na(Recr))

RecrModA <- glm(Recr ~ z + a,
                data = RecrData,
                family = binomial())

RecrModZ <- glm(Recr ~ z,
                data = RecrData,
                family = binomial())

# Recruit size ~ parent size
RcszData <- filter(data, Recr == 1)

RcszModA <- lm(Rcsz ~ z + a,
               data = RcszData)

RcszModZ <- lm(Rcsz ~ z,
               data = RcszData)

# Model Testing
Anova(GrowModA)
Anova(SurvModA)
Anova(RecrModA)
Anova(ReproModA)
Anova(RcszModA)

# Most of these look like they're worth using both terms. Except 
# Recruitment. Looks like thats only a function of size, so we'll
# use a different Model
RecrMod <- glm(Recr ~ a, data = RecrData, family = binomial)

# Now extract coefficients and rename some of them
m_par <- c(surv    = coef(SurvModA),
           grow    = coef(GrowModA),
           grow_sd = summary(GrowModA)$sigma,
           repr    = coef(ReproModA),
           recr    = coef(RecrMod),
           rcsz    = coef(RcszModA),
           rcsz_sd = summary(RcszModA)$sigma)

names(m_par) <-  sub(pattern = "(Intercept)", 
                     replace = "int", 
                     x = names(m_par),
                     fixed = TRUE)

names(m_par) <-  sub(pattern = ".", 
                     replace = "_", 
                     x = names(m_par),
                     fixed = TRUE)

# Define functions for the IPM
g_z1z <- function(z, z1, a, m_par) {
  mu <- m_par['grow_int'] + 
    z * m_par['grow_z'] + 
    a * m_par['grow_a']
  out <- dnorm(z1, mean = mu, sd = m_par['grow_sd'])
  return(out) 
}

s_z <- function(z, a, m_par) {
  linear_mu <- m_par['surv_int'] + 
    z * m_par['surv_z'] + 
    a * m_par['surv_a']
  out <- 1 / (1 + exp(-linear_mu))
  return(out)
}

pb_z <- function(z, a, m_par) {
  if(a == 0) {
    out <- 0
  } else {
    linear_mu <- m_par['repr_int'] + 
      z * m_par['repr_z'] + 
      a * m_par['repr_a']
    out <- 1 / (1 + exp(-linear_mu))
  }
  return(out)
}

pr_z <- function(a, m_par) { 
  linear_mu <- m_par['recr_int'] + 
    a * m_par['recr_a']
  out <- 1 / (1 + exp(-linear_mu))
  return(out)
}  

c_z1z <- function(z, z1, a, m_par) { 
  mu <- m_par['rcsz_int'] + 
    z * m_par['rcsz_z'] +
    a * m_par['rcsz_a']
  out <- dnorm(z1, mean = mu, sd = m_par['rcsz_sd'])
  return(out)
}



# Now, define out kernels!
# We need functions to define them though. 
P_z1z <- function(z, z1, a, m_par) {
  out <- s_z(z, a, m_par) * g_z1z(z, z1, a, m_par)
  return(out)
}

F_z1z <- function(z, z1, a, m_par) {
  out <- s_z(z, a, m_par) *
    pb_z(z, a, m_par) *
    (1/2) * 
    pr_z(a, m_par) * 
    c_z1z(z, z1, a, m_par)
  return(out)
}

# Create parameters describing each one. Below is my attempt at describing this,
# but Dylan has an alternative way. I'll construct those further down and try 
# to compare the two approaches to see how bad I am at doing this on my own

#system.time({
# nMeshP <- 200
# L <- 1.6 #min(data$z) * 0.9
# U <- 3.7 #max(data$z) * 1.1
# MeshPoints <- seq(L, U, length.out = nMeshP)
# BinWidth <- MeshPoints[2] - MeshPoints[1]
# Ages <- seq(0, 20, by = 1)
# 
# PMatList <- list()
# FMatList <- list()
# for(Age in Ages) {
#   GMat <- BinWidth * outer(MeshPoints,
#                            MeshPoints,
#                            FUN = g_z1z,
#                            a = Age,
#                            m_par = m_par)
#   GMat <- GMat / matrix(as.vector(apply(GMat,
#                                         2,
#                                         sum)),
#                         nrow = nMeshP,
#                         ncol = nMeshP,
#                         byrow = TRUE)
#   
#   PMatList[[Age + 1]] <- GMat %*% diag(s_z(MeshPoints,
#                                            a = Age, 
#                                            m_par = m_par))
#   
#   FMatList[[Age + 1]] <- F_z1z(MeshPoints,
#                                MeshPoints,
#                                a = Age,
#                                m_par = m_par)
#   
# }
#})
# Now, Dylan's implementation!
# Calculate the mesh points, mesh width and store with upper/lower 
# bounds and max age

#system.time({
mk_intpar <- function(m, L, U, M) {
  h <- (U - L) / m
  meshpts  <-  L + ((1:m) - 1/2) * h
  na <- M + 2
  return( list(meshpts = meshpts, M = M, na = na, h = h, m = m) )
}


# make initial state distribution (uniform size, age 0 only)
mk_init_nt0 <- function(i_par) {
  nt <- lapply(seq_len(i_par$na), function(x) rep(0, i_par$m))
  nt[[1]] <- rep(1, i_par$m) / (i_par$h * i_par$m)
  return(nt)
}

# Build the list of age/process specific kernels + store the 
# integration parameters in the same list
mk_age_IPM <- function(i_par, m_par) {
  within(i_par, {
    F <- P <- list()
    for (ia in seq_len(na)) {
      F[[ia]] <- outer(meshpts, meshpts, F_z1z, a = ia-1, m_par = m_par) * h
      P[[ia]] <- outer(meshpts, meshpts, P_z1z, a = ia-1, m_par = m_par) * h
    }
    rm(ia)
  })
}

# Create structural parameters
StructPars <- mk_intpar(200, # Number of meshpoints
                        1.6, # smallest sizes
                        3.7, # largest sizes
                        M = 20) # maximum allowable age

init_n0 <- mk_init_nt0(i_par = StructPars)

IpmKernels <- mk_age_IPM(i_par = StructPars,
                         m_par = m_par)


#})
# System.time results: 
# Dylan: 0.53 second
# Sam: 0.27 seconds
# Haha!

# Now right and left iteration. Right iteration creates the right
# eigenvector, left iteration creates the left eigenvector
r_iter <- function(x, F, P) {
  na <- length(x)
  xnew <- list(0)
  # multiply 
  for (ia in seq_len(na)) {
    xnew[[1]] <- (xnew[[1]] + F[[ia]] %*% x[[ia]])[,,drop=TRUE]
  }
  for (ia in seq_len(na-1)) {
    xnew[[ia+1]] <- (P[[ia]] %*% x[[ia]])[,,drop=TRUE]
  }
  xnew[[na]] <- xnew[[na]] + (P[[na]] %*% x[[na]])[,,drop=TRUE]
  return(xnew)
}

# See what the hell this thing is doing
FirstIt <- r_iter(x = init_n0, F = IpmKernels$F, P = IpmKernels$P)

# It's moved all individuals from our initial vector (which is just a uniform
# distribution of individuals at Age 0), on to the next year. So I guess
# we just run this for a lot of years and see what happens

YearSeq <- seq(1, 100, 1)
Output <- list()
Lambdas <- list()

for(i in YearSeq) {
  if(i == 1) {
    New_n0 <- r_iter(init_n0, 
                     F = IPM_sys$F, 
                     P = IPM_sys$P)
    Lambda <- sum(unlist(New_n0)) / sum(unlist(init_n0))
  } else {
    New_n0 <- r_iter(New_n0, 
                     F = IPM_sys$F, 
                     P = IPM_sys$P)
    Lambda <- sum(unlist(New_n0)) / sum(unlist(Output[[i-1]]))
  }  
  Output[[i]] <- New_n0
  Lambdas[[i]] <- Lambda

}

for(i in YearSeq) {
  if(i == 1) {
    New_n0 <- r_iter(init_n0, 
                     F = IpmKernels$F, 
                     P = IpmKernels$P)
    Lambda <- sum(unlist(New_n0)) / sum(unlist(init_n0))
  } else {
    New_n0 <- r_iter(New_n0, 
                     F = IpmKernels$F, 
                     P = IpmKernels$P)
    Lambda <- sum(unlist(New_n0)) / sum(unlist(Output[[i-1]]))
  }  
  Output[[i]] <- New_n0
  Lambdas[[i]] <- Lambda
  
}

Lambdas[[100]]

