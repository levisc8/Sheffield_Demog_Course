# Age x Size IPM

# Sam Levin
# 1/25/18

rm(list = ls())
library(car)
library(dplyr)

source('R/WD_Paths.R')

# Load data, reset working directory
setwd(data.wd)
data <- read.csv('25.1.18/Soay_Age_Size/soay_demog_data.csv',
                 stringsAsFactors = FALSE) %>%
  as_tibble()

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
  mu <- m_par['grow_int'] + z * m_par['grow_z'] + a * m_par['grow_a']
  out <- dnorm(z1, mean = mu, sd = m_par['grow_sd'])
  return(out) 
}

s_z <- function(z, a, m_par) {
  linear_mu <- m_par['surv_int'] + z * m_par['surv_z'] + a * m_par['grow_a']
  out <- 1 / (1 + exp(-linear_mu))
  return(out)
}

pb_z <- function(z, a, m_par) {
  if(a == 0) {
    out <- 0
  } else {
    linear_mu <- m_par['repr_int'] + z * m_par['repr_z'] + a * m_par['repr_a']
    out <- 1 / (1 + exp(-linear_mu))
  }
  return(out)
}

pr_z <- function(a, m_par) { 
  linear_mu <- m_par['recr_int'] + a * m_par['recr_a']
  out <- 1 / (1 + exp(-linear_mu))
  return(out)
}  

c_z1z <- function(z, z1, a, m_par) { 
  mu <- m_par['rcsz_int'] + z * m_par['rcsz_z'] + a * m_par['rcsz_a']
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
  out <- s_z(z, a, m_par) * pb_z(z, a, m_par) *
    pr_z(a, m_par) * c_z1z(z, z1, a, m_par)
  return(out/2)
}

L <- min(data$z) * 0.9
U <- max(data$z) * 1.1
MeshPoints <- seq(L, U, length.out = 250)
BinWidth <- MeshPoints[2] - MeshPoints[1]
Ages <- seq(0, max(data$a, na.rm = TRUE), by = 1)

library(Matrix)
GMat <- BinWidth * outer(MeshPoints,
                         MeshPoints,
                         FUN = g_z1z)

