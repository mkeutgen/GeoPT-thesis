#### New Simulation GeoPT Data Missing Values ####
#### Gert Jan Weltje & Maxime Keutgen De Greef ###
#### January 17 2023 #############################

library(tidyverse)
library(compositions)
library(MASS)
library(xtable)
set.seed(10)

# Logit and BLR Mean Function
logit <- function(x){log((x)/(1-x))}
logitInv <- function(x){1/(1+exp(-x))}
# BLR Mean function
blr.mean <- function(x){
  # Map a vector from [0,1] to the real space. Takes the mean and return as its output the logit-inverse mean.
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}


# Let a composition of major elements, trace elements, oxygen and 
# the rest.


## 1. Create the composition
# Major elements 
name.oxmajors <- c("SiO2", "TiO2", "Al2O3", "Fe2O3T", "MnO",
               "MgO", 
               "CaO", "Na2O", "K2O", "P2O5")

name.majors <- c("Si","Ti","Al","Fe","Mn","Mg","Ca","Na","K","P")

oxmajors <- rexp(length(name.majors),3) %>% clo()*0.99
# Fraction of element in the major oxides

frac_el <- 1/100*c(46.75, 47.87, 52.93, 69.94, 77.44, 60.31, 71.47, 74.18, 83.01, 43.64)

# Fraction of oxygen in the major oxides  
frac_ox <- rep(1,times=length(frac_el))-frac_el
# Sum of element majors : 0.65116
element.major <- oxmajors * frac_el
names(element.major) <- name.majors

oxygen.fractions <- oxmajors * frac_ox
# Sum of oxygen : 0.339 
oxygen <- sum(oxygen.fractions)
names(oxygen) <- "O"
# 10 Traces elements
trace.elements <- rexp(length(name.oxmajors),3) %>% clo()*0.008
names(trace.elements) <- c("trace1",
  "trace2",
  "trace3",
  "trace4",
  "trace5",
  "trace6",
  "trace7",
  "trace8",
  "trace9",
  "trace10",
  "trace11"
)
mu <- c(element.major,oxygen,trace.elements)
mu["rest"] <- 1-sum(mu)
# Generating variance covariance matrices for element major, oxygen,
# trace elements.

# Variance matrix for major elements : 

# Sample from a MV normal model with mu = alr(element.major), 
# sigma = diag(abs(alr(element.major))) (additive logistic
# normal model).
samp <- mvrnorm(100,mu=alr(element.major),Sigma =
          diag(abs(alr(element.major))))

# Major elements of the first element
el1 <- samp[sample(nrow(samp),size = 1)] %>% alrInv()
el1 <- el1 %>% unclass()*sum(element.major)
# Perturb the major elements with a perturbation vector 
# which follows some distribution around 0
perturbation.vector <- mvrnorm(n=1,mu=rep(0,length(element.major)),Sigma = diag(abs(element.major)))
perturbation.vector
logit(el1) 
el1
rnorm(n=1,mu=0,sd=element.major)
# Trace elements of the first element

sum(element.major)


sum(el1*sum(element.major))


mvrnorm(1,mu = alr(element.major),
        Sigma = diag(abs(alr(element.major))))

mvrnorm(1,alr(element.major),Sigma=alr(diag(element.major))

