# Logit Normality in Data Acquisition Process in Geochemistry, 
# comparing with other assumptions (ilr normality, clr normality, R normality)

# Several Parts, first is the data acquisition process algorithm, presented hereunder : 

########################
# DATA ACQUISITION ALGO
#######################

library(MASS)
library(tidyverse)
library(compositions)
library(xtable)
library(rootSolve)
library(EnvStats)
library(ggpubr)
library(quadprog)

# Clear Workspace


# Proportion of Oxygen in Major Oxides 
frac_el <- 1/100*c(46.75, 47.87, 52.93, 69.94, 77.44, 60.31, 71.47, 74.18, 83.01, 43.64)
frac_ox <- rep(1,times=length(frac_el))-frac_el


# Let a true 22 parts composition roughly inspired from GEOPT48 : 
true <- c(Si = 0.26501011, Ti = 0.00504876, Al = 0.10090108, Fe = 0.05051534, 
          Mn = 0.00114497, Mg = 0.00636078, Ca = 0.02521995, Na = 0.04770775, 
          K = 0.03117402, P = 0.00256124, O = 0.45967706, Ba = 0.00128169, 
          Sr = 0.00095118, Zr = 0.0004889, Ce = 0.0001817, F = 0.00013441, 
          Cr = 0.0001311, Nb = 9.66e-05, Zn = 7.195e-05, S = 6.746e-05, 
          Rb = 6.407e-05, U = 0.00120988)

name.majors <- names(true)[1:10]
name.traces <- names(true)[12:21]
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

blr.min <- function(mu){
  D <- length(mu)
  Dmat <- 2 * diag(D)
  dvec <- 2 * mu
  
  Amat <- cbind(rep(1, D), diag(D))
  bvec <- c(1, rep(0, D))
  
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  x <- sol$solution
  return(x)
}

opti.blr.m <- function(x) {
  blr.m <- apply(x,2,blr.mean)
  opti.m <- blr.min(blr.m)
  names(opti.m) <- names(blr.m)
  return(opti.m)
}

mu <- true[c(c("Si", "Ti", "Al", "Fe", "Mn", "Mg", "Ca", "Na", "K", "P", 
               "Ba", "Sr", "Zr", "Ce", "F", "Cr", "Nb", "Zn", "S", "Rb", "U"
))] %>% logit() %>% unclass()

# mu <- mu.parameters[[27]]
# Sigma <- variance.R.space[[27]]

GJ.sim <- function(true,Sigma){
  mu <- logit(true)[1:21]
  sample.inR <- MASS::mvrnorm(n=100,mu = mu,Sigma = Sigma)
  sample.inH <- sample.inR %>% as_tibble() %>% apply(2,logitInv) %>% as_tibble()
  # Keep only the major elements :
  sample.inH.onlymajor <- sample.inH[1:10]
  
  for (i in 1:nrow(sample.inR)){
    sample.inH.onlymajor[i,] <- clo(sample.inH.onlymajor[i,])*sum(true[1:10])
  }
  
  
  # Now one should perturb this sample by perturbating each observation
  # with a measurement error term. That's each row should deviate from 
  # the exact proportion (0.535644 %) by some random error term...
  
  # 1 Map each column of the matrix of major elements in the simplex
  logit.transf.sample <- logit(sample.inH.onlymajor)
  # Let a multiplicative error term which is lognormal distributed with
  # mean 0 and standard deviation 0.005.
  
  multiplicative.error.term <- exp(rnorm(100,0,2.5E-2))
  # Perturb the major elements 
  logit.transf.sample.perturbed <- logit.transf.sample*multiplicative.error.term
  # Then perform logitINV transf
  perturbed.sample.major <- logit.transf.sample.perturbed %>% logitInv()
  
  #  Step 2 : Compute the oxygen for each composition (deterministic quantity, not random !)
  oxygen.computed <- c()
  for  (i in 1:nrow(perturbed.sample.major)){
    oxygen.computed[i] <- sum(frac_el^{-1}*perturbed.sample.major[i,]-perturbed.sample.major[i,])
  }
  names(perturbed.sample.major) <- name.majors
  
  perturbed.sample.major$O <- oxygen.computed
  
  # Step 3 draw randomly from traces, now Sigma is 1E-4*diag(abs(mu))
  # Keep only the trace elements :
  sample.inH.onlytraces <- sample.inH[12:21]
  # Close it
  for (i in 1:nrow(sample.inR)){
    sample.inH.onlytraces[i,] <- clo(sample.inH.onlytraces[i,])*sum(true[12:21])
  }
  # Perturb the sample with a lognormal error term : 
  # 1 Map each column of the matrix of major elements in the simplex
  logit.transf.sample <- logit(sample.inH.onlytraces)
  # Let a multiplicative error term which is lognormal distributed with
  # mean 0 and standard deviation 0.005.
  
  multiplicative.error.term <- exp(rnorm(100,0,2.5E-2))
  # Perturb the major elements 
  logit.transf.sample.perturbed <- logit.transf.sample*multiplicative.error.term
  # Then perform logitINV transf
  perturbed.sample.traces <- logit.transf.sample.perturbed %>% logitInv()
  
  
  df <- data.frame(perturbed.sample.major,perturbed.sample.traces)
  
  virtual.U <- 1-rowSums(df)
  sample.U <- ifelse(virtual.U>=0,virtual.U,NA)
  mean.nonnegu <- sample.U %>% log() %>% mean(na.rm=T)
  sd.nonnegu <- sample.U %>% log() %>% sd(na.rm=T)
  nb.missingval <- length(sample.U[is.na(sample.U)] )
  replacement.val <- rnorm(nb.missingval,mean.nonnegu,sd.nonnegu) %>% exp()
  df$U <- sample.U
  df$U <- replace(df$U,which(is.na(df$U)),replacement.val)
  df
  #df.w.U <- clo(df)
  #df <- clo(df[-22])
  return(df)
}  


# Applying data acqusition process algorithm to generate 500 datasets :


 #list.df <- list()
 #m <- 500
 #for (i in 1:500){
#   list.df[[i]] <- GJ.sim(mu)
# }

#saveRDS(list.df,file="src/generated_geopt48_data.RData" )

 

