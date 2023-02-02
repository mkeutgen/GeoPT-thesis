# New Sim Gert Jan 25 January 2023 #
library(MASS)
library(tidyverse)
library(compositions)
library(xtable)

frac_el <- 1/100*c(46.75, 47.87, 52.93, 69.94, 77.44, 60.31, 71.47, 74.18, 83.01, 43.64)
frac_ox <- rep(1,times=length(frac_el))-frac_el


# Let a true 22 parts composition rougly inspired from GEOPT48 : 
true <- c(Si = 0.26501011, Ti = 0.00504876, Al = 0.10090108, Fe = 0.05051534, 
          Mn = 0.00114497, Mg = 0.00636078, Ca = 0.02521995, Na = 0.04770775, 
          K = 0.03117402, P = 0.00256124, O = 0.45967706, Ba = 0.00128169, 
          Sr = 0.00095118, Zr = 0.0004889, Ce = 0.0001817, F = 0.00013441, 
          Cr = 0.0001311, Nb = 9.66e-05, Zn = 7.195e-05, S = 6.746e-05, 
          Rb = 6.407e-05, Rest = 0.00120988)

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



set.seed(1)


mu <- true[c(c("Si", "Ti", "Al", "Fe", "Mn", "Mg", "Ca", "Na", "K", "P", 
               "Ba", "Sr", "Zr", "Ce", "F", "Cr", "Nb", "Zn", "S", "Rb", "Rest"
))] %>% ilr() %>% unclass()
#n <- length(mu)
#A <- matrix(runif(n^2)*2-1, ncol=n) 
#Sigma <- t(A) %*% A
# 100 observations of major and oxygen elements 
sample.inR <- MASS::mvrnorm(n=100,mu = mu,Sigma = 1E-2*diag(abs(mu)))
sample.inS <- sample.inR %>% ilrInv() %>% as_tibble()
# Keep only the major elements :
sample.inS.onlymajor <- sample.inS[1:10]

for (i in 1:nrow(sample.inR)){
  sample.inS.onlymajor[i,] <- clo(sample.inS.onlymajor[i,])*sum(true[1:10])
  }




# Now one should perturb this sample by perturbating each observation
# with a measurement error term. That's each row should deviate from 
# the exact proportion (0.535644 %) by some random error term...

# 1 Map each column of the matrix of major elements in the simplex
logit.transf.sample <- logit(sample.inS.onlymajor)
# Let a multiplicative error term which is lognormal distributed with
# mean 0 and standard deviation 0.005.

multiplicative.error.term <- exp(rnorm(100,0,0.5E-2))
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
sample.inR <- MASS::mvrnorm(n=100,mu = mu,Sigma = 1E-2*diag(abs(mu)))
sample.inS <- sample.inR %>% ilrInv() %>% as_tibble()
# Keep only the trace elements :
sample.inS.onlytraces <- sample.inS[12:21]
# Close it
for (i in 1:nrow(sample.inR)){
  sample.inS.onlytraces[i,] <- clo(sample.inS.onlytraces[i,])*sum(true[12:21])
}
# Perturb the sample with a lognormal error term : 
# 1 Map each column of the matrix of major elements in the simplex
logit.transf.sample <- logit(sample.inS.onlytraces)
# Let a multiplicative error term which is lognormal distributed with
# mean 0 and standard deviation 0.005.

multiplicative.error.term <- exp(rnorm(100,0,0.5E-2))
# Perturb the major elements 
logit.transf.sample.perturbed <- logit.transf.sample*multiplicative.error.term
# Then perform logitINV transf
perturbed.sample.traces <- logit.transf.sample.perturbed %>% logitInv()

# Step 4 Rest :

df <- data.frame(perturbed.sample.major,perturbed.sample.traces)
# rest follows a lognormal 
df
Rest <- true["Rest"]*exp(rnorm(100,0,0.5E-2))
df$Rest <- Rest
names(df) <- names(true)
df %>% colMeans()

blr.mean.vector <- apply(df,2,blr.mean)
blr.mean.vector
true
