# New Sim 24 February, U is a missing value if sum of the parts  is above 1  
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

frac_el <- 1/100*c(46.75, 47.87, 52.93, 69.94, 77.44, 60.31, 71.47, 74.18, 83.01, 43.64)
frac_ox <- rep(1,times=length(frac_el))-frac_el


# Let a true 22 parts composition rougly inspired from GEOPT48 : 
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
))] %>% ilr() %>% unclass()


GJ.sim <- function(mu){
  sample.inR <- MASS::mvrnorm(n=200,mu = mu,Sigma = 1E-2*diag(abs(mu)))
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
  
  multiplicative.error.term <- exp(rnorm(200,0,0.5E-2))
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
  sample.inR <- MASS::mvrnorm(n=200,mu = mu,Sigma = 1E-2*diag(abs(mu)))
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
  
  multiplicative.error.term <- exp(rnorm(200,0,0.5E-2))
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
  names(df) <- names(true)  
  df.w.U <- clo(df)
  df <- clo(df[-22])
  return(df)
}  


#list.df <- list()
#m <- 500
#for (i in 1:500){
#  list.df[[i]] <- GJ.sim(mu)
#}

#saveRDS(list.df,file="src/twentyfourfebsimdata.RData")

##################
# ESTIMATION ALGO
##################
  
# import data

list.df <- readRDS("src/twentyfourfebsimdata.RData")



# Modify blr.mean, objective function 

blr.mu <- apply(list.df[[1]],2,blr.mean)
blr.mu.opti <- opti.blr.m(list.df[[1]])





list.matrices <- list()
m <- 500
for (i in 1:4){
  list.matrices[[i]] <- matrix(ncol=21,nrow=m)
  
}

for (i in 1:m){
  list.matrices[[1]][i,] <- list.df[[i]] %>% colMeans()
  list.matrices[[2]][i,] <- apply(list.df[[i]],2,blr.mean)
  list.matrices[[3]][i,]  <- ilr(list.df[[i]]) %>% colMeans() %>% 
    ilrInv() %>% unclass()
  list.matrices[[4]][i,] <- opti.blr.m(list.df[[i]]) 
  
}


df1 <- list.matrices[[1]] %>% as_tibble()
df2 <- list.matrices[[2]] %>% as_tibble() 
df3 <- list.matrices[[3]] %>% as_tibble() 
df4 <- list.matrices[[4]] %>% as_tibble() 

colnames(df1) <- names(true) 
colnames(df2) <- names(true)
colnames(df3) <- names(true)
colnames(df4) <- names(true)

df1 <- df1 %>% mutate(type="arithmetic mean")
df2 <- df2 %>% mutate(type="blr mean")
df3 <- df3 %>% mutate(type="ilr mean")
df4 <- df4 %>% mutate(type="optimized blr mean")


df <- bind_rows(df1,df2,df3,df4)

df.output <- df %>% pivot_longer(cols = !type,names_to = "component",values_to = "value")


True.Df <- as_tibble(rownames="component",true)
True.Df <- True.Df[1:21,] 


ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free",ncol=3,shrink = T)+theme_bw()+
  scale_x_continuous(labels=scales::scientific)+
  scale_y_continuous(labels=scales::scientific)+
  geom_vline(aes(xintercept=value),data=True.Df,show.legend = T,size=.5,color="darkblue")+
  theme(legend.position = "bottom")



True.Df <- as_tibble(rownames="component",true)
True.Df

# True.Df$component <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14")


ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free",ncol = 3)+theme_bw()+
  geom_vline(aes(xintercept=value),data=True.Df,show.legend = T,size=.5,color="darkblue")+
  theme(legend.position = "bottom")




blr.mean.vector <- apply(GJ.sim(mu)
                         ,2,blr.mean)
blr.mean.vector

true <- true[1:21] 
# Computing biases
arithm.bias <- data.frame(matrix(ncol=21,nrow=500))
blr.bias <-    data.frame(matrix(ncol=21,nrow=500))
ilr.bias <-    data.frame(matrix(ncol=21,nrow=500))
optim.blr.bias <- data.frame(matrix(ncol=21,nrow=500))

for (i in 1:m){
  arithm.bias[i,] <- df1[1:21][i,]-true
  blr.bias[i,] <- df2[1:21][i,]-true
  ilr.bias[i,] <- df3[1:21][i,]-true
  optim.blr.bias[i,] <- df4[1:21][i,]-true
}



estim.bias.arithm <- arithm.bias %>% colMeans() 
estim.bias.blr <- blr.bias %>% colMeans() 
estim.bias.ilr <- ilr.bias %>% colMeans() 
estim.bias.optim.blr <- optim.blr.bias %>% colMeans() 

df.bias <- data.frame(estim.bias.arithm,estim.bias.blr,estim.bias.ilr,estim.bias.optim.blr,
                      chem.element=names(true))
df.bias %>%
  pivot_longer(cols=!chem.element) %>%
  ggplot(aes(x=chem.element,y=value,color=name))+
  geom_point(size=2,alpha=.5)+theme_bw()


df.bias
colnames(arithm.bias) <- names(true)
colnames(blr.bias) <- names(true)
colnames(ilr.bias) <- names(true)

arithm.bias$type <- "arithm.mean"
blr.bias$type <- "blr.mean"
ilr.bias$type <- "ilr.mean"

df.boxplot <- data.frame(rbind(df1,df2,df3))
df.boxplot %>% pivot_longer(cols=!type) %>% 
  ggplot(aes(x=name,y=value,color=type)) + geom_boxplot()+
  theme_bw()


#df.bias <- data.frame(arithm.bias,blr.bias,ilr.bias)
# Showing the bias for each chemical element :
df.bias$chem_element <- rownames(df.bias)
#df.bias %>% pivot_longer(cols = !chem_element) %>%
#ggplot(aes(x=chem_element,y=value,color=name),alpha=.4)+geom_point()+
#theme_bw()

