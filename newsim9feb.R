# New Sim 9 February, U is a missing value if sum of the parts  is above 1  
library(MASS)
library(tidyverse)
library(compositions)
library(xtable)
library(rootSolve)
library(EnvStats)
library(ggpubr)
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






mu <- true[c(c("Si", "Ti", "Al", "Fe", "Mn", "Mg", "Ca", "Na", "K", "P", 
               "Ba", "Sr", "Zr", "Ce", "F", "Cr", "Nb", "Zn", "S", "Rb", "U"
))] %>% ilr() %>% unclass()


GJ.sim <- function(mu){
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
  
  # Step 4 Undetected :
  
  df <- data.frame(perturbed.sample.major,perturbed.sample.traces)
  # Undetected follows a lognormal 
  U <- logitInv(logit(true["U"])*exp(rnorm(100,0,0.5E-2)))
  df$U <- U
  names(df) <- names(true)
  df %>% colMeans()
  return(df)
}
list.df <- list()
m <- 500
for (i in 1:500){
list.df[[i]] <- GJ.sim(mu)
}


#saveRDS(list.df,file="twentyfourthfebGJsim.RData")

list.df <- readRDS("src/firstfebGJsim.RData")

list.matrices <- list()
m <- 500
for (i in 1:3){
  list.matrices[[i]] <- matrix(ncol=length(true),nrow=m)
  
}

for (i in 1:m){
  list.matrices[[1]][i,] <- list.df[[i]] %>% colMeans()
  list.matrices[[2]][i,] <- apply(list.df[[i]],2,blr.mean)
  list.matrices[[3]][i,]  <- ilr(list.df[[i]]) %>% colMeans() %>% 
    ilrInv() %>% unclass()
}


df1 <- list.matrices[[1]] %>% as_tibble()
df2 <- list.matrices[[2]] %>% as_tibble() 
df3 <- list.matrices[[3]] %>% as_tibble() 

colnames(df1) <- names(true) 
colnames(df2) <- names(true)
colnames(df3) <- names(true)

df1 <- df1 %>% mutate(type="arithmetic mean")
df2 <- df2 %>% mutate(type="blr mean")
df3 <- df3 %>% mutate(type="ilr mean")


df <- bind_rows(df1,df2,df3)

df.output <- df %>% pivot_longer(cols = !type,names_to = "component",values_to = "value")


ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free",ncol=3,shrink = T)+theme_bw()+
  scale_x_continuous(labels=scales::scientific)+
  scale_y_continuous(labels=scales::scientific)+
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

# Computing biases
arithm.bias <- data.frame(matrix(ncol=22,nrow=500))
blr.bias <-    data.frame(matrix(ncol=22,nrow=500))
ilr.bias <-    data.frame(matrix(ncol=22,nrow=500))

for (i in 1:m){
  arithm.bias[i,] <- df1[1:22][i,]-true
  blr.bias[i,] <- df2[1:22][i,]-true
  ilr.bias[i,] <- df3[1:22][i,]-true
}



estim.bias.arithm <- arithm.bias %>% colMeans() 
estim.bias.blr <- blr.bias %>% colMeans() 
estim.bias.ilr <- ilr.bias %>% colMeans() 

df.bias <- data.frame(estim.bias.arithm,estim.bias.blr,estim.bias.ilr,
                      chem.element=names(true))
df.bias %>%
  pivot_longer(cols=!chem.element) %>%
  ggplot(aes(x=chem.element,y=value,color=name))+
  geom_point(size=2,alpha=.5)+theme_bw()



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


# About the distribution of 1-U, which is a lognormal quantity which can be above 1, we have that :

# true compo.wo.U :
true.compo.wo.U <- true[1:21]
# Lognormal distribution on the real line : \mathcal{N}_(+)(log(true.compo.wo.U),1E-4*log(true.compo.wo.U))

# sd for each component :

log.mu <- log(true.compo.wo.U)
log.sigmasq <- 1E-4*log(true.compo.wo.U)

N <- (sum(exp(2*log.mu+log.sigmasq)*(exp(log.sigmasq-1))))
D <- sum(exp(log.mu+log.sigmasq/2))^2
sigmasq_z <- log(N/D + 1)

mu_z <- log(sum(exp(log.mu+sigmasq_z/2)))-sigmasq_z/2



# Lognormal with parameters mu_z, sigmasq_z :
sample <- rnorm(1E6,mu_z,sigmasq_z) %>% exp() %>% as_tibble()
mean(sample$value)
median(sample$value)

# Theoretical U mu, exp of the mu parameter
U.mu <- median(1-sample$value)
U.virtual <- 1-sample$value
tibble(U.virtual) %>% ggplot(aes(x=U.virtual))+
  geom_histogram(fill="lightblue",colour="black")+theme_bw()+
  geom_vline(xintercept=0,colour="red",size=1)


# Actually mu is known when sum of known parts is under 1,
# missing when it is above or equal to 1 :
sample$U <- ifelse(sample$value<=1,1-sample$value,NA)
median(sample$U,na.rm = T) %>% log()
hist(sample$U)

xbar <- mean(sample$U,na.rm = T)

non.neg.u <- sample$U[!is.na(sample$U)]
# MLE for var(U) and E(U) : 
var.u <- sd(log(non.neg.u))^2
mean.u <- mean(log(sample$U),na.rm = T)

# How well is non neg U described by a lognormal model ? Plot and test
# Generate random sample
lognormal.law <- rnorm(length(non.neg.u),mean = mean.u,sd = sqrt(var.u)) %>% exp()
lognormal.law.trunc <- rnormTrunc(length(non.neg.u),mean = mean.u,sd = sqrt(var.u),min=log(0),max = log(1)) %>% exp()
normal.law.trunc <- rnormTrunc(length(non.neg.u),mean=mean(non.neg.u),
                               sd = sd(non.neg.u),min = min(non.neg.u),max=max(non.neg.u))


#  max(non.neg.u)
#[1] 0.3808868  -> log : -0.965253
# min(non.neg.u)
# 3.715005e-07 -> log : -14.80572
lognormal.law.trunc <- rnormTrunc(length(non.neg.u),
                                  mean = mean.u,sd = sqrt(var.u),min=-14.80572,max = -0.96525) %>% exp()


exp(mean(log(lognormal.law.trunc)))
exp(mean(log(lognormal.law)))
exp(mean(log(non.neg.u)))
mean(normal.law.trunc)
mean(non.neg.u)


t <- rexp(n = length(non.neg.u),1/mean(non.neg.u))

ks.test(non.neg.u,t)
rexp(n = length(non.neg.u),1/mean(non.neg.u)) %>% as_tibble() %>% ggplot(aes(x=value))+geom_histogram()


rnormTrunc(1E2,mu_z,sigmasq_z,min = -Inf) %>% exp() 

df <- tibble(non.neg.u,lognormal.law,lognormal.law.trunc) 

ggplot(df,aes(x=lognormal.law))+geom_histogram(fill="lightblue",colour="black")+
  theme_bw()+annotate("text",x = 0.3,y=4E4,label=paste("D statistic =",round(kslognorm$statistic,digits = 5)),size=6)+
  annotate("text",x = 0.3,y=2E4,label="pvalue < 0.01",size=6)+
  scale_x_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4),limits = c(0,0.4))

# Ks test

kslognorm <- ks.test(non.neg.u,lognormal.law) # D = 0.102
kslognorm.trunc <- ks.test(non.neg.u,lognormal.law.trunc) # D = 0.12013
ksnorm.trunc <- ks.test(non.neg.u,normal.law.trunc) # D = 0.13091


# Histogram plots

h1 <- ggplot(df,aes(x=non.neg.u))+geom_histogram(fill="lightblue",colour="black")+theme_bw()
h2 <- ggplot(df,aes(x=lognormal.law))+geom_histogram(fill="lightblue",colour="black")+
  theme_bw()+annotate("text",x = 0.3,y=4E4,label=paste("D statistic =",round(kslognorm$statistic,digits = 5)),size=3)+
  annotate("text",x = 0.3,y=2E4,label="pvalue < 0.01",size=3)+
  scale_x_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4),limits = c(0,0.4))
h3 <- ggplot(df,aes(x=lognormal.law.trunc))+geom_histogram(fill="lightblue",colour="black")+
  theme_bw()+annotate("text",x = 0.3,y=4E4,label=paste("D statistic =",round(kslognorm.trunc$statistic,digits = 5)),size=3)+
  annotate("text",x = 0.3,y=2E4,label="pvalue < 0.01",size=3)+
  scale_x_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4),limits = c(0,0.4))
h4 <- ggplot(df,aes(x=normal.law.trunc))+geom_histogram(fill="lightblue",colour="black")+
  theme_bw()+annotate("text",x = 0.3,y=4E4,label=paste("D statistic =",round(ksnorm.trunc$statistic,digits = 5)),size=3)+
  annotate("text",x = 0.3,y=2E4,label="pvalue < 0.01",size=3)+
  scale_x_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4),limits = c(0,0.4))

ggarrange(h1,h2,h3,h4)

df %>% select(non.neg.u,lognormal.law.trunc) %>%
  pivot_longer(cols=everything()) %>% ggplot(aes(x=value,fill=name)) + geom_histogram(alpha=.5)+theme_bw()

hist(non.neg.u)
hist(lognormal.law.trunc)

lognormal.law.trunc %>% max()
  
exp(sd(log(non.neg.u)))
exp(sd(log(lognormal.law.trunc)))

ks.test(lognormal.law,non.neg.u)







