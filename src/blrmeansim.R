
# Properties of the blr mean estimator
# A
# prominent case is the normal distribution on the simplex (Mateu-Figueras and
#                                                           Pawlowsky-Glahn 2008) that is followed by a D-part random composition x, if
# the random vector of its orthonormal coordinates follows a multivariate normal
# distribution on R D−1 . Accordingly, the density for a given coordinate representation
# z of x is obtained as


# Sample from a MV normal : 

# Justification for the use of the blr mean instead of arithmetic mean of the compo data, geometric mean,
# mean of the ilr transform data
library(MASS)
library(tidyverse)
library(compositions)
# Sampling from a multivariate normal distribution in the simplex is equivalent to sampling from a MV normal distribution
# in R^{D-1} then doing the ilrInverse transformation which maps to the D-dimensional unit simplex :
set.seed(636)
# number of observations in each sample
nobs <- 100
# number of simulations 
nsim <- 10000

TrueMean <- ilrInv(c(1,2))
list.matrices <- list()
for (i in 1:length(TrueMean) ){
  list.matrices[[i]] <- matrix(ncol=nsim,nrow=nobs)
}

list.compositions <- list()
for (i in 1:nsim){
  list.compositions[[i]] <- ilrInv(mvrnorm(100, c(1,2), matrix(c(1,0.5,0.5,1),nrow = 2) ) )
}

for (j in 1:length(TrueMean)){
  for (i in 1:nsim){
  list.matrices[[j]][,i] <- list.compositions[[i]][,j]  
  }
}

# Logit mapping
logit <- function(x){log((x)/(1-x))}



# Matrix of simulations, with 100 rows, 10000 columns.

blr.mean <- function(x){
  # Map a vector from [0,1] to the real space. Takes the mean and return as its output the logit-inverse mean.
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}


true.mean <- ilrInv(c(1,2))
as_tibble(seq(1:length(true.mean)),as_vector(true.mean))

true.mean <- as_tibble(cbind(true.mean,seq(1:length(true.mean))))
colnames(true.mean) <- c("value","component")

list.arithm.mean <- lapply(list.matrices,colMeans)

list.blr.mean <- list()
for (i in 1:length(TrueMean)){
  list.blr.mean[[i]] <- apply(list.matrices[[i]],2,blr.mean)
}

res <- lapply(list.compositions,ilr)

yay <- lapply(res,colMeans)
ilr.mean <- lapply(yay,ilrInv)

list.ilr.mean <- list()
for (i in 1:length(TrueMean)){
  list.ilr.mean[[i]] <- lapply(ilr.mean,function(x) x[i]) %>% unlist()
}

# We have three estimators of mean for 3 components.
# List.ilr.mean
# List.blr.mean
# List.arithm.mean

# We want a dataframe with the following structure :
# A column indicating whether it is first, second or third component
# A column indicating whether it is blr, arith or ilr mean
# A column indicating the value
# Total rows should be 3*nsim -> 30k
dim(list.ilr.mean[[1]])
str(list.ilr.mean)
names(list.ilr.mean) <- seq(1:length(list.ilr.mean))
names(list.blr.mean) <- seq(1:length(list.ilr.mean))
names(list.arithm.mean) <- seq(1:length(list.ilr.mean))

df1 <- list.ilr.mean %>% as_tibble() %>% mutate(type="ilr mean")
df2 <- list.blr.mean %>% as_tibble() %>% mutate(type="blr mean")
df3 <- list.arithm.mean %>% as_tibble() %>% mutate(type="arithm mean")

df <- bind_rows(df1,df2,df3)
df.output <- df %>% pivot_longer(cols = !type,names_to = "component",values_to = "value")



ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free")+theme_bw()+
  geom_vline(aes(xintercept=value),data=true.mean,show.legend = T,size=.5)


ggplot(df.output,aes(y=value))+geom_boxplot(aes(fill=type),alpha=.7)+
  facet_wrap(component~.,scales = "free")+theme_bw()



# Functionalize it
Mu <- c(1,2)
SigmaSquare <- matrix(c(1,0.5,0.5,1),nrow=2)
set.seed(636)


blr.distrib <- function(Mu,SigmaSquare,nobs=100,nsim=1000){
  TrueMean <- ilrInv(Mu)
  list.matrices <- list()
  for (i in 1:length(TrueMean) ){
    list.matrices[[i]] <- matrix(ncol=nsim,nrow=nobs)
  }
  
  list.compositions <- list()
  for (i in 1:nsim){
    list.compositions[[i]] <- ilrInv(mvrnorm(nobs, Mu,SigmaSquare))
  }
  
  for (j in 1:length(TrueMean)){
    for (i in 1:nsim){
      list.matrices[[j]][,i] <- list.compositions[[i]][,j]  
    }
  }
  
  # Matrix of simulations, with 100 rows, 10000 columns.
  
  blr.mean <- function(x){
    # Map a vector from [0,1] to the real space. Takes the mean and return as its output the logit-inverse mean.
    logit.t <- sapply(x,logit)
    m <- mean(logit.t,na.rm=T)
    result <- (exp(m))/( 1 + exp(m))
    return(result)
  }
  
  

  true.mean <- as_tibble(cbind(TrueMean,seq(1:length(TrueMean))))
  colnames(true.mean) <- c("value","component")
  
  list.arithm.mean <- lapply(list.matrices,colMeans)
  
  list.blr.mean <- list()
  for (i in 1:length(TrueMean)){
    list.blr.mean[[i]] <- apply(list.matrices[[i]],2,blr.mean)
  }
  
  res <- lapply(list.compositions,ilr)
  
  yay <- lapply(res,colMeans)
  ilr.mean <- lapply(yay,ilrInv)
  
  list.ilr.mean <- list()
  for (i in 1:length(TrueMean)){
    list.ilr.mean[[i]] <- lapply(ilr.mean,function(x) x[i]) %>% unlist()
  }
  
  # We have three estimators of mean for 3 components.
  # List.ilr.mean
  # List.blr.mean
  # List.arithm.mean
  
  # We want a dataframe with the following structure :
  # A column indicating whether it is first, second or third component
  # A column indicating whether it is blr, arith or ilr mean
  # A column indicating the value
  # Total rows should be 3*nsim -> 30k
  dim(list.ilr.mean[[1]])
  str(list.ilr.mean)
  names(list.ilr.mean) <- seq(1:length(list.ilr.mean))
  names(list.blr.mean) <- seq(1:length(list.ilr.mean))
  names(list.arithm.mean) <- seq(1:length(list.ilr.mean))
  
  df1 <- list.ilr.mean %>% as_tibble() %>% mutate(type="ilr mean")
  df2 <- list.blr.mean %>% as_tibble() %>% mutate(type="blr mean")
  df3 <- list.arithm.mean %>% as_tibble() %>% mutate(type="arithm mean")
  
  df <- bind_rows(df1,df2,df3)
  df.output <- df %>% pivot_longer(cols = !type,names_to = "component",values_to = "value")
  result <- list(df.output,true.mean)
  return(result)
}  

small.dimension.sim <- blr.distrib(Mu,SigmaSquare)
ggplot(small.dimension.sim[[1]],aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free")+theme_bw()+
  geom_vline(aes(xintercept=value),data=small.dimension.sim[[2]],show.legend = T,size=.5)

# Large dimensional vector :
set.seed(100)
Mu <- (rexp(10,rate = 1/3))
round(ilrInv(Mu),3)
# n <- length(Mu) 
# length(Mu)
# A <- matrix(runif(n^2)*2-1, ncol=n) 
# 
# SigmaSquare <- t(A) %*% A
SigmaSquare <- diag(abs(1/Mu),10)
Mu
ilrInv(Mu)
SigmaSquare
large.dimension.sim <- blr.distrib(Mu,SigmaSquare,nobs = 100,nsim=1000)
large.dimension.sim[[1]]$component <- factor(large.dimension.sim[[1]]$component,levels=c(1:11))
large.dimension.sim[[1]]$component
ggplot(large.dimension.sim[[1]],aes(x=value))+geom_density(aes(color=type,fill=type,bins=100),alpha=.7)+
  facet_wrap(~as.factor(component),scales = "free",nrow = 4)+theme_bw()+
  geom_vline(aes(xintercept=value),data=large.dimension.sim[[2]],show.legend = T,size=.5)

ggplot(large.dimension.sim[[1]],aes(x=value))+geom_histogram(aes(fill=type),alpha=.7,bins = )+
  facet_wrap(~as.factor(component),scales = "free",nrow = 4)+theme_bw()+
  geom_vline(aes(xintercept=value),data=large.dimension.sim[[2]],show.legend = T,size=.5)


#########################
## Testing the ILR Mean #
#########################

# A multivariate normal composition has mv normal coordinates : 
set.seed(10)
Mu <- c(-1/2,4)
SigmaSquare <- matrix(c(1,0.8,0.8,1),nrow = 2 )
nobs <- 100
nsim <- 1000
ilrInv(Mu)
TrueMean <- ilrInv(Mu)
list.matrices <- list()
for (i in 1:length(TrueMean) ){
  list.matrices[[i]] <- matrix(ncol=nsim,nrow=nobs)
}

list.compositions <- list()
for (i in 1:nsim){
  list.compositions[[i]] <- ilrInv(mvrnorm(nobs, Mu,SigmaSquare))
}

for (j in 1:length(TrueMean)){
  for (i in 1:nsim){
    list.matrices[[j]][,i] <- list.compositions[[i]][,j]  
  }
}
# List matrices with 3 elements, one for each parts
# For each parts, 100 rows corresponding to the sample, 1000 columns corresponding to a different simulation
# We compute in a univariate way blr.mean :
logit <- function(x){log((x)/(1-x))}
blr.mean <- function(x){
  # Map a vector from [0,1] to the real space. Takes the mean and return as its output the logit-inverse mean.
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}
blr.median <- function(x){
  # Map a vector from [0,1] to the real space. Takes the mean and return as its output the logit-inverse mean.
  logit.t <- sapply(x,logit)
  m <- median(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}


list.blr.mean <- list()
for (i in 1:length(TrueMean)){
  list.blr.mean[[i]] <- apply(list.matrices[[i]],2,blr.mean)
}
list.blr.median <- list()
for (i in 1:length(TrueMean)){
  list.blr.median[[i]] <- apply(list.matrices[[i]],2,blr.median)
}
list.arithm.mean <- lapply(list.matrices,colMeans)
res <- lapply(list.compositions,ilr)

yay <- lapply(res,colMeans)
ilr.mean <- lapply(yay,ilrInv)

list.ilr.mean <- list()
for (i in 1:length(TrueMean)){
  list.ilr.mean[[i]] <- lapply(ilr.mean,function(x) x[i]) %>% unlist()
}
names(list.ilr.mean) <- seq(1:length(list.ilr.mean))
names(list.blr.mean) <- seq(1:length(list.ilr.mean))
names(list.arithm.mean) <- seq(1:length(list.ilr.mean))
names(list.blr.median) <- seq(1:length(list.ilr.mean))
df1 <- list.ilr.mean %>% as_tibble() %>% mutate(type="ilr mean")
df2 <- list.blr.mean %>% as_tibble() %>% mutate(type="blr mean")
df3 <- list.arithm.mean %>% as_tibble() %>% mutate(type="arithm mean")
df4 <- list.blr.median %>% as_tibble %>% mutate(type="blr median")
df <- bind_rows(df1,df2,df3,df4)
df.output <- df %>% pivot_longer(cols = !type,names_to = "component",values_to = "value")
ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free")+theme_bw()

library(ggtern)
# We can plot the sample of 100 observations from the 1st simulation :
firstsim.df <- tibble(list.matrices[[1]][,1],list.matrices[[2]][,1],list.matrices[[3]][,1])

names(firstsim.df) <- c("FirstPart","SecondPart","ThirdPart")
ggtern(firstsim.df)+geom_point(aes(x=FirstPart,y=SecondPart,z=ThirdPart),alpha=.05,color="blue")+theme_bw()

# distribution of estimators of the mean of a mv-normal in the simplex 
