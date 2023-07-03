#### Mean Compo Dataset
library(MASS)
library(tidyverse)
library(compositions)
library(xtable)
library(rootSolve)
library(EnvStats)
library(ggpubr)
library(quadprog)
library(missMethods)

frac_el <- 1/100*c(46.75, 47.87, 52.93, 69.94, 77.44, 60.31, 71.47, 74.18, 83.01, 43.64)
frac_ox <- rep(1,times=length(frac_el))-frac_el


# Code to import all processed dataframes, returns a named list of the processed dataframes.
source("src/importdata.R")
# Blr (really a logit transformation) 
logit <- function(x){log((x)/(1-x))}
logitInv <- function(x){1/(1+exp(-x))}

blr <- function(dataframe) apply(dataframe,2,logit)
blrInv <- function(dataframe) apply(dataframe,2,logitInv)
# 10 Major Elements
names.full.wO <- c("SiO2","TiO2","Al2O3","Fe2O3T","MnO","MgO","CaO","Na2O","K2O","P2O5",
                   "Ag","Ba","Cu","Rb","Zr","Sr","Li","Ce","Zn","V")
names.full <- c("Si","Ti","Al","Fe","Mn","Mg","Ca","Na","K","P","O",
                "Ag","Ba","Cu","Rb","Zr","Sr","Li","Ce","Zn","V")

majors.el <- c("SiO2","TiO2","Al2O3","Fe2O3T","MnO","MgO","CaO","Na2O","K2O","P2O5")
frac_el <- 1/100*c(46.75, 47.87, 52.93, 69.94, 77.44, 60.31, 71.47, 74.18, 83.01, 43.64)
frac_ox <- rep(1,times=length(frac_el))-frac_el



# Df of means
list.c <- list()
for (i in 1:29){
list.c[[i]] <- list.df[[i]] %>% blr() %>% colMeans() %>% logitInv() %>% sort(,decreasing=TRUE)
list.c[[i]] <- list.c[[i]][1:21]
}

proper.format <- list.df
for (i in 1:29){
proper.format[[i]]$Fe2O3T <- proper.format[[i]]$Fe2O3T + proper.format[[i]]$`Fe(II)O`
proper.format[[i]]$`Fe(II)O` <- NULL
}

major.el <- proper.format


for (i in 1:29){
  major.el[[i]] <- major.el[[i]][1:10]
  # multiply major elements mass fraction with oxygen
  major.el[[i]] <- major.el[[i]] %>% as.matrix() %*% diag(frac_el) %>% as_tibble()
}

for (i in 1:29){
#proper.format[[i]] <- proper.format[[i]][names.full.wO]
t <- proper.format[[i]][majors.el] %>% as.matrix() %*% diag(frac_ox) %>% rowSums()
major.el[[i]]$O <- t
proper.format[[i]] <- as_tibble(cbind(major.el[[i]],proper.format[[i]][11:20]))
names(proper.format[[i]])[1:11] <- names.full[1:11]
}



blr.df <- lapply(proper.format,blr) 
lapply(blr.df,var)

lapply(proper.format,ncol)
var(blr.df[[1]])
var(blr.df[[2]])
# Process Infinity values :
blr.df[[2]]
list.imputed.blr.df <- list()
list.var.df <- list()
for (i in 1:length(blr.df)){
list.imputed.blr.df[[i]] <- replace(blr.df[[i]],blr.df[[i]]==-Inf,NA)
list.imputed.blr.df[[i]] <- impute_EM(list.imputed.blr.df[[i]])
list.var.df[[i]] <- var(list.imputed.blr.df[[i]]) 
}
names(list.var.df) <- names(list.df)

saveRDS(list.var.df,file="src/list_var_df.RData")

list.var.df[[2]] %>% colnames()



meancompo <- lapply(lapply(blr.df,colMeans),logitInv)

for (i in 1:length(meancompo)){
  meancompo[[i]]["U"] <- 1 - sum(meancompo[[i]]) 
}

saveRDS(meancompo,file = "src/meancompo.RData")

View(df)
df
name.majors <- names(df)[1:10]
name.traces <- names(df)[12:21]


# Logit and BLR Mean Function
# BLR Mean function
blr.mean <- function(x){
  # Map a vector from [0,1] to the real space. Takes the mean and return as its output the logit-inverse mean.
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}

# What is the order of magnitude of the multiplicative error term ?

j <- 26

sim.df1 <- as_tibble(matrix(nrow=nrow(proper.format[[j]] ) ,ncol=ncol(proper.format[[j]])))
naive.norm <- sim.df1
naive.lognorm <- sim.df1
logit.df1 <- apply(proper.format[[j]],2,logit)
log.df1 <- apply(proper.format[[j]],2,log)
colmean.df1 <- colMeans(logit.df1,na.rm = T)
sd.df1 <- apply(logit.df1,2,sd,na.rm=T)

naive.norm.m <- colMeans(proper.format[[j]],na.rm=T)
naive.norm.sd <- apply(proper.format[[j]],2,sd,na.rm=T)

naive.lognorm.m <- colMeans(log.df1,na.rm=T)
naive.lognorm.sd <- apply(log.df1,2,sd,na.rm=T)

rnorm(23,mean = colmean.df1[1],sd=sd.df1[1])
for (i in 1:ncol(proper.format[[j]])){
sim.df1[i] <- rnorm(nrow(proper.format[[j]]),colmean.df1[i],sd.df1[i]) %>% logitInv()
naive.norm[i] <- rnorm(nrow(proper.format[[j]]),naive.norm.m[i],naive.norm.sd[i])
naive.lognorm[i] <- rnorm(nrow(proper.format[[j]]),naive.lognorm.m[i],naive.lognorm.sd[i]) %>% exp()
}

names(naive.norm) <- names(naive.lognorm) <- names(sim.df1) <- names(proper.format[[1]])

real.data <- proper.format[[j]] 
real.data$type <- "real"
sim.df1$type <- "logit-normal simulated"
naive.norm$type <- "normal simulated"
naive.lognorm$type <- "lognormal simulated"

df <- rbind(real.data,sim.df1,naive.norm,naive.lognorm) %>% as_tibble() %>% pivot_longer(cols=!type)
df %>% ggplot(aes(x=value,colour=type))+facet_wrap(.~name,scales = "free")+geom_density()
df
list.df.logit <- proper.format 
for (i in 1:29){
  list.df.logit[[i]] <- apply(proper.format[[i]],2,logit)
}
l <- list() 
for (i in 1:29){
l[[i]] <- list.df[[i]] %>% rowSums()
}
l <- unlist(l)
l
sd(unlist(l))
