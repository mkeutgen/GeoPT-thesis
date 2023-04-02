# Evaluating the BLR Mean Estimator


library(MASS)
library(tidyverse)
library(compositions)
library(xtable)
library(rootSolve)
library(EnvStats)
library(ggpubr)
library(quadprog)

setwd("src/sim/")

file.list <- list.files(pattern='*.RData')

list.sim <- lapply(file.list,readRDS)

names(list.sim) <- file.list

list.sim$sim_list_GeoPT1.RData
mean.compo <- readRDS()
mean.compo$X <- NULL
rownames(mean.compo) <- mean.compo$name
mean.compo$name <- NULL
mean.compo$U <- mean.compo$NA.
mean.compo$NA. <- NULL

rownames(mean.compo) <- gsub(pattern = " ",replacement = "",rownames(mean.compo))





# Computing the bias for several estimators of the mean for one dataset


true <- mean.compo.simp[[26]]

mu.parameters[26] %>% as_vector() %>% logitInv()
list.df <- list.sim$sim_list_GeoPT48.RData

mean.compo




# Custom Functions

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

blr.mean.fun <- function(dataset){
  output <-  dataset %>% as_tibble() %>% apply(MARGIN = 2,logit) %>% colMeans() %>% logitInv()
}
optimizing.blr.mean <- function(dataset){
  # mu_star is the initial mean estimate from the blr transformed data 
  mu_star <- dataset %>% as_tibble() %>% apply(MARGIN = 2,logit) %>% colMeans()
  # sigma 
  sigma <- dataset %>% apply(MARGIN=2,logit) %>% apply(MARGIN=2,sd)
  # Optimization : sum(logitInv(mu_star+alpha*sigma)) = 1
  f1 <- function(alpha)(1-sum(logitInv(mu_star+alpha*sigma)))
  # alpha : 
  alpha <- uniroot(f1, c(-1,1))$root
  mu_optimized <- logitInv(mu_star+alpha*sigma)
  return(mu_optimized)
} 



#alpha.v <- c()
#for (i in 1:length(list.df)){
#  alpha.v[i] <- optimizing.blr.mean(list.df[[i]])[2]
#}

# Comparing the performance of different estimators of the mean composition

opti.blr.df <- lapply(list.df,optimizing.blr.mean) %>% bind_rows() %>% mutate(type="opti.blr")

blr.df <- lapply(list.df,blr.mean.fun) %>% bind_rows() %>% mutate(type="blr.df")

col.mean.df <- lapply(list.df,colMeans) %>% bind_rows() %>% mutate(type="arithm.mean")

ilr.mean.df <- lapply(lapply(lapply(list.df,ilr),colMeans),ilrInv) %>% bind_rows() %>% mutate(type="ilr mean") 
names(ilr.mean.df) <- names(blr.df)

combined.df <- bind_rows(opti.blr.df,blr.df,col.mean.df,ilr.mean.df)


df.output <- combined.df %>% pivot_longer(cols = !type,names_to = "component",values_to = "value")

True.Df <- as_tibble(true) %>% pivot_longer(cols=everything())
names(True.Df) <- c("component","value")
 

ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.4)+
  facet_wrap(~component,scales = "free",ncol=3,shrink = T)+theme_bw()+
  scale_x_continuous(labels=scales::scientific)+
  scale_y_continuous(labels=scales::scientific)+
  theme(legend.position = "bottom")+
  geom_vline(aes(xintercept=value),data=True.Df,show.legend = T,size=.5,color="darkblue")

View(True.Df)

True.Df <- as_tibble(rownames="component",true)
True.Df
names(True.Df) <- c("component","value")



# Computing biases
arithm.bias <- data.frame(matrix(ncol=22,nrow=500))
blr.bias <-    data.frame(matrix(ncol=22,nrow=500))
ilr.bias <-    data.frame(matrix(ncol=22,nrow=500))
optim.blr.bias <- data.frame(matrix(ncol=22,nrow=500))
closed.blr.bias <- data.frame(matrix(ncol=22,nrow=500))

sweep(combined.df,true,"/")
sweep(combined.df, MARGIN = 2, unlist(true), FUN = `/`)


for (i in 1:m){
  arithm.bias[i,] <- df1[1:22][i,]/true
  blr.bias[i,] <- df2[1:22][i,]/true
  ilr.bias[i,] <- df3[1:22][i,]/true
  optim.blr.bias[i,] <- df4[1:22][i,]/true
  
}


estim.bias.arithm     <- apply(arithm.bias,MARGIN = 2,geometricmean)
estim.bias.blr        <- apply(blr.bias,MARGIN = 2,geometricmean)     
estim.bias.ilr        <- apply(ilr.bias,MARGIN = 2,geometricmean)
estim.bias.optim.blr  <- apply(optim.blr.bias,MARGIN = 2,geometricmean)
estim.bias.closed.blr <- apply(closed.blr.bias,MARGIN = 2,geometricmean)

estim.bias.arithm %>% geometricmean() #   0.397432
estim.bias.blr %>% geometricmean      #   0.3841558
estim.bias.ilr %>% geometricmean() #      0.3836476
estim.bias.optim.blr %>% geometricmean() # 0.382191


df.bias <- data.frame(estim.bias.arithm,estim.bias.blr,estim.bias.ilr,
                      estim.bias.optim.blr,estim.bias.closed.blr,
                      chem.element=names(true[1:21]))
df.bias %>%
  pivot_longer(cols=!chem.element) %>%
  ggplot(aes(x=chem.element,y=value,color=name))+
  geom_point(size=2,alpha=.4)+theme_bw()+scale_y_continuous(limits = c(0,2.5))
df.bias %>% pivot_longer(cols = !c("chem.element")) %>% 
  ggplot(aes(x=chem.element,y=value))+geom_point()+
  facet_wrap(.~name,nrow=3)+
  geom_hline(yintercept = 1,colour="red")+theme_bw()+scale_y_log10()









