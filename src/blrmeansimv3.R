#######################
### BLR MEAN SIM V3 ###

# Note on this version :
# blr mean siv v1 had the drawback of always favoring ilr.mean.
# Now let's test the performance of the blr estimator on a real
# situation where X and T are sampled independently such that their
# sum exceed unity. 

library(MASS)
library(tidyverse)
library(compositions)

# Such a vector should be a mixture, a closed composition of major
# elements and some trace elements
set.seed(4)
major <- rexp(4,6)
major <- clo(major)
# Order of magnitude between major elements and trace is often 4 (ppm
# vs %. Therefore :

traces <- rexp(10,6)
traces <- clo(traces)/1E4

# parts 
sum(traces,major)

true <- c(major,traces)

# We cannot use ilrInv to generate compositions with these parameters
# because it would artificially close compositions since ilrInv is a 
# mapping from the Euclidean Space to the SIMPLEX. WE can however cbind
# MV from major the and MV from traces
set.seed(4)
major <- rexp(4,6)
major <- clo(major)
# Order of magnitude between major elements and trace is often 4 (ppm
# vs %. Therefore :

traces <- rexp(10,6)
traces <- clo(traces)/1E8

true <- c(major,traces)
true
# number of individuals in each simulation
n <- 100
# m number of simulations
m <- 100
majorcomp <- mvrnorm(n,unclass(ilr(major)),diag(abs(ilr(major))) ) %>%
  ilrInv() %>% as_tibble()

tracescomp <- mvrnorm(n,unclass(ilr(traces)),diag(abs(ilr(traces)))) %>%
  ilrInv() %>% as_tibble()/1E4

names(majorcomp) <- seq(1:4) %>% as.character()
names(tracescomp) <- seq(from=5,to=9) %>% as.character()
df <- bind_cols(majorcomp,tracescomp)
names(df) <- seq(1:length(df)) %>% as.character()

arithm.mean.vec <- df %>% colMeans()
blr.mean.vec <- apply(df,2,blr.mean)
ilr.mean.vec <- ilr(df) %>% colMeans() %>% ilrInv() %>% unclass()
c(major,traces)
arithm.mean.vec
ilr.mean.vec
# Simulation, list of 100 dataframes : 

list.majorcomp <- list()
list.tracescomp <- list()
list.df <- list()

for (i in 1:m){
  list.majorcomp[[i]] <- mvrnorm(n,unclass(ilr(major)),diag(abs(ilr(major))) ) %>%
    ilrInv() %>% as_tibble()
  names(list.majorcomp[[i]]) <- seq(1:4) %>% as.character()
  
  list.tracescomp[[i]] <- mvrnorm(n,unclass(ilr(traces)),diag(abs(ilr(traces)))) %>%
    ilrInv() %>% as_tibble()/1E4
  
  names(list.tracescomp[[i]]) <- seq(from=5,to=14) %>% as.character()
  
  list.df[[i]] <- bind_cols(list.majorcomp[[i]],list.tracescomp[[i]])
  names(list.df[[i]]) <- seq(1:length(list.df[[i]])) %>% as.character()  
}

# Now for each of these dataframes, compute the 3 estimators for centrality.
# Gather the results in dataframes of dim number of parts x number of sim

list.matrices <- list()
for (i in 1:3){
  list.matrices[[i]] <- matrix(ncol=length(true),nrow=m)
  
}
for (i in 1:m){
  list.matrices[[1]][i,] <- list.df[[i]] %>% colMeans()
  list.matrices[[2]][i,] <- apply(list.df[[i]],2,blr.mean)
  list.matrices[[3]][i,]  <- ilr(list.df[[i]]) %>% colMeans() %>% 
    ilrInv() %>% unclass()
}

df1 <- list.matrices[[1]] %>% as_tibble() %>% mutate(type="arithmetic mean")
df2 <- list.matrices[[2]] %>% as_tibble() %>% mutate(type="blr mean")
df3 <- list.matrices[[3]] %>% as_tibble() %>% mutate(type="ilr mean")

df <- bind_rows(df1,df2,df3)
df.output <- df %>% pivot_longer(cols = !type,names_to = "component",values_to = "value")

ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free")+theme_bw()
True.Df <- as_tibble(rownames="component",true)

# True.Df$component <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14")
df.output$component <- substr(df.output$component,2,length(df.output$component))
df.output$component <- factor(df.output$component,levels=as.character(seq(1:14)) )

True.Df$component <- factor(True.Df$component,levels=as.character(seq(1:14)) )

ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free")+theme_bw()+
  geom_vline(aes(xintercept=value),data=True.Df,show.legend = T,size=.5,color="darkred")+
  theme(legend.position = "bottom")
