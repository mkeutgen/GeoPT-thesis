#######################
### BLR MEAN SIM V4 ###

# Note on this version :
# blr mean siv v1 had the drawback of always favoring ilr.mean.
# Now let's test the performance of the blr estimator on a real
# situation where X and T are sampled independently such that their
# sum exceed unity. 



# Logit and BLR Mean Function
logit <- function(x){log((x)/(1-x))}
# BLR Mean function
blr.mean <- function(x){
  # Map a vector from [0,1] to the real space. Takes the mean and return as its output the logit-inverse mean.
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}


library(MASS)
library(tidyverse)
library(compositions)

# Such a vector should be a mixture, a closed composition of major
# elements and some trace elements
set.seed(4)
major <- rexp(3,6)
major <- clo(major)
# Order of magnitude between major elements and trace is often 4 (ppm
# vs %. Therefore :
#excess <- 1/100
 excess <- 5/100
n.traces <- 6
# excess <- 10/100 
traces <- rexp(n.traces,6)
traces <- clo(traces)*excess
# Such that traces amount to one percent of total
sum(traces,major)
round(c(major,traces),3)
true <- c(major,traces) # %>% clo()


# We cannot use ilrInv to generate compositions with these parameters
# because it would artificially close compositions since ilrInv is a 
# mapping from the Euclidean Space to the SIMPLEX. WE can however cbind
# MV from major the and MV from traces

# number of individuals in each simulation
n <- 100
# m number of simulations
m <- 100


list.majorcomp <- list()
list.tracescomp <- list()
list.df <- list()

for (i in 1:m){
  list.majorcomp[[i]] <- mvrnorm(n,unclass(ilr(major)),diag(abs(ilr(major))) ) %>%
    ilrInv() %>% as_tibble()
  names(list.majorcomp[[i]]) <- seq(1:4) %>% as.character()
  
  list.tracescomp[[i]] <- mvrnorm(n,unclass(ilr(traces)),diag(abs(ilr(traces)))) %>%
    ilrInv() %>% as_tibble()*excess
  
  names(list.tracescomp[[i]]) <- seq(from=5,to=length(traces)) %>% as.character()
  
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
true
# True.Df$component <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14")
df.output$component <- substr(df.output$component,2,length(df.output$component))
df.output$component <- factor(df.output$component,levels=as.character(seq(1:14)) )

True.Df$component <- factor(True.Df$component,levels=as.character(seq(1:14)) )

ggplot(df.output,aes(x=value))+geom_density(aes(fill=type),alpha=.7)+
  facet_wrap(~component,scales = "free")+theme_bw()+
  geom_vline(aes(xintercept=value),data=True.Df,show.legend = T,size=.5,color="darkblue")+
  theme(legend.position = "bottom")


