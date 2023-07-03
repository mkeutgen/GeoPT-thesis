# Evaluating the BLR Mean Estimator


library(MASS)
library(tidyverse)
library(compositions)
library(xtable)
library(rootSolve)
library(EnvStats)
library(ggpubr)
library(quadprog)



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
geometric.mean.fun <- function(dataset){
  output <- dataset %>% apply(MARGIN = 2,geometricmean)
}


optimizing.blr.mean <- function(dataset, min = -10, max = 10) {
  # mu_star is the initial mean estimate from the blr transformed data 
  mu_star <- dataset %>% as_tibble() %>% apply(MARGIN = 2, logit) %>% colMeans()
  
  # sigma 
  sigma <- dataset %>% apply(MARGIN = 2, logit) %>% apply(MARGIN = 2, sd)
  
  # Optimization: sum(logitInv(mu_star + alpha * sigma)) = 1
  f2 <- function(alpha) 1-sum(logitInv(mu_star + alpha * sigma)) 
  
  alpha <- uniroot(f2, c(min, max))$root  # alpha
  
  mu_optimized <- logitInv(mu_star + alpha * sigma)
  
  return(mu_optimized)
}

dataset <- list.sim$sim_list_GeoPT1.RData[[11]]

for (i in 1:length(list.sim$sim_list_GeoPT1.RData)){
  optimizing.blr.mean(list.sim$sim_list_GeoPT1.RData[[i]], min = -0.9, max = 1)
}

optimizing.blr.mean <- function(dataset, min = -10, max = 10) {
  # mu_star is the initial mean estimate from the blr transformed data 
  mu_star <- dataset %>% as_tibble() %>% apply(MARGIN = 2, logit) %>% colMeans()
  
  # sigma 
  sigma <- dataset %>% apply(MARGIN = 2, logit) %>% apply(MARGIN = 2, sd)
  
  # Optimization: sum(logitInv(mu_star + alpha * sigma)) = 1
  f2 <- function(alpha) 1-sum(logitInv(mu_star + alpha * sigma)) 
  
  alpha <- tryCatch({
    uniroot(f2, c(min, max))$root
  }, error = function(e) {
    if (e$message == "f() values at end points not of opposite sign") {
      mid <- (min + max) / 2
      if (f2(mid) > 0) {
        max <- mid
      } else {
        min <- mid
      }
      uniroot(f2, c(min, max))$root
    } else {
      stop(e)
    }
  })
  
  mu_optimized <- logitInv(mu_star + alpha * sigma)
  
  return(mu_optimized)
}

# Test the function with different datasets




optimized.blr.mean.l <- list()


list.df <- list.sim[[1]]
list.optiblr.df <- list()

# For each of the 29 rocks
for (j in 1:29){
  # For each of the 200 simulated datasets
for (i in 1 :length(list.sim[[j]])){
  optimized.blr.mean.l[[i]] <- optimizing.blr.mean(list.sim[[j]][[i]])
}
list.optiblr.df[[j]] <- bind_rows(optimized.blr.mean.l[[i]])
}



geombias.mse.fun(mean.compo$GeoPT1.csv,sim_list_GeoPT1.csv)

alpha_fun(list.sim[[21]][[2]])

sapply(list.optiblr.df,colMeans) %>% colSums(is.na(df)) # 2, 3, 4, 9, 11, 18, 24, 29

names(list.optiblr.df) <- names(list.sim)
list.blr.opti.df <- list()
for (i in 1:29){
list.blr.opti.df[[i]] <- lapply(list.sim[[i]],blr.mean.fun) %>% bind_rows() %>% mutate(type="blr mean")
}


lapply(list.df,optimizing.blr.mean)

# READING DATA


setwd("src/sim/")

file.list <- list.files(pattern='*.RData')
mean.compo <- readRDS("~/Documents/GeoPTManuscript/src/meancompo.RData")
list.sim <- lapply(file.list,readRDS)

# loop through the outer list
for (i in seq_along(list.sim)) {
  
  # loop through the inner list
  for (j in seq_along(list.sim[[i]])) {
    
    # remove rows with elements above 1
    list.sim[[i]][[j]] <- list.sim[[i]][[j]][apply(list.sim[[i]][[j]], 1, function(x) all(x <= 1)), ]
    
  }
  
}


names(list.sim) <- file.list








# Computing the bias for several estimators of the mean for one dataset, here GeoPT48 :

true <- mean.compo$GeoPT2.csv
list.df <- list.sim$sim_list_GeoPT2.RData


geombias.mse.fun <- function(true,list.df,index){  
  for (i in 1:length(list.df)){
    names(list.df[[i]]) <- names(true)
  }

  
  #alpha.v <- c()
  #for (i in 1:length(list.df)){
  #  alpha.v[i] <- optimizing.blr.mean(list.df[[i]])[2]
  #}
  
  # Comparing the performance of different estimators of the mean composition
  
  opti.blr.df <- lapply(list.df,optimizing.blr.mean) %>% bind_rows() %>% mutate(type="optimized blr mean")
  
  
  blr.df <- lapply(list.df,blr.mean.fun) %>% bind_rows() %>% mutate(type="blr mean")
  
  l2.blr.df <- sapply(X = blr.df[-23], FUN = blr.min) %>% 
    as_tibble() %>% mutate(type ="l2 blr mean")
  
  col.mean.df <- lapply(list.df,colMeans) %>% bind_rows() %>% mutate(type="arithmetic mean")
  
  ilr.mean.df <- lapply(lapply(lapply(list.df,ilr),colMeans),ilrInv) %>% bind_rows() %>% mutate(type="ilr mean") 
  names(ilr.mean.df) <- names(blr.df)
  
  closed.blr <- lapply(list.df,blr.mean.fun) %>% bind_rows() %>% clo() %>% as_tibble() %>% mutate(type="closed blr mean")  
  
  

  geometric.mean <- lapply(list.df,geometric.mean.fun) %>% bind_rows() %>% mutate(type="geometric mean")
    
  
  combined.df <- bind_rows(opti.blr.df,blr.df,col.mean.df,ilr.mean.df,closed.blr,geometric.mean,l2.blr.df)
  

  # Bias Dataframe
  
  bias.df <- sweep(combined.df[1:22], MARGIN = 2, true, FUN = `/`)
  
  # MSE Dataframe
  
  mse.df <- exp(sweep(combined.df[1:22], 2, true, "-")^2)
  

  bias.df$type <- combined.df$type
  mse.df$type <- combined.df$type 
  
  
  bias.per.element <- bias.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
    group_by(type,chem.el) %>% summarize(geomean = geometricmean(value))
  
  
  bias.plot <- bias.per.element %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_point()+
    labs(y="Geometric Mean Bias",x="Chemical Element")+theme_bw()+
    theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,3,5))
  
  total.bias <- bias.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
    group_by(type) %>% summarize(geomean = geometricmean(value)) 
  
  
  mse.per.element <- mse.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
    group_by(type,chem.el) %>% summarize(geomean = geometricmean(value))
  
  mse.plot <- mse.per.element %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter()+
    labs(y="Geometric MSE",x="Chemical Element")+theme_bw()+
    theme(legend.position="bottom") + scale_y_log10()
  
  total.mse <- mse.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
    group_by(type) %>% summarize(geomean = geometricmean(value)) 
  

  output <- list(bias.df,total.bias,bias.per.element,bias.plot,mse.df,total.mse,mse.plot)
}


list.geombias <- list()
for (i in 1:29){
list.geombias[[i]] <- geombias.mse.fun(mean.compo[[i]],list.sim[[i]])[[2]] %>%
  mutate(name = names(mean.compo)[i] )
}

geombias.mse.fun(mean.compo[[i]],list.sim[[i]])[[2]] 









tbias.sim_list_GeoPT1.csv   <- geombias.mse.fun(mean.compo[[1]],list.sim[[1]])[[2]]
tbias.sim_list_GeoPT2.csv   <- geombias.mse.fun(mean.compo[[5]],list.sim[[5]])[[2]]
tbias.sim_list_GeoPT20.csv  <- geombias.mse.fun(mean.compo[[6]],list.sim[[6]])[[2]]
tbias.sim_list_GeoPT21.csv  <- geombias.mse.fun(mean.compo[[7]],list.sim[[7]])[[2]]
tbias.sim_list_GeoPT22.csv  <- geombias.mse.fun(mean.compo[[8]],list.sim[[8]])[[2]]
tbias.sim_list_GeoPT25.csv  <- geombias.mse.fun(mean.compo[[10]],list.sim[[10]])[[2]]
tbias.sim_list_GeoPT3.csv   <- geombias.mse.fun(mean.compo[[12]],list.sim[[12]])[[2]]
tbias.sim_list_GeoPT32.csv  <- geombias.mse.fun(mean.compo[[13]],list.sim[[13]])[[2]]
tbias.sim_list_GeoPT34.csv  <- geombias.mse.fun(mean.compo[[14]],list.sim[[14]])[[2]]
tbias.sim_list_GeoPT35.csv  <- geombias.mse.fun(mean.compo[[15]],list.sim[[15]])[[2]]
tbias.sim_list_GeoPT36.csv  <- geombias.mse.fun(mean.compo[[16]],list.sim[[16]])[[2]]
tbias.sim_list_GeoPT38A.csv<- geombias.mse.fun( mean.compo[[19]],list.sim[[19]])[[2]] 
tbias.sim_list_GeoPT39.csv  <- geombias.mse.fun(mean.compo[[20]],list.sim[[20]])[[2]]
tbias.sim_list_GeoPT39A.csv<- geombias.mse.fun( mean.compo[[21]],list.sim[[21]])[[2]] 
tbias.sim_list_GeoPT4.csv   <- geombias.mse.fun(mean.compo[[22]],list.sim[[22]])[[2]]
tbias.sim_list_GeoPT41.csv  <- geombias.mse.fun(mean.compo[[23]],list.sim[[23]])[[2]]
tbias.sim_list_GeoPT46.csv  <- geombias.mse.fun(mean.compo[[25]],list.sim[[25]])[[2]]
tbias.sim_list_GeoPT5.csv   <- geombias.mse.fun(mean.compo[[27]],list.sim[[27]])[[2]]

list.total.bias <- list(tbias.sim_list_GeoPT1.csv ,
     tbias.sim_list_GeoPT2.csv ,
     tbias.sim_list_GeoPT20.csv,
     tbias.sim_list_GeoPT21.csv,
     tbias.sim_list_GeoPT22.csv,
     tbias.sim_list_GeoPT25.csv,
     tbias.sim_list_GeoPT3.csv ,
     tbias.sim_list_GeoPT32.csv,
     tbias.sim_list_GeoPT34.csv,
     tbias.sim_list_GeoPT35.csv,
     tbias.sim_list_GeoPT36.csv,
     tbias.sim_list_GeoPT38A.csv,
     tbias.sim_list_GeoPT39.csv,
     tbias.sim_list_GeoPT39A.csv,
     tbias.sim_list_GeoPT4.csv ,
     tbias.sim_list_GeoPT41.csv,
     tbias.sim_list_GeoPT46.csv,
     tbias.sim_list_GeoPT5.csv  )


name_totalbias <- c("GeoPT1.csv", 
                "GeoPT2.csv", 
                "GeoPT20.csv", 
                "GeoPT21.csv", 
                "GeoPT22.csv", 
                "GeoPT25.csv", 
                "GeoPT3.csv",  
                "GeoPT32.csv", 
                "GeoPT34.csv", 
                "GeoPT35.csv", 
                "GeoPT36.csv", 
                "GeoPT38A.csv",
                "GeoPT39.csv", 
                "GeoPT39A.csv",
                "GeoPT4.csv",  
                "GeoPT41.csv", 
                "GeoPT46.csv", 
                "GeoPT5.csv")

for (i in 1:length(list.total.bias)){
list.total.bias[[i]] <- list.total.bias[[i]] %>% mutate(name = name_totalbias[i])
}

geombias.df <- do.call(rbind, list.total.bias)
unique(geombias.df$name)


geombias.df %>% 
  filter(!type %in% c("l2 blr mean")) %>%
  ggplot(aes(x = name, y = geomean, color = type)) + 
  geom_point() + scale_y_log10(breaks=c(1,1.1,1.2,1.3,1.4))+
  theme_bw()



##########################
# MSE ####################
##########################

tmse.sim_list_GeoPT1.csv   <- geombias.mse.fun(mean.compo[[1]],list.sim[[1]])[[6]]
tmse.sim_list_GeoPT2.csv   <- geombias.mse.fun(mean.compo[[5]],list.sim[[5]])[[6]]
tmse.sim_list_GeoPT20.csv  <- geombias.mse.fun(mean.compo[[6]],list.sim[[6]])[[6]]
tmse.sim_list_GeoPT21.csv  <- geombias.mse.fun(mean.compo[[7]],list.sim[[7]])[[6]]
tmse.sim_list_GeoPT22.csv  <- geombias.mse.fun(mean.compo[[8]],list.sim[[8]])[[6]]
tmse.sim_list_GeoPT25.csv  <- geombias.mse.fun(mean.compo[[10]],list.sim[[10]])[[6]]
tmse.sim_list_GeoPT3.csv   <- geombias.mse.fun(mean.compo[[12]],list.sim[[12]])[[6]]
tmse.sim_list_GeoPT32.csv  <- geombias.mse.fun(mean.compo[[13]],list.sim[[13]])[[6]]
tmse.sim_list_GeoPT34.csv  <- geombias.mse.fun(mean.compo[[14]],list.sim[[14]])[[6]]
tmse.sim_list_GeoPT35.csv  <- geombias.mse.fun(mean.compo[[15]],list.sim[[15]])[[6]]
tmse.sim_list_GeoPT36.csv  <- geombias.mse.fun(mean.compo[[16]],list.sim[[16]])[[6]]
tmse.sim_list_GeoPT38A.csv<- geombias.mse.fun( mean.compo[[19]],list.sim[[19]])[[6]] 
tmse.sim_list_GeoPT39.csv  <- geombias.mse.fun(mean.compo[[20]],list.sim[[20]])[[6]]
tmse.sim_list_GeoPT39A.csv<- geombias.mse.fun( mean.compo[[21]],list.sim[[21]])[[6]] 
tmse.sim_list_GeoPT4.csv   <- geombias.mse.fun(mean.compo[[22]],list.sim[[22]])[[6]]
tmse.sim_list_GeoPT41.csv  <- geombias.mse.fun(mean.compo[[23]],list.sim[[23]])[[6]]
tmse.sim_list_GeoPT46.csv  <- geombias.mse.fun(mean.compo[[25]],list.sim[[25]])[[6]]
tmse.sim_list_GeoPT5.csv   <- geombias.mse.fun(mean.compo[[27]],list.sim[[27]])[[6]]

list.total.mse <- list(tmse.sim_list_GeoPT1.csv ,
                        tmse.sim_list_GeoPT2.csv ,
                        tmse.sim_list_GeoPT20.csv,
                        tmse.sim_list_GeoPT21.csv,
                        tmse.sim_list_GeoPT22.csv,
                        tmse.sim_list_GeoPT25.csv,
                        tmse.sim_list_GeoPT3.csv ,
                        tmse.sim_list_GeoPT32.csv,
                        tmse.sim_list_GeoPT34.csv,
                        tmse.sim_list_GeoPT35.csv,
                        tmse.sim_list_GeoPT36.csv,
                        tmse.sim_list_GeoPT38A.csv,
                        tmse.sim_list_GeoPT39.csv,
                        tmse.sim_list_GeoPT39A.csv,
                        tmse.sim_list_GeoPT4.csv ,
                        tmse.sim_list_GeoPT41.csv,
                        tmse.sim_list_GeoPT46.csv,
                        tmse.sim_list_GeoPT5.csv  )


name_totalmse <- c("GeoPT1.csv", 
                    "GeoPT2.csv", 
                    "GeoPT20.csv", 
                    "GeoPT21.csv", 
                    "GeoPT22.csv", 
                    "GeoPT25.csv", 
                    "GeoPT3.csv",  
                    "GeoPT32.csv", 
                    "GeoPT34.csv", 
                    "GeoPT35.csv", 
                    "GeoPT36.csv", 
                    "GeoPT38A.csv",
                    "GeoPT39.csv", 
                    "GeoPT39A.csv",
                    "GeoPT4.csv",  
                    "GeoPT41.csv", 
                    "GeoPT46.csv", 
                    "GeoPT5.csv")

for (i in 1:length(list.total.mse)){
  list.total.mse[[i]] <- list.total.mse[[i]] %>% mutate(name = name_totalmse[i])
}

geommse.df <- do.call(rbind, list.total.mse)


geommse.df %>% 
  filter(!type %in% c("l2 blr mean")) %>%
  ggplot(aes(x = name, y = geomean, color = type)) + 
  geom_point() + scale_y_log10(breaks=c(1,1.1,1.2,1.3,1.4))+labs(y="Geometric MSE")+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))




geombias.mse.fun(mean.compo[[4]],list.sim[[4]])[[1]][-23] 
row.geom.mean <- apply(geombias.mse.fun(mean.compo[[4]],list.sim[[4]])[[1]][-23], 1, geometricmean)


geombias.df <- do.call(rbind, list.geombias)
unique(geombias.df$name)


geombias.df %>% 
   filter(!type %in% c("l2 blr mean")) %>%
   ggplot(aes(x = name, y = geomean, color = type)) + 
   geom_point(alpha=.8,size=2) + scale_y_log10(breaks=c(1,1.1,1.2,1.3,1.4))+
   theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))

tapply(geombias.df$geomean,INDEX = geombias.df$type,geometricmean





geombias.mse.fun(mean.compo[[1]],list.sim[[1]])[[2]]

geopt1.bias.o <- geombias.mse.fun(mean.compo$GeoPT1.csv,list.sim$sim_list_GeoPT1.RData)

geopt1.bias.o[[2]]



geopt2.bias.o <- geombias.mse.fun(mean.compo$GeoPT2.csv,list.sim$sim_list_GeoPT2.RData)

geopt2.bias.o

geopt2.bias.o



geopt48.bias.o <- geombias.mse.fun(mean.compo$`GeoPT48 .csv`,list.sim$sim_list_GeoPT48.RData)

# 3
geopt16.bias.o <- geombias.mse.fun(mean.compo$`GeoPT16 .csv`,list.sim$sim_list_GeoPT16.RData)
# 4
geopt19.bias.o <- geombias.mse.fun(mean.compo$`GeoPT19 .csv`,list.sim$sim_list_GeoPT19.RData)
# 5
geopt2.bias.o <- geombias.mse.fun(mean.compo$GeoPT2.csv,list.sim$sim_list_GeoPT2.RData)


ggarrange(p1,p2,p3,common.legend = T,nrow = 3,legend = "bottom")

# 6
geopt20.bias.o <- geombias.mse.fun(mean.compo$`GeoPT20 .csv`,list.sim$sim_list_GeoPT20.RData)
# 8
geopt22.bias.o <- geombias.mse.fun(mean.compo$`GeoPT22 .csv`,list.sim$sim_list_GeoPT22.RData)

geopt32.bias.o <- geombias.mse.fun(mean.compo$`GeoPT32 .csv`,list.sim$sim_list_GeoPT32.RData)

geopt46.bias.o <- geombias.mse.fun(mean.compo$`GeoPT46 .csv`,list.sim$sim_list_GeoPT46.RData)

geopt38a.bias.o <- geombias.mse.fun(mean.compo$GeoPT38A.csv,list.sim$sim_list_GeoPT38A.RData)

geopt22.bias.o[[2]]
geopt20.bias.o[[2]]
geopt32.bias.o[[2]]
geopt46.bias.o[[2]]
geopt2.bias.o[[2]]

p1 <- geopt2.bias.o[[3]] %>% filter(type %in% c("arithmetic mean","geometric mean","blr mean","ilr mean")) %>%
  ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT2 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))

p2 <- geopt20.bias.o[[3]] %>% filter(type %in% c("arithmetic mean","geometric mean","blr mean","ilr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT20 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))


geopt20.bias.o[[3]] %>% filter(type %in% c("arithmetic mean","geometric mean","blr mean","ilr mean")) %>% 
  ggplot(aes(x=chem.el,y=geomean,fill=type))+geom_col(position = "jitter")+
  labs(y="GeoPT20 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))


p3 <- geopt22.bias.o[[3]] %>% filter(type %in% c("arithmetic mean","geometric mean","blr mean","ilr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT 22 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))

p4 <- geopt32.bias.o[[3]] %>% filter(type %in% c("arithmetic mean","geometric mean","blr mean","ilr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT 32 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))

p5 <- geopt46.bias.o[[3]] %>% filter(type %in% c("arithmetic mean","geometric mean","blr mean","ilr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT 46 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))

p6 <- geopt38a.bias.o[[3]] %>% filter(type %in% c("arithmetic mean","geometric mean","blr mean","ilr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT 38A : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))


ggarrange(p1,p2,p3,p4,p5,p6,common.legend = T,legend="bottom")

list.bias.output <- list()
for (i in 1:length(mean.compo)){
  list.bias.output[[i]] <- geombias.fun(mean.compo[[i]],list.sim[[i]])
}

o <- geombias.fun(mean.compo[[1]],sim_list_GeoPT1.csv)


sim_list_GeoPT1.csv[[1]] %>% apply(,MARGIN=2,FUN=geometricmean)
list.bias.output[[1]] <- geombias.fun(mean.compo[[1]],list.sim[[1]])
list.bias.output[[2]] <- geombias.fun(mean.compo[[2]],list.sim[[2]])
list.bias.output[[3]] <- geombias.fun(mean.compo[[3]],list.sim[[3]])
list.bias.output[[4]] <- geombias.fun(mean.compo[[4]],list.sim[[4]])
list.bias.output[[5]] <- geombias.fun(mean.compo[[5]],list.sim[[5]])
list.bias.output[[6]] <- geombias.fun(mean.compo[[6]],list.sim[[6]])
list.bias.output[[7]] <- geombias.fun(mean.compo[[7]],list.sim[[7]])
list.bias.output[[8]] <- geombias.fun(mean.compo[[8]],list.sim[[8]])
list.bias.output[[9]] <- geombias.fun(mean.compo[[9]],list.sim[[9]])
list.bias.output[[10]] <- geombias.fun(mean.compo[[10]],list.sim[[10]])
list.bias.output[[21]] <- geombias.fun(mean.compo[[21]],list.sim[[21]])

list.bias.output[[5]]

list.bias.output[[26]] <- geombias.fun(mean.compo[[26]],list.sim[[26]])
list.bias.output[[27]] <- geombias.fun(mean.compo[[27]],list.sim[[27]])

list.bias.output$
list.bias.output$`GeoPT41 .csv`

length(list.bias.output) <- length(mean.compo)
names(list.bias.output) <- names(mean.compo)

list.bias.output[[8]]
list.bias.output[[21]]

i <- 8
list.bias.output[[i]] <- geombias.fun(mean.compo[[i]],list.sim[[i]])[[3]]

geopt19.bias.o
geopt19.bias.o[[3]]


geopt2.bias.o[[3]]$type %>% unique()

# 
p1 <- geopt2.bias.o[[3]] %>% filter(type %in% c("closed blr mean","l2 blr mean","ilr mean","optimized blr mean")) %>%
  ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT2 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))

p2 <- geopt20.bias.o[[3]] %>% filter(type %in% c("closed blr mean","l2 blr mean","ilr mean","optimized blr mean"))) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT20 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))


p3 <- geopt22.bias.o[[3]] %>% filter(type %in% c("closed blr mean","l2 blr mean","ilr mean","optimized blr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT 22 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))

p4 <- geopt32.bias.o[[3]] %>% filter(type %in% c("closed blr mean","l2 blr mean","ilr mean","optimized blr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT 32 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))

p5 <- geopt46.bias.o[[3]] %>% filter(type %in% c("closed blr mean","l2 blr mean","ilr mean","optimized blr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT 46 : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))

p6 <- geopt38a.bias.o[[3]] %>% filter(type %in% c("closed blr mean","l2 blr mean","ilr mean","optimized blr mean")) %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter(size=2)+
  labs(y="GeoPT 38A : Geometric Mean Bias",x="Chemical Element")+theme_bw()+
  theme(legend.position="bottom") + scale_y_log10(limits = c(3E-1, 5),breaks=c(0.3,0.3,1,2,3,5))


ggarrange(p1,p2,p3,p4,p5,p6,common.legend = T,legend="bottom")




