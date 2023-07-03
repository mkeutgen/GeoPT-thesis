 #############################
 ## ESTIMATION OF THE BIAS ###
 #############################
 
 setwd("MathGeosciPaper/Code/")
 
 library(MASS)
 library(tidyverse)
 library(compositions)
 library(xtable)
 library(rootSolve)
 library(EnvStats)
 library(ggpubr)
 library(quadprog)
 library(missMethods)

 #### Import Simulated Datasets ####

file.list <- list.files("sim/")
file.list
mean.compo <- readRDS("~/Documents/GeoPTManuscript/MathGeosciPaper/Code/meancompo.RData")
list.sim <- lapply(paste("sim/",file.list,sep = ""),readRDS)

names(list.sim) <- file.list

### Custom Functions needed to compute several estimators of the mean #####
 # Custom logit & logitINV functions 
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
 
 
 optimizing.blr.mean <- function(dataset, min = -10, max = 10) {
   # mu_star is the initial mean estimate from the blr transformed data 
   logit.df <- logit(dataset)
   
   mu_star <- logit.df %>% colMeans()
   
   # sigma 
   sigma <- logit.df %>% apply(MARGIN = 2, sd)
   
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
 
 
 
 # Blr mean (Hypercube)
 blr.mean <- function(imputed.dataset){
   logit(imputed.dataset) %>% colMeans() %>% logitInv()  
 } 
 
 
 # Ilr Mean (Simplex)
 ilr.mean <- function(dataset){
   ilr(dataset) %>% colMeans() %>% ilrInv()
 }
 
 # Geometric mean (R+) 
 geometric.mean <- function(imputed.dataset){
   imputed.dataset %>% log() %>% colMeans() %>% exp()
 }
 # Arithmetic mean (R)
 arithmetic.mean <- function(imputed.dataset){
   imputed.dataset  %>% colMeans()
 }
 
 
 ## 
 #####
 replace_negatives <- function(df) {
   df[!rowSums(df>1,na.rm = T),]
   df[df < 0] <- 1E-9
   df <- df[!apply(df, 1, function(x) any(x > 1)),]
   return(df)
 }
 

 
 
#### GeomBias MSE function #####

 
bias.fun <- function(true,list.df){
  list.df <- lapply(list.df,replace_negatives)
  
  opti.blr.df <- lapply(list.df,optimizing.blr.mean) %>% bind_rows() %>% mutate(type="optimized blr mean")
  
  
  blr.df <- lapply(list.df,blr.mean) %>% bind_rows() %>% mutate(type="blr mean")
  
  col.mean.df <- lapply(list.df,arithmetic.mean) %>% bind_rows() %>% mutate(type="arithmetic mean")
  
  closed.blr <- blr.df[,-ncol(blr.df)] %>% clo() %>% as_tibble() %>% mutate(type="closed blr mean")  
  
  ilr.mean.df <- lapply(lapply(lapply(list.df,ilr),colMeans),ilrInv) %>% bind_rows() %>% mutate(type="ilr mean") 
  names(ilr.mean.df) <- names(blr.df)
  
  
  geometric.mean <- lapply(list.df,geometric.mean) %>% bind_rows() %>% mutate(type="geometric mean")
  
  
  combined.df <- bind_rows(opti.blr.df,blr.df,col.mean.df,closed.blr,ilr.mean.df,geometric.mean)
  bias.df <- sweep(combined.df[1:22], MARGIN = 2, true[1:22], FUN = `-`)
  mse.df <- (sweep(combined.df[1:21], 2, log(true[1:21]), "-")^2)
  bias.df$type <- combined.df$type
  mse.df$type <- combined.df$type
  output <- list(combined.df,bias.df,mse.df)
  return(output)
}
 


 geombias.mse.fun <- function(true,list.df){  
   
   #alpha.v <- c()
   #for (i in 1:length(list.df)){
   #  alpha.v[i] <- optimizing.blr.mean(list.df[[i]])[2]
   #}
   
   # Comparing the performance of different estimators of the mean composition
  list.df <- lapply(list.df,replace_negatives)
  
  
   opti.blr.df <- lapply(list.df,optimizing.blr.mean) %>% bind_rows() %>% mutate(type="optimized blr mean")
   
   
   blr.df <- lapply(list.df,blr.mean) %>% bind_rows() %>% mutate(type="blr mean")
   
   col.mean.df <- lapply(list.df,arithmetic.mean) %>% bind_rows() %>% mutate(type="arithmetic mean")
   
   closed.blr <- blr.df[,-ncol(blr.df)] %>% clo() %>% as_tibble() %>% mutate(type="closed blr mean")  
   
   ilr.mean.df <- lapply(lapply(lapply(list.df,ilr),colMeans),ilrInv) %>% bind_rows() %>% mutate(type="ilr mean") 
   names(ilr.mean.df) <- names(blr.df)
   
   
   geometric.mean <- lapply(list.df,geometric.mean) %>% bind_rows() %>% mutate(type="geometric mean")
   
   
   combined.df <- bind_rows(opti.blr.df,blr.df,col.mean.df,closed.blr,ilr.mean.df,geometric.mean)
   
   # Bias Dataframe
   # In log coordinates :
   #bias.df <- sweep(log(combined.df[1:22]), MARGIN = 2, log(true), FUN = '-')
   # In R : 
   bias.df <- sweep(logit(combined.df[1:22]), MARGIN = 2, logit(true[1:22]), FUN = `-`)
   
   
   # MSE Dataframe
   
   mse.df <- (sweep(logit(combined.df[1:22]), 2, logit(true[1:22]), "-")^2)
   
   
   bias.df$type <- combined.df$type
   mse.df$type <- combined.df$type 
   
   bias.df %>% pivot_longer(cols=!type,names_to="chem.el") %>% group_by(type,chem.el)  
   
   
   bias.per.element <- bias.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
     group_by(type,chem.el) %>% summarize(geomean = mean(value))
   
   
   bias.per.element
   
   bias.plot <- bias.per.element %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter()+
     labs(y="Geometric Mean Bias",x="Chemical Element")+theme_bw()+
     theme(legend.position="bottom") + scale_y_log10()
   
   total.bias <- bias.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
     group_by(type) %>% summarize(geomean = mean(value)) 
   
   
   mse.per.element <- mse.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
     group_by(type,chem.el) %>% summarize(geomean = mean(value))
   
   mse.plot <- mse.per.element %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter()+
     labs(y="Geometric MSE",x="Chemical Element")+theme_bw()+
     theme(legend.position="bottom") + scale_y_log10()
   
   total.mse <- mse.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
     group_by(type) %>% summarize(geomean = mean(value)) 
   
   
   output <- list(bias.df,total.bias,bias.per.element,bias.plot,mse.df,total.mse,mse.plot)
 }


 
 
 
 
geopt1.o <- geombias.mse.fun(true = mean.compo$GeoPT1.csv,list.df = list.sim$sim_list_GeoPT1.RData)
geopt2.o <- geombias.mse.fun(true=mean.compo$GeoPT2.csv,list.df = list.sim$sim_list_GeoPT2.RData)
geopt20.o <- geombias.mse.fun(true=mean.compo$`GeoPT20 .csv`,list.df = list.sim$sim_list_GeoPT20.RData)
geopt21.o <- geombias.mse.fun(true=mean.compo$`GeoPT21 .csv`,list.df = list.sim$sim_list_GeoPT21.RData)
geopt22.o <- geombias.mse.fun(true=mean.compo$`GeoPT22 .csv`,list.df = list.sim$sim_list_GeoPT22.RData)
geopt25.o <- geombias.mse.fun(true=mean.compo$`GeoPT25 .csv`,list.df = list.sim$sim_list_GeoPT25.RData)
geopt3.o <- geombias.mse.fun(true=mean.compo$GeoPT3.csv,list.df = list.sim$sim_list_GeoPT3.RData)
geopt32.o <- geombias.mse.fun(true=mean.compo$`GeoPT32 .csv`,list.df = list.sim$sim_list_GeoPT32.RData)
geopt34.o <- geombias.mse.fun(true=mean.compo$`GeoPT34 .csv`,list.df = list.sim$sim_list_GeoPT34.RData)
geopt35.o <- geombias.mse.fun(true=mean.compo$`GeoPT35 .csv`,list.df = list.sim$sim_list_GeoPT35.RData)
geopt36.o <- geombias.mse.fun(true=mean.compo$`GeoPT36 .csv`,list.df = list.sim$sim_list_GeoPT36.RData)



meandf <- function(dataset){
   dataset <- as_tibble(dataset)
   df <- dataset[-ncol(dataset)]%>%
      rowMeans() %>% as_tibble()
   df[ncol(df)+1] <- dataset[ncol(dataset)]
   return(df)
}



list.dataframe <- list(geopt1.o [[1]],
                       geopt2.o [[1]],
                       geopt20.o[[1]],
                       geopt21.o[[1]],
                       geopt22.o[[1]],
                       geopt25.o[[1]],
                       geopt3.o [[1]],
                       geopt32.o[[1]],
                       geopt34.o[[1]],
                       geopt35.o[[1]],
                       geopt36.o[[1]])
                       
                       
            
result <- lapply(list.dataframe,meandf)

names(result) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25",
                   "geopt3","geopt32","geopt34","geopt35","geopt36")                      



# combine the tibbles using bind_rows(), and add a column to indicate the original tibble
combined_tib <- bind_rows(result, .id = "dataset")

combined_tib %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
   theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean Bias")

mean.df <- combined_tib %>% group_by(dataset,type) %>% summarize(mean=mean(value)) 
mean.df %>% ggplot(aes(x=dataset,y=mean,color=type))+geom_point()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Mean Bias")
############"
## MSE DF
############

list.dataframe <- list(geopt1.o [[5]],
                       geopt2.o [[5]],
                       geopt20.o[[5]],
                       geopt21.o[[5]],
                       geopt22.o[[5]],
                       geopt25.o[[5]],
                       geopt3.o [[5]],
                       geopt32.o[[5]],
                       geopt34.o[[5]],
                       geopt35.o[[5]],
                       geopt36.o[[5]])


result <- lapply(list.dataframe,meandf)

names(result) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25",
                   "geopt3","geopt32","geopt34","geopt35","geopt36")                      



# combine the tibbles using bind_rows(), and add a column to indicate the original tibble
combined_tib <- bind_rows(result, .id = "dataset")

combined_tib %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean MSE")+scale_y_log10()

mean.df <- combined_tib %>% group_by(dataset,type) %>% summarize(mean=mean(value)) 
mean.df %>% ggplot(aes(x=dataset,y=mean,color=type))+geom_point()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+scale_y_log10()+labs(y="Mean MSE")
