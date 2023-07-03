###########################
### ESTIMATING BIAS #######
### FINAL FILE ############


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
blr.mean <- function(dataset){
  logit(dataset) %>% colMeans() %>% logitInv()  
} 


# Ilr Mean (Simplex)
ilr.mean <- function(dataset){
  ilr(dataset) %>% colMeans() %>% ilrInv()
}

# Geometric mean (R+) 
geometric.mean <- function(dataset){
  dataset %>% log() %>% colMeans() %>% exp()
}
# Arithmetic mean (R)
arithmetic.mean <- function(dataset){
  dataset  %>% colMeans()
}


# Replace negative functions
replace_negatives <- function(df) {
  df[!rowSums(df>1,na.rm = T),]
  df[df < 0] <- 1E-9
  df <- df[!apply(df, 1, function(x) any(x > 1)),]
  return(df)
}

## A. NO MISSING VALUES ###


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
   bias.df <- sweep(logit(combined.df[1:22]), MARGIN = 2, logit(true[1:22]), FUN = `-`)
   
   
   # MSE Dataframe
   mse.df <- (sweep(logit(combined.df[1:22]), 2, logit(true[1:22]), "-")^2)
   
   
   bias.df$type <- combined.df$type
   mse.df$type <- combined.df$type 
   
   bias.df %>% pivot_longer(cols=!type,names_to="chem.el") %>% group_by(type,chem.el)  

   output <- list(bias.df,mse.df)
 }


list.11.sim <- list(list.sim$sim_list_GeoPT1.RData,list.sim$sim_list_GeoPT2.RData,
list.sim$sim_list_GeoPT20.RData,list.sim$sim_list_GeoPT21.RData,
list.sim$sim_list_GeoPT22.RData,list.sim$sim_list_GeoPT25.RData,
list.sim$sim_list_GeoPT3.RData,list.sim$sim_list_GeoPT32.RData,
list.sim$sim_list_GeoPT34.RData,list.sim$sim_list_GeoPT35.RData,
list.sim$sim_list_GeoPT36.RData)
 
 list.11.true.comp <- list(mean.compo$GeoPT1.csv,   mean.compo$GeoPT2.csv,  
                           mean.compo$`GeoPT20 .csv`,mean.compo$`GeoPT21 .csv`,
                           mean.compo$`GeoPT22 .csv`,mean.compo$`GeoPT25 .csv`,mean.compo$GeoPT3.csv,
                           mean.compo$`GeoPT32 .csv`,mean.compo$`GeoPT34 .csv`,mean.compo$`GeoPT35 .csv`,
                           mean.compo$`GeoPT36 .csv`)
 
 output.nomval <- list()
 output.nomval.bias <- list()
 output.nomval.mse <- list()
 
 for (i in 1:11){
   output.nomval[[i]] <- geombias.mse.fun(true = list.11.true.comp[[i]],list.11.sim[[i]])
 }
 
 for (i in 1:11){
  output.nomval.bias[[i]] <- output.nomval[[i]][[1]]
 }
 
 for (i in 1:11){
   output.nomval.mse[[i]] <- output.nomval[[i]][[2]]
 }
 
 
 names(output.nomval.bias) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25",
                          "geopt3","geopt32","geopt34","geopt35","geopt36")                      
 names(output.nomval.mse) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25",
                                "geopt3","geopt32","geopt34","geopt35","geopt36")                      
 

meandf <- function(dataset){
   dataset <- as_tibble(dataset)
   df <- dataset[-ncol(dataset)]%>%
     rowMeans() %>% as_tibble()
   df[ncol(df)+1] <- dataset[ncol(dataset)]
   return(df)
 }
 
mergingdf.fun <- function(list.df){
  result <- lapply(list.df,meandf)
  names(result) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25","geopt3",
                     "geopt32","geopt34","geopt35","geopt36")                      
  combined_tib <- bind_rows(result, .id = "dataset")
}
 
nomval.bias <- mergingdf.fun(output.nomval.bias)
bias.plot.nomval <- nomval.bias %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean Bias")

## MSE 
nomval.mse <- mergingdf.fun(output.nomval.mse)
mse.plot.nomval <- nomval.mse %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean MSE")

# B. 10 % Missing Values 

# We need a function to generate missing values

fun.mval <- function(tibble,propmissing){
  output <- gen.mval.fun(tibble[-ncol(tibble)],prop = propmissing)
  output$U <- output$sum.vec + tibble$U
  output$sum.vec <- NULL
  return(output)
}



# And a function to impute those missing values 
imputed.df.fun <- function(dataset){
  # Imputed df in logit space
  dataset <- dataset[!rowSums(dataset>1,na.rm = T),]
  imputed.df.logit <- logit(dataset)  %>% impute_EM() 
  imputed.df <- logitInv(imputed.df.logit)
  sum.imp <- c()
  for (i in 1:nrow(dataset)){
    # For each row, sum the missing elements which are imputed and substract them from U
    sum.imp[i] <- sum(imputed.df[i,which(is.na(dataset[i,]))] )
  }
  # Yields a set of realistic compositions, from compositions with missing values 
  imputed.df$U <- abs(imputed.df$U - sum.imp)
  return(imputed.df)
}

gen.mval.fun <- function(df,prop){
  # inner function, select random elements and take their sum
  random_row_sum <- function(row) {
    # select 2 val
    selected_indices <- sample(length(row), floor(prop*length(row)))
    selected_values <- row[selected_indices]
    sum(selected_values)
    return(list(sum(selected_values),selected_indices))
  }
  
  list.sim$sim_list_GeoPT2.RData[[4]] %>% names() 
  
  
  # create a sum vector 
  sum.vec <- c()
  for (i in 1:nrow(df)){
    random.row.o <- random_row_sum(df[i,])
    sum.vec[i] <- random.row.o[[1]]
    df[i,c(random.row.o[[2]])] <- NA
  }
  output <- cbind(df,sum.vec) %>% as_tibble()
  return(output)
}




# 10 % missing values 
prop10.na.geopt1.list <- list()
prop10.na.geopt2.list <- list()
prop10.na.geopt20.list <- list()
prop10.na.geopt21.list <- list()
prop10.na.geopt22.list <- list()
prop10.na.geopt25.list <- list()
prop10.na.geopt3.list <- list()
prop10.na.geopt32.list <- list()
prop10.na.geopt34.list <- list()
prop10.na.geopt35.list <- list()
prop10.na.geopt36.list <- list()

# Generate the datasets with 10 % missing values
for (i in 1:200) { prop10.na.geopt1.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT1.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt2.list [[i]]   <- fun.mval(list.sim$sim_list_GeoPT2.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt20.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT20.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt21.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT21.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt22.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT22.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt25.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT25.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt3.list [[i]]  <- fun.mval(list.sim$sim_list_GeoPT3.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt32.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT32.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt34.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT34.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt35.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT35.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt36.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT36.RData[[i]],0.1)                                     }

list.11.sim.10mval <- list(prop10.na.geopt1.list,prop10.na.geopt2.list ,prop10.na.geopt20.list,prop10.na.geopt21.list,prop10.na.geopt22.list,prop10.na.geopt25.list,prop10.na.geopt3.list ,
                           prop10.na.geopt32.list,prop10.na.geopt34.list,
                           prop10.na.geopt35.list,prop10.na.geopt36.list)


list.11.sim.10mval <- lapply(list.11.sim.10mval, function(x) lapply(x, imputed.df.fun))

list.11.sim.10mval[[1]][[1]] %>% colMeans()
list.11.sim[[1]][[1]] %>% colMeans()


output.10mval <- list()
output.10mval.bias <- list()
output.10mval.mse <- list()

for (i in 1:11){
  output.10mval[[i]] <- geombias.mse.fun(true = list.11.true.comp[[i]],list.11.sim.10mval[[i]])
}

for (i in 1:11){
  output.10mval.bias[[i]] <- output.10mval[[i]][[1]]
}

for (i in 1:11){
  output.10mval.mse[[i]] <- output.10mval[[i]][[2]]
}


names(output.10mval.bias) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25",
                               "geopt3","geopt32","geopt34","geopt35","geopt36")                      

names(output.10mval.mse) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25",
                              "geopt3","geopt32","geopt34","geopt35","geopt36")                      


mval10.bias <- mergingdf.fun(output.10mval.bias)
bias.plot.10mval <- mval10.bias %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean Bias")

## MSE 
mval10.mse <- mergingdf.fun(output.10mval.mse)
mse.plot.10mval <- mval10.mse %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean MSE")


# 20 % missing values

prop20.na.geopt1.list <- list()
prop20.na.geopt2.list <- list()
prop20.na.geopt20.list <- list()
prop20.na.geopt21.list <- list()
prop20.na.geopt22.list <- list()
prop20.na.geopt25.list <- list()
prop20.na.geopt3.list <- list()
prop20.na.geopt32.list <- list()
prop20.na.geopt34.list <- list()
prop20.na.geopt35.list <- list()
prop20.na.geopt36.list <- list()

# Generate the datasets with 10 % missing values
for (i in 1:200) { prop20.na.geopt1.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT1.RData[[i]],  0.2)                                     }
for (i in 1:200) { prop20.na.geopt2.list [[i]]   <- fun.mval(list.sim$sim_list_GeoPT2.RData[[i]],0.2)                                     }
for (i in 1:200) { prop20.na.geopt20.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT20.RData[[i]],0.2)                                     }
for (i in 1:200) { prop20.na.geopt21.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT21.RData[[i]],0.2)                                     }
for (i in 1:200) { prop20.na.geopt22.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT22.RData[[i]],0.2)                                     }
for (i in 1:200) { prop20.na.geopt25.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT25.RData[[i]],0.2)                                     }
for (i in 1:200) { prop20.na.geopt3.list [[i]]  <- fun.mval(list.sim$sim_list_GeoPT3.RData[[i]], 0.2)                                     }
for (i in 1:200) { prop20.na.geopt32.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT32.RData[[i]],0.2)                                     }
for (i in 1:200) { prop20.na.geopt34.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT34.RData[[i]],0.2)                                     }
for (i in 1:200) { prop20.na.geopt35.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT35.RData[[i]],0.2)                                     }
for (i in 1:200) { prop20.na.geopt36.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT36.RData[[i]],0.2)                                     }

list.11.sim.20mval <- list(prop20.na.geopt1.list,prop20.na.geopt2.list ,prop20.na.geopt20.list,prop20.na.geopt21.list,prop20.na.geopt22.list,prop20.na.geopt25.list,prop20.na.geopt3.list ,
                           prop20.na.geopt32.list,prop20.na.geopt34.list,
                           prop20.na.geopt35.list,prop20.na.geopt36.list)


list.11.sim.20mval <- lapply(list.11.sim.20mval, function(x) lapply(x, imputed.df.fun))



output.20mval <- list()
output.20mval.bias <- list()
output.20mval.mse <- list()

for (i in 1:11){
  output.20mval[[i]] <- geombias.mse.fun(true = list.11.true.comp[[i]],list.11.sim.20mval[[i]])
}

for (i in 1:11){
  output.20mval.bias[[i]] <- output.20mval[[i]][[1]]
}

for (i in 1:11){
  output.20mval.mse[[i]] <- output.20mval[[i]][[2]]
}


names(output.20mval.bias) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25",
                               "geopt3","geopt32","geopt34","geopt35","geopt36")                      

names(output.20mval.mse) <- c("geopt1" ,"geopt2","geopt20","geopt21","geopt22","geopt25",
                              "geopt3","geopt32","geopt34","geopt35","geopt36")                      


mval20.bias <- mergingdf.fun(output.20mval.bias)
bias.plot.20mval <- mval20.bias %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean Bias")

## MSE 
mval20.mse <- mergingdf.fun(output.20mval.mse)
mse.plot.20mval <- mval20.mse %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean MSE")


# 50 % missing values

prop50.na.geopt1.list <- list()
prop50.na.geopt2.list <- list()
prop50.na.geopt20.list <- list()
prop50.na.geopt21.list <- list()
prop50.na.geopt22.list <- list()
prop50.na.geopt25.list <- list()
prop50.na.geopt3.list <- list()
prop50.na.geopt32.list <- list()
prop50.na.geopt34.list <- list()
prop50.na.geopt35.list <- list()
prop50.na.geopt36.list <- list()

# Generate the datasets with 10 % missing values
for (i in 1:200) { prop50.na.geopt1.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT1.RData[[i]],  0.5)                                     }
for (i in 1:200) { prop50.na.geopt2.list [[i]]   <- fun.mval(list.sim$sim_list_GeoPT2.RData[[i]],0.5)                                     }
for (i in 1:200) { prop50.na.geopt20.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT20.RData[[i]],0.5)                                     }
for (i in 1:200) { prop50.na.geopt21.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT21.RData[[i]],0.5)                                     }
for (i in 1:200) { prop50.na.geopt22.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT22.RData[[i]],0.5)                                     }
for (i in 1:200) { prop50.na.geopt25.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT25.RData[[i]],0.5)                                     }
for (i in 1:200) { prop50.na.geopt3.list [[i]]  <- fun.mval(list.sim$sim_list_GeoPT3.RData[[i]], 0.5)                                     }
for (i in 1:200) { prop50.na.geopt32.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT32.RData[[i]],0.5)                                     }
for (i in 1:200) { prop50.na.geopt34.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT34.RData[[i]],0.5)                                     }
for (i in 1:200) { prop50.na.geopt35.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT35.RData[[i]],0.5)                                     }
for (i in 1:200) { prop50.na.geopt36.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT36.RData[[i]],0.5)                                     }

list.11.sim.50mval <- list(prop50.na.geopt1.list,prop50.na.geopt2.list ,prop50.na.geopt20.list,prop50.na.geopt21.list,prop50.na.geopt22.list,prop50.na.geopt25.list,prop50.na.geopt3.list ,
                           prop50.na.geopt32.list,prop50.na.geopt34.list,
                           prop50.na.geopt35.list,prop50.na.geopt36.list)


list.11.sim.50mval <- lapply(list.11.sim.50mval, function(x) lapply(x, imputed.df.fun))



output.50mval <- list()
output.50mval.bias <- list()
output.50mval.mse <- list()

for (i in 1:11){
  output.50mval[[i]] <- geombias.mse.fun(true = list.11.true.comp[[i]],list.11.sim.50mval[[i]])
}

for (i in 1:11){
  output.50mval.bias[[i]] <- output.50mval[[i]][[1]]
}

for (i in 1:11){
  output.50mval.mse[[i]] <- output.50mval[[i]][[2]]
}


names(output.50mval.bias) <- c("geopt1" ,"geopt2","geopt50","geopt21","geopt22","geopt25",
                               "geopt3","geopt32","geopt34","geopt35","geopt36")                      

names(output.50mval.mse) <- c("geopt1" ,"geopt2","geopt50","geopt21","geopt22","geopt25",
                              "geopt3","geopt32","geopt34","geopt35","geopt36")                      


mval50.bias <- mergingdf.fun(output.50mval.bias)
bias.plot.50mval <- mval50.bias %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean Bias")

## MSE 
mval50.mse <- mergingdf.fun(output.50mval.mse)
mse.plot.50mval <- mval50.mse %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean MSE")



# 70 % missing values

prop70.na.geopt1.list <- list()
prop70.na.geopt2.list <- list()
prop70.na.geopt20.list <- list()
prop70.na.geopt21.list <- list()
prop70.na.geopt22.list <- list()
prop70.na.geopt25.list <- list()
prop70.na.geopt3.list <- list()
prop70.na.geopt32.list <- list()
prop70.na.geopt34.list <- list()
prop70.na.geopt35.list <- list()
prop70.na.geopt36.list <- list()

# Generate the datasets with 10 % missing values
for (i in 1:200) { prop70.na.geopt1.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT1.RData[[i]],  0.5)                                     }
for (i in 1:200) { prop70.na.geopt2.list [[i]]   <- fun.mval(list.sim$sim_list_GeoPT2.RData[[i]],0.5)                                     }
for (i in 1:200) { prop70.na.geopt20.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT20.RData[[i]],0.5)                                     }
for (i in 1:200) { prop70.na.geopt21.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT21.RData[[i]],0.5)                                     }
for (i in 1:200) { prop70.na.geopt22.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT22.RData[[i]],0.5)                                     }
for (i in 1:200) { prop70.na.geopt25.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT25.RData[[i]],0.5)                                     }
for (i in 1:200) { prop70.na.geopt3.list [[i]]  <- fun.mval(list.sim$sim_list_GeoPT3.RData[[i]], 0.5)                                     }
for (i in 1:200) { prop70.na.geopt32.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT32.RData[[i]],0.5)                                     }
for (i in 1:200) { prop70.na.geopt34.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT34.RData[[i]],0.5)                                     }
for (i in 1:200) { prop70.na.geopt35.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT35.RData[[i]],0.5)                                     }
for (i in 1:200) { prop70.na.geopt36.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT36.RData[[i]],0.5)                                     }

list.11.sim.70mval <- list(prop70.na.geopt1.list,prop70.na.geopt2.list ,prop70.na.geopt20.list,prop70.na.geopt21.list,prop70.na.geopt22.list,prop70.na.geopt25.list,prop70.na.geopt3.list ,
                           prop70.na.geopt32.list,prop70.na.geopt34.list,
                           prop70.na.geopt35.list,prop70.na.geopt36.list)


list.11.sim.70mval <- lapply(list.11.sim.70mval, function(x) lapply(x, imputed.df.fun))



output.70mval <- list()
output.70mval.bias <- list()
output.70mval.mse <- list()

for (i in 1:11){
  output.70mval[[i]] <- geombias.mse.fun(true = list.11.true.comp[[i]],list.11.sim.70mval[[i]])
}

for (i in 1:11){
  output.70mval.bias[[i]] <- output.70mval[[i]][[1]]
}

for (i in 1:11){
  output.70mval.mse[[i]] <- output.70mval[[i]][[2]]
}


names(output.70mval.bias) <- c("geopt1" ,"geopt2","geopt70","geopt21","geopt22","geopt25",
                               "geopt3","geopt32","geopt34","geopt35","geopt36")                      

names(output.70mval.mse) <- c("geopt1" ,"geopt2","geopt70","geopt21","geopt22","geopt25",
                              "geopt3","geopt32","geopt34","geopt35","geopt36")                      


mval70.bias <- mergingdf.fun(output.70mval.bias)
bias.plot.70mval <- mval70.bias %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean Bias")

mval70.bias %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
    theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean Bias")

## MSE 
mval70.mse <- mergingdf.fun(output.70mval.mse)
mse.plot.70mval <- mval70.mse %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+labs(y="Boxplot of Mean MSE")



# Saving plots

# Bias plots
ggsave(plot = bias.plot.nomval,filename = "plots/bias_nomval.jpeg",dpi = "retina",width = 8,height = 6)
ggsave(plot = bias.plot.10mval,filename = "plots/bias_10mval.jpeg",dpi = "retina",width = 8,height = 6)
ggsave(plot = bias.plot.20mval,filename = "plots/bias_20mval.jpeg",dpi = "retina",width = 8,height = 6)
ggsave(plot = bias.plot.50mval,filename = "plots/bias_50mval.jpeg",dpi = "retina",width = 8,height = 6)
ggsave(plot = bias.plot.70mval,filename = "plots/bias_70mval.jpeg",dpi = "retina",width = 8,height = 6)

# MSE plots
ggsave(plot = mse.plot.nomval,filename = "plots/mse_nomval.jpeg",dpi = "retina",width = 8,height = 6)
ggsave(plot = mse.plot.10mval,filename = "plots/mse_10mval.jpeg",dpi = "retina",width = 8,height = 6)
ggsave(plot = mse.plot.20mval,filename = "plots/mse_20mval.jpeg",dpi = "retina",width = 8,height = 6)
ggsave(plot = mse.plot.50mval,filename = "plots/mse_50mval.jpeg",dpi = "retina",width = 8,height = 6)
ggsave(plot = mse.plot.70mval,filename = "plots/mse_70mval.jpeg",dpi = "retina",width = 8,height = 6)



