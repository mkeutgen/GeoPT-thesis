
###############################
# Missing Values 10 % Mval2 ###
###############################

# Load libraries

library(MASS)
library(tidyverse)
library(compositions)
library(xtable)
library(rootSolve)
library(EnvStats)
library(ggpubr)
library(quadprog)
library(missMethods)
library(Amelia)

mean.compo <- readRDS("~/Documents/GeoPTManuscript/MathGeosciPaper/Code/meancompo.RData")

setwd("/home/maxime/Documents/GeoPTManuscript/MathGeosciPaper/Code/sim/")

# Read the original data files (without missing values)

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
geometricmean.fun <- function(imputed.dataset){
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



file.list <- list.files(pattern='*.RData')
mean.compo <- readRDS("~/Documents/GeoPTManuscript/MathGeosciPaper/Code/meancompo.RData")
list.sim <- lapply(file.list,readRDS)
names(list.sim) <- file.list


# Define a custom function that selects random cells and adds them together
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


# Create a function to randomly replace a fraction of the non unkown elements in each column with missing values
# Then, this function add the content of the values now set to missing to the unknown mass fraction
# This should affect the bias of the unknown part, but not of other parts
tibble <- list.sim$sim_list_GeoPT2.RData[[1]]

fun.mval <- function(tibble,propmissing=0.2){
  output <- gen.mval.fun(tibble[-ncol(tibble)],prop = propmissing)
  output$U <- output$sum.vec + tibble$U
  output$sum.vec <- NULL
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

for (i in 1:200) { prop10.na.geopt1.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT1.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt2.list[[i]]   <- fun.mval(list.sim$sim_list_GeoPT2.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt20.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT20.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt21.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT21.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt22.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT22.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt25.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT25.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt3.list [[i]]  <- fun.mval(list.sim$sim_list_GeoPT3.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt32.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT32.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt34.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT34.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt35.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT35.RData[[i]],0.1)                                     }
for (i in 1:200) { prop10.na.geopt36.list[[i]]  <- fun.mval(list.sim$sim_list_GeoPT36.RData[[i]],0.1)                                     }


setwd("~/Documents/GeoPTManuscript/MathGeosciPaper/Code/")

saveRDS(object = prop10.na.geopt1.list  , file = "missingvalues/prop10/prop10_geopt1.RData")
saveRDS(object = prop10.na.geopt2.list  , file = "missingvalues/prop10/prop10_geopt2.RData")
saveRDS(object = prop10.na.geopt20.list , file = "missingvalues/prop10/prop10_geopt20.RData")
saveRDS(object = prop10.na.geopt21.list , file = "missingvalues/prop10/prop10_geopt21.RData")
saveRDS(object = prop10.na.geopt22.list , file = "missingvalues/prop10/prop10_geopt22.RData")
saveRDS(object = prop10.na.geopt25.list , file = "missingvalues/prop10/prop10_geopt25.RData")
saveRDS(object = prop10.na.geopt3.list  , file = "missingvalues/prop10/prop10_geopt3.RData")
saveRDS(object = prop10.na.geopt32.list , file = "missingvalues/prop10/prop10_geopt32.RData")
saveRDS(object = prop10.na.geopt34.list , file = "missingvalues/prop10/prop10_geopt34.RData")
saveRDS(object = prop10.na.geopt35.list , file = "missingvalues/prop10/prop10_geopt35.RData")
saveRDS(object = prop10.na.geopt36.list , file = "missingvalues/prop10/prop10_geopt36.RData")

geombias.mse.fun <- function(true,list.df){  
  
  #alpha.v <- c()
  #for (i in 1:length(list.df)){
  #  alpha.v[i] <- optimizing.blr.mean(list.df[[i]])[2]
  #}
  
  # Comparing the performance of different estimators of the mean composition
  list.df <- lapply(list.df,imputed.df.fun)
  list.df <- lapply(list.df,replace_negatives)
  
  
  opti.blr.df <- lapply(list.df,optimizing.blr.mean) %>% bind_rows() %>% mutate(type="optimized blr mean")
  
  
  blr.df <- lapply(list.df,blr.mean) %>% bind_rows() %>% mutate(type="blr mean")
  
  col.mean.df <- lapply(list.df,arithmetic.mean) %>% bind_rows() %>% mutate(type="arithmetic mean")
  
  closed.blr <- blr.df[,-ncol(blr.df)] %>% clo() %>% as_tibble() %>% mutate(type="closed blr mean")  
  
  ilr.mean.df <- lapply(lapply(lapply(list.df,ilr),colMeans),ilrInv) %>% bind_rows() %>% mutate(type="ilr mean") 
  names(ilr.mean.df) <- names(blr.df)
  
  geometric.mean <- lapply(list.df,geometricmean.fun) %>% bind_rows() %>% mutate(type="geometric mean")
  
  
  combined.df <- bind_rows(opti.blr.df,blr.df,col.mean.df,closed.blr,ilr.mean.df,geometric.mean)
  
  # Bias Dataframe
  # In log coordinates :
  #bias.df <- sweep(log(combined.df[1:22]), MARGIN = 2, log(true), FUN = '-')
  # In R : 
  bias.df <- sweep(logit(combined.df[1:22] ), MARGIN = 2, logit(true[1:22]), FUN = `-`)
  
  
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




## 
#####
replace_negatives <- function(df) {
  df[!rowSums(df>1,na.rm = T),]
  df[df < 0] <- 1E-9
  df <- df[!apply(df, 1, function(x) any(x > 1)),]
  return(df)
}


geopt1.o <- geombias.mse.fun(true = mean.compo$GeoPT1.csv,  list.df  = prop10.na.geopt1.list)
geopt2.o <- geombias.mse.fun(true=mean.compo$GeoPT2.csv,    list.df  = prop10.na.geopt2.list)
geopt20.o <- geombias.mse.fun(true=mean.compo$`GeoPT20 .csv`,list.df = prop10.na.geopt20.list)
geopt21.o <- geombias.mse.fun(true=mean.compo$`GeoPT21 .csv`,list.df = prop10.na.geopt21.list)
geopt22.o <- geombias.mse.fun(true=mean.compo$`GeoPT22 .csv`,list.df = prop10.na.geopt22.list)
geopt25.o <- geombias.mse.fun(true=mean.compo$`GeoPT25 .csv`,list.df = prop10.na.geopt25.list)
geopt3.o <- geombias.mse.fun(true=mean.compo$GeoPT3.csv,     list.df = prop10.na.geopt3.list)
geopt32.o <- geombias.mse.fun(true=mean.compo$`GeoPT32 .csv`,list.df = prop10.na.geopt32.list)
geopt34.o <- geombias.mse.fun(true=mean.compo$`GeoPT34 .csv`,list.df = prop10.na.geopt34.list)
geopt35.o <- geombias.mse.fun(true=mean.compo$`GeoPT35 .csv`,list.df = prop10.na.geopt35.list)
geopt36.o <- geombias.mse.fun(true=mean.compo$`GeoPT36 .csv`,list.df = prop10.na.geopt36.list)



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
