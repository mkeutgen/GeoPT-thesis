
###############################
## Generating Missing Values ##
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



file.list <- list.files(pattern='*.RData')
mean.compo <- readRDS("~/Documents/GeoPTManuscript/MathGeosciPaper/Code/meancompo.RData")
list.sim <- lapply(file.list,readRDS)
names(list.sim) <- file.list


# Define a custom function that selects random cells and adds them together
gen.mval.fun <- function(df,prop=0.2){
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

prop10.na.geopt2.list <- list()
prop10.na.geopt20.list <- list()
prop10.na.geopt22.list <- list()
prop10.na.geopt25.list <- list()
prop10.na.geopt32.list <- list()
prop10.na.geopt35.list <- list()
prop10.na.geopt36.list <- list()
prop10.na.geopt38.list <- list()
prop10.na.geopt41.list <- list()
prop10.na.geopt46.list <- list()

geopt2.prop10.corrected <- list()
for (i in 1:200){prop10.na.geopt2.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT2.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt20.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT20.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt22.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT22.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt25.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT25.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt32.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT32.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt35.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT35.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt36.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT36.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt38.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT38.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt41.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT41.RData[[i]],0.1)}
for (i in 1:200){prop10.na.geopt46.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT46.RData[[i]],0.1)}

setwd("~/Documents/GeoPTManuscript/MathGeosciPaper/Code/")
saveRDS(object = prop10.na.geopt2.list,file =  "missingvalues/prop10/prop10_geopt2.RData")
saveRDS(object = prop10.na.geopt20.list,file = "missingvalues/prop10/prop10_geopt20.RData")
saveRDS(object = prop10.na.geopt22.list,file = "missingvalues/prop10/prop10_geopt22.RData")
saveRDS(object = prop10.na.geopt25.list,file = "missingvalues/prop10/prop10_geopt25.RData")
saveRDS(object = prop10.na.geopt32.list,file = "missingvalues/prop10/prop10_geopt32.RData")
saveRDS(object = prop10.na.geopt35.list,file = "missingvalues/prop10/prop10_geopt35.RData")
saveRDS(object = prop10.na.geopt36.list,file = "missingvalues/prop10/prop10_geopt36.RData")
saveRDS(object = prop10.na.geopt38.list,file = "missingvalues/prop10/prop10_geopt38.RData")
saveRDS(object = prop10.na.geopt41.list,file = "missingvalues/prop10/prop10_geopt41.RData")
saveRDS(object = prop10.na.geopt46.list,file = "missingvalues/prop10/prop10_geopt46.RData")





## 
#####
replace_negatives <- function(df) {
  df[!rowSums(df>1,na.rm = T),]
  df[df < 0] <- 1E-9
  df <- df[!apply(df, 1, function(x) any(x > 1)),]
  return(df)
}




#### GeomBias MSE function #####


optimizing.blr.mean <- function(dataset, min = -1, max = 1) {
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


blr.mean <- function(imputed.dataset){
  logit(imputed.dataset) %>% colMeans() %>% logitInv()  
} 


geombias.mse.fun <- function(true,list.df){  
  
  #alpha.v <- c()
  #for (i in 1:length(list.df)){
  #  alpha.v[i] <- optimizing.blr.mean(list.df[[i]])[2]
  #}
  
  # Comparing the performance of different estimators of the mean composition

  list.df <- lapply(list.df,imputed.df.fun)
  
  
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
  #bias.df <- sweep(log(combined.df[1:22]), MARGIN = 2, log(true), FUN = `-`)
  # In R : 
  bias.df <- sweep(combined.df[1:22], MARGIN = 2, true, FUN = `-`)
  
  
  # MSE Dataframe
  
  mse.df <- (sweep(log(combined.df[1:22]), 2, log(true), "-")^2)
  
  
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


prop10.na.geopt2.list.o <-  geombias.mse.fun(true = mean.compo$GeoPT2.csv,list.df = prop10.na.geopt2.list )
prop10.na.geopt20.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT20 .csv`,list.df = prop10.na.geopt20.list)
prop10.na.geopt22.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT22 .csv`,list.df = prop10.na.geopt22.list)
prop10.na.geopt25.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT25 .csv`,list.df = prop10.na.geopt25.list)
prop10.na.geopt32.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT32 .csv`,list.df = prop10.na.geopt32.list)
prop10.na.geopt35.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT35 .csv`,list.df = prop10.na.geopt35.list)
prop10.na.geopt36.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT35 .csv`,list.df = prop10.na.geopt36.list)
#prop10.na.geopt38.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT38 .csv`,list.df = prop10.na.geopt38.list)
prop10.na.geopt41.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT41 .csv`,list.df = prop10.na.geopt41.list)
prop10.na.geopt46.list.o <- geombias.mse.fun(true = mean.compo$`GeoPT46 .csv`,list.df = prop10.na.geopt46.list)


geommeandf <- function(dataset){
  dataset <- as_tibble(dataset)
  df <- dataset[-ncol(dataset)]%>% log %>%
    rowMeans() %>% exp() %>% as_tibble()
  df[ncol(df)+1] <- dataset[ncol(dataset)]
  return(df)
}

list.output <- list(prop10.na.geopt2.list.o[[1]] ,
                    prop10.na.geopt20.list.o[[1]],
                    prop10.na.geopt22.list.o[[1]],
                    prop10.na.geopt25.list.o[[1]],
                    prop10.na.geopt32.list.o[[1]],
                    prop10.na.geopt35.list.o[[1]],
                    prop10.na.geopt36.list.o[[1]],
                    prop10.na.geopt41.list.o[[1]],
                    prop10.na.geopt46.list.o[[1]])
result <- lapply(list.dataframe,geommeandf)

names(result) <- c( "geopt2",
                   "geopt20",
                   "geopt22",
                   "geopt25",
                   "geopt32",
                   "geopt35",
                   "geopt36",
                   "geopt41",
                   "geopt46" )

combined_tib <- bind_rows(result, .id = "dataset")

combined_tib %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle=45,hjust=1))+
  labs(y="Empirical Mean Bias")



