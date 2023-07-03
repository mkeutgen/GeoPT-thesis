
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




# 20 % missing values 

prop20.na.geopt2.list <- list()
prop20.na.geopt20.list <- list()
prop20.na.geopt22.list <- list()
prop20.na.geopt25.list <- list()
prop20.na.geopt32.list <- list()
prop20.na.geopt35.list <- list()
prop20.na.geopt36.list <- list()
prop20.na.geopt38.list <- list()
prop20.na.geopt41.list <- list()
prop20.na.geopt46.list <- list()

geopt2.prop20.corrected <- list()
for (i in 1:200){
prop20.na.geopt2.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT2.RData[[i]],0.2)
}


for (i in 1:length(sim_list_GeoPT2.csv)){ 
prop20.na.geopt2.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT2.RData[[i]],0.2)
}

for (i in 1:200){prop20.na.geopt20.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT20.RData[[i]],0.2)}
for (i in 1:200){prop20.na.geopt22.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT22.RData[[i]],0.2)}
for (i in 1:200){prop20.na.geopt25.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT25.RData[[i]],0.2)}
for (i in 1:200){prop20.na.geopt32.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT32.RData[[i]],0.2)}
for (i in 1:200){prop20.na.geopt35.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT35.RData[[i]],0.2)}
for (i in 1:200){prop20.na.geopt36.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT36.RData[[i]],0.2)}
for (i in 1:200){prop20.na.geopt38.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT38.RData[[i]],0.2)}
for (i in 1:200){prop20.na.geopt41.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT41.RData[[i]],0.2)}
for (i in 1:200){prop20.na.geopt46.list[[i]] <- fun.mval(list.sim$sim_list_GeoPT46.RData[[i]],0.2)}

setwd("~/Documents/GeoPTManuscript/MathGeosciPaper/Code/")
saveRDS(object = prop20.na.geopt2.list,file = "missingvalues/prop20/prop20_geopt2.RData")
saveRDS(object = prop20.na.geopt20.list,file = "missingvalues/prop20/prop20_geopt20.RData")
saveRDS(object = prop20.na.geopt22.list,file = "missingvalues/prop20/prop20_geopt22.RData")
saveRDS(object = prop20.na.geopt25.list,file = "missingvalues/prop20/prop20_geopt25.RData")
saveRDS(object = prop20.na.geopt32.list,file = "missingvalues/prop20/prop20_geopt32.RData")
saveRDS(object = prop20.na.geopt35.list,file = "missingvalues/prop20/prop20_geopt35.RData")
saveRDS(object = prop20.na.geopt36.list,file = "missingvalues/prop20/prop20_geopt36.RData")
saveRDS(object = prop20.na.geopt38.list,file = "missingvalues/prop20/prop20_geopt38.RData")
saveRDS(object = prop20.na.geopt41.list,file = "missingvalues/prop20/prop20_geopt41.RData")
saveRDS(object = prop20.na.geopt46.list,file = "missingvalues/prop20/prop20_geopt46.RData")

######
## Part 2 : estimating the bias when unknown elements are present
#####
replace_negatives <- function(df) {
  df[!rowSums(df>1,na.rm = T),]
  df[df < 0] <- 1E-9
  df <- df[!apply(df, 1, function(x) any(x > 1)),]
  return(df)
}


# Apply the function to each data frame in the list



# Given a dataset with missing values :
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





optimizing.blr.mean.na <- function(imputed.dataset, min = -10, max = 10) {
  # mu_star is the initial mean estimate from the blr transformed data 
  logit.imputed.df <- logit(imputed.dataset)
  
  mu_star <- logit.imputed.df %>% colMeans()
  
  # sigma 
  sigma <- logit.imputed.df %>% apply(MARGIN = 2, sd)
  
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
blr.mean.na <- function(imputed.dataset){
  logit(imputed.dataset) %>% colMeans() %>% logitInv()  
} 


# Ilr Mean (Simplex)
ilr.mean <- function(dataset){
  ilr(dataset) %>% colMeans() %>% ilrInv()
}

# Geometric mean (R+) 
geometric.mean.na <- function(imputed.dataset){
  imputed.dataset %>% log() %>% colMeans() %>% exp()
}
# Arithmetic mean (R)
arithmetic.mean.na <- function(imputed.dataset){
  imputed.dataset  %>% colMeans()
}



geombias.mse.na.fun <- function(true,list.df){  
 
  #alpha.v <- c()
  #for (i in 1:length(list.df)){
  #  alpha.v[i] <- optimizing.blr.mean(list.df[[i]])[2]
  #}
  
  # Comparing the performance of different estimators of the mean composition

  imputed.df <- lapply(list.df,imputed.df.fun)
  
  opti.blr.df <- lapply(imputed.df,optimizing.blr.mean.na) %>% bind_rows() %>% mutate(type="optimized blr mean")
  
  
  blr.df <- lapply(imputed.df,blr.mean.na) %>% bind_rows() %>% mutate(type="blr mean")
  
  col.mean.df <- lapply(imputed.df,arithmetic.mean.na) %>% bind_rows() %>% mutate(type="arithmetic mean")

  
  closed.blr <- blr.df[,-ncol(blr.df)] %>% clo() %>% as_tibble() %>% mutate(type="closed blr mean")  
  
  ilr.mean.df <- lapply(lapply(lapply(list.df,ilr),colMeans),ilrInv) %>% bind_rows() %>% mutate(type="ilr mean") 
  names(ilr.mean.df) <- names(blr.df)
  
  
  geometric.mean <- lapply(imputed.df,geometric.mean.na) %>% bind_rows() %>% mutate(type="geometric mean")
  
  
  combined.df <- bind_rows(opti.blr.df,blr.df,col.mean.df,closed.blr,ilr.mean.df,geometric.mean)

  # Bias Dataframe
  
  bias.df <- sweep(combined.df[1:22], MARGIN = 2, true, FUN = `/`)
  
  
  # MSE Dataframe
  
  mse.df <- exp(sweep(combined.df[1:22], 2, true, "-")^2)
  
  
  bias.df$type <- combined.df$type
  mse.df$type <- combined.df$type 
  
  bias.df %>% pivot_longer(cols=!type,names_to="chem.el") %>% group_by(type,chem.el)  
  
  
  bias.per.element <- bias.df %>% pivot_longer(cols =!type,names_to = "chem.el") %>%
    group_by(type,chem.el) %>% summarize(geomean = geometricmean(value))
  
  
  bias.per.element
  
  bias.plot <- bias.per.element %>% ggplot(aes(x=chem.el,y=geomean,color=type))+geom_jitter()+
    labs(y="Geometric Mean Bias",x="Chemical Element")+theme_bw()+
    theme(legend.position="bottom") + scale_y_log10()
  
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

mean.compo$GeoPT2.csv
prop20.na.geopt2.list[[1]] %>% colMeans(na.rm=T)
output.geopt2 <- geombias.mse.na.fun(true = mean.compo$GeoPT2.csv,list.df = prop20.na.geopt2.list)
output.geopt20 <- geombias.mse.na.fun(true = mean.compo$`GeoPT20 .csv`,list.df = prop20.na.geopt20.list)
output.geopt22 <- geombias.mse.na.fun(true = mean.compo$`GeoPT22 .csv`,list.df = prop20.na.geopt22.list)
output.geopt25 <- geombias.mse.na.fun(true = mean.compo$`GeoPT25 .csv` ,list.df = prop20.na.geopt25.list)
output.geopt32 <- geombias.mse.na.fun(true = mean.compo$`GeoPT32 .csv` ,list.df = prop20.na.geopt32.list)
output.geopt35 <- geombias.mse.na.fun(true = mean.compo$`GeoPT35 .csv` ,list.df = prop20.na.geopt35.list)
output.geopt36 <- geombias.mse.na.fun(true = mean.compo$`GeoPT36 .csv` ,list.df = prop20.na.geopt36.list)
output.geopt41 <- geombias.mse.na.fun(true = mean.compo$`GeoPT41 .csv` ,list.df = prop20.na.geopt41.list)
output.geopt46 <- geombias.mse.na.fun(true = mean.compo$`GeoPT46 .csv` ,list.df = prop20.na.geopt46.list)



geommeandf <- function(dataset){
  dataset <- as_tibble(dataset)
  df <- dataset[-ncol(dataset)]%>% log %>%
  rowMeans() %>% exp() %>% as_tibble()
  df[ncol(df)+1] <- dataset[ncol(dataset)]
  return(df)
}


list.dataframe <- list(output.geopt2[[1]],
output.geopt20[[1]],
output.geopt22[[1]],
output.geopt25[[1]],
output.geopt32[[1]],
output.geopt35[[1]],
output.geopt36[[1]],
output.geopt41[[1]],
output.geopt46[[1]])



result <- lapply(list.dataframe,geommeandf)

names(result) <- c("geopt2","geopt20","geopt22","geopt25",
                   "geopt32","geopt35","geopt36","geopt41","geopt46")
# combine the tibbles using bind_rows(), and add a column to indicate the original tibble
combined_tib <- bind_rows(result, .id = "dataset")


combined_tib %>% ggplot(aes(x=dataset,y=value,color=type))+geom_boxplot()+
  theme_bw()+theme(legend.position = "bottom")+labs(y="Empirical Geometric Mean Bias")








##### Total Bias Plot
total.bias.geopt2 <- output.geopt2[[2]]
total.bias.geopt20 <- output.geopt20[[2]]
total.bias.geopt22 <- output.geopt22[[2]]
total.bias.geopt25 <- output.geopt25[[2]]
total.bias.geopt32 <- output.geopt32[[2]]
total.bias.geopt35 <- output.geopt35[[2]]
total.bias.geopt36 <- output.geopt36[[2]]
total.bias.geopt41 <- output.geopt41[[2]]
total.bias.geopt46 <- output.geopt46[[2]]

list.total.bias <- list(total.bias.geopt2,
                        total.bias.geopt20,
                        total.bias.geopt22,
                        total.bias.geopt25,
                        total.bias.geopt32,
                        total.bias.geopt35,
                        total.bias.geopt36,
                        total.bias.geopt41,
                        total.bias.geopt46)

names(list.total.bias) <- c("geopt2","geopt20","geopt22","geopt25",
                                              "geopt32","geopt35","geopt36","geopt41","geopt46")


name_totalbias <- c("geopt2","geopt20","geopt22","geopt25",
                    "geopt32","geopt35","geopt36","geopt41","geopt46")

for (i in 1:length(list.total.bias)){
  list.total.bias[[i]] <- list.total.bias[[i]] %>% mutate(name = name_totalbias[i])
}

geombias.df <- do.call(rbind, list.total.bias)
unique(geombias.df$name)


geombias.df %>% 
  ggplot(aes(x = name, y = geomean, color = type)) + 
  geom_point() + scale_y_log10()+scale_y_continuous(breaks = c(1,1.5,1.75,2,2.5,3))+
  theme_bw()+theme(legend.position = "bottom")+labs(y="Geometric Mean Bias")



## MSE Bias

total.mse.geopt2  <-  output.geopt2[[6]]
total.mse.geopt20 <- output.geopt20[[6]]
total.mse.geopt22 <- output.geopt22[[6]]
total.mse.geopt25 <- output.geopt25[[6]]
total.mse.geopt32 <- output.geopt32[[6]]
total.mse.geopt35 <- output.geopt35[[6]]
total.mse.geopt36 <- output.geopt36[[6]]
total.mse.geopt41 <- output.geopt41[[6]]
total.mse.geopt46 <- output.geopt46[[6]]

list.total.mse <- list(total.mse.geopt2 ,
                       total.mse.geopt20,
                       total.mse.geopt22,
                       total.mse.geopt25,
                       total.mse.geopt32,
                       total.mse.geopt35,
                       total.mse.geopt36,
                       total.mse.geopt41,
                       total.mse.geopt46)

name_totalmse <- c("geopt2","geopt20","geopt22","geopt25",
                    "geopt32","geopt35","geopt36","geopt41","geopt46")


for (i in 1:length(list.total.mse)){
  list.total.mse[[i]] <- list.total.mse[[i]] %>% mutate(name = name_totalmse[i])
}

geommse.df <- do.call(rbind, list.total.mse)

geommse.df %>% 
  filter(!type %in% c("l2 blr mean")) %>%
  ggplot(aes(x = name, y = geomean, color = type)) + 
  geom_point() + scale_y_log10(breaks=c(1,1.1,1.2,1.3,1.4))+labs(y="Geometric MSE")+
  theme_bw()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))


bias.list <- list(total.bias.geopt2,
                  total.bias.geopt20,
                  total.bias.geopt22,
                  total.bias.geopt25,
                  total.bias.geopt32,
                  total.bias.geopt35,
                  total.bias.geopt36,
                  total.bias.geopt41,
                  total.bias.geopt46) 

names(bias.list) <- c("geopt2","geopt20","geopt22","geopt25",
                      "geopt32","geopt35","geopt36","geopt41","geopt46")



total.bias <- bind_rows(bias.list, .id = "dataset")

saveRDS(total.bias,file = "missingvalues/prop20/totalbias.RData")

tbias20 <- readRDS(file = "missingvalues/prop20/totalbias.RData")
tbias50 <- readRDS(file = "missingvalues/prop50/totalbias.RData")
tbias20


ilr(x=c(0.5,0.4,0.2)) %>% ilrInv()
