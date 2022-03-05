# Outlier detection based on reconstruction error of the data matrix.
library(MASS)
library(readr)
library(tidyverse)
library(compositions)

data <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/GeoPT46 .csv")
RMSEfunction <- function(data){
  # Takes as input a processed dataframe without missing values. Transform it using the ilr transform from the simplex 
  # to the Euclidean D-1 real space.
  X <- data %>% select(-Rest) %>% clo() %>% ilr() %>% as_tibble()
  # Store the means of each columns.
  mu = colMeans(X)
  # Perform a PCA
  Xpca <- prcomp(X)
  
  # Reconstruction error, mean squared error 
  mse<-function(x_hat,x) rowMeans((x_hat-x)^2)
  
  # How reconstruction error decreases with the number of components retained ? 
  Components <- seq(1:length(Xpca$sdev))
  
  Xhat.list <- list()
  root.mse.list <- list()
  root.mse <- c()
  for (i in 1:length(Components)){
    Xhat.list[[i]] = Xpca$x[,1:Components[i]] %*% t(Xpca$rotation[,1:Components[i]])
    Xhat.list[[i]] = scale(Xhat.list[[i]], center = -mu, scale = FALSE)
    root.mse.list[[i]] <- sqrt(mse(Xhat.list[[i]],X))
    root.mse[i] <- mse(Xhat.list[[i]],X) %>% sum() %>% sqrt()
  }
  # The function returns a list with the MSE for each rows (mse.list) and the root mse summed for all rows,
  # This RMSE should decreases with the number of Principal Component Increasing.
  rootmsedf <- cbind(root.mse,Components) %>% as_tibble
  output <- list(root.mse.list,root.mse,rootmsedf)
  names(output) <- c("Root MSE for each rows","Summed RMSE","Dataframe with rootmse and number of PC retained")
  return(output)
}
names(list.df)
# GeoPT 1
RMSEfunction(list.df[[1]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()
# GeoPT 9
RMSEfunction(list.df[[9]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()
# GeoPT 29
RMSEfunction(list.df[[11]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()
# GeoPT 46
RMSEfunction(list.df[[25]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()

list.rmseplot <- list()

for (i in 1:length(list.df) ){
  list.rmseplot[[i]] <- RMSEfunction(list.df[[i]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()
}
names(list.rmseplot) <- names(list.df)
strsubst(t,' ','_')
str_replace(t," ","_")
filenamesplot <- paste(substr(names(list.rmseplot),1,stop = 7),".jpeg")
for (i in 1:length(list.rmseplot)){
  ggsave(filename = paste("rmseplot/",filenamesplot[i]),list.rmseplot[[i]])  
}
# RMSE is much higher for projection of 1-2 PC for all rocks older than GeoPT22. Valuable information.

# Flagging outliers in the sub PC space which explains more than 75 % of the variability of the data
# GeoPT 1 

RMSEfunction(list.df[[1]])

X <- list.df[[10]] %>% select(-Rest) %>% clo() %>% ilr() %>% as_tibble()
# Store the means of each columns.
mu = colMeans(X)
# Perform a PCA
Xpca <- prcomp(X)
index <- min(which(cumsum(Xpca$sdev^2)/sum(Xpca$sdev^2) >= share.variance))

# Reconstruction error, mean squared error 
mse<-function(x_hat,x) rowMeans((x_hat-x)^2)

# How reconstruction error decreases with the number of components retained ? 
Components <- index
Xhat = Xpca$x[,1:Components] %*% t(Xpca$rotation[,1:Components])
Xhat = scale(Xhat, center = -mu, scale = FALSE)
mse.vector <- mse(Xhat,X) %>% as_tibble()

ggplot(mse.vector,aes(x=value))+geom_histogram(bins = 100)+theme_bw()
