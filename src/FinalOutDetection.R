#### Final Outlier Detection ####

# Libraries & importing all the datasets
library(MASS)
library(readr)
library(tidyverse)
library(compositions)

# Import all data 
source("importdata.R")

ordered.names <- c("GeoPT1.csv","GeoPT2.csv","GeoPT3.csv","GeoPT4.csv","GeoPT5.csv","GeoPT6.csv",
                   "GeoPT8.csv","GeoPT14 .csv","GeoPT16 .csv","GeoPT19 .csv","GeoPT20 .csv",
                   "GeoPT21 .csv","GeoPT22 .csv","GeoPT23 .csv","GeoPT25 .csv", "GeoPT29 .csv",
                   "GeoPT32 .csv","GeoPT34 .csv","GeoPT35 .csv","GeoPT36 .csv","GeoPT37 .csv","GeoPT38 .csv",
                   "GeoPT38A.csv","GeoPT39 .csv","GeoPT39A.csv","GeoPT41 .csv","GeoPT43.csv",
                   "GeoPT46 .csv","GeoPT48 .csv")
list.df <- list.df[ordered.names]


#### RMSE Function and RMSE plot to determine optimal k, subspace of the PCA ####
RMSEfunction <- function(data){
  # Takes as input a processed dataframe without missing values. Transform it using the ilr transform from the simplex 
  # to the Euclidean D-1 real space.
  X <- data %>% dplyr::select(-Rest) %>% clo() %>% ilr() %>% as_tibble()
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
  
  euclideanlist <- list()
  length(euclideanlist) <- length(Components)
  
  for (i in 1:length(Components)){
    Xhat.list[[i]] = Xpca$x[,1:Components[i]] %*% t(Xpca$rotation[,1:Components[i]])
    Xhat.list[[i]] = scale(Xhat.list[[i]], center = -mu, scale = FALSE)
    for (j in 1:nrow(X)){
      euclideanlist[[i]][j] <- sqrt(sum((Xhat.list[[i]][j,]- X[j,])^2))
    }
    root.mse[i] <- sqrt(mse(Xhat.list[[i]],X))
  }
  # The function returns a list with the MSE for each rows (mse.list) and the root mse summed for all rows,
  # This RMSE should decrease with the number of Principal Component Increasing.
  rootmsedf <- cbind(root.mse,Components) %>% as_tibble
  output <- list(euclideanlist,root.mse,rootmsedf)
  names(output) <- c("Euclidean distance for each rows","Summed RMSE","Dataframe with rootmse and number of PC retained")
  return(output)
}

# Store the RMSE function outputs and RMSE plots.
list.rmseplot <- list()
RMSE.outputs <- list()

for (i in 1:length(list.df)){
  RMSE.outputs[[i]] <- RMSEfunction(list.df[[i]])
}

# How does RMSE decrease with number of components retained
for (i in 1:length(list.df) ){
  list.rmseplot[[i]] <- RMSE.outputs[[i]][[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()
}
names(list.rmseplot) <- names(list.df)
filenamesplot <- paste(substr(names(list.rmseplot),1,stop = 7),".jpeg")
# Export RMSEplots 
for (i in 1:length(list.rmseplot)){
  ggsave(filename = paste("rmseplot/",filenamesplot[i]),list.rmseplot[[i]])  
}

#### Second method, simulate from MV Normal ####
ncomponent <- function(dataset,cutoff=.95){
  # Computes the ilr coordinates of each datapoints in the dataset
  dataset <- dataset %>% dplyr::select(-Rest) %>% clo() %>% ilr() %>% as_tibble()
  # Estimate the Mean & VarCov Matrix of the dataset
  xbar <- colMeans(dataset)
  sigmahat <- var(dataset)
  
  # Create a empty simulation dataset
  sim.dataset <- list()
  pca.simdataset <- list()
  ncompo <- c()
  for (i in 1:100){
    sim.dataset[[i]] <- mvrnorm(nrow(dataset),mu = xbar,Sigma = sigmahat) %>% as_tibble()
    pca.simdataset[[i]] <- prcomp(sim.dataset[[i]])
    ncompo[i] <- which(cumsum(pca.simdataset[[i]]$sdev^2)/sum(pca.simdataset[[i]]$sdev^2)>=cutoff)  %>% min()
  }
  # Return the number of components which should be kept
  ncompo <- ceiling(mean(ncompo))
  return(ncompo)
}

# Number of meaningful components under the null

ncomponents <- lapply(list.df,ncomponent)

#### TO BE CONTINUED ####


euclideandist.kapprox <- function(data,k=1){
  # Takes as input a processed dataframe without missing values. Transform it using the ilr transform from the simplex 
  # to the Euclidean D-1 real space.
  X <- data %>% dplyr::select(-Rest) %>% clo() %>% ilr() %>% as_tibble()
  # Store the means of each columns.
  mu = colMeans(X)
  # Perform a PCA
  Xpca <- prcomp(X,center = TRUE)
  
  # Reconstruction error, mean squared error 
  mse<-function(x_hat,x) rowMeans((x_hat-x)^2)
  
  # How reconstruction error decreases with the number of components retained ? 
  Components <- k
  
  Xhat <- list()
  root.mse.list <- list()
  root.mse <- c()
  
  euclideandist <- c()

  Xhat = Xpca$x[,1:Components] %*% t(Xpca$rotation[,1:Components])
  Xhat = scale(Xhat, center = -mu, scale = FALSE)
  for (j in 1:nrow(X)){
    euclideandist[j] <- sqrt(sum((Xhat[j,]- X[j,])^2))
  }
  
  return(euclideandist)
}

euclideandist.list <- list()

for (i in 1:length(list.df)){
  euclideandist.list[[i]] <- euclideandist.kapprox(list.df[[i]],ncomponents[[i]])
}


names(euclideandist.list) <- names(list.df)

which(euclideandist.list[["GeoPT19 .csv"]]>3)

histograms.l <- list()
for (i in 1:length(list.df)){
  histograms.l[[i]] <- euclideandist.list[[i]] %>% 
    as_tibble() %>% ggplot(aes(x=value))+geom_histogram(bins = 100)+theme_bw()
}

# Filenames repair, remove csv extension
filenamesplot <- paste(str_extract(names(list.df), '.*(?=\\.csv)'),".jpeg")
# Remove empty space
filenamesplot <- gsub(" ", "", filenamesplot, fixed = TRUE)

for (i in 1:length(list.rmseplot)){
  ggsave(filename = paste("histplot/",filenamesplot[i]),histograms.l[[i]])  
}

histograms.l[[5]]
names(histograms.l) <- names(list.df)
names(ncomponents) <- names(list.df)

pca.simdataset <- prcomp(sim.dataset[[1]])
"%>% as_tibble  %>% ggplot() + geom_histogram(bins = 100)"


pca.simdataset
# Number of principal components under normality assumption

# What should be the distribution of the euclidean distances of observation ?
# Generate a dataframe

find.cutoff <- function(dataframe,nComp=2,quantile=.975){
# INPUT : DF
# OUTPUT : 97.5 % estimated quantile of the distribution of the euclid distances between 
  # a MV-normal sample and its k-dimensional reconstruction from a sampe.
dataframe <- dataframe %>% dplyr::select(-Rest) %>% clo() %>% ilr() %>% as_tibble()

Xhat <- colMeans(dataframe)
Sigmahat <- var(dataframe)
X <- mvrnorm(4000,mu = Xhat,Sigma = Sigmahat)
# Perform PCA
mu = colMeans(X)
Xpca = prcomp(X)

Xhat = Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
Xhat = scale(Xhat, center = -mu, scale = FALSE)
euclideandist <- c()
for (j in 1:nrow(X)){
  euclideandist[j] <- sqrt(sum((Xhat[j,]- X[j,])^2))
}

return(quantile(euclideandist,probs = quantile))

}
quantiles <- c()
for (i in 1:length(list.df)){
  quantiles[i] <- find.cutoff(list.df[[i]],ncomponents[[i]])
}
list.outliers <- list()

for (i in length(euclideandist.list)){
  list.outliers[[i]] <- which(euclideandist.list[[i]]<quantiles[[i]])
}
list.outliers <- list()
for (i in 1:29){
list.outliers[[i]] <- which(euclideandist.list[[i]] > quantiles[[i]])
}
