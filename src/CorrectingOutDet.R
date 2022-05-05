#### Final Outlier Detection BASED ON SCALED PCA ####

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
RMSEfun <- function(data,k=1){
  X <- data %>% dplyr::select(-Rest) %>% clo() %>% ilr() %>% as_tibble()
  
  Xpca = prcomp(X,scale. = T)
  
  mu = colMeans(X)
  
  col.sd <- apply(X,2,sd)
  
  nComp <- k
  
  Xhat = Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
  
  Xhat <- Xhat %*% diag(col.sd)
  
  Xhat = scale(Xhat, center = -mu, scale = FALSE)
  # Euclidean distance between each observations X and Xhat : 
  ED <- sqrt(rowSums((X-Xhat)^2))
  RMSE <- sqrt(mean(ED^2))
  output <- list(RMSE,ED)
  names(output) <- c("RMSE","L2 distance between X and Xhat")
  return(output)
}
rmse.list <- list()
length(rmse.list) <- length(list.df)
for (j in 1:length(list.df) ){
  t[j] <- length(prcomp(list.df[[j]],scale. = T,center = T)$sdev)
  for (i in 1:(t[j]-4) ){
    rmse.list[[j]][i] <- RMSEfun(list.df[[j]],i)[[1]]
  }
}

list.scaledRMSEplots <- list()
length(list.scaledRMSEplots) <- length(list.df)

for (i in 1:length(list.df)){
df <- tibble(rmse.list[[i]],1:length(rmse.list[[i]]))
names(df) <- c("RMSE","Components")
list.scaledRMSEplots[[i]] <- df %>% ggplot(aes(x=Components,y=RMSE))+geom_point()+theme_bw()
}

# Filenames repair, remove csv extension
filenamesplot <- paste(str_extract(names(list.df), '.*(?=\\.csv)'),".jpeg")
# Remove empty space
filenamesplot <- gsub(" ", "", filenamesRMSEplot, fixed = TRUE)

# Export RMSEplots 
for (i in 1:length(list.scaledRMSEplots)){
  ggsave(filename = paste("rmseplotSCALED/",filenamesplot[i],sep = ""),list.scaledRMSEplots[[i]])  
}

#### RMSE PLOT DONE GOOD JOB. ####

#### Second method, simulate from MV Normal ####
ncomponent <- function(dataset,cutoff=.95){
  # Computes the ilr coordinates of each datapoints in the dataset
  dataset <- dataset %>% dplyr::select(-Rest) %>% clo() %>% ilr() %>% as_tibble()
  # Estimate the Mean & VarCov Matrix of the dataset
  xbar <- colMeans(dataset)
  sigmahat <- var(dataset)
  
  # Create a empty simulation dataset and scaled PCA
  sim.dataset <- list()
  pca.simdataset <- list()
  ncompo <- c()
  for (i in 1:100){
    sim.dataset[[i]] <- mvrnorm(nrow(dataset),mu = xbar,Sigma = sigmahat) %>% as_tibble()
    pca.simdataset[[i]] <- prcomp(sim.dataset[[i]],scale. = T)
    ncompo[i] <- which(cumsum(pca.simdataset[[i]]$sdev^2)/sum(pca.simdataset[[i]]$sdev^2)>=cutoff)  %>% min()
  }
  # Return the number of components which should be kept
  ncompo <- ceiling(mean(ncompo))
  return(ncompo)
}

# Number of meaningful components under the null

ncomponents <- lapply(list.df,ncomponent)

#### Euclidean Distance between original data matrix  ####

eucliddist.l <- list()
length(eucliddist.l) <- length(list.df)
names(eucliddist.l) <- names(list.df)

length(histograms.l) <- length(list.df)
names(histograms.l) <- names(list.df)

for (i in 1:length(list.df)){
eucliddist.l[[i]] <-  RMSEfun(list.df[[i]],k=ncomponents[[i]])[[2]]
histograms.l[[i]] <-   eucliddist.l[[i]] %>% as_tibble() %>% 
  ggplot(aes(x=value))+geom_histogram(bins = 50)+theme_bw()

}



for (i in 1:length(histograms.l)){
ggsave(filename = paste("histplotSCALED/",filenamesplot[i]),histograms.l[[i]] )  
}

#HISTPLOTS DONE

pca.simdataset
# Number of principal components under normality assumption

# What should be the distribution of the euclidean distances of observation ?
# Generate a dataframe

find.cutoff <- function(dataframe,nComp,quantile=.975){
  # INPUT : DF
  # OUTPUT : 97.5 % estimated quantile of the distribution of the euclid distances between 
  # a MV-normal sample and its k-dimensional reconstruction from a sampe.
  dataframe <- dataframe %>% dplyr::select(-Rest) %>% clo() %>% ilr() %>% as_tibble()
  Xhat <- colMeans(dataframe)
  Sigmahat <- var(dataframe)
  X <- mvrnorm(5000,mu = Xhat,Sigma = Sigmahat)
  
  Xpca = prcomp(X,scale. = T)
  mu = colMeans(X)
  
  col.sd <- apply(X,2,sd)
  
  
  Xhat = Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
  
  Xhat <- Xhat %*% diag(col.sd)
  
  Xhat = scale(Xhat, center = -mu, scale = FALSE)
  # Euclidean distance between each observations X and Xhat : 
  ED <- sqrt(rowSums((X-Xhat)^2))
  quant <- quantile(ED,probs = quantile)
  output <- list(quant,ED)
  names(output) <- c("1-alpha sample quantile","theoretical euclidean distances")
  return(output)
  
}
# Test the function
find.cutoff(list.df[[2]] , ncomponents[[2]], quantile = .975)

quantiles <- c()
for (i in 1:length(list.df)){
  quantiles[i] <- find.cutoff(list.df[[i]],ncomponents[[i]])[[1]]
}
list.outliers <- list()

list.outliers <- list()
for (i in 1:29){
  list.outliers[[i]] <- which(eucliddist.l[[i]] > quantiles[[i]])
}
list.outliers


# DONE DONE DONE :D 