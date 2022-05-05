# Outlier detection based on reconstruction error of the data matrix.
library(MASS)
library(readr)
library(tidyverse)
library(compositions)
setwd("~/Documents/MStatistics/MA2/Thesis/Repository/src/")
data <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/GeoPT46 .csv")
View(data)
names(list.df)
ordered.names <- c("GeoPT1.csv","GeoPT2.csv","GeoPT3.csv","GeoPT4.csv","GeoPT5.csv","GeoPT6.csv",
                   "GeoPT8.csv","GeoPT14 .csv","GeoPT16 .csv","GeoPT19 .csv","GeoPT20 .csv",
                   "GeoPT21 .csv","GeoPT22 .csv","GeoPT23 .csv","GeoPT25 .csv", "GeoPT29 .csv",
                   "GeoPT32 .csv","GeoPT34 .csv","GeoPT35 .csv","GeoPT36 .csv","GeoPT37 .csv","GeoPT38 .csv",
                   "GeoPT38A.csv","GeoPT39 .csv","GeoPT39A.csv","GeoPT41 .csv","GeoPT43 .csv",
                   "GeoPT46 .csv","GeoPT48 .csv")

ncomponents[ordered.names]
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
# GeoPT 1
RMSEfunction(list.df[[1]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()
# GeoPT 9
RMSEfunction(list.df[[9]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()
# GeoPT 29
RMSEfunction(list.df[[11]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()
# GeoPT 46
RMSEfunction(list.df[[25]])[[3]] %>% ggplot(aes(x=Components,y=root.mse))+geom_point()+theme_bw()

list.rmseplot <- list()

# Store the RMSE function outputs.

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

for (i in 1:length(list.rmseplot)){
  ggsave(filename = paste("rmseplot/",filenamesplot[i]),list.rmseplot[[i]])  
}
# RMSE is much higher for projection of 1-2 PC for all rocks older than GeoPT22. Valuable information.

# Flagging outliers in the sub PC space which explains more than 95 % of the variability of the data
# Number of PC under assumption of multivariate normality.
dataset <- list.df[[10]]  

ncomponent <- function(dataset,cutoff=.95){
  # Computes the ilr coordinates of each datapoints in the dataset
  dataset  <- dataset %>% ilr() %>% as_tibble()
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
  length(euclideanlist) <- length(Components)
  
    Xhat = Xpca$x[,1:Components] %*% t(Xpca$rotation[,1:Components])
    Xhat = scale(Xhat, center = -mu, scale = FALSE)
    for (j in 1:nrow(X)){
      euclideandist[j] <- sqrt(sum((Xhat[j,]- X[j,])^2))
    }

  return(euclideandist)
}

euclideandist.list <- list()
names(euclideandist.list) <- names(list.df)
for (i in 1:length(list.df)){
euclideandist.list[[i]] <- euclideandist.kapprox(list.df[[i]],ncomponents[[i]])
}
which(euclideandist.list[["GeoPT19 .csv"]]>3)
histograms.l <- list()
for (i in 1:length(list.df)){
histograms.l[[i]] <- euclideandist.list[[i]] %>% 
  as_tibble() %>% ggplot(aes(x=value))+geom_histogram(bins = 100)+theme_bw()
}

filenamesplot <- paste(substr(names(list.df),1,stop = 7),".jpeg")
for (i in 1:length(list.rmseplot)){
  ggsave(filename = paste("histplot/",filenamesplot[i]),histograms.l[[i]])  
}

histograms.l[[5]] %>% as_tibble() %>% ggplot(aes(x=value))+geom_histogram(bins = 100)+theme_bw()
names(histograms.l) <- names(list.df)
names(ncomponents) <- names(list.df)

pca.simdataset <- prcomp(sim.dataset[[1]])
"%>% as_tibble  %>% ggplot() + geom_histogram(bins = 100)"
RMSEfunction(list.df[[i]])[[1]][[ncomponents[[i]]]]


pca.simdataset
# Number of principal components under normality assumption
which(cumsum(pca.simdataset$sdev^2)/sum(pca.simdataset$sdev^2)>=.95)  %>% min()

summary(prcomp(dataset))


xbar <- rnorm(66)
summary(prcomp(mvrnorm(60,mu = xbar,Sigma = varcov)))

# What should be the distribution of the euclidean distances of observation ?
# Generate a dataframe
# dataframe <- list.df[[1]]
 Xhat <- colMeans(dataframe)
Sigmahat <- var(dataframe)
X <- mvrnorm(1000,mu = Xhat,Sigma = Sigmahat)
norm.euclideandist.kapprox <- function(X,k=2){
#   # Takes as input a processed dataframe without missing values. Transform it using the ilr transform from the simplex 
#   # to the Euclidean D-1 real space.
#   # Store the means of each columns.
X <- X %>% ilr()
mu = colMeans(X)
#   # Perform a PCA
Xpca <- prcomp(X,center = TRUE)

euclideandist <- c()
length(euclideanlist) <- length(Components)
Xhat <- Xpca$x[,1:k] %*% t(Xpca$rotation[,1:k])
Xhat <- scale(Xhat, center = -mu, scale = FALSE)
for (j in 1:nrow(X)){
euclideandist[j] <- sqrt(sum((Xhat[j,]- X[j,])^2))
}
   
return(euclideandist)
}

norm.euclideandist.kapprox(X,k = 3) %>% as_tibble() %>% ggplot(aes(x=value))+geom_histogram(bins = 1000)+theme_bw()

out.obs <- which(out>=0.3)



# Testing Multivariate Normality Assumption on the reduced dataset :
  # Without outliers 
  # On only the n meaningful number of PCs
# For dataframe 
ID <- "GeoPT48 .csv"
list.df[[ID]]
histograms.l[[ID]]
out <- RMSEfunction(list.df[[ID]])[[1]][[ncomponents[[ID]]]]

# Cutoff above the 0.3 mark. 
out.obs <- which(out>=0.3)

RMSEfunction(list.df[[ID]])
# Output without outliers
list.df[[ID]][-out.obs,]
# Outliers are : 9 30 33 49 55
# PCA without these 5 outliers who contribute the most to the reconstruction error
# Perform the PCA on the contaminated dataset such that the PCs will be influenced by these outliers and then test for normality

X <- list.df[[ID]] %>% ilr() %>% as_tibble()
mu = colMeans(X)
PCA <- prcomp(X,scale. = T)

Xhat = PCA$x[,1:ncomponents[[ID]]] %*% t(PCA$rotation[,1:ncomponents[[ID]]])

Xhat = scale(Xhat, center = -mu, scale = FALSE)

X %>% head()

Xhat <- Xhat %>% as_tibble() 

Xhat %>% ggplot() + geom_histogram(aes(x=V30))

apply(Xhat,2,shapiro.pvalue)

apply(Xhat[-out.obs,],2,shapiro.pvalue)
MVN::mvn(Xhat[-out.obs,])

summary(prcomp(Xhat[-out.obs,]))

df <- PCA$x[,1:ncomponents[[ID]]] %>% as_tibble()
df <- df[-out.obs,]
apply(df,2,shapiro.pvalue)

# Perform the PCA on the uncontaminated dataset
PCA <- prcomp(list.df[[ID]][-out.obs,],scale. = T)




df <- PCA$x[,1:ncomponents[[ID]]] %>% as_tibble()
df %>% ggplot() + geom_histogram(aes(x=PC9))

apply(df,2,shapiro.pvalue)

# Outlier detection seems pretty inefficient with the PCA-method at adressing normality assumption. May be because we are in this case looking at orthogonal outliers which influence
# on normality of the marginals is limited when compared to score outliers


# Number of components for each dataframe and proportion of components to column 
names.df <- ncomponents[ordered.names] %>% names()
ncomp <- ncomponents[ordered.names] %>% as.numeric()
ncol <- lapply(list.df[ordered.names],ncol) %>% as.numeric()
df <- tibble(names.df,ncomp,ncol) %>% mutate(prop = ncomp/ncol)
write.csv(df,file="summaryPCA.csv")
