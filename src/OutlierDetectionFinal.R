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
names(list.df) <- gsub(" ","",names(list.df))


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
  ggsave(filename = paste("rmseplotSCALED",filenamesplot[i],sep = "/"),list.scaledRMSEplots[[i]])  
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
  for (i in 1:500){
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
plot.based.optimal <- c()
for (i in 1:length(rmse.list)){
plot.based.optimal[i] <- max(which(rmse.list[[i]]>1))+1
}
plot.based.optimal
sim.based.optimal <- ncomponents %>% as.numeric()
# export csv
write.csv(tibble(names(list.df),plot.based.optimal,sim.based.optimal),file = "sim_or_rmseplot.csv")

#### Euclidean Distance between original data matrix  ####

eucliddist.l <- list()
length(eucliddist.l) <- length(list.df)
names(eucliddist.l) <- names(list.df)

histograms.l <- list()
length(histograms.l) <- length(list.df)
names(histograms.l) <- names(list.df)

for (i in 1:length(list.df)){
eucliddist.l[[i]] <-  RMSEfun(list.df[[i]],k=ncomponents[[i]])[[2]]
histograms.l[[i]] <-   eucliddist.l[[i]] %>% as_tibble() %>% 
  ggplot(aes(x=value))+geom_histogram(bins = 50)+theme_bw()

}



for (i in 1:length(histograms.l)){
ggsave(filename = paste("histplotSCALED",filenamesplot[i],sep = "/"),histograms.l[[i]] )  
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
  X <- mvrnorm(50000,mu = Xhat,Sigma = Sigmahat)
  
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
# Test the function, 29th dataset, 4 outliers detected

df.29 <- find.cutoff(list.df[[29]] , ncomponents[[29]], quantile = .975)[[2]]
df.29 %>% as_tibble %>% ggplot(aes(x=value))+geom_histogram(bins=60)+theme_bw()
eucliddist.l[[29]] <-  RMSEfun(list.df[[29]],k=ncomponents[[29]])[[2]]
df <- as.tibble(eucliddist.l[[29]])
df$type <- "Real dataset"
df.sim <- df.29 %>% as_tibble
df.sim$type <- "Simulated dataset"
df.tot <- rbind(df.sim,df)
ggplot(aes(x=value),data=df.tot)+geom_density(aes(fill=type),alpha=.5)+geom_vline(xintercept = quantiles[29],show.legend = T)+
  scale_x_continuous(breaks=c(0,0.5,1,1.25,2,2.5,round(quantiles[29],digits = 3)))+theme_bw()


# DF[[10]] (GeoPT19), no outliers detected
df.10 <- find.cutoff(list.df[[10]] , ncomponents[[10]], quantile = .975)[[2]]
df.10 %>% as_tibble %>% ggplot(aes(x=value))+geom_histogram(bins=60)+theme_bw()
eucliddist.l[[10]] <-  RMSEfun(list.df[[10]],k=ncomponents[[10]])[[2]]
df <- as.tibble(eucliddist.l[[10]])
df$type <- "Real dataset"
df.sim <- df.10 %>% as_tibble
df.sim$type <- "Simulated dataset"
df.tot <- rbind(df.sim,df)
ggplot(aes(x=value),data=df.tot)+geom_density(aes(fill=type),alpha=.5)+geom_vline(xintercept = quantiles[10],show.legend = T)+
  scale_x_continuous(breaks=c(0,0.5,1,2,2.5,3,round(quantiles[10],digits = 3),4.5,5,6))+theme_bw()


df.39 <- list.df[["GeoPT39.csv"]]
df39.mean <- df.39 %>% ilr() %>% colMeans() %>% ilrInv() %>%
  unclass() %>% as_tibble() %>% t()
df39.mean
colnames(df39.mean) <- colnames(df.39)

df39.mean

quantiles <- c()
for (i in 1:length(list.df)){
  quantiles[i] <- find.cutoff(list.df[[i]],ncomponents[[i]])[[1]]
}

list.outliers <- list()

histograms.l[["GeoPT1.csv"]]+geom_vline(xintercept = quantiles[[1]],color="red",show_guide=T)+
  scale_x_continuous(breaks=c(0,0.5,1,1.5,2,2.5,1,round(quantiles[1],digits = 3)))

histograms.l[["GeoPT16 .csv"]]+geom_vline(xintercept = quantiles[[9]],color="red",show_guide=T)+
  +   scale_x_continuous(breaks=c(0,0.5,1,1.5,2,2.5,1,round(quantiles[9],digits = 3)))

histograms.l[["GeoPT36 .csv"]]+geom_vline(xintercept = quantiles[[20]],color="red",show_guide=T)+
  scale_x_continuous(breaks=c(0,0.5,1,1.5,2,2.5,1,round(quantiles[20],digits = 3)))

histograms.l[["GeoPT48 .csv"]]+geom_vline(xintercept = quantiles[[29]],color="red",show_guide=T)+
  scale_x_continuous(breaks=c(0,0.5,1,1.25,2,2.5,1,round(quantiles[29],digits = 3)))



names(histograms.l)

list.outliers <- list()
for (i in 1:29){
  list.outliers[[i]] <- which(eucliddist.l[[i]] > quantiles[[i]])
}
list.outliers

summarytable <- tibble(names(list.df),quantiles,as.character(list.outliers))
write.csv(summarytable,file = "summarytableoutdet.csv")

# DONE DONE DONE :D 

# Export dataframes without outliers
list.outfree <- list()
for (i in 1:length(list.df)){
list.outfree[[i]] <- list.df[[i]][-list.outliers[[i]],]
}

for (i in 1:length(list.outfree) ){
  write.csv(list.outfree[[i]],file = paste0("~/Documents/MA2/Thesis/Repository/data/outfree/",names(list.df)[[i]]),row.names = FALSE)
}

