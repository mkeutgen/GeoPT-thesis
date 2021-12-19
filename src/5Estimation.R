library(tidyverse)
library(readr)
library(ggtern)
library(patchwork)
library(compositions)
library(robCompositions)
library(ggfortify)
library(ggrepel)
library(robustbase)
library(rrcov)

geomean <- function(x){
  # compute the geometric mean of a vector
  exp(mean(log(x),na.rm=TRUE)) 
}

setwd("/home/max/Documents/MStatistics/MA2/Thesis/Repository/")
data <- read_csv("data/raw/GeoPT48 -84Ra.csv")

sel <-c("SiO2","TiO2","Al2O3","Fe2O3T","MnO","MgO","CaO","Na2O",
        "K2O","P2O5")

df.majors.raw <- select(data,all_of(sel))

# data cleaning we remove all columns filled with NA
df.majors <- df.majors.raw[rowSums(is.na(df.majors.raw)) != ncol(df.majors.raw),]

# closure operation
df <- data.frame(clo(df.majors))
dim(df.majors.raw) # 97 observations
dim(df) # 86 observations

# remove columns containing even only one missing value
df.na<- df.majors.raw[rowSums(is.na(df.majors.raw)) == 0,] # 62 observations
geomean.v <- sapply(rbind(df.majors),geomean)
for (i in 1:ncol(df.majors)){
  df[,i][is.na(df.majors[,i])] <- geomean.v[i] 
}
# Now there are 2 dataframes. df.na where missing values were removed. df where missing values were imputed
# Their centered log ratio transform are respectively :
clr.df <- clr(df)
clr.df.na <- clr(df.na)
ilr.df <- ilr(df)
ilr.df.na <- ilr(df.na)
# Flagging outliers : compute mcd

ilrdf.mcd <- covMcd(ilr.df)
ilrdf.mcd.na <- covMcd(ilr.df.na)

# Diagnostic plot 

plot(ilrdf.mcd,which = c("distance"),classic = TRUE) # 12 classical outliers when removing rows with missing values
plot(ilrdf.mcd.na,which = c("distance"),classic = TRUE) # 11 classical outliers when imputing missing values by geometric mean

# Conducting PCA without outliers
# Filzmozer & Gregorich propose to use mcd : https://link.springer.com/article/10.1007/s11004-020-09861-6#Sec4
threshold <- sqrt(qchisq(p = 0.975, df = ncol(ilr.df)))
outliers <-  which(sqrt(ilrdf.mcd$mah) >= threshold) # gives the 
outliers
round(sqrt(ilrdf.mcd$mah),3)
outliers <-  which(sqrt(ilrdf.mcd$mah) >= 50)

# PCA without outliers 

# restricted dataset
pca.df <- prcomp(clr.df[-outliers,],scale. = T)
       
pca.df.na <- prcomp(clr.df.na[-outliers,],scale. = T)

# Biplots 
# wo outliers :
autoplot(pca.df,loadings=T,loadings.label=T)+theme_bw()

# Projecting outliers in the first 2 PC space :
# What is happening : clr.df[outliers] %*% t(pca.df$rotation)
proj.outliers <- predict(pca.df,clr.df[outliers])[,1:2]

# Biplot with outliers and normal points

df.notout <- data.frame(pca.df$x[,1:2]) %>% mutate(outlier="no")

df.out <- data.frame(proj.outliers) %>% mutate(outlier="yes")

df.tot <- data.frame(rbind(df.notout,df.out))

ggplot(df.tot,aes(x=PC1,y=PC2,colour=outlier))+geom_point()
# with outliers
plot <- autoplot(pca.df,loadings=F,loadings.label=F)+theme_bw()
plot + geom_point(data = df.out,aes(x=PC1,y=PC2),colour="purple")+theme_bw()
autoplot(pca.df.na,loadings=T,loadings.label=T)+theme_bw()

# We can confidently remove these 13 outlying points out of 97 obs. 

# Provide mean and covariance estimate : 
clr.mean <- colMeans(clr.df[-outliers,])
clr.cov <- Cov(data.frame(clr.df[-outliers,]))


mean.estimate <- as.vector(clrInv(clr.mean))
cov.estimate <- clrInv(as.matrix(clr.cov$cov))
mean.estimate
cov.estimate # matrix of constant

# Repeat for 32 rocks ? 