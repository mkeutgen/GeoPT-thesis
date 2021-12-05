#########################
## Motivating Example.R #
#########################


# Load libraries

library(tidyverse)
library(readr)
library(ggtern)
library(patchwork)
library(compositions)
library(ggfortify)
library(ggrepel)
library(robCompositions)
set.seed(1)
# Start with most recent dataset

data <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/GeoPT48 -84Ra.csv")
# Geometric mean of this sample of compositions : 
geomean <- function(x){
  # compute the geometric mean of a vector
  exp(mean(log(x),na.rm=TRUE)) 
}


# Principal Component Analysis
sel <-c("SiO2","TiO2","Al2O3","Fe2O3T","MnO","MgO","CaO","Na2O","K2O","P2O5")
df.majors <- select(data,all_of(sel))
df.majors <- data.frame(clo(df.majors))

# replace missing values of each column by its geometric mean
geomean.v <- sapply(rbind(df.majors),geomean)
for (i in 1:ncol(df.majors)){
  df.majors[,i][is.na(df.majors[,i])] <- geomean.v[i] 
}
# Now the dataframe does not have any missing values. 
# We use the centered log-ratio transformation.
# Initialize empties dataframe with the centered ratio and log centered ratio
# dataframes
cr.df <- data.frame()
clr.df <- data.frame()
# cr.df, divide each entries in a column by the geometric mean of this column
cr.df <- sweep(df.majors,MARGIN = 2,FUN="/",STATS = geomean.v)
# clr.df is the natural logarithm of cr.df
clr.df <- log(cr.df)
# Perform Singular Value Decomposition on clr.df : 
pca.clr <- prcomp(clr.df,scale = T,rank. = ncol(clr.df)-1 )
summary(pca.clr)
biplot(pca.clr,choices = c(4,5))
# Observation 53 seems outlying.
biplot(pca.clr)
# Let's check that.
# Extract the score matrix
X <- pca.clr$x
df <- data.frame(X[,c(1,2)])
ggplot(aes(x=PC1,y=PC2),data=df)+geom_point()
df.fmcd <- covMcd(df,alpha = 0.95)

plot(df.fmcd,which = "tolEllipse",classic=TRUE)

covPlot(df,which="dd",m.cov=animals.fmcd)
## compositional PCA, non-robust
p_comp <- pcaCoDa(df.majors, method = "classical")
## compositional PCA, robust
set.seed(234) # to reproduce the result exactly
p_comp_rob <- pcaCoDa(df.majors, method = "robust")
summary(p_comp_rob)
biplot(p_comp)
biplot(p_comp_rob)

p_comp_rob
# Outlier detection with the outCoDa function which uses Fast-MCD :
# References : "Outlier detection for compositional data using robust methods.
# Math. Geosciences""

df.outl <- df.majors[od$outlierIndex,]
pca.res <- pcaCoDa(df.outl)
biplot(pca.res)
df.outl.clr <- clr(df.major.wo.outliers)
pca.res <- prcomp(df.outl.clr,scale = T,center =T)

autoplot(pca.res)
