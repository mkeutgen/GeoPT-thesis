# The Modeling
library(readr)
library(tidyverse)
library(compositions)
library(moments)
library(nortest)
library(mvtnorm)
library(MVN)
setwd("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/")
# CLR Transform, assessing the normality assumption
# import dataset
dataframe <- read_csv("GeoPT48 .csv")

# Proportion of unique values in the dataframe
prop.unique <- function(x) length((unique(x)))/length(x)
apply(dataframe,2,prop.unique)
# One sees that in major elements, no values were imputed by the blr mean. This shows that the cutoff in the df cleaning
# function were successful in eliminating dubious analyses w.r.t major elements.
# However, some elements like Indium, Germanium or Bore have less than 30 % of values which are not the blr.mean... 





dataframe.clr <- dataframe %>% select(!c(X1,Rest)) %>%
  clo() %>% clr() %>%
  as_tibble()
# Blr (really a logit transformation) 
logit.transf <- function(x) log(x/(1-x))
blr <- function(dataframe) apply(dataframe,2,logit.transf)

# Goodness of fit test, we do Shapiro Wilk and not Kolmogorov-Smirnov. 

shapiro.pvalue <- function(x) shapiro.test(x)$p.value
pvalues <- apply(dataframe, 2, shapiro.pvalue)
# Proportion of columns where normality assumption is not rejected 
sum(ifelse(pvalues>0.05,1,0))/length(pvalues)
# 0.03076923 with CLR
# 0.046875 with ILR
summary(prcomp(dataframe))

# Which transformation of the data is the most likely to lead to an acceptable normality assumption ?

# Import all processed datasets at once :

# Get the files path :
filePaths <- list.files("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/", "\\.csv$", full.names = TRUE)
# Get the files name :
fileNames <- list.files("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/", "\\.csv$", full.names = FALSE)


list.df <- list()

for (i in 1:length(filePaths)){
  list.df[[i]] <- read_csv(file = filePaths[[i]])
}

names(list.df) <- fileNames
# Remove NA.df, why does it appear in the first place ? 
list.df <- list.df[names(list.df)!=c("NA.csv")]

list.df.ilr <- list() 
list.df.clr <- list()  
list.df.alr <- list()  
list.df.blr <- list()

for (i in 1:length(list.df)){
  list.df.ilr[[i]] <- list.df[[i]] %>% ilr() %>% as_tibble()
}

for (i in 1:length(list.df)){
  list.df.clr[[i]] <- list.df[[i]] %>% clr() %>% as_tibble()
}

for (i in 1:length(list.df)){
  list.df.alr[[i]] <- list.df[[i]] %>% alr() %>% as_tibble()
}

for (i in 1:length(list.df)){
  list.df.blr[[i]] <- list.df[[i]] %>% blr() %>% as_tibble()
}


# Different data transformations :
# ALR additive log ratios 
# CLR centred log ratios
# ILR isometric log ratios

# Let's start with a simple dataframe, then replicate with the list of 29 dataframes previously computed. 
dataframe <- read_csv("GeoPT48 .csv")
dataframe <- dataframe %>% select(!c(Rest)) %>%
  clo() %>% as_tibble()

# Dataframes ALR CLR ILR BLR 
dataframe.alr <- dataframe %>% alr(ivar = 1) %>% as_tibble()
dataframe.clr <- dataframe %>% clr() %>% as_tibble()
dataframe.ilr <- dataframe %>% ilr() %>% as_tibble()
colnames(dataframe.ilr) <- colnames(dataframe.clr)
dataframe.blr <- dataframe %>% blr() %>% as_tibble()

## Testing for Multivariate Normality. https://doi.org/10.1111/j.2517-6161.1982.tb01195.x
# Aitchinson : "For d-dimensional compositional data sets we have applied Kolmogorov-Smirnov and
# Cramer-von Mises tests in their Stephens (1974) versions to all d marginal distributions, to all
# -!d(d - I) bivariate angle distributions, and to the distribution of d-dimensional radii. For the
# Skyelava compositions we have tested in this way both the additive aN9 and the multiplicative
# mN910gistic-normal models. For the additive version not a single one of the battery of92 tests
# gives a significant indication of non-normality at the 5 per cent significance level; for the
# multiplicative version only one of the marginal tests gives evidence of any departure from
# normality, at the 1 per cent significance level. Application of the battery of tests to another 20
# data sets of different geological types similarly encourages the view that transformed-normal
# distributions may have an important practical role to play in the analysis of compositional
# data."
# We argue against using Kolmogorov Smirnov test after having done simulations that showed pvalues
# will be artificially too high.
dataframe.ilr[1:10]


normality.test <- mvn(dataframe.ilr, subset = NULL, mvnTest = c("mardia"), covariance = TRUE, tol = 1e-25, alpha = 0.5,
    scale = FALSE, desc = TRUE, transform = "square",
    univariateTest = c("AD"),#"SW", CVM", "Lillie", "SF", "AD"
    univariatePlot = "none", multivariatePlot = "none",
    multivariateOutlierMethod = "none", bc = FALSE, bcType = "rounded",
    showOutliers = TRUE, showNewData = FALSE)


scaled.ilr <- scale(dataframe.ilr,center=T,scale = F) %>% as_tibble()
# Under the assumption that the standardize dataset is sampled from a standard normal distribution, we remove observations 
# whose values are above or under a cutoff which is the 0.995 quantile of the std normal distribution 

cutoff.out <- qnorm(0.995) # 0.5 % observations under (-cutoff) and 0.5 % observations above this cutoff.

# Let a vector
x <- c(rnorm(10),c(3,-3))
x[which(!(x>(-1.96) & x<1.96))]

out <- function(x)which((x>(-cutoff.out) & x<cutoff.out))
scaled.ilr[apply(scaled.ilr,2,out)]


subset(x,x> -1.96)

View(scaled.ilr)
apply(scaled.ilr,2,min)




pvalue.ad.t <- function(x) ad.test(x)$p.value
pvalue.cvm.t <- function(x) cvm.test(x)$p.values
pvalue.shap.t <- function(x) shapiro.test(x)$p.value

alr.pvalue.ad   <- apply(dataframe.alr, 2, pvalue.ad.t)
alr.pvalue.cvm  <- apply(dataframe.alr,2, pvalue.cvm.t)
alr.pvalue.shap <- apply(dataframe.alr,2,pvalue.shap.t)

clr.pvalue.ad   <- apply(dataframe.clr, 2, pvalue.ad.t)
clr.pvalue.cvm  <- apply(dataframe.clr,2, pvalue.cvm.t)
clr.pvalue.shap <- apply(dataframe.clr,2,pvalue.shap.t)

ilr.pvalue.ad   <- apply(dataframe.ilr, 2, pvalue.ad.t)
ilr.pvalue.cvm  <- apply(dataframe.ilr,2, pvalue.cvm.t)
ilr.pvalue.shap <- apply(dataframe.ilr,2,pvalue.shap.t)

ilr.pvalues.df <- cbind(ilr.pvalue.ad,  
                            ilr.pvalue.cvm, 
                            ilr.pvalue.shap)

View(ilr.pvalues.df)


 # Proportion of columns where normality assumption is not rejected 
sum(ifelse(pvalue.ad>0.05,1,0))/length(pvalues)
sum(ifelse(pvalue.ad>0.05,1,0))/length(pvalues)

sum(ifelse(pvalues>0.05,1,0))/length(pvalues)


# Hypothesis Tests for Normality : we need to test if the sample Y_1,..,Y_n which are the transformed composition comes from a normal distribution with unknown parameter
# mu^2 sigma^2
nrow(dataframe.alr)
View(dataframe.alr)


par(mfrow=c(2,2))
hist(dataframe.alr$SiO2,breaks=50)
hist(dataframe.clr$SiO2,breaks=50)
hist(dataframe.ilr$SiO2,breaks=50)
hist(dataframe.blr$SiO2,breaks=50)

par(mfrow=c(2,2))
hist(dataframe.alr$TiO2,breaks=50)
hist(dataframe.clr$TiO2,breaks=50)
hist(dataframe.ilr$TiO2,breaks=50)
hist(dataframe.blr$TiO2,breaks=50)


par(mfrow=c(2,2))
hist(dataframe.alr$Al2O3,breaks=50)
hist(dataframe.clr$Al2O3,breaks=50)
hist(dataframe.ilr$Al2O3,breaks=50)
hist(dataframe.blr$Al2O3,breaks=50)


par(mfrow=c(2,2))
hist(dataframe.alr$Cr,breaks=50)
hist(dataframe.clr$Cr,breaks=50)
hist(dataframe.ilr$Cr,breaks=50)
hist(dataframe.blr$Cr,breaks=50)

# How reasonable is the Normality Assumption ? 
shapiro.pvalue <- function(x) shapiro.test(x)$p.value

pvalues.alr <- apply(dataframe.alr, 2, shapiro.pvalue)
pvalues.clr <- apply(dataframe.clr, 2, shapiro.pvalue)
pvalues.ilr <- apply(dataframe.ilr, 2, shapiro.pvalue)
pvalues.blr <- apply(dataframe.blr, 2, shapiro.pvalue)

# Number of times normality assumption is rejected : 
propH0 <- function(x) sum(ifelse(x>0.01,1,0))/length(x)
propH0.alr <- propH0(pvalues.alr)
propH0.clr <- propH0(pvalues.clr)
propH0.ilr <- propH0(pvalues.ilr)
propH0.blr <- propH0(pvalues.blr)

# Normality assumption unreasonable. Transforming the data ? Reason for this is high kurtosis. Lepdokurtosis (fat right tail)
apply(dataframe.clr,2,skewness)

# Removing outliers ? https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/ 
# First we need to decide which columns should get their outliers removed.
# Outlier detection should be pursued ONLY if the number of cells replaced by their blr.mean is not too large. 
# What is too large ? We propose a cutoff of 25 %.

cut.off <- 0.75
result <- apply(dataframe,2,prop.unique)>=cut.off
df <- dataframe[result]
df.clr <- df %>% clr() %>% as_tibble()
df.clr[10]
element<- element[-which(df.clr[1] %in% out[[1]]),]
out <- list()
for (i in 1:ncol(df)){
  out[[i]] <- boxplot(df.clr[i], plot=FALSE)$out
}

# Weird bug going out here
df.clr[[15]][-which(df.clr[[15]] %in% out[[15]]  )]
# To be continued here....



bx <- boxplot(scale(dataframe.clr$SiO2,center = T,scale=T))
hist(dataframe.clr$SiO2)

outliers <- boxplot(dataframe.clr$SiO2, plot=FALSE)$out

x<-dataframe.clr

x<- x[-which(dataframe.clr$SiO2 %in% outliers),]
outliers
shapiro.test(x$SiO2)
kurtosis(dataframe.ilr)
hist(dataframe.clr$Th)
hist(dataframe.ilr$Th)


outliers.list <- list()
for (i in 1:ncol(dataframe.clr)){
  outliers.list[[i]] <- boxplot(dataframe.clr[[i]], plot=FALSE)$out
  
}

hist(x$SiO2)
