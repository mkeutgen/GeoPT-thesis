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

logit.transf <- function(x) log(x/(1-x))
blr <- function(dataframe) apply(dataframe,2,logit.transf)

dataframe <- read_csv("GeoPT48 .csv")

# Proportion of unique values in the data frame
prop.unique <- function(x) length((unique(x)))/length(x)
apply(dataframe,2,prop.unique)
# One sees that in major elements, no values were imputed by the blr mean. This shows that the cutoff in the df cleaning
# function were successful in eliminating dubious analyses w.r.t major elements.
# However, some elements like Indium, Germanium or Bore have less than 30 % of values which are not the blr.mean... 


# Let's start with a simple dataframe.

dataframe.pr <- dataframe %>% select(!c(Rest)) %>%
  clo() %>% as_tibble()

# Dataframes ALR CLR ILR BLR 
dataframe.alr <- dataframe.pr %>% alr() %>% as_tibble()
dataframe.clr <- dataframe.pr %>% clr() %>% as_tibble()
dataframe.ilr <- dataframe.pr %>% ilr() %>% as_tibble()
colnames(dataframe.ilr) <- colnames(dataframe.clr)
dataframe.blr <- dataframe.pr %>% blr() %>% as_tibble()



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



 # Proportion of columns where normality assumption is not rejected 
sum(ifelse(ilr.pvalue.ad>0.05,1,0))/length(ilr.pvalue.ad)
sum(ifelse(clr.pvalue.ad>0.05,1,0))/length(clr.pvalues.ad)

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

# Chrome, 94 % unique values
par(mfrow=c(2,2))
hist(dataframe.alr$Cr,breaks=50)
hist(dataframe.clr$Cr,breaks=50)
hist(dataframe.ilr$Cr,breaks=50)
hist(dataframe.blr$Cr,breaks=50)

# Nd 89 % unique values
par(mfrow=c(2,2))
hist(dataframe.alr$Nd,breaks=50)
hist(dataframe.clr$Nd,breaks=50)
hist(dataframe.ilr$Nd,breaks=50)
hist(dataframe.blr$Nd,breaks=50)

# Th 91 % unique values
par(mfrow=c(2,2))
hist(dataframe.alr$Th,breaks=50)
hist(dataframe.clr$Th,breaks=50)
hist(dataframe.ilr$Th,breaks=50)
hist(dataframe.blr$Th,breaks=50)


shapiro.pvalue(dataframe.blr$Th)
skewness(dataframe.blr$Th)
skewness(dataframe.blr$SiO2)
skewness(dataframe.blr$Al2O3)
skewness(dataframe.blr$Ag)


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
# What is too large ? We propose a cutoff of 25 %. If more than 25 % of values in the column were imputed, we shall not
# do outlier detection on these rows.

# How many dimensions are necessary for an observation to be considered as an outlier ? 
prop.unique <- function(x) length((unique(x)))/length(x)
apply(dataframe,2,prop.unique)

outlier.detection <- function(dataframe,cut.off=0.75,dimcutoff=5){
  # outliers detection is done only on columns with less than 25 % (cut.off) of replicated values   
  result <- apply(dataframe,2,prop.unique)>=cut.off 
  dataframe <- dataframe[result]
  df.ilr <- dataframe %>% ilr() %>% as_tibble()
  out <- list()
  for (i in 1:ncol(df.clr)){
    out[[i]] <- boxplot(df.ilr[i], plot=FALSE,range=1.5)$out #threshold should be 3.81387 if we want to have 99 \% coverage probability
  }
  # List the outliers
  list.out <- list()
  
  for (i in 1:length(df.clr)){
    list.out[[i]] <- which(df.ilr[[i]] %in% out[[i]])
  }
  vec.out <- unlist(list.out)
  table <- table(vec.out)
  outliers <- c()
  ndim <- floor(ncol(df.clr)/dimcutoff)
  for (i in 1:length(table)){
    outliers[i] <- ifelse(table[i]>=ndim,names(table)[i],NA)
  }
  outliers <- na.omit(outliers)
  outliers
  df.wo.outliers <- df.ilr[!rownames(df.clr)%in% c(outliers),]
  result <- list(df.wo.outliers,outliers)
  names(result) <- c("Dataframe without outliers","Outliers")
  return(result)
}

# Normality on dataset with outliers
df.with.outliers.shap.pval <- apply(dataframe,2,pvalue.shap.t)
sum(ifelse(df.with.outliers.shap.pval>=0.01,1,0))/ncol(dataframe) # 1 % of non rejection of the null.




df.wo.out <- outlier.detection(dataframe)
df.wo.out[[1]] %>% clr() %>% prcomp() %>% biplot()

# Normality on dataset without the 7 outliers (outside the 95 % probability range)
df.wo.outliers.shap.pval <- apply(df.wo.out[[1]],2,pvalue.shap.t)
# 36 % of non-rejection of the null H0 at alpha=0.01
sum(ifelse(df.wo.outliers.shap.pval>=0.01,1,0))/ncol(df.wo.out[[1]])
summary(prcomp(dataframe.clr))

