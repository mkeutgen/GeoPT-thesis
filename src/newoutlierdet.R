# The Modeling
library(readr)
library(tidyverse)
library(compositions)
library(moments)
library(nortest)
library(mvtnorm)
library(MVN)
library(robCompositions)
setwd("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/")
# CLR Transform, assessing the normality assumption
# import dataset
# OUTLIER DETECTION : ! FAT DATASETS
library(readr)
library(tidyverse)
library(compositions)
library(moments)
library(nortest)
library(mvtnorm)
library(MVN)
library(robCompositions)

library(ggfortify)
library(ggrepel)
library(robustbase)
library(rrcov)
logit.transf <- function(x) log(x/(1-x))
blr <- function(dataframe) apply(dataframe,2,logit.transf)
setwd("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/")

dataframe <- read_csv("GeoPT48 .csv")

acomp(dataframe)

outlierplot(acomp(dataframe))
#outlierplot(acomp(dataframe))
#Error in covMcd(c(-0.0100033736503616, -0.0289052328168178, -0.104994139214679,  : 
#                    n <= p -- you can't be serious!


# Set negative values to 0
dataframe[dataframe < 0] <- 0

ilr.df <- dataframe %>% ilr()
ilrdf.mcd <- covMcd(ilr.df)
plot(ilrdf.mcd,which = c("distance"),classic = TRUE) 
#robCompositions::isomLR(dataframe)
# ilr(dataframe)
#outCoDa(dataframe,method = "robust")
#?outCoDa
#oD <- outCoDa(dataframe)
#outCoDa

x <- dataframe
?outlierplot()
acomp(dataframe)


View(dataframe)

# Proportion of unique values in the data frame
prop.unique <- function(x) length((unique(x)))/length(x)
apply(dataframe,2,prop.unique)
# One sees that in major elements, no values were imputed by the blr mean. This shows that the cutoff in the df cleaning
# function were successful in eliminating dubious analyses w.r.t major elements.
# However, some elements like Indium, Germanium or Bore have less than 30 % of values which are not the blr.mean... 

# Let's start with a simple dataframe.

dataframe.pr <- dataframe %>% select(!c(Rest)) %>%
  clo() %>% as_tibble()

# Multivariate normality is too hard to request. We shall request instead normality of the marginals.
hist(unique(dataframe$Yb))


# 9 observations from ~ N(0,1) and 1 outlier with value 10.
x <- c(rnorm(9),1E2,1E4)
z <- (x-mean(x))/sqrt(var(x))

