### Script One v0.1 ########
## Cleaning the data #######
## Processing the data #####
## Handling missing values #
############################


# Load libraries

library(tidyverse)
library(readr)

# Start with most recent dataset

data <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/GeoPT48 -84Ra.csv")

# Vectors with major and minor elements
majors <- c("SiO2","TiO2","Al2O3","Fe2O3T","Fe(II)O","MnO","MgO","CaO","Na2O","K2O","P2O5","H2O+","CO2","LOI")
minors <- c("Ag","As","B","Ba","Be","Bi","Br","C(org)","C(tot)","Cd","Ce","Cl","Co","Cr",
            "Cs","Cu","Dy","Er","Eu","F","Ga","Gd","Ge","Hf","Hg","Ho","I","In","La","Li",
            "Lu","Mo","Nb","Nd","Ni","Pb","Pd","Pr","Rb","S","Sb","Sc","Se","Sm","Sn","Sr",
            "Ta","Tb","Te","Th","Tl","Tm","U","V","W","Y","Yb","Zn","Zr")

# Homogeneity of units : ppm -> %
minors.df <- 10^(-4)*select(data,one_of(minors))


true.majors <- c("SiO2","TiO2","Al2O3","Fe2O3T",
                 "MnO","MgO","CaO","Na2O","K2O","P2O5")
amounts <- data[,true.majors]

Other <- rowSums(minors.df,na.rm=TRUE)
# df
df <- data.frame(cbind(amounts,Other))

View(df)
library(compositions)


# We further remove FeIIO, H2O+, CO2 and LOI
data <- data %>% select(!c(LOI,CO2,`H2O+`,`Fe(II)O`)) 

# We now almost have a complete dataset (without missing values). Compute the proportion of missing values
# per columns 
prop.na <- function (x){
  # compute the proportion of NA in a vector x
  sum(length(which(is.na(x))))/length(x)
}


na_count <-sapply(df, prop.na)
na_count
# The proportion of missing values for the major elements is between 12.3 and 18.5 %. 
# This relatively low proportion of missing values allow us to replace it by a naive estimate, like 
# the geometric mean of its non-missing observation 

geomean <- function(x){
  # compute the geometric mean of a vector
  exp(mean(log(x),na.rm=TRUE)) 
}

geom_mean <- sapply(df,geomean)
geom_mean
# 
# if is.na(data) == TRUE
# for i in 1:length(data)
# {
# data[i] = ifelse(test = is.na(data[i]),yes = geom_mean[i],no = data[i]) 
# }

