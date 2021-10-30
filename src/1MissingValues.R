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


# We first restrict ourselves to the major elements 
data <- data %>% select(majors)
# We further remove FeIIO, H2O+, CO2 and LOI
data <- data %>% select(!c(LOI,CO2,`H2O+`,`Fe(II)O`)) 

# We now almost have a complete dataset (without missing values)



