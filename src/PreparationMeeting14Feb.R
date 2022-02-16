# Meeting of 14th February : Preparation Document
# Last time, 3 topics were discussed :
# Missing Values Algorithm => implemented in the blr function.
# PCA => now can be started
# Computing Missing Values row-wise and column-wise => implemented in the function. 
library(readr)
library(tidyverse)
library(compositions)
df <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/ GeoPT48 -84Ra.csv")
df.raw <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/GeoPT48 -84Ra.csv")

df <- df[-1]
# Remove useless first column
df
df.raw
df.raw[traces] <- 10E-4 * df.raw[traces] 
# Absurdity : 
rowSums(df.raw[-1],na.rm = T)

# Our method largely diminishes proportion of major elements 
# because we now hypothetize about the unknown elements.
# Cutoffs ?

# PCA on clr transformed dataframe without missing values. To be continued. 

summary(prcomp(as_tibble(clr(df)),scale. = T))
#########################################################################################
# Structuring the thesis
# Introduction
# Chapter 1 Motivating Examples
# Chapter 2 Fundamentals of Compositional Data
# Chapter 3 Missing Values Problem
# Chapter 4 Outlier Problem
# Chapter 5 Estimation of a rock's composition, summarizing each rock by its mean
# Chapter 6 PCA with the 29 different rocks, summarized by their means, 
# to find the associations between elements (example Mg is associated with Fe,Na, etc..).
# Conclusion...
########################################################################################