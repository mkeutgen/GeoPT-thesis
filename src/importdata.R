library(readr)
library(tidyverse)
library(compositions)
library(moments)
library(nortest)
library(mvtnorm)


# Code to import all processed dataframes, returns a named list of the processed dataframes.
# Blr (really a logit transformation) 
logit.transf <- function(x) log(x/(1-x))
blr <- function(dataframe) apply(dataframe,2,logit.transf)


filePaths <- list.files("~/Documents/MA2/Thesis/Repository/data/processed/", "\\.csv$", full.names = TRUE)
# Get the files name :
fileNames <- list.files("~/Documents/MA2/Thesis/Repository/data/processed/", "\\.csv$", full.names = FALSE)


list.df <- list()

for (i in 1:length(filePaths)){
  list.df[[i]] <- read_csv(file = filePaths[[i]])
}

names(list.df) <- fileNames
# Remove NA.df, why does it appear in the first place ? 
list.df <- list.df[names(list.df)!=c("NA.csv")]

## Different data transformations :
## ALR additive log ratios 
## CLR centred log ratios
## ILR isometric log ratios
#
#
#list.df.ilr <- list() 
#list.df.clr <- list()  
#list.df.alr <- list()  
#list.df.blr <- list()
#
#for (i in 1:length(list.df)){
#  list.df.ilr[[i]] <- list.df[[i]] %>% ilr() %>% as_tibble()
#}
#
#for (i in 1:length(list.df)){
#  list.df.clr[[i]] <- list.df[[i]] %>% clr() %>% as_tibble()
#}
#
#for (i in 1:length(list.df)){
#  list.df.alr[[i]] <- list.df[[i]] %>% alr() %>% as_tibble()
#}
#
#for (i in 1:length(list.df)){
#  list.df.blr[[i]] <- list.df[[i]] %>% blr() %>% as_tibble()
#}
#
#
