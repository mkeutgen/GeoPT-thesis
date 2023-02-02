# Loading Libraries
library(readr)
library(tidyverse)
library(compositions) 
library(moments)
library(nortest)
library(mvtnorm)
library(MVN)

source("importdata.R")

names(list.df)
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


var(list.df.ilr[[10]])*1000
