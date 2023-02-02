#### Apply Weltje Algorithm to create the "processed" missing.values-free dataframes. ####
#########################################################################################
library("tidyverse")
library("readr")
library("compositions")

#####
source(file="WeltjeAlgorithm.R")


# Get the full dataframe
filePaths <- list.files("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/",
                        "\\.csv$", full.names = TRUE)[-1]
filePaths
names.df <- c()
for (i in 1:length(filePaths)){
  names.df[i] <- substr(filePaths[i],65,nchar(filePaths[i]))
}

# Import all dataframes at once
list.df <- list()
for (i in 1:length(filePaths)){
  list.df[[i]] <- read_csv(filePaths[i])}

names(list.df) <- names.df


# Column Repair
list.df[[3]]<- list.df[[3]] %>% rename(Fe2O3T = Fe2O3)
list.df[[4]]<- list.df[[4]] %>% rename(Fe2O3T = Fe2O3)
list.df[[6]]<- list.df[[6]] %>% rename(Fe2O3T = Fe2O3)
list.df[[7]]<- list.df[[7]] %>% rename(Fe2O3T = Fe2O3)
list.df[[8]]<- list.df[[8]] %>% rename(Fe2O3T = Fe2O3)
list.df[[9]]<- list.df[[9]] %>% rename(Fe2O3T = Fe2O3)
list.df[[10]]<- list.df[[10]] %>% rename(Fe2O3T = Fe2O3)

# Initiate list of imputed.df
imputed.df <- list()

# Applying Algo.


imputed.df[[1]] <- GWAlgo(list.df[[1]] )
imputed.df[[2]] <- GWAlgo(list.df[[2]] )
imputed.df[[3]] <- GWAlgo(list.df[[3]] )
imputed.df[[4]] <- GWAlgo(list.df[[4]] )
imputed.df[[5]] <- GWAlgo(list.df[[5]] )
imputed.df[[6]] <- GWAlgo(list.df[[6]] )
imputed.df[[7]] <- GWAlgo(list.df[[7]] )
imputed.df[[8]] <- GWAlgo(list.df[[8]] )
imputed.df[[9]] <- GWAlgo(list.df[[9]] )
imputed.df[[10]]<- GWAlgo.mod(list.df[[10]])
imputed.df[[11]]<- GWAlgo.mod(list.df[[11]])
imputed.df[[12]]<- GWAlgo.mod(list.df[[12]])
imputed.df[[13]]<- GWAlgo.mod(list.df[[13]])
imputed.df[[14]]<- GWAlgo.mod(list.df[[14]])
imputed.df[[15]]<- GWAlgo.mod(list.df[[15]])
imputed.df[[16]]<- GWAlgo.mod(list.df[[16]])
imputed.df[[17]]<- GWAlgo.mod(list.df[[17]])


imputed.df[[18]]<- GWAlgo(list.df[[18]])
imputed.df[[19]]<- GWAlgo(list.df[[19]])
imputed.df[[20]]<- GWAlgo(list.df[[20]])
imputed.df[[21]]<- GWAlgo(list.df[[21]])
imputed.df[[22]]<- GWAlgo(list.df[[22]])
imputed.df[[23]]<- GWAlgo(list.df[[23]])
imputed.df[[24]]<- GWAlgo(list.df[[24]])
imputed.df[[25]]<- GWAlgo(list.df[[25]])
imputed.df[[26]]<- GWAlgo(list.df[[26]])
imputed.df[[27]]<- GWAlgo.mod(list.df[[27]])
imputed.df[[28]]<- GWAlgo(list.df[[28]])
imputed.df[[29]]<- GWAlgo(list.df[[29]])

names(imputed.df)

# Export the 29 dataframes to csv
names(imputed.df) <- names(list.df)
names(imputed.df) <- substr(names(imputed.df),1,8)
imputed.df["GeoPT1"] <-  imputed.df["GeoPT1 -"] 
imputed.df["GeoPT2"] <-  imputed.df["GeoPT2 -"] 
imputed.df["GeoPT3"] <-  imputed.df["GeoPT3 -"] 
imputed.df["GeoPT4"] <-  imputed.df["GeoPT4 -"] 
imputed.df["GeoPT5"] <-  imputed.df["GeoPT5 -"] 
imputed.df["GeoPT6"] <-  imputed.df["GeoPT6 -"] 
imputed.df["GeoPT8"] <-  imputed.df["GeoPT8 -"] 



for (i in 1:length(imputed.df)){
  imputed.df[[i]][imputed.df[[i]] < 0] <- 0
}
for (i in 1:length(imputed.df) ){
  write.csv(imputed.df[[i]],file = paste0("~/Documents/MStatistics/MA2/
                                          Thesis/Repository/data/processed/",
                                          names(imputed.df[i]),".csv"),row.names = FALSE)
}



