#### Testing Variance-Covariance is equal across datasets ####

library(tidyverse)

fileNames <- list.files("~/Documents/MA2/Thesis/Repository/data/outfree/", "\\.csv$", full.names = FALSE)

filePaths <- list.files("~/Documents/MA2/Thesis/Repository/data/outfree", "\\.csv$", full.names = TRUE)

fileNames
filePaths
list.df <- list()

for (i in 1:length(filePaths)){
  list.df[[i]] <- read_csv(file = filePaths[[i]])
}

names(list.df) <- fileNames
# Remove NA.df, why does it appear in the first place ? 
list.df <- list.df[names(list.df)!=c("NA.csv")]

# Different data transformations :
# ALR additive log ratios 
# CLR centred log ratios
# ILR isometric log ratios


list.df.ilr <- list() 
