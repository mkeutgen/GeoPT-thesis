
# Get the full dataframe
filePaths <- list.files("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/", "\\.csv$", full.names = TRUE)[-1]
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







imputed.df <- list()

imputed.df[[1]] <- blr_imputation(list.df[[1]])[[1]]
# imputed.df[[2]] <- blr_imputation(list.df[[2]])[[1]]

imputed.df[[3]] <- blr_imputation(list.df[[3]])[[1]]
imputed.df[[4]] <- blr_imputation(list.df[[4]])[[1]]
imputed.df[[5]] <- blr_imputation(list.df[[5]])[[1]]
imputed.df[[6]] <- blr_imputation(list.df[[6]])[[1]]

imputed.df[[7]] <- blr_imputation(list.df[[7]])[[1]]
imputed.df[[8]] <- blr_imputation(list.df[[8]])[[1]]
# Do not work with df 9 
imputed.df[[9]] <- blr_imputation(list.df[[9]])[[1]]


imputed.df[[10]] <- blr_imputation(list.df[[10]])[[1]]
# Issue with 11 and 12
imputed.df[[11]] <- blr_imputation.mod(list.df[[11]])[[1]]
imputed.df[[12]] <- blr_imputation.mod(list.df[[12]])[[1]]
imputed.df[[13]] <- blr_imputation.mod(list.df[[13]])[[1]]
imputed.df[[14]] <- blr_imputation.mod(list.df[[14]])[[1]]
imputed.df[[15]] <- blr_imputation.mod(list.df[[15]])[[1]]
imputed.df[[16]] <- blr_imputation.mod(list.df[[16]])[[1]]
imputed.df[[17]] <- blr_imputation.mod(list.df[[17]])[[1]]
imputed.df[[18]] <- blr_imputation.mod(list.df[[18]])[[1]]
imputed.df[[19]] <- blr_imputation.mod(list.df[[19]])[[1]]
imputed.df[[20]] <- blr_imputation.mod(list.df[[20]])[[1]]
imputed.df[[21]] <- blr_imputation.mod(list.df[[21]])[[1]]
imputed.df[[22]] <- blr_imputation.mod(list.df[[22]])[[1]]
imputed.df[[23]] <- blr_imputation.mod(list.df[[23]])[[1]]
imputed.df[[24]] <- blr_imputation(list.df[[24]])[[1]]
imputed.df[[25]] <- blr_imputation(list.df[[25]])[[1]]
imputed.df[[26]] <- blr_imputation(list.df[[26]])[[1]]
imputed.df[[27]] <- blr_imputation.mod(list.df[[27]])[[1]]
imputed.df[[28]] <- blr_imputation(list.df[[28]])[[1]]
imputed.df[[29]] <- blr_imputation(list.df[[29]])[[1]]

# Export the 29 dataframes to csv
names(imputed.df) <- names(list.df)
for(i in 1:length(imputed.df) ){
  write.csv(imputed.df[[i]],file = paste("~/Documents/MStatistics/MA2/Thesis/Repository/data/processed/",names(imputed.df[i])))
}
