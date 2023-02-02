### Script Zero v1 ########
## Cleaning the data #####
##########################

#### 1. Loading libraries ####

library(readxl)
library(tidyverse)

#### 2. Importing all the late excel files ####

# extract names of the datasets
sheetnames <- excel_sheets("GeoPT Raw data late Table 1.xlsx")
# extract the dataframes and gather them in a list 
mylist <- lapply(excel_sheets("GeoPT Raw data late Table 1.xlsx"),
                 read_excel,path="GeoPT Raw data late Table 1.xlsx")
# name the dataframes
names(mylist) <- sheetnames

# initialize a series of list :
# list of the identifiers of the rocks
list.ID <- list()
# list of long.format data frames
list.lf <- list()
# list of wide.format data frames 
list.wf <- list()
# name the list
names(list.wf) <- names(mylist)

# Transform the 13 datasets 

for (i in 1:13){
  list.ID[i] <- substring(colnames(mylist[[i]][10]),1,1)
  
}
for (i in 1:13){
# delete first two rows
mylist[[i]] <-  mylist[[i]][-c(1,2),]
# clean colnames
colnames(mylist[[i]]) <- make.unique(str_remove(colnames( mylist[[i]])
                                                , "\\.+.*"))
# set the name of the 2nd column to "Unit"
# delete the unit column
mylist[[i]] <- mylist[[i]][-2]

}
for (i in 1:13){
list.lf[[i]] <- mylist[[i]] %>% pivot_longer(col=starts_with(
  list.ID[[i]]),names_to="Laboratory",values_to="Concentration",
  names_repair = "universal")
}

# clean the names of the 'Lab' column, set it to 'Lab.Code'
names(list.lf[[2]])[1]<- "Lab.Code"
names(list.lf[[12]])[1]<- "Lab.Code"
names(list.lf[[13]])[1]<- "Lab.Code"


# transform to wide format
for (i in 1:13){
list.wf[[i]] <- list.lf[[i]] %>% 
  pivot_wider(names_from = Lab.Code, 
  values_from = Concentration)
}

#### Transform manually the 3rd dataframe ####

mylist[[3]] <- mylist[[3]][-75,]
names(mylist[[3]])[1]<- "Lab.Code"
# transform to long format
list.lf[[3]] <- mylist[[3]] %>% 
  pivot_longer(col=starts_with(list.ID[[3]]),names_to="Laboratory",
               values_to="Concentration",names_repair = "universal")
# transform to wide format
list.wf[[3]] <- list.lf[[3]] %>% pivot_wider(names_from = Lab.Code, 
                                             values_from = Concentration)


#### Export dataframes to .csv ####
# set the names of the list.wf to list of df
names(list.wf) <- names(mylist)
# export dataframes in .csv
for(i in 1:length(list.wf) ){
  write.csv(list.wf[[i]],
            file = paste("/home/max/Documents/MStatistics/Thesis/Repository/data/raw/",
                         names(list.wf[i]),'.csv',sep=""))
}

#### 3. Importing all the early excel files ####

# extract names of the datasets
sheetnames <- excel_sheets("GeoPT Raw data early Table 1.xlsx")
# extract the dataframes and gather them in a list 
mylist <- lapply(excel_sheets("GeoPT Raw data early Table 1.xlsx"),
                 read_excel,path="GeoPT Raw data early Table 1.xlsx")
# name the dataframes
names(mylist) <- sheetnames

# initialize a series of list :
# list of the identifiers of the rocks
list.ID <- list()
# list of long.format data frames
list.lf <- list()
# list of wide.format data frames 
list.wf <- list()
# name the list
names(list.wf) <- names(mylist)

# Transform the 19 datasets 

for (i in 1:19){
  list.ID[i] <- substring(colnames(mylist[[i]][10]),1,1)
  
}
for (i in 1:19){
  # delete first two rows
  mylist[[i]] <-  mylist[[i]][-c(1,2),]
  # clean colnames
  colnames(mylist[[i]]) <- make.unique(
    str_remove(colnames( mylist[[i]]), "\\.+.*"))
  # set the name of the 2nd column to "Unit"
  # delete the unit column
  mylist[[i]] <- mylist[[i]][-2]
  
}
for (i in 1:19){
  list.lf[[i]] <- mylist[[i]] %>% 
    pivot_longer(col=starts_with(list.ID[[i]]),
                 names_to="Laboratory",values_to="Concentration",
                 names_repair = "universal")
}

names(list.lf) <- sheetnames


# Troubleshoot, dataframe GeoPT11 -83R (8), GeoPT27 -78R (17), 'GeoPT30 -75R (19)"
# These dataframes should be cleaned manually. For now, we will simply drop these 3.

list.lf <- list.lf[-c(8,17,19)]

# Transform to wide.format
for (i in 1:16){
  list.wf[[i]] <- list.lf[[i]] %>% pivot_wider(names_from = Lab.Code,
                                               values_from = Concentration)
}
#### Export dataframes to .csv ####
# set the names of the list.wf to list of df
names(list.wf) <- names(list.lf)
# export dataframes in .csv
for(i in 1:length(list.wf) ){
  write.csv(list.wf[[i]],
            file = paste("/home/max/Documents/MStatistics/Thesis/Repository/data/raw/",
                         names(list.wf[i]),'.csv',sep=""))
}
