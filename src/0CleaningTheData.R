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
mylist <- lapply(excel_sheets("GeoPT Raw data late Table 1.xlsx"),read_excel,path="GeoPT Raw data late Table 1.xlsx")
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
colnames(mylist[[i]]) <- make.unique(str_remove(colnames( mylist[[i]]), "\\.+.*"))
# set the name of the 2nd column to "Unit"
# delete the unit column
mylist[[i]] <- mylist[[i]][-2]

}
for (i in 1:13){
list.lf[[i]] <- mylist[[i]] %>% pivot_longer(col=starts_with(list.ID[[i]]),names_to="Laboratory",values_to="Concentration",names_repair = "universal")
}

# clean the names of the 'Lab' column, set it to 'Lab.Code'
names(list.lf[[2]])[1]<- "Lab.Code"
names(list.lf[[12]])[1]<- "Lab.Code"
names(list.lf[[13]])[1]<- "Lab.Code"


# transform to wide format
for (i in 1:13){
list.wf[[i]] <- list.lf[[i]] %>% pivot_wider(names_from = Lab.Code, values_from = Concentration)
}

#### Transform manually the 3rd dataframe ####

mylist[[3]] <- mylist[[3]][-75,]
names(mylist[[3]])[1]<- "Lab.Code"
# transform to long format
list.lf[[3]] <- mylist[[3]] %>% pivot_longer(col=starts_with(list.ID[[3]]),names_to="Laboratory",values_to="Concentration",names_repair = "universal")
# transform to wide format
list.wf[[3]] <- list.lf[[3]] %>% pivot_wider(names_from = Lab.Code, values_from = Concentration)


#### Export dataframes to .csv ####
# set the names of the list.wf to list of df
names(list.wf) <- names(mylist)
# export dataframes in .csv
for(i in 1:length(list.wf) ){
  write.csv(list.wf[[i]],file = paste("/home/max/Documents/MStatistics/Thesis/Repository/data/raw/",names(list.wf[i]),'.csv',sep=""))
}


#### Create a dataframe containing key infos on the sample analyzed ####
dput(summary$Name.of.the.rock)
ID <- c("GeoPT48", "GeoPT46", "GeoPT43", "GeoPT41", "GeoPT39", "GeoPT39A", 
          "GeoPT38", "GeoPT38A", "GeoPT37", "GeoPT36", "GeoPT35", "GeoPT34", 
          "GeoPT32")
name <- c("Monzonite", "Granodiorite", "Dolerite", "Andesite", "Syenite", 
          "Nepheline Syenite", "Ardnamurchan gabbro", "Modified harzburgite", 
          "Rhyolite", "Gabbro", "Tonalite", "Granite", "Woodstock basalt")
year <- c(2021L, 2020L, 2018L, 2017L, 2016L, 2016L, 2016L, 2016L, 2015L, 
          2015L, 2014L, 2014L, 2013L)
samplename <- c("MzBP-1", "HG-1", "ADS-1", "ORA-1", "SyMP-1", "MNS-1", "OU-7", 
                "HARZ01", "ORPT-1", "GSM-1", "TLM-1", "GRI-1", "WG-1")
description <- tibble(ID,name,year,samplename)

# export to .csv 

write.csv(description,"/home/max/Documents/MStatistics/Thesis/Repository/data/raw/description.csv")
