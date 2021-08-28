S#######################
## Script Zero ########
## Cleaning the data ##
#######################

# 1. Attempt 

library(readxl)
library(tidyverse)
library(readxl)
library(tidyverse)

# Import all datasets

path <- read_excel("GeoPT Raw data late Table 1.xlsx")
sheetnames <- excel_sheets("GeoPT Raw data late Table 1.xlsx")
mylist <- lapply(excel_sheets("GeoPT Raw data late Table 1.xlsx"),read_excel,path="GeoPT Raw data late Table 1.xlsx")
# name the dataframes
names(mylist) <- sheetnames
# names(mylist) <- substr(names(mylist),1,7)
names.df <- names(mylist)
names(list.wf) <- names(mylist)
# mylist
list.identifier <- list()
list.ID <- list()
list.lf <- list()
list.wf <- list()
# For LOOP of DOOM

for (i in 1:13){
  list.ID[i] <- substring(colnames(mylist[[i]][10]),1,1)
  
}
for (i in 1:13){
# delete first two rows
mylist[[i]] <-  mylist[[i]][-c(1,2),]
# clean colnames
colnames(mylist[[i]]) <- make.unique(str_remove(colnames( mylist[[i]]), "\\.+.*"))
# set the name of the 2nd column to "Unit"
# names(mylist[[i]])[2]<- "Unit"
# delete the unit column
mylist[[i]] <- mylist[[i]][-2]

}
for (i in 1:13){
list.lf[[i]] <- mylist[[i]] %>% pivot_longer(col=starts_with(list.ID[[i]]),names_to="Laboratory",values_to="Concentration",names_repair = "universal")
}

# Lab.Code
names(list.lf[[2]])[1]<- "Lab.Code"
names(list.lf[[12]])[1]<- "Lab.Code"
names(list.lf[[13]])[1]<- "Lab.Code"


# transform to wide format
for (i in 1:13){
list.wf[[i]] <- list.lf[[i]] %>% pivot_wider(names_from = Lab.Code, values_from = Concentration)
}

# Do manually the 3rd dataframe
mylist[[3]] <- mylist[[3]][-75,]
names(mylist[[3]])[1]<- "Lab.Code"
# transform to long format
list.lf[[3]] <- mylist[[3]] %>% pivot_longer(col=starts_with(list.ID[[3]]),names_to="Laboratory",values_to="Concentration",names_repair = "universal")
# transform to wide format
list.wf[[3]] <- list.lf[[3]] %>% pivot_wider(names_from = Lab.Code, values_from = Concentration)
#### Ended Here, To resume tomorrow, 28th august 2021
# write to csv


for(i in 1:length(list.wf) ){
  write.csv(list.wf[[i]],file = paste("/home/max/Documents/MStatistics/Thesis/Repository/Data/",names(list.wf[i]),'.csv',sep=""))
}
GeoPT36_78Ra <- read_csv("Data/GeoPT36 -78Ra.csv")
View(GeoPT36_78Ra)
print("CONGRATS U BOI") 
