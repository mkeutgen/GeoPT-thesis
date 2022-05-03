# Weltje Algorithm Illustration
library(xtable)
library(readr)
library(tidyverse)
chocolate <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/src/illustrationalgo/chocolate.csv")
row.names(chocolate) <- chocolate[1]
chocolate <- chocolate[-1]
rownames(chocolate) <- NULL
xtable(chocolate[-1]/100)
chocolate
chocolate0 <- 1/100*chocolate
# Check for feasibility constraint :
rowSums(chocolate0,na.rm = T)
# Only rows 1 & 2 will be affected by the algorithm

imputed.chocolate0 <- impute_na(chocolate0)
rowsum(imputed.chocolate0)
imputed.chocolate0
rowSums(imputed.chocolate0)
xtable(apply(chocolate0,2,blr.mean) %>% as_tibble())

xtable(round(imputed.chocolate0,3))
# Divide the entries in rows which do not satisfy the feasibility constraint by the rowsum 

dataframe <- chocolate0
imputed.dataframe <- imputed.chocolate0
scaling <- function(dataframe){
  imputed.dataframe <- impute_na(dataframe)
  for (i in 1:nrow(dataframe)){
    # if statement check if composition is feasible, if not feasible, scale it
    dataframe[i,] <- if (rowSums(imputed.dataframe)[i]>1){
      sweep(dataframe[i,],2,rowSums(imputed.dataframe,na.rm=T)[i],"/")
    } # else no need to scale it 
    else {
      dataframe[i,]
    }
  }
return(dataframe)    
}
View(chocolate0)
chocolate1 <- scaling(chocolate0)
ad1 <- abs(sum(rowSums(chocolate1,na.rm=T))-sum(rowSums(chocolate0,na.rm=T)))

chocolate2 <- scaling(chocolate1)
ad2 <- abs(sum(rowSums(chocolate2,na.rm=T))-sum(rowSums(chocolate1,na.rm=T)))

chocolate3 <- scaling(chocolate2)
ad3 <- abs(sum(rowSums(chocolate3,na.rm=T))-sum(rowSums(chocolate2,na.rm=T)))

chocolate4 <- scaling(chocolate3)
ad4 <- abs(sum(rowSums(chocolate4,na.rm=T))-sum(rowSums(chocolate3,na.rm=T)))

conv.df <- cbind(c(1,2,3,4),c(ad1,ad2,ad3,ad4)) %>% as_tibble() 
colnames(conv.df) <- c("Iterations","DCS")

conv.df
ggplot(conv.df,aes(y=DCS,x=Iterations))+geom_line()+theme_bw()+scale_y_log10()
chocolate4
# Quadratic convergence rate
toLatex(chocolate0)

rowSums(imputed.chocolate0,na.rm = T)
choco.blrmean <- apply(chocolate[-1],2,blr.mean)

