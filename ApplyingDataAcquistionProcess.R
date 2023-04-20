# Apply the Data Acquisition Process Algorithm to the 29 different datasets
setwd("~/Documents/GeoPTManuscript/")

library(MASS)
library(tidyverse)
library(compositions)
library(xtable)
library(rootSolve)
library(EnvStats)
library(ggpubr)
library(quadprog)

# Blr (really a logit transformation) 
logit <- function(x){log((x)/(1-x))}
logitInv <- function(x){1/(1+exp(-x))}

blr <- function(dataframe) apply(dataframe,2,logit)
blrInv <- function(dataframe) apply(dataframe,2,logitInv)

source("DataAcquisitionProcess.R")

# Import the simplified (only first 20 elements compositions of 29 datasets in the GeoPT dataset) : 

mean.compo.simp <- readRDS(file = "src/meancompo.RData")
variance.R.space <- readRDS("src/list_var_df.RData")



mu.parameters <-  list()
for (i in 1:29){
  mu.parameters[[i]] <- mean.compo.simp[[i]]
}


names(mu.parameters) <- names(mean.compo.simp)


sim_list_GeoPT1.csv   <- list()
sim_list_GeoPT14.csv  <- list()
sim_list_GeoPT16.csv  <- list()
sim_list_GeoPT19.csv  <- list()
sim_list_GeoPT2.csv   <- list()
sim_list_GeoPT20.csv  <- list()
sim_list_GeoPT21.csv  <- list()
sim_list_GeoPT22.csv  <- list()
sim_list_GeoPT23.csv  <- list()
sim_list_GeoPT25.csv  <- list()
sim_list_GeoPT29.csv  <- list()
sim_list_GeoPT3.csv   <- list()
sim_list_GeoPT32.csv  <- list()
sim_list_GeoPT34.csv  <- list()
sim_list_GeoPT35.csv  <- list()
sim_list_GeoPT36.csv  <- list()
sim_list_GeoPT37.csv  <- list()
sim_list_GeoPT38.csv  <- list()
sim_list_GeoPT38A.csv <- list()
sim_list_GeoPT39.csv  <- list()
sim_list_GeoPT39A.csv <- list()
sim_list_GeoPT4.csv   <- list()
sim_list_GeoPT41.csv  <- list()
sim_list_GeoPT43.csv  <- list()
sim_list_GeoPT46.csv  <- list()
sim_list_GeoPT48.csv  <- list()
sim_list_GeoPT5.csv   <- list()
sim_list_GeoPT6.csv   <- list()
sim_list_GeoPT8.csv   <- list()


GJ.sim(mu.parameters[[4]],Sigma = variance.R.space[[4]]) 

View(mean.compo.full)

GJ.sim(mu.parameters[[3]],Sigma = variance.R.space[[3]]) %>% rowSums()





for (i in 1:200){sim_list_GeoPT1.csv[[i]]  <- GJ.sim(mu.parameters[[1]],Sigma = variance.R.space[[1]])   }
for (i in 1:200){sim_list_GeoPT14.csv[[i]] <- GJ.sim(mu.parameters[[2]],Sigma = variance.R.space[[2]])   }
for (i in 1:200){sim_list_GeoPT16.csv[[i]] <- GJ.sim(mu.parameters[[3]],Sigma = variance.R.space[[3]])   }
for (i in 1:200){sim_list_GeoPT19.csv[[i]] <- GJ.sim(mu.parameters[[4]],Sigma = variance.R.space[[4]])   }
for (i in 1:200){sim_list_GeoPT2.csv[[i]]   <-GJ.sim(mu.parameters[[5]],Sigma = variance.R.space[[5]])   }
for (i in 1:200){sim_list_GeoPT20.csv[[i]]  <-GJ.sim(mu.parameters[[6]],Sigma = variance.R.space[[6]])   }
for (i in 1:200){sim_list_GeoPT21.csv[[i]]  <-GJ.sim(mu.parameters[[7]],Sigma = variance.R.space[[7]])   }
for (i in 1:200){sim_list_GeoPT22.csv[[i]] <- GJ.sim(mu.parameters[[8]],Sigma = variance.R.space[[8]])   }
for (i in 1:200){sim_list_GeoPT23.csv[[i]]  <-GJ.sim(mu.parameters[[9]],Sigma = variance.R.space[[9]])   }
for (i in 1:200){sim_list_GeoPT25.csv[[i]] <-GJ.sim(mu.parameters[[10]],Sigma= variance.R.space[[10]])  }
for (i in 1:200){sim_list_GeoPT29.csv[[i]] <-GJ.sim(mu.parameters[[11]],Sigma= variance.R.space[[11]])  }
for (i in 1:200){sim_list_GeoPT3.csv[[i]]  <-GJ.sim(mu.parameters[[12]],Sigma= variance.R.space[[12]])  }
for (i in 1:200){sim_list_GeoPT32.csv[[i]]<- GJ.sim(mu.parameters[[13]],Sigma= variance.R.space[[13]])   }
for (i in 1:200){sim_list_GeoPT34.csv[[i]] <-GJ.sim(mu.parameters[[14]],Sigma= variance.R.space[[14]])   }
for (i in 1:200){sim_list_GeoPT35.csv[[i]] <-GJ.sim(mu.parameters[[15]],Sigma= variance.R.space[[15]])   }
for (i in 1:200){sim_list_GeoPT36.csv[[i]] <-GJ.sim(mu.parameters[[16]],Sigma= variance.R.space[[16]])   }
for (i in 1:200){sim_list_GeoPT37.csv[[i]] <-GJ.sim(mu.parameters[[17]],Sigma= variance.R.space[[17]])   }
for (i in 1:200){sim_list_GeoPT38.csv[[i]] <-GJ.sim(mu.parameters[[18]],Sigma= variance.R.space[[18]])   }
for (i in 1:200){sim_list_GeoPT38A.csv[[i]]<-GJ.sim(mu.parameters[[19]],Sigma= variance.R.space[[19]])   }
for (i in 1:200){sim_list_GeoPT39.csv[[i]] <-GJ.sim(mu.parameters[[20]],Sigma= variance.R.space[[20]])   }
for (i in 1:200){sim_list_GeoPT39A.csv[[i]]<-GJ.sim(mu.parameters[[21]],Sigma= variance.R.space[[21]])   }
for (i in 1:200){sim_list_GeoPT4.csv[[i]] <- GJ.sim(mu.parameters[[22]],Sigma= variance.R.space[[22]])   }
for (i in 1:200){sim_list_GeoPT41.csv[[i]]<- GJ.sim(mu.parameters[[23]],Sigma= variance.R.space[[23]])   }
for (i in 1:200){sim_list_GeoPT43.csv[[i]]<- GJ.sim(mu.parameters[[24]],Sigma= variance.R.space[[24]])   }
for (i in 1:200){sim_list_GeoPT46.csv[[i]]<- GJ.sim(mu.parameters[[25]],Sigma= variance.R.space[[25]])   }
for (i in 1:200){sim_list_GeoPT48.csv[[i]]<- GJ.sim(mu.parameters[[26]],Sigma= variance.R.space[[26]])   }
for (i in 1:200){sim_list_GeoPT5.csv[[i]] <- GJ.sim(mu.parameters[[27]],Sigma= variance.R.space[[27]])   }
for (i in 1:200){sim_list_GeoPT6.csv[[i]] <- GJ.sim(mu.parameters[[28]],Sigma= variance.R.space[[28]])   }
for (i in 1:200){sim_list_GeoPT8.csv[[i]] <- GJ.sim(mu.parameters[[29]],Sigma= variance.R.space[[29]])   }

# Write it to RDATA BECAUSE THESE ARE LISTS and we want hUNDREDS OF simulated DATASETS


#sim_list_GeoPT1.csv   <- sim_list_GeoPT1.csv  %>% bind_rows()
#sim_list_GeoPT14.csv  <- sim_list_GeoPT14.csv %>% bind_rows()
#sim_list_GeoPT16.csv  <- sim_list_GeoPT16.csv %>% bind_rows()
#sim_list_GeoPT19.csv  <- sim_list_GeoPT19.csv %>% bind_rows()
#sim_list_GeoPT2.csv   <- sim_list_GeoPT2.csv  %>% bind_rows()
#sim_list_GeoPT20.csv  <- sim_list_GeoPT20.csv %>% bind_rows()
#sim_list_GeoPT21.csv  <- sim_list_GeoPT21.csv %>% bind_rows()
#sim_list_GeoPT22.csv  <- sim_list_GeoPT22.csv %>% bind_rows()
#sim_list_GeoPT23.csv  <- sim_list_GeoPT23.csv %>% bind_rows()
#sim_list_GeoPT25.csv  <- sim_list_GeoPT25.csv %>% bind_rows()
#sim_list_GeoPT29.csv  <- sim_list_GeoPT29.csv %>% bind_rows()
#sim_list_GeoPT3.csv   <- sim_list_GeoPT3.csv  %>% bind_rows()
#sim_list_GeoPT32.csv  <- sim_list_GeoPT32.csv %>% bind_rows()
#sim_list_GeoPT34.csv  <- sim_list_GeoPT34.csv %>% bind_rows()
#sim_list_GeoPT35.csv  <- sim_list_GeoPT35.csv %>% bind_rows()
#sim_list_GeoPT36.csv  <- sim_list_GeoPT36.csv %>% bind_rows()
#sim_list_GeoPT37.csv  <- sim_list_GeoPT37.csv %>% bind_rows()
#sim_list_GeoPT38.csv  <- sim_list_GeoPT38.csv %>% bind_rows()
#sim_list_GeoPT38A.csv <- sim_list_GeoPT38A.csv%>% bind_rows()
#sim_list_GeoPT39.csv  <- sim_list_GeoPT39.csv %>% bind_rows()
#sim_list_GeoPT39A.csv <- sim_list_GeoPT39A.csv%>% bind_rows()
#sim_list_GeoPT4.csv   <- sim_list_GeoPT4.csv  %>% bind_rows()
#sim_list_GeoPT41.csv  <- sim_list_GeoPT41.csv %>% bind_rows()
#sim_list_GeoPT43.csv  <- sim_list_GeoPT43.csv %>% bind_rows()
#sim_list_GeoPT46.csv  <- sim_list_GeoPT46.csv %>% bind_rows()
#sim_list_GeoPT48.csv  <- sim_list_GeoPT48.csv %>% bind_rows()
#sim_list_GeoPT5.csv   <- sim_list_GeoPT5.csv  %>% bind_rows()
#sim_list_GeoPT6.csv   <- sim_list_GeoPT6.csv  %>% bind_rows()
#sim_list_GeoPT8.csv   <- sim_list_GeoPT8.csv  %>% bind_rows()
#
#
#write.csv(x = sim_list_GeoPT1.csv  ,file = "src/sim/sim_list_GeoPT1.csv")
#write.csv(x = sim_list_GeoPT14.csv ,file = "src/sim/sim_list_GeoPT14.csv")
#write.csv(x = sim_list_GeoPT16.csv ,file = "src/sim/sim_list_GeoPT16.csv")
#write.csv(x = sim_list_GeoPT19.csv ,file = "src/sim/sim_list_GeoPT19.csv")
#write.csv(x = sim_list_GeoPT2.csv  ,file = "src/sim/sim_list_GeoPT2.csv")
#write.csv(x = sim_list_GeoPT20.csv ,file = "src/sim/sim_list_GeoPT20.csv")
#write.csv(x = sim_list_GeoPT21.csv ,file = "src/sim/sim_list_GeoPT21.csv")
#write.csv(x = sim_list_GeoPT22.csv ,file = "src/sim/sim_list_GeoPT22.csv")
#write.csv(x = sim_list_GeoPT23.csv ,file = "src/sim/sim_list_GeoPT23.csv")
#write.csv(x = sim_list_GeoPT25.csv ,file = "src/sim/sim_list_GeoPT25.csv")
#write.csv(x = sim_list_GeoPT29.csv ,file = "src/sim/sim_list_GeoPT29.csv")
#write.csv(x = sim_list_GeoPT3.csv  ,file = "src/sim/sim_list_GeoPT3.csv")
#write.csv(x = sim_list_GeoPT32.csv ,file = "src/sim/sim_list_GeoPT32.csv")
#write.csv(x = sim_list_GeoPT34.csv ,file = "src/sim/sim_list_GeoPT34.csv")
#write.csv(x = sim_list_GeoPT35.csv ,file = "src/sim/sim_list_GeoPT35.csv")
#write.csv(x = sim_list_GeoPT36.csv ,file = "src/sim/sim_list_GeoPT36.csv")
#write.csv(x = sim_list_GeoPT37.csv ,file = "src/sim/sim_list_GeoPT37.csv")
#write.csv(x = sim_list_GeoPT38.csv ,file = "src/sim/sim_list_GeoPT38.csv")
#write.csv(x = sim_list_GeoPT38A.csv,file = "src/sim/sim_list_GeoPT38A.csv")
#write.csv(x = sim_list_GeoPT39.csv ,file = "src/sim/sim_list_GeoPT39.csv")
#write.csv(x = sim_list_GeoPT39A.csv,file = "src/sim/sim_list_GeoPT39A.csv")
#write.csv(x = sim_list_GeoPT4.csv  ,file = "src/sim/sim_list_GeoPT4.csv")
#write.csv(x = sim_list_GeoPT41.csv ,file = "src/sim/sim_list_GeoPT41.csv")
#write.csv(x = sim_list_GeoPT43.csv ,file = "src/sim/sim_list_GeoPT43.csv")
#write.csv(x = sim_list_GeoPT46.csv ,file = "src/sim/sim_list_GeoPT46.csv")
#write.csv(x = sim_list_GeoPT48.csv ,file = "src/sim/sim_list_GeoPT48.csv")
#write.csv(x = sim_list_GeoPT5.csv  ,file = "src/sim/sim_list_GeoPT5.csv")
#write.csv(x = sim_list_GeoPT6.csv  ,file = "src/sim/sim_list_GeoPT6.csv")
#write.csv(x = sim_list_GeoPT8.csv  ,file = "src/sim/sim_list_GeoPT8.csv")

saveRDS(object = sim_list_GeoPT1.csv  ,file = "src/sim/sim_list_GeoPT1.RData")
saveRDS(object = sim_list_GeoPT14.csv ,file = "src/sim/sim_list_GeoPT14.RData") 
saveRDS(object = sim_list_GeoPT16.csv ,file = "src/sim/sim_list_GeoPT16.RData") 
saveRDS(object = sim_list_GeoPT19.csv ,file = "src/sim/sim_list_GeoPT19.RData") 
saveRDS(object = sim_list_GeoPT2.csv  ,file = "src/sim/sim_list_GeoPT2.RData")
saveRDS(object = sim_list_GeoPT20.csv ,file = "src/sim/sim_list_GeoPT20.RData") 
saveRDS(object = sim_list_GeoPT21.csv ,file = "src/sim/sim_list_GeoPT21.RData") 
saveRDS(object = sim_list_GeoPT22.csv ,file = "src/sim/sim_list_GeoPT22.RData") 
saveRDS(object = sim_list_GeoPT23.csv ,file = "src/sim/sim_list_GeoPT23.RData") 
saveRDS(object = sim_list_GeoPT25.csv ,file = "src/sim/sim_list_GeoPT25.RData") 
saveRDS(object = sim_list_GeoPT29.csv ,file = "src/sim/sim_list_GeoPT29.RData") 
saveRDS(object = sim_list_GeoPT3.csv  ,file = "src/sim/sim_list_GeoPT3.RData")
saveRDS(object = sim_list_GeoPT32.csv ,file = "src/sim/sim_list_GeoPT32.RData") 
saveRDS(object = sim_list_GeoPT34.csv ,file = "src/sim/sim_list_GeoPT34.RData") 
saveRDS(object = sim_list_GeoPT35.csv ,file = "src/sim/sim_list_GeoPT35.RData") 
saveRDS(object = sim_list_GeoPT36.csv ,file = "src/sim/sim_list_GeoPT36.RData") 
saveRDS(object = sim_list_GeoPT37.csv ,file = "src/sim/sim_list_GeoPT37.RData") 
saveRDS(object = sim_list_GeoPT38.csv ,file = "src/sim/sim_list_GeoPT38.RData") 
saveRDS(object = sim_list_GeoPT38A.csv,file = "src/sim/sim_list_GeoPT38A.RData")
saveRDS(object = sim_list_GeoPT39.csv ,file = "src/sim/sim_list_GeoPT39.RData") 
saveRDS(object = sim_list_GeoPT39A.csv,file = "src/sim/sim_list_GeoPT39A.RData")
saveRDS(object = sim_list_GeoPT4.csv  ,file = "src/sim/sim_list_GeoPT4.RData")
saveRDS(object = sim_list_GeoPT41.csv ,file = "src/sim/sim_list_GeoPT41.RData") 
saveRDS(object = sim_list_GeoPT43.csv ,file = "src/sim/sim_list_GeoPT43.RData") 
saveRDS(object = sim_list_GeoPT46.csv ,file = "src/sim/sim_list_GeoPT46.RData") 
saveRDS(object = sim_list_GeoPT48.csv ,file = "src/sim/sim_list_GeoPT48.RData") 
saveRDS(object = sim_list_GeoPT5.csv  ,file = "src/sim/sim_list_GeoPT5.RData")
saveRDS(object = sim_list_GeoPT6.csv  ,file = "src/sim/sim_list_GeoPT6.RData")
saveRDS(object = sim_list_GeoPT8.csv  ,file = "src/sim/sim_list_GeoPT8.RData")


































