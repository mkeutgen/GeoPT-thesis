#### Testing Variance-Covariance is equal across datasets ####
#### A. IMPORTING DATASETS

library(tidyverse)
library(compositions)
library(MVTests)
library(moments)
library(mvnormtest)

fileNames <- list.files("~/Documents/MA2/Thesis/Repository/data/outfree/", "\\.csv$", full.names = FALSE)
filePaths <- list.files("~/Documents/MA2/Thesis/Repository/data/outfree", "\\.csv$", full.names = TRUE)
list.df <- list()


ordered.names <- c("GeoPT1.csv","GeoPT2.csv","GeoPT3.csv","GeoPT4.csv","GeoPT5.csv","GeoPT6.csv",
                   "GeoPT8.csv","GeoPT14 .csv","GeoPT16 .csv","GeoPT19 .csv","GeoPT20 .csv",
                   "GeoPT21 .csv","GeoPT22 .csv","GeoPT23 .csv","GeoPT25 .csv", "GeoPT29 .csv",
                   "GeoPT32 .csv","GeoPT34 .csv","GeoPT35 .csv","GeoPT36 .csv","GeoPT37 .csv","GeoPT38 .csv",
                   "GeoPT38A.csv","GeoPT39 .csv","GeoPT39A.csv","GeoPT41 .csv","GeoPT43.csv",
                   "GeoPT46 .csv","GeoPT48 .csv")

for (i in 1:length(filePaths)){
  list.df[[i]] <- read_csv(file = filePaths[[i]])
}

names(list.df) <- fileNames


# Amalgamating and finding a common composition.


list.df.ilr <- list() 

tab <- table(unlist(lapply(list.df, names)))
common.names <- names(tab[tab == length(list.df)])

#  As a first step, we perform a very conservative amalgamation,
# all traces elements are amalgamated.

majors <- c("SiO2","TiO2","Al2O3","Fe(II)O","Fe2O3T","MnO","MgO","CaO","Na2O",
            "K2O","P2O5")

#minors <- setdiff(common.names,majors)
minors <- c("Ag", "As", "Ba", "Be", "Cd", "Ce", "Co", "Cr", "Cs", "Cu", 
  "Dy", "Er", "Eu", "Ga", "Gd", "Ge", "Hf", "Ho", "La", "Li", "Lu", 
  "Mo", "Nb", "Nd", "Ni", "Pb", "Pr", "Rb", "Sb", "Sc", "Sm", "Sn", 
  "Sr", "Ta", "Tb", "Th", "Tl", "Tm", "U", "V", "W", "Y", "Yb", 
  "Zn", "Zr")
#ree <- c("La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu")
#chalco <- c("Sn","As","Pb","Zn","Ag")

traces <- list.df[[1]] %>% select(all_of(minors)) %>% rowSums() 
reduced.df <- list.df[[1]] %>% select(all_of(majors))

reduced.df <- cbind(reduced.df,traces) %>% clo() %>% as_tibble() 

Xbar <- reduced.df %>% ilr() %>% colSums()
Sigma <- reduced.df %>% ilr() %>% var()

amalgam_fn <- function(dataframe){
  traces <- dataframe %>% select(minors) %>% rowSums() 
  reduced.df <- dataframe %>% select(majors)
  reduced.df <- cbind(reduced.df,traces) %>% clo() %>% as_tibble() 
  return(reduced.df)
}

reduced.df.l <- list()

for (i in 1:length(list.df)){
reduced.df.l[[i]] <- amalgam_fn(list.df[[i]])
  }



# ILR coordinates of reduced.df.l :
ilr.red.l <- list()
for (i in 1:length(list.df)){
  ilr.red.l[[i]] <- reduced.df.l[[i]] %>% ilr() %>% as_tibble()
}

# B. PRIOR TO TESTING EQUALITY OF COVARIANCE MATRICES, WE SHOULD FIRST
# CHECK HOW REASONABLE
# THE NORMALITY ASSUMPTION ACTUALLY IS. WE PERFORM SHAPIRO TEST
# ON THE MARGINALS
shap.pval <- function(x) shapiro.test(x)$p.value


propH0rejnorm <- function(dataframe){
  count.H0rej <- ifelse(round(apply(dataframe,2,shap.pval),10)<=0.01,1,0) %>% sum()
propH0rej <- count.H0rej/ncol(dataframe)
return(propH0rej)
}

# Testing null hypothesis (Shapiro Wilk test) on marginals of ilr
# transformed datasets
round(lapply(ilr.red.l,propH0rejnorm) %>% as_vector(),digits=3)
# Transform into matrix and test equality of cov matrices 
ilr.red.l <- lapply(ilr.red.l,as.matrix)
Schott2001(ilr.red.l) # pval significant
Schott2007(ilr.red.l) # pval significant
Srivastava2014(ilr.red.l) # pval significant



#### C. Testing for homogeneity of covariance matrices
library(covTestR)

names(iris_ls) <- irisSpecies
iris_ls
for (i in 1:length(ilr.red.l)){
  ilr.red.l[[i]] <- ilr.red.l[[i]] %>% as.matrix()
}
# So now we've a list of the ilr coordinates of the 29 dataframes. 

Ahmad2017(ilr.red.l[c(1,2)]) # not significant
Schott2001(ilr.red.l[c(1,2)]) # significant
Schott2007(ilr.red.l[c(1,2)]) # significant
Srivastava2005(ilr.red.l[c(1,2)]) 
Srivastava2007(ilr.red.l[c(1,2)]) # not significant
Srivastava2014(ilr.red.l[c(1,2)]) # significant
SrivastavaYanagihara2010(ilr.red.l[c(1,2)]) # not significant
Ishii2016(ilr.red.l[c(1,2)]) # not significant
Chaipitak2013(ilr.red.l[c(1,2)]) # not significant




list.kurt <- list()
for (i in 1:length(ilr.red.l)){
  list.kurt[[i]] <- kurtosis(ilr.red.l[[i]])
}
df <- do.call("rbind",list.kurt) 
df
Ishii

Schott2001(ilr.red.l,parameter=28) # p-value < 2.2e-16
Schott2007(ilr.red.l,parameter=28) # p-value < 2.2e-16
Srivastava2014(ilr.red.l) # p-value < 2.2e-16

# All tests are significant, null hypothesis of equality of covariance
# matrices is rejected at all levels.

# D. creating the "meancomp" dataset
par(mfrow=c(2,2))
hist(ilr.red.l[[1]][,1])
hist(ilr.red.l[[1]][,2])
hist(ilr.red.l[[1]][,3])
hist(ilr.red.l[[1]][,4])

lapply(ilr.red.l,nrow)


ilr.list.df <- lapply(list.df,select(all_of(common.names)))

list.df$GeoPT1.csv %>% dplyr::select(all_of(common.names))
lapply(list.df,select(all_of(common.names)))
co.list.df <- list()
co.list.df[[1]]  <- list.df[[1]] %>% select(all_of(common.names))
co.list.df[[2]]  <- list.df[[2]] %>% select(all_of(common.names))
co.list.df[[3]]  <- list.df[[3]] %>% select(all_of(common.names))
co.list.df[[4]]  <- list.df[[4]] %>% select(all_of(common.names))
co.list.df[[5]]  <- list.df[[5]] %>% select(all_of(common.names))
co.list.df[[6]]  <- list.df[[6]] %>% select(all_of(common.names))
co.list.df[[7]]  <- list.df[[7]] %>% select(all_of(common.names))
co.list.df[[8]]  <- list.df[[8]] %>% select(all_of(common.names))
co.list.df[[9]]  <- list.df[[9]] %>% select(all_of(common.names))
co.list.df[[10]] <- list.df[[10]] %>% select(all_of(common.names))
co.list.df[[11]] <- list.df[[11]] %>% select(all_of(common.names))
co.list.df[[12]] <- list.df[[12]] %>% select(all_of(common.names))
co.list.df[[13]] <- list.df[[13]] %>% select(all_of(common.names))
co.list.df[[14]] <- list.df[[14]] %>% select(all_of(common.names))
co.list.df[[15]] <- list.df[[15]] %>% select(all_of(common.names))
co.list.df[[16]] <- list.df[[16]] %>% select(all_of(common.names))
co.list.df[[17]] <- list.df[[17]] %>% select(all_of(common.names))
co.list.df[[18]] <- list.df[[18]] %>% select(all_of(common.names))
co.list.df[[19]] <- list.df[[19]] %>% select(all_of(common.names))
co.list.df[[20]] <- list.df[[20]] %>% select(all_of(common.names))
co.list.df[[21]] <- list.df[[21]] %>% select(all_of(common.names))
co.list.df[[22]] <- list.df[[22]] %>% select(all_of(common.names))
co.list.df[[23]] <- list.df[[23]] %>% select(all_of(common.names))
co.list.df[[24]] <- list.df[[24]] %>% select(all_of(common.names))
co.list.df[[25]] <- list.df[[25]] %>% select(all_of(common.names))
co.list.df[[26]] <- list.df[[26]] %>% select(all_of(common.names))
co.list.df[[27]] <- list.df[[27]] %>% select(all_of(common.names))
co.list.df[[28]] <- list.df[[28]] %>% select(all_of(common.names))
co.list.df[[29]] <- list.df[[29]] %>% select(all_of(common.names))


ilr.co.listdf <- lapply(co.list.df,ilr)
ilr.df <- lapply(co.list.df,ilr)
t <- lapply(ilr.df,mean)
df <- lapply(t,ilrInv) 


df.o <- data.frame(matrix(unlist(df), nrow=length(df), byrow=TRUE))

colnames(df.o) <- colnames(co.list.df[[1]])

df.o$name <- names(list.df)
df.o
write.csv(x=df.o,file="/home/maxime/Documents/MA2/Thesis/Repository/
          data/meancomp.csv",row.names = F)
bind_rows()
