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

for (i in 1:length(filePaths)){
  list.df[[i]] <- read_csv(file = filePaths[[i]])
}

names(list.df) <- fileNames

# Amalgamating and finding a common composition.


list.df.ilr <- list() 

tab <- table(unlist(lapply(list.df, names)))
common.names <- names(tab[tab == length(list.df)])

# Amalgamating elements based on Goldschmidt ? As a first step, we perform a very conservative amalgamation,
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
common.names %>% as_tibble() %>% View()
reduced.df %>% colnames()

reduced.df.l <- list()

for (i in 1:length(list.df)){
reduced.df.l[[i]] <- amalgam_fn(list.df[[i]])
  }

# ILR coordinates of reduced.df.l :
ilr.red.l <- list()
for (i in 1:length(list.df)){
  ilr.red.l[[i]] <- reduced.df.l[[i]] %>% ilr() %>% as_tibble()
}

# B. PRIOR TO TESTING EQUALITY OF COVARIANCE MATRICES, WE SHOULD FIRST CHECK HOW REASONABLE
# THE NORMALITY ASSUMPTION ACTUALLY IS. WE PERFORM SHAPIRO TEST ON THE MARGINALS
propH0rejnorm <- function(dataframe){
  count.H0rej <- ifelse(round(apply(dataframe,2,shap.pval),10)<=0.01,1,0) %>% sum()
propH0rej <- count.H0rej/ncol(dataframe)
return(propH0rej)
}
round(lapply(ilr.red.l,propH0rejnorm) %>% as_vector(),digits=3)

mu_sigma_fun <- function(dataframe){
  Xbar <- reduced.df %>% ilr() %>% colMeans()
  Sigma <- reduced.df %>% ilr() %>% var()
  output <- list(Xbar,Sigma)
  names(output) <- c("Xbar","VarCov Mat")
  return(output)
}
mu_sigma_fun(list.df[[2]])


mu_sigma_fun(list.df[[1]])[[2]] %>% eigen()

mu_sigma_fun(list.df[[2]])[[2]]
mu_sigma_fun(list.df[[3]])[[2]]
mu_sigma_fun(list.df[[4]])[[2]]
mu_sigma_fun(list.df[[5]])[[2]]
mu_sigma_fun(list.df[[6]])[[2]]



library(mvnormtest)

mshapiro.test(as.matrix(list.df[[1]]))
mshapiro.test(as.matrix(list.df[[2]]))
mshapiro.test(as.matrix(list.df[[3]]))
mshapiro.test(as.matrix(list.df[[4]]))
mshapiro.test(as.matrix(list.df[[5]]))
mshapiro.test(as.matrix(list.df[[6]]))
mshapiro.test(as.matrix(list.df[[7]]))
mshapiro.test(as.matrix(list.df[[8]]))
mshapiro.test(as.matrix(list.df[[9]]))
apply(list.df[[1]],2,kurtosis)
apply(list.df[[2]],2,kurtosis)


# Positive kurtosis everywhere.



apply(list.df[[1]],2,kurtosis)
apply(list.df[[2]],2,kurtosis)
apply(list.df[[3]],2,kurtosis)
apply(list.df[[4]],2,kurtosis)
apply(list.df[[5]],2,kurtosis)
apply(list.df[[6]],2,kurtosis)
apply(list.df[[7]],2,kurtosis)
apply(list.df[[8]],2,kurtosis)

str(list.df[[1]])


reduced.df.l[[29]] %>% ilr() %>% cov()

reduced.df.l[[28]] %>% ilr() %>% cov()



df29 <- reduced.df.l[[29]] %>% ilr() %>% as_tibble()
df28 <- reduced.df.l[[28]] %>% ilr() %>% as_tibble()

df29$group <- 29
df28$group <- 28
df28n29 <- rbind(df29,df28)
BoxM(df28n29[1:11],group = df28n29$group)


#### C. Testing for homogeneity of covariance matrices
library(covTestR)

irisSpecies <- unique(iris$Species)

iris_ls <- lapply(irisSpecies, 
                  function(x){as.matrix(iris[iris$Species == x, 1:4])}
)

names(iris_ls) <- irisSpecies
iris_ls
for (i in 1:length(ilr.red.l)){
  ilr.red.l[[i]] <- ilr.red.l[[i]] %>% as.matrix()
}
# So now we've a list of the ilr coordinates of the 29 dataframes. 
# There are a bunch of tests available within 
# https://rdrr.io/cran/covTestR/man/homogeneityStatistics.html
#Ahmad2017(x, ...)
#
#BoxesM(x, ...)
#
#Chaipitak2013(x, ...)
#
#Ishii2016(x, ...)
#
#Schott2001(x, ...)
#
#Schott2007(x, ...)
#
#Srivastava2007(x, ...)
#
#Srivastava2014(x, ...)
#
#SrivastavaYanagihara2010(x, ...)

# Which test to choose ???
#Ahmad, R. (2017). Location-invariant test of homogeneity of large-dimensional covariance matrices. Journal of Statistical Theory and Practice, 11(4):731-745. 10.1080/15598608.2017.1308895
#
#Chaipitak, S. and Chongcharoen, S. (2013). A test for testing the equality of two covariance matrices for high-dimensional data. Journal of Applied Sciences, 13(2):270-277. 10.3923/jas.2013.270.277
#
#Ishii, A., Yata, K., and Aoshima, M. (2016). Asymptotic properties of the first pricipal component and equality tests of covariance matrices in high-dimesion, low-sample-size context. Journal of Statistical Planning and Inference, 170:186-199. 10.1016/j.jspi.2015.10.007
#
#Schott, J (2001). Some Tests for the Equality of Covariance Matrices. Journal of Statistical Planniing and Inference. 94(1), 25-36. 10.1016/S0378-3758(00)00209-3
#
#Schott, J. (2007). A test for the equality of covariance matrices when the dimension is large relative to the sample sizes. Computational Statistics & Data Analysis, 51(12):6535-6542. 10.1016/j.csda.2007.03.004
#
#Srivastava, M. S. (2007). Testing the equality of two covariance matrices and independence of two sub-vectors with fewer observations than the dimension. InInternational Conference on Advances in InterdisciplinaryStistics and Combinatorics, University of North Carolina at Greensboro, NC, USA.
#
#Srivastava, M., Yanagihara, H., and Kubokawa T. (2014). Tests for covariance matrices in high dimension with less sample size. Journal of Multivariate Analysis, 130:289-309. 10.1016/j.jmva.2014.06.003
#
#Srivastava, M. and Yanagihara, H. (2010). Testing the equality of several covariance matrices with fewer observation that the dimension. Journal of Multivariate Analysis, 101(6):1319-1329. 10.1016/j.jmva.2009.12.010Ahmad2017(ilr.red.l[c(1,2)])

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

Ahmad2017(ilr.red.l,parameter=28) # p-value < 2.2e-16
Schott2001(ilr.red.l,parameter=28) # p-value < 2.2e-16
Schott2007(ilr.red.l,parameter=28) # p-value < 2.2e-16
Srivastava2007(ilr.red.l) # p-value = 1
Srivastava2014(ilr.red.l) # p-value < 2.2e-16
SrivastavaYanagihara2010(ilr.red.l) # p-value = 1
Ishii2016(ilr.red.l) # p-value = 0.05926
Chaipitak2013(ilr.red.l) # p-value < 2.2e-16

# SrivastavaYanagihara2010 not significant and Ishii2016 not significant neither and 
# Srivastava2007 not significant.

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



Ahmad2017(ilr.co.listdf) 
Schott2001(ilr.co.listdf)
Schott2007(ilr.co.listdf) 
Srivastava2007(ilr.co.listdf) 
Srivastava2014(ilr.co.listdf) 
SrivastavaYanagihara2010(ilr.co.listdf)
Ishii2016(ilr.co.listdf) 
Chaipitak2013(ilr.co.listdf) 

ilr.df <- lapply(co.list.df,ilr)
t <- lapply(ilr.df,mean)
df <- lapply(t,ilrInv) 
df

df.o <- data.frame(matrix(unlist(df), nrow=length(df), byrow=TRUE))
colnames(df.o) <- colnames(co.list.df[[1]])

df.o$name <- names(list.df)
write.csv(x=df.o,file="/home/maxime/Documents/MA2/Thesis/Repository/data/meancomp.csv")

            