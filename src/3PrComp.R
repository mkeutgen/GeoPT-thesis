# Loading libraries

# Load libraries
#########################
## Principal Components.R #
#########################

library(tidyverse)
library(readr)
library(ggtern)
library(patchwork)
library(compositions)
library(ggfortify)
library(ggrepel)
set.seed(1)
# Start with most recent dataset

data <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/GeoPT48 -84Ra.csv")
# Geometric mean of this sample of compositions : 
geomean <- function(x){
  # compute the geometric mean of a vector
  exp(mean(log(x),na.rm=TRUE)) 
}


# Principal Component Analysis
sel <-c("SiO2","TiO2","Al2O3","Fe2O3T","MnO","MgO","CaO","Na2O","K2O","P2O5")
df.majors <- select(data,all_of(sel))
df.majors <- data.frame(clo(df.majors))

# replace missing values of each column by its geometric mean
geomean.v <- sapply(rbind(df.majors),geomean)
for (i in 1:nrow(df.majors)){
  df.majors[,i][is.na(df.majors[,i])] <- geomean.v[i] 
}

pairwise.log.df <- df.majors %>% transmute(
  "log(TiO2_SiO2)" = log(TiO2/SiO2),
  "log(Al2O3_TiO2)" = log(Al2O3/TiO2),
  "log(Fe2O3T_Al2O3)" = log(Fe2O3T/Al2O3),
  "log(MnO_Fe2O3T)" = log(MnO/Fe2O3T),
  "log(MgO_MnO)" = log(MgO/MnO),
  "log(CaO_MgO)" = log(CaO/MgO),
  "log(Na2O_CaO)" = log(Na2O/CaO),
  "log(K2O_Na2O)" = log(K2O/Na2O),
  "log(P2O5_K2O)" = log(P2O5/K2O)
)
pca.result <- prcomp(na.omit(pairwise.log.df),scale = T)
autoplot(pca.result,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)+theme_bw()

# w.o lower.right outlier
pairwise.log.df <- na.omit(pairwise.log.df)
which.max(pairwise.log.df$`log(MnO_Fe2O3T)`)
which.max(pairwise.log.df[-40,]$`log(MnO_Fe2O3T)`)
pairwise.log.df <- pairwise.log.df[-40,]
pairwise.log.df <- pairwise.log.df[-54,]
pca.result <- prcomp(pairwise.log.df,scale = T)
autoplot(pca.result,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)+theme_bw()

# Robust ?
library(robCompositions)
robpca <- pcaCoDa(df.majors,method = "robust")
plot(robpca)
biplot(robpca)

# w.o upper right outlier (extreme negative MgO/MnO concentration)
pca.result <- prcomp(na.omit(pairwise.log.df[-96,]),scale. = T)
autoplot(pca.result,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)+theme_bw()
# w.o upper right outlier (extreme MnO/Fe2O3 concentration)
df.majors <- df.majors %>% mutate(name = "monzonite")
# We import a second data frame GeoPT46, a grano-diorite. 
df2 <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/GeoPT46 -84Ra.csv")	
sel <-c("SiO2","TiO2","Al2O3","Fe2O3T","MnO","MgO","CaO","Na2O","K2O","P2O5")
df2.majors <- select(data,all_of(sel))
geomean.v <- sapply(rbind(df2.majors),geomean)
for (i in 1:ncol(df2.majors)){
  df2.majors[,i][is.na(df2.majors[,i])] <- geomean.v[i] 
}

df2.majors <- df2.majors %>% mutate(name="granite")

df.merged <- rbind(df.majors,df2.majors)
# ALR Transformation
# We divide all elements by SiO2 
vector <- df.majors$SiO2

mat <- as.matrix(df.merged[-ncol(df.merged)])
mat <- sweep(mat,MARGIN=1,FUN="/",STATS=vector)
mat <- log(mat)

log.df.major <- data.frame(mat)[-1]
name <- df.merged$name
log.df.major.merged <- cbind(log.df.major,name)

pca.result <- prcomp(log.df.major,scale = T)

autoplot(pca.result, data = log.df.major.merged,colour="name",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)+theme_bw()
summary(log.df.major)

# Outlying MnO/Fe2O3 observation
