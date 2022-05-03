# Geometry of the sample space. The need for ILR transformation.
library(compositions)
library(tidyverse)
library(ggtern)
# Let a compositional dataset hydrochem :
set.seed(123)
data("Hydrochem")
data.tibble <- Hydrochem[sample(nrow(Hydrochem),10),] %>% select(c(Na,Ca,Mg)) %>% clo() %>% as_tibble()
Hydrochem
x <- data[1,]
y <- data[2,]
data.comp <- acomp(data)

data.tibble[nrow(data.tibble)+1,] 
# Plot a ternary diagram with the sum of the two vectors in the simplex

data.tibble[nrow(data.tibble)+1,] <- perturbe(data.comp[1],data.comp[2]) %>% as.numeric() %>% t()
result <- data.tibble[c(1,2,11),]
result
ggtern(result,mapping = aes(x=Na,y=Ca,z=Mg))+geom_point()

ilr.result <- result %>% ilr() %>% as_tibble()
ilr.result
dist(ilr.result)
dist(result)
compo <- c(0.2,0.4,0.4)
geomean <- geometricmean(compo)
D <- length(compo)
j <- 1
sqrt((D-j)/(D-j+1))*log(compo[1]/sqrt(compo[-1]))
ilr(compo)

ggplot(ilr.result,aes(x=V1,y=V2))+geom_point()
