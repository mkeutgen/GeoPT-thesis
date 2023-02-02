#########################
## Motivating Example.R #
#########################


# Load libraries

library(tidyverse)
library(readr)
library(ggtern)
library(patchwork)
library(compositions)
library(ggfortify)
library(ggrepel)
set.seed(1)
# Start with most recent dataset

data <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/
                 GeoPT48 -84Ra.csv")
# Geometric mean of this sample of compositions : 
geomean <- function(x){
  # compute the geometric mean of a vector
  exp(mean(log(x),na.rm=TRUE)) 
}


three.parts <- c("Cr","Ce","Ba")
three.parts.df <- select(data,all_of(three.parts))
three.parts.df <- data.frame(clo(three.parts.df))
# Close it
colMeans(three.parts.df,na.rm = T)
three.parts.df <- data.frame(clo(three.parts.df))
geom <- sapply(rbind(three.parts.df),geomean)
arithm <- colMeans(three.parts.df,na.rm = T)
arithm <- clo(arithm)
geom <- clo(geom)
print(rbind(round(geom,3),round(arithm,3)))

ggtern(data=three.parts.df,aes(Cr,Ce,Ba)) + 
  geom_mask() +
  geom_point(fill="black",alpha=.5) + 
  geom_point(aes(arithm[1],arithm[2],arithm[3],color="arithmetic mean"))+
  geom_point(aes(geom[1],geom[2],geom[3],color="geometric mean"),alpha=1)+
  theme_bw()
 
three.parts.df %>% ilr() %>% colMeans() %>% ilrInv()
sapply(three.parts.df,geometricmean)
sum(geom)
#################################
######## Missing value problem ##
#################################
three.parts <- c("Lu","Ce","Ba")
three.parts.df <- select(data,all_of(three.parts))
summary(three.parts.df)

ggtern(data=three.parts.df,aes(Lu,Ce,Ba)) + 
  geom_mask() +
  geom_point(fill="black",alpha=.5) + 
  geom_point(aes(arithm[1],arithm[2],arithm[3],color="arithmetic mean"))+
  geom_point(aes(geom[1],geom[2],geom[3],color="geometric mean"),alpha=1)+
  theme_bw()
geom <- sapply(rbind(three.parts.df),geomean)
arithm <- colMeans(rbind(three.parts.df),na.rm = T)
arithm

three.parts.df.cena <- three.parts.df[!is.na(three.parts.df$Ce),]
three.parts.df.bana <- three.parts.df[!is.na(three.parts.df$Ba),]
# now df has only Lu N.A : 
geom <- sapply(rbind(three.parts.df.bana),geomean)
arithm <- colMeans(rbind(three.parts.bana),na.rm = T)
geom
arithm
ggtern(data=three.parts.df.bana,aes(Lu,Ce,Ba)) + 
  geom_mask() +
  geom_point(fill="black",alpha=.5) + 
  geom_point(aes(arithm[1],arithm[2],arithm[3],color="arithmetic mean"))+
  geom_point(aes(geom[1],geom[2],geom[3],color="geometric mean"),alpha=1)+
  theme_bw()

# replacing with detection limit and removing outlier

df.outlierfree <- three.parts.df.bana[-16,]
df.outlierfree[is.na(df.outlierfree$Lu),] <- t(c(0.001,geom[2],geom[3]))
geom <- sapply(rbind(df.outlierfree),geomean)
arithm <- colMeans(rbind(df.outlierfree),na.rm = T)
geom
arithm

ggtern(data=df.outlierfree,aes(Lu,Ce,Ba)) + 
  geom_mask() +
  geom_point(fill="black",alpha=.5) + 
  geom_point(aes(arithm[1],arithm[2],arithm[3],color="arithmetic mean"))+
  geom_point(aes(geom[1],geom[2],geom[3],color="geometric mean"),alpha=1)+
  theme_bw()

# 4 parts composition
sel <- c("SiO2","Na2O","K2O","CaO")
df.4 <- select(data,all_of(sel))
df.4.cl <- data.frame(clo(df.4))
df4.plot <- ggplot(df.4.cl,aes(x=Na2O,y=K2O))+
  geom_point(color="midnightblue",size = 2)+theme_bw()+
  geom_smooth(method="lm",size=1,se=FALSE,color="mediumpurple2")



# 3 parts composition
sel <- c("Na2O","CaO","K2O")
df.3 <- select(data,all_of(sel))
df.3.cl <- data.frame(clo(df.3))
df3.plot <- ggplot(df.3.cl,aes(x=Na2O,y=K2O))+
  geom_point(color="midnightblue",size=2)+
  theme_bw()+geom_smooth(method="lm",size=1,se = FALSE,
                         color="mediumpurple2")

df3.plot + df4.plot

# logratio : 
# 4 parts composition
df.4.cl <- df.4.cl %>% mutate(
  log.Na2O = log(Na2O/CaO),
  log.K2O = log(K2O/CaO)
)
df.3.cl <- df.3.cl %>% mutate(
  log.Na2O = log(Na2O/CaO),
  log.K2O = log(K2O/CaO)
)
logratio4plot <- ggplot(df.4.cl,aes(x=log.Na2O,y=log.K2O))+
  geom_point(color="blue",size=2)+
  labs(y="log(K2O/CaO)",x="log(Na2O/CaO)")+theme_bw()

logratio3plot <- ggplot(df.3.cl,aes(x=log.Na2O,y=log.K2O))+
  geom_point(color="blue",size=2)+
  labs(y="log(K2O/CaO)",x="log(Na2O/CaO)")+theme_bw()

logratio3plot + logratio4plot
# 3 parts composition
sel <- c("Na2O","K2O","CaO")
df.3 <- select(data,all_of(sel))
df.3.cl <- data.frame(clo(df.3))
ggplot(df.3.cl,aes(x=Na2O,y=CaO))+geom_point()

