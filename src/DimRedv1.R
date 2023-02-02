# DimRed version 2

# Loading libraries
library(readr)
library(tidyverse)
library(compositions)
library(MVTests)
library(moments)
library(mvnormtest)
library(ggfortify)
library(ggpubr)

# Loading datasets and data cleaning
meancomp <- read_csv("~/Documents/MA2/Thesis/Repository/data/meancomp.csv")
meancomp.df <- meancomp %>% select(name, everything())
rownames(meancomp) <- NULL
descri <- read_csv("~/Documents/MA2/Thesis/Repository/data/raw/description.csv")
descri$Classification

descri$ID <- paste(descri$ID,".csv",sep = "") 
descri$name <- descri$ID
meancomp
ordered.names <- c("GeoPT1.csv","GeoPT2.csv","GeoPT3.csv","GeoPT4.csv"
                   ,"GeoPT5.csv","GeoPT6.csv",
                   "GeoPT8.csv","GeoPT14.csv","GeoPT16.csv",
                   "GeoPT19.csv","GeoPT20.csv",
                   "GeoPT21.csv","GeoPT22.csv","GeoPT23.csv",
                   "GeoPT25.csv", "GeoPT29.csv",
                   "GeoPT32.csv","GeoPT34.csv","GeoPT35.csv",
                   "GeoPT36.csv","GeoPT37.csv","GeoPT38.csv",
                   "GeoPT38A.csv","GeoPT39.csv","GeoPT39A.csv",
                   "GeoPT41.csv","GeoPT43.csv",
                   "GeoPT46.csv","GeoPT48.csv")

df <- merge(meancomp,descri) 
# Order the df
df <- df %>%
  slice(match(ordered.names, ID))


meancomp.mat <- meancomp[-1] %>% clr() 
var(meancomp.mat)

df$acidbase <- ifelse(df$SiO2>.66,"Acidic",
                      ifelse(df$SiO2>.52,"Intermediate",
                             ifelse(df$SiO2>.45,"Basic","Ultrabasic")))

clr.df <- df[2:57] %>% clr()


pca <- prcomp(clr.df,scale. = F,center = T)

df$ID <- factor(df$ID,levels=ordered.names)
loadplot <- autoplot(pca, loadings = T, loadings.label = T,
         data = df,label=F,size=1)+theme_bw()

autoplot(pca, loadings = F, loadings.label = F,
         data = df, colour = 'ID',label=TRUE,size=0.0001)+theme_bw()
pointplot <- autoplot(pca, loadings = F, loadings.label = F,
         data = df, colour = 'acidbase')+theme_bw()+theme(legend.position = "bottom")

pointplot
amal.class
notamal <- ggarrange(loadplot,pointplot,nrow = 2)
ggsave("notamal.png",notamal,width=7.5,height=15)

# AMALGAMATIONS 
lile.v <- c("Sr","Rb","Ba")
hfse.v <- c("Ta","Nb","Ce","Hf","Sm","Y")
lree.v <- c("La","Ce","Pr","Nd","Sm")
hree.v <- c("Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu")
#chalco <- c("Sn","As","Pb","Zn","Ag")

meancomp.amal <- meancomp.df %>% mutate(lile = Sr+Rb+Ba,
                                        hfse = Ta+Nb+Ce+Hf+Y,
                                        lree = La+Ce+Pr+Nd+Sm,
                                        hree = Gd+Tb+Dy+Ho+Er+Tm+Yb+Lu)

amal.meancomp <- meancomp.amal %>% select(-lile.v,-hfse.v,-lree.v,-hree.v) 
amal.meancomp.mat <- amal.meancomp[-1] %>% clo() %>% clr() %>% as_tibble()
amal.pca <- prcomp(amal.meancomp.mat,scale=F,retx = T)
dim(amal.meancomp.mat)
# First two PC -> 63.4 % of total variance 
amal.pca %>% summary()
amal.load <- autoplot(amal.pca, loadings = T, loadings.label = T,
         data = df)+theme_bw()+scale_x_reverse()
amal.class <- autoplot(amal.pca, loadings = F, loadings.label = F,
         data = df, colour = 'Classification')+theme_bw()+
  scale_x_reverse()+theme(legend.position = "bottom")
amal <- ggarrange(amal.load,amal.class,nrow = 2)
ggsave("amal.png",amal,width=7.5,height=15)


autoplot(amal.pca, loadings = F, loadings.label = F,
         data = df, colour = 'ID',label=TRUE,size=0.0001)+theme_bw()+
  scale_x_reverse()

parts10 <- amalgam::amalgam(clr.df,10)
weights10 <- t(parts10$weights)
weights10 <- apply(weights10, 1, \(x) names(x)[as.logical(x)])

parts11 <- amalgam::amalgam(clr.df,11)
weights11 <- t(parts11$weights)
weights11 <- apply(weights11, 1, \(x) names(x)[as.logical(x)])


weight.fn <- function(n){
parts12 <- amalgam::amalgam(clr.df,n)
weights12 <- t(parts12$weights)
weights12 <- apply(weights12, 1, \(x) names(x)[as.logical(x)])
return(weights12)}
weights2 <- weight.fn(2)
weights3 <- weight.fn(3)
weights4 <- weight.fn(4)
weights5 <- weight.fn(5)
weights6 <- weight.fn(6)
weights7 <- weight.fn(7)
weights8 <- weight.fn(8)
weights9 <- weight.fn(9)
weights10 <- weight.fn(10)

weights13 <- weight.fn(13)
weights14 <- weight.fn(14)
weights15 <- weight.fn(15)
weights16 <- weight.fn(16)
weights17 <- weight.fn(17)
weights18 <- weight.fn(18)
weights19 <- weight.fn(19)
weights20 <- weight.fn(20)

weights50 <- weight.fn(50)
weights51 <- weight.fn(51)
weights52 <- weight.fn(52)
weights53 <- weight.fn(53)



ilrBase(D=5)
