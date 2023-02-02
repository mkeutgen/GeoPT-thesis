library(readr)
library(tidyverse)
library(compositions)
library(MVTests)
library(moments)
library(mvnormtest)
library(ggfortify)


# Importing the meancomp dataset
meancomp <- read_csv("~/Documents/MA2/Thesis/Repository/data/meancomp.csv")
meancomp <- meancomp %>% select(name, everything())
rownames(meancomp) <- NULL


descri <- read_csv("~/Documents/MA2/Thesis/Repository/data/raw/description.csv")
descri$Classification
descri$ID <- paste(descri$ID,".csv",sep = "") 
descri$name <- descri$ID
meancomp
descri

# COME ON 
ordered.names <- c("GeoPT1.csv","GeoPT2.csv","GeoPT3.csv","GeoPT4.csv","GeoPT5.csv","GeoPT6.csv",
                   "GeoPT8.csv","GeoPT14.csv","GeoPT16.csv","GeoPT19.csv","GeoPT20.csv",
                   "GeoPT21.csv","GeoPT22.csv","GeoPT23.csv","GeoPT25.csv", "GeoPT29.csv",
                   "GeoPT32.csv","GeoPT34.csv","GeoPT35.csv","GeoPT36.csv","GeoPT37.csv","GeoPT38.csv",
                   "GeoPT38A.csv","GeoPT39.csv","GeoPT39A.csv","GeoPT41.csv","GeoPT43.csv",
                   "GeoPT46.csv","GeoPT48.csv")

df <- merge(meancomp,descri) 
# Order the df
df <- df %>%
  slice(match(ordered.names, ID))

meancomp.mat <- meancomp[-1] %>% clr() 

df$acidbase <- ifelse(df$SiO2>.66,"Acidic",ifelse(df$SiO2>.52,"Intermediate",ifelse(df$SiO2>.45,"Basic","Ultrabasic")))

clr.df <- df[2:57] %>% clr()


pca <- prcomp(clr.df,scale. = F,center = T)

df$ID <- factor(df$ID,levels=ordered.names)
autoplot(pca, loadings = T, loadings.label = T,
         data = df, colour = 'ID',label=TRUE,size=0.0001)+theme_bw()

autoplot(pca, loadings = F, loadings.label = F,
         data = df, colour = 'ID',label=TRUE,size=0.0001)+theme_bw()
autoplot(pca, loadings = F, loadings.label = F,
         data = df, colour = 'acidbase')+theme_bw()

biplot(pca,choices = c(2:3))
# 18, 13, 26


df.wo.3 <- df[-c(18,13,26)]
colnames(df.wo.3)
df.wo.3.clr <- df.wo.3[2:55] %>% clr()
pca.mod <- prcomp(df.wo.3.clr,scale=T,center = T)
autoplot(pca.mod, loadings = T, loadings.label = T,
         data = df.wo.3, colour = 'Classification')+theme_bw()





df %>% filter(ID %in% ordered.names)

descri %>% select(ID,Classification)
descri

dftot$column_label

dftot <- read_csv("~/Documents/MA2/Thesis/Repository/data/totaldataset.csv")
dftot$column_label
ilr.df <- dftot[c(-1,-2)] %>% ilr()

pca <- prcomp(ilr.df,scale. = T)

autoplot(pca, loadings = T, loadings.label = F,
         data = dftot, colour = 'column_label')

merge(iris,df)

dftot$ID <- dftot$column_label

names(dftot)
names(descri)
descri.df <- descri %>% select(ID,Classification)
output <- merge(descri.df,dftot)
dftot
names(descri.df) <- c("column_label","classification") 
output <- merge(descri.df,dftot)
# OUTPUT 
output %>% head()

clr.df <- output[c(-1,-2,-3)] %>% clr()



#pca <- prcomp(clr.df,scale. = T)
#
#autoplot(pca, loadings = F, loadings.label = T,
#         data = output, colour = 'classification')+theme_bw()
#
#unique(output$classification)
#
#Mafic.df        <- output %>% filter(classification=="Mafic")
#Felsic.df       <- output %>% filter(classification=="Felsic")
#Ultramafic.df   <- output %>% filter(classification=="Ultramafic")
#Intermediate.df <- output %>% filter(classification=="Intermediate")
#
#clr.Mafic.df       <- clr(Mafic.df[c(-1,-2,-3)]) %>% as_tibble()      
#clr.Felsic.df      <- clr(Felsic.df[c(-1,-2,-3)]) %>% as_tibble()
#clr.Ultramafic.df  <- clr(Ultramafic.df[c(-1,-2,-3)]) %>% as_tibble()
#clr.Intermediate.df<- clr(Intermediate.df[c(-1,-2,-3)]) %>% as_tibble()
#
#pcaMafic        <- prcomp(clr.Mafic.df,scale. = T)              
#pcaFelsic       <- prcomp(clr.Felsic.df,scale. = T)     
#pcaUltraMaf     <- prcomp(clr.Ultramafic.df,scale. = T)  
#pcaIntermediate <- prcomp(clr.Intermediate.df,scale. = T)
#
#pcaMafic       %>% biplot()  
#pcaFelsic      %>% biplot()
#pcaUltraMaf    %>% biplot()
#pcaIntermediate%>% biplot()
#
# https://www.r-bloggers.com/2021/05/linear-discriminant-analysis-in-r/
clr.meancomp <- meancomp[-57] %>% clr() %>% as_tibble()
meancomp$name
descri$Classification
meancomp$name
meancomptot <- merge(meancomp,descri)
meancomptot
pca <- prcomp(clr.meancomp,center = T,scale. = T)
summary(pca)
autoplot(pca, loadings = F, loadings.label = F,
         data = meancomptot, colour = 'Classification')

View(descri)
round(eigen(var(clr.meancomp))$values/cumsum(eigen(var(clr.meancomp))$values),digits = 3)
summary(pca)
cor(clr.meancomp)

# compute total variance
variance = pca$sdev^2 / sum(pca$sdev^2)

# Scree plot
qplot(c(1:29), variance) +
  geom_line() +
  geom_point(size=1)+
  xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 0.5)+theme_bw()


### AMALGAMATIONS
library(amalgam)
meancomp.df <- meancomp[2:56]
amalgam3 <- amalgam(meancomp.df,n.amalgams = 3,objective=objective.keepDist,
                    # # objective = objective.keepWADIST # another distance
                    # objective = objective.keepSKL # another distance
                    # objective = objective.maxRDA, # if maximizing RDA
                    #z = iris[,5], # only needed if maximizing RDA
                    #asSLR = FALSE, # if TRUE, n.amalgams must be even
                    #shrink = FALSE) # toggles James-Stein type shrinkage
# visualize results)
)
amalgam3 <- amalgam(meancomp.df,n.amalgams = 3,objective=objective.keepDist)
t(amalgam3$weights)
amalgam4 <- amalgam(meancomp.df,n.amalgams = 4,objective=objective.keepDist)
amalgam5 <- amalgam(meancomp.df,n.amalgams = 5,objective=objective.keepDist)
amalgam6 <- amalgam(meancomp.df,n.amalgams = 6,objective=objective.keepDist)
amalgam7 <- amalgam(meancomp.df,n.amalgams = 7,objective=objective.keepDist)
amalgam8 <- amalgam(meancomp.df,n.amalgams = 8,objective=objective.keepDist)
amalgam9 <- amalgam(meancomp.df,n.amalgams = 9,objective=objective.keepDist)
amalgam10 <- amalgam(meancomp.df,n.amalgams = 10,objective=objective.keepDist)

amal <- function(df,x=10){
one <- amalgam(df,n=x)
out <- apply(t(one$weights), 1, \(x) names(x)[as.logical(x)])
return(out)
}

groups5 <- amal(meancomp.df,x=5)
groups6 <- amal(meancomp.df,x=6)
groups7 <- amal(meancomp.df,x=7)
groups8 <- amal(meancomp.df,x=8)
groups9 <- amal(meancomp.df,x=9)


groups6
groups5

tengroups <- amal(meancomp.df,n=10)
twentygroups <- amal(meancomp.df,n=20)
fortygroups <- amal(meancomp.df,n=40)
fortygroups


amalgam3 <- amalgam(meancomp.df,n.amalgams = 3,objective=objective.keepSKL)

View(meancomp)


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
amal.pca <- prcomp(amal.meancomp.mat,scale=T,retx = T)

autoplot(amal.pca, loadings = T, loadings.label = T,
         data = df, colour = 'Classification')+theme_bw()+scale_x_reverse()


autoplot(amal.pca, loadings = F, loadings.label = F,
         data = df, colour = 'ID',label=TRUE,size=0.0001)+theme_bw()+scale_x_reverse()

# anomaly of GeoPT39A, ultrabasic rock yet surronded by acidic rock ???


# AMALGAM MAX
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
amal.pca <- prcomp(amal.meancomp.mat,scale=T,retx = T)

autoplot(amal.pca, loadings = T, loadings.label = T,
         data = df, colour = 'acidbase')+theme_bw()+scale_x_reverse()

