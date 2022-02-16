library("tidyverse")
library("readr")
library("compositions")

dataframe <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/GeoPT48 -84Ra.csv")


blr_imputation <- function(dataframe,cutoff.major=1/2,cutoff.trace=3/4,cutoff.col=0.95){
## BLR IMPUTATION FUNCTION ##
# Blr Imputation Function replace iteratively missing values with the blr mean
# Input : Dataframe, rock's composition.
# 3 hyperparameters : "cutoff.major", laboratories with more than 1/2 (default) missing values in
# the major elements are omitted. "cutoff.trace", laboratories with less than 3/4 (defaut) missing values
# are omitted. 
# "cutoff.col" analytes (chemical elements) with more than 95 % (default) missing values are omitted.   

## PART ONE : Removing NA's 

# Vector of major elements
majors <- c("Si", "Ti", "Al", "Fe3", "Fe2", "Mn", "Mg", 
            "Ca", "Na", "K", "P")

# Vector of trace elements
traces <- c("Ag", "As", "Au", "B", "Ba", "Be", "Bi", "Br", 
            "Cd", "Ce", "Cl", "Co", "Cr", "Cs", "Cu", "Dy", "Er", "Eu", "F", 
            "Ga", "Gd", "Ge", "Hf", "Hg", "Ho","In", "Ir", "La", "Li", 
            "Lu", "Mo", "Nb", "Nd", "Ni", "Os", "Pb", "Pd", "Pr", "Pt", 
            "Rb", "Re", "Rh", "Ru", "S", "Sb", "Sc", "Se", "Sm", "Sn", "Sr", 
            "Ta", "Tb", "Te", "Th", "Tl", "Tm", "U", "V", "W", "Y", "Yb", 
            "Zn", "Zr")





# Vector of major oxides
ox.majors <- c("SiO2", "TiO2", "Al2O3", "Fe2O3T", "Fe(II)O", "MnO", "MgO", 
               "CaO", "Na2O", "K2O", "P2O5")


# Fraction of element in the major oxides
frac_el <- 1/100*c(46.75, 47.87, 52.93, 69.94, 77.74, 77.44, 60.31, 71.47, 74.18, 83.01, 43.64)

# Fraction of oxygen in the major oxides  
frac_ox <- rep(1,times=length(frac_el))-frac_el

# Some dataframes have * instead of NA
dataframe[dataframe=="*"]=NA


# Some dataframes don't have numeric columns yet 
dataframe[c(traces,ox.majors)] <- lapply(dataframe[c(traces,ox.majors)], function(x) {
  if(is.character(x)) as.numeric(as.character(x)) else x
})


# Unit Conformity 
dataframe[traces] <- dataframe[traces] *10E-4
# Sum to one
dataframe <- 1/100*dataframe %>% select(ox.majors,traces)



# Elementary cleaning, removing columns with no values in the major elements

df.majors <- select(dataframe,all_of(ox.majors))
df.traces <- select(dataframe,all_of(traces))

# Cut off rows with more than 50 % of major analytes are dropped this step is 
# important to avoid infinity values in the logit space
m <- which(rowSums(is.na(df.majors)) < cutoff.major*ncol(df.majors))

# Cut off rows with more than 75 % trace elements are dropped
n <- which(rowSums(is.na(df.traces)) < cutoff.trace *ncol(df.traces))

rowindex <- intersect(m,n)
colindex <- which(colSums(is.na(dataframe)) < cutoff.col*nrow(dataframe))



dataframe <- dataframe[rowindex,colindex]



# Part 2 : Compute feasibility matrix
## Feasible Function 
feasible <- function(data.frame){
  # Feasibility matrix
  
  f.mat <- data.frame %>% summarize(X = rowSums(data.frame %>% select(ox.majors)*frac_el,na.rm = T),
                                    OX = rowSums(data.frame %>% select(ox.majors)*frac_ox,na.rm = T),
                                    Tr = rowSums(data.frame %>% select(intersect(names(colindex),traces)),na.rm = T))
  # Scaling this matrix to unit sum
  f.mat <- as_tibble(clo(f.mat))
  f.mat <- f.mat %>% mutate(R = 1-X-OX-Tr,
                            blrX = log(X/(1-X)),
                            blrTr = log(Tr/(1-Tr)),
                            blrOX = log(OX/(1-OX))
                            
  )
  f.mat.final <- f.mat %>% select(X,OX,Tr,R)
  
  blrX.m <- mean(f.mat$blrX)
  varX.m <- var(f.mat$blrX,na.rm = T)
  blrTr.m <- mean(f.mat$blrTr)
  blrOx.m <- mean(f.mat$blrOX)
  
  Inv.X <- (exp(blrX.m) )/ ( 1 + exp(blrX.m))
  Inv.Tr <- (exp(blrTr.m) )/ ( 1 + exp(blrTr.m))
  Inv.Ox <- (exp(blrOx.m) )/ ( 1 + exp(blrOx.m))
  # R is :
  Inv.R <- 100-(Inv.X + Inv.Tr + Inv.Ox)
  result <- list(Inv.X,Inv.Tr,Inv.R,f.mat.final)
  names(result) <- c("Estimate X","Estimate Tr","Estimate R","Feasibility Matrix")
  return(result)
  
}
# Next we impute missing values by their means after logit transform. 
# Convergence criterion : Euclidean distance of successive estimates X' and T'


# Replace iteratively NAs values row-wise

#for i in 1:ncol(df)

logit <- function(x){log((x)/(1-x))}

blr.mean <- function(x){
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}
df <- dataframe

blr.mean.l <- as.list(sapply(df,blr.mean))
names(blr.mean.l) <- names(df)

# imputing NA with the mean calculated 
# Replace Na's by mean of the logit transformed observed values column after column :


list.df <- list()
length(list.df) <- nrow(df)
list.row <- list()
length(list.row) <- nrow(df)
blr.mean.list <- list()
length(blr.mean.list) <- nrow(df)


# Initialize it 

for (i in 2:(nrow(df))) {
  list.df[[1]] <- df
  blr.mean.list[[i]] <- as.list(sapply(list.df[[i-1]],blr.mean))
  names(blr.mean.list[[i]]) <- names(df)
  list.row[[i]] <- replace_na(df[1:i,],blr.mean.list[[i]]) 
  list.df[[i]] <- as_tibble(rbind(list.row[[i]],df[(i+1):nrow(df),]))
}

# Does replacing the i-th rows affect the feasibility matrix ?

feas.mat <- sapply(list.df,feasible)
feas.mat.full <- feas.mat[,nrow(df)][[4]]
# Preliminary conclusions :
# X at 1st iteration 0.6134728
# X at 30th iteration 0.611168
# X at 58th iteration 
# 

full.df <- list.df[[nrow(df)]]
full.df <- full.df[1:nrow(df),]
full.df <- as_tibble(clo(full.df))

result <- list(full.df,feas.mat.full,feas.mat)
names(result) <- c("Imputed Dataframe","Feasibility Matrix full df","Estimated Values")
return(result)
}

# TEST THE FUNCTION
dataset <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/GeoPT48 -84Ra.csv")
dataset_full <- blr_imputation(dataset)
head(dataset_full)
