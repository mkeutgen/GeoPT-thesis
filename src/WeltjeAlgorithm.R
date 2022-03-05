######################
## Weltje Algorithm ##
######################

library("tidyverse")
library("readr")
library("compositions")

# Vector of major oxides
ox.majors <- c("SiO2", "TiO2", "Al2O3", "Fe2O3T", "Fe(II)O", "MnO", "MgO", 
               "CaO", "Na2O", "K2O", "P2O5")


# Fraction of element in the major oxides
frac_el <- 1/100*c(46.75, 47.87, 52.93, 69.94, 77.74, 77.44, 60.31, 71.47, 74.18, 83.01, 43.64)

# Fraction of oxygen in the major oxides  
frac_ox <- rep(1,times=length(frac_el))-frac_el

# Logit and BLR Mean Function
logit <- function(x){log((x)/(1-x))}
# BLR Mean function
blr.mean <- function(x){
  # Map a vector from [0,1] to the real space. Takes the mean and return as its output the logit-inverse mean.
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}

# Cleaning function
df_cleaning <- function(dataframe,cutoff.major=1/2,cutoff.trace=3/4,cutoff.col=.95,carbon=FALSE,platinoids=TRUE){
  # Input : Dataframe, rock's composition.
  # 5 hyperparameters : "cutoff.major", laboratories with more than 1/2 (default) missing values in
  # the major elements are omitted. "cutoff.trace", laboratories with less than 3/4 (defaut) missing values
  # are omitted. 
  # "cutoff.col" analytes (chemical elements) with more than 95 % (default) missing values are omitted.   
  # "carbon" should the columns "Carbon(tot)" and "Carbon(org)" be included in the cleaned dataframe, BOOLEAN.
  # "platinoids" should the platinoids metals be included in the cleaned dataframe, BOOLEAN. 
  
  
  # Vector of major elements
  majors <- c("Si", "Ti", "Al", "Fe3", "Fe2", "Mn", "Mg", 
              "Ca", "Na", "K", "P")
  
  # Vector of trace elements
  traces.carbon <- c("Ag", "As", "Au", "B", "Ba", "Be", "Bi", "Br", "C(org)", "C(tot)", 
              "Cd", "Ce", "Cl", "Co", "Cr", "Cs", "Cu", "Dy", "Er", "Eu", "F", 
              "Ga", "Gd", "Ge", "Hf", "Hg", "Ho", "I", "In", "Ir", "La", "Li", 
              "Lu", "Mo", "N", "Nb", "Nd", "Ni", "Os", "Pb", "Pd", "Pr", "Pt", 
              "Rb", "Re", "Rh", "Ru", "S", "Sb", "Sc", "Se", "Sm", "Sn", "Sr", 
              "Ta", "Tb", "Te", "Th", "Tl", "Tm", "U", "V", "W", "Y", "Yb", 
              "Zn", "Zr")
  traces.wocarbon <- c("Ag", "As", "Au", "B", "Ba", "Be", "Bi", "Br", 
              "Cd", "Ce", "Cl", "Co", "Cr", "Cs", "Cu", "Dy", "Er", "Eu", "F", 
              "Ga", "Gd", "Ge", "Hf", "Hg", "Ho", "I", "In", "Ir", "La", "Li", 
              "Lu", "Mo", "N", "Nb", "Nd", "Ni", "Os", "Pb", "Pd", "Pr", "Pt", 
              "Rb", "Re", "Rh", "Ru", "S", "Sb", "Sc", "Se", "Sm", "Sn", "Sr", 
              "Ta", "Tb", "Te", "Th", "Tl", "Tm", "U", "V", "W", "Y", "Yb", 
              "Zn", "Zr")
 

traces <-    if(carbon==T){
     traces.carbon
      }  else {traces.wocarbon} 
  
 traces <- if (platinoids == F){
     c("Ag", "As","B", "Ba", "Be", "Bi", "Br", 
                         "Cd", "Ce", "Cl", "Co", "Cr", "Cs", "Cu", "Dy", "Er", "Eu", "F", 
                         "Ga", "Gd", "Ge", "Hf", "Hg", "Ho", "La", "Li", 
                         "Lu", "Mo", "Nb", "Nd", "Ni", "Pb", "Pr", 
                         "Rb", "S", "Sb", "Sc", "Se", "Sm", "Sn", "Sr", 
                         "Ta", "Tb", "Th", "Tl", "Tm", "U", "V", "W", "Y", "Yb", 
                         "Zn", "Zr")} else {traces.wocarbon}
  
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
  dataframe[traces] <- dataframe[traces] *1e-04
  # Sum to one
  dataframe <- 1/100*dataframe %>% select(ox.majors,traces)
  
  
  
  # Elementary cleaning, removing columns with no values in the major elements
  
  df.majors <- select(dataframe,all_of(ox.majors))
  df.traces <- select(dataframe,all_of(traces))
  
  # Cut off rows with more than 50 % of major analytes are dropped this step is 
  # important to avoid infinity values in the logit space (values too close to zero in the simplex space)
  m <- which(rowSums(is.na(df.majors)) < cutoff.major*ncol(df.majors))
  
  # Cut off rows with more than 75 % trace elements are dropped
  n <- which(rowSums(is.na(df.traces)) < cutoff.trace *ncol(df.traces))
  
  rowindex <- intersect(m,n)
  # Cut off columns with more than 95 % missing values are dropped
  # The Rest column, which variables are under the cutoff ?
  remainal <- which(colSums(is.na(dataframe)) >= cutoff.col*nrow(dataframe))
  # Compute the Rest column
  dataframe["Rest"] <- rowSums(dataframe[remainal],na.rm = T)
  
  colindex <- which(colSums(is.na(dataframe)) < cutoff.col*nrow(dataframe) )

  
  dataframe <- as_tibble(dataframe[rowindex,colindex])
  return(dataframe)
}

# Feasibility function
feasibility <- function(dataframe){


  
f.mat <- dataframe %>% summarize(X = rowSums(dataframe %>% select(ox.majors)*frac_el,na.rm = T),
                                  OX = rowSums(dataframe %>% select(ox.majors)*frac_ox,na.rm = T),
                                  Tr = rowSums(dataframe %>% select(!ox.majors),na.rm = T))
return(f.mat)
}
# Impute NA function
impute_na <- function(dataframe){
  blr.mean.l <- as.list(sapply(dataframe,blr.mean))
  names(blr.mean.l) <- names(dataframe)
  replace_na(dataframe,blr.mean.l)
}

# Scaling function 
scalingfunction <- function(dataframe,imputed.dataframe){
  # scale every rows whose sum is above 1
  feasibility.matrix <- feasibility(imputed.dataframe)  
  for (i in 1:nrow(dataframe)){
    # if statement check if composition is feasible, if not feasible, scale it
    dataframe[i,] <- if (rowSums(feasibility.matrix)[i]>1){
      sweep(dataframe[i,],2,rowSums(feasibility.matrix,na.rm=T)[i],"/")
    } # else no need to scale it 
    else {
      dataframe[i,]
    }
    
  }
  return(dataframe)
}

# Iterative mean function
iterative.mean <- function(dataframe,epsilon=1E-8){
list.df <- list()
ADS.list <- list()

list.df[[1]] <- dataframe
list.df[[2]] <- scalingfunction(dataframe,impute_na(dataframe))
i <- 2
epsilon <- 10E-8
# Do while the difference in mean is above a "physical" threshold, argue why this is better than conventional MSE.
while (abs(sum(rowSums(list.df[[i-1]],na.rm=T))-sum(rowSums(list.df[[i]],na.rm=T))) > epsilon){
  list.df[[i+1]] <- scalingfunction(list.df[[i]],impute_na(list.df[[i]]))
  ADS.list[i] <- abs(sum(rowSums(list.df[[i-1]],na.rm=T))-sum(rowSums(list.df[[i]],na.rm=T)))
  i <- i + 1
}
ADS <- abs(mean(rowSums(list.df[[length(list.df)]],na.rm=T))-mean(rowSums(list.df[[length(list.df)-1]],na.rm=T)))
# Impute NA of the dataframe with converged values :
dataframe <- impute_na(list.df[[length(list.df)]])
# Final Output
output <- scalingfunction(dataframe,dataframe)
# Rest column 
final.output <- output %>% mutate(Rest=1-rowSums(output)) 

result <- list(final.output,paste(i,"iterations needeed until convergence"),paste("At final iteration, ADS is",ADS),ADS.list)
names(result) <- c("Output","Iterations","ADS","ADS.list")
return(result)
}

# GW Algorithm, takes a messy dataframe and do :
# 1 Cleaning, drop arbitrary rows & columns.
# 2. Iterative Mean Estimation until convergence. 
GWAlgo <- function(dataframe,carbon=FALSE,platinoids=TRUE){
# 2 Clean the data
dataframe <- df_cleaning(dataframe,platinoids=TRUE)
# 3 Call the iterative mean function
result <- iterative.mean(dataframe)[[1]] 
return(result)
}

# More restrictive version of GW Algorithm, where platinoids and carbon columns are not included.
# This is required for early files where these chemical elements were often not analyzed. This function should
# be considered as obsolete for the latest analysis from the GeoPT Project.
GWAlgo.mod <- function(dataframe,carbon=FALSE,platinoids=FALSE){
  # 2 Clean the data
  dataframe <- df_cleaning(dataframe,platinoids=FALSE)
  # 3 Call the iterative mean function
  result <- iterative.mean(dataframe)[[1]]
  return(result)
}

# # ## Testing the function
# # # Let a raw dataframe.
setwd("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/")
# # 
dataframe <- read_csv(file = "GeoPT48 -84Ra.csv")
dataframe

# # 
df <- iterative.mean(df_cleaning(dataframe,carbon = F,platinoids = F))
conv.data <- cbind(seq(1:length(unlist(df$ADS.list))),unlist(df$ADS.list)) %>% as_tibble()
# Convergence rate for GeoPT48 
output <- ggplot(conv.data,aes(x=V1,y=V2))+geom_line() + scale_y_log10()+theme_bw()+labs(y="ADS",x="Iterations") # Quadratic !
output

names(list.df[25])

convergenceplot(list.df[[26]])
