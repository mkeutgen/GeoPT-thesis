library("tidyverse")
library("readr")
library("compositions")
# Loading the dataset to test the function

setwd("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/")
dataframe <- read_csv("GeoPT48 -84Ra.csv")

# Logit and BLR Mean Function
logit <- function(x){log((x)/(1-x))}

blr.mean <- function(x){
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}



# Cleaning function
df_cleaning <- function(dataframe,cutoff.major=1/2,cutoff.trace=3/4,cutoff.col=.95){
  # Vector of major elements
  majors <- c("Si", "Ti", "Al", "Fe3", "Fe2", "Mn", "Mg", 
              "Ca", "Na", "K", "P")
  
  # Vector of trace elements
  traces <- c("Ag", "As", "Au", "B", "Ba", "Be", "Bi", "Br", "C(org)", "C(tot)", 
              "Cd", "Ce", "Cl", "Co", "Cr", "Cs", "Cu", "Dy", "Er", "Eu", "F", 
              "Ga", "Gd", "Ge", "Hf", "Hg", "Ho", "I", "In", "Ir", "La", "Li", 
              "Lu", "Mo", "N", "Nb", "Nd", "Ni", "Os", "Pb", "Pd", "Pr", "Pt", 
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


# Apply data cleaning function on dataframe
dataframe <- df_cleaning(dataframe)
View(dataframe %>% select(!ox.majors))
# Now that we cleaned the data and deleted absurd laboratories and columns, we move on to compute the feasibility for each row.
feasibility <- function(dataframe){
f.mat <- dataframe %>% summarize(X = rowSums(dataframe %>% select(ox.majors)*frac_el,na.rm = T),
                                  OX = rowSums(dataframe %>% select(ox.majors)*frac_ox,na.rm = T),
                                  Tr = rowSums(dataframe %>% select(!ox.majors),na.rm = T))
return(f.mat)
}

rowSums(feasibility(dataframe))
ifelse(rowSums(f.mat)[1]>1,dataframe[1])

# ITERATION ONE :
# If row-wise sum of the feasibility matrix is above one, scale the row, else do nothing


# scale every rows whose sum is above 1
feasibility.matrix <- feasibility(dataframe)  
for (i in 1:nrow(dataframe)){
  # if statement check if composition is feasible, if not feasible, scale it
  dataframe[i,] <- if (rowSums(feasibility.matrix)[i]>1){
    sweep(dataframe[i,],2,rowSums(dataframe,na.rm=T)[i],"/")
  } # else no need to scale it 
  else {
    dataframe[i,]
  }
}
return(dataframe)
}
# Scale observations: 
dataframe <- scalingfunction(dataframe = dataframe)  
# Impute missing values
dataframe.imp <- 
rowSums(dataframe,na.rm=T)
# Then scale the original dataframe :
dataframe/rowSums(dataframe.prime,na.rm=T)[i]
# Then compute the missing values 

impute_na <- function(dataframe){
  blr.mean.l <- as.list(sapply(dataframe,blr.mean))
  names(blr.mean.l) <- names(dataframe)
  replace_na(dataframe,blr.mean.l)
}
# Impute NA :
imputed.df <- impute_na(dataframe)
# We need to modify the scaling function, it should scale the dataframe WITHOUT imputed values on the basis of the dataframe
# with imputed values 
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
# Scale the df
scalingfunction(dataframe)
# The algorithm :


# 2. Scale it
dataframe <- read_csv("GeoPT48 -84Ra.csv")
# 1. Cleaning step
dataframe <- df_cleaning(dataframe)
# Create a list of n dataframe :
df.iteration1 <- scalingfunction(dataframe,impute_na(dataframe))
df.iteration2 <- scalingfunction(df.iteration1,impute_na(df.iteration1))
df.iteration3 <- scalingfunction(df.iteration2,impute_na(df.iteration2))

mean(rowSums(df.iteration2,na.rm=T)) 
mean(rowSums(df.iteration1,na.rm=T)) 


rowSums(scalingfunction(impute_na(dataframe)))
# Set n to 10
n <- 10
# Fill this list
for (i in 1:n){
  list.df[[i]] <- dataframe
}
list.df <- list()
list.df[[1]] <- dataframe
list.df[[2]] <- scalingfunction(list.df[[1]],impute_na(list.df[[1]]))
i <- 2
# Do while the difference in mean is above a threshold

while (mean(rowSums(list.df[[i-1]],na.rm=T))-mean(rowSums(list.df[[i]],na.rm=T)) > 10E-8){
  list.df[[i+1]] <- scalingfunction(list.df[[i]],impute_na(list.df[[i]]))
  i <- i + 1

}
list.df
list.df[[10]] == list.df[[9]]

# Convergence, debug the while loop.


