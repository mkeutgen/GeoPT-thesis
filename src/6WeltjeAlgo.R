
###########################################
## GJ Weltje Algorithm
## Iterative Estimation of Mean Composition
## GeoPT Thesis
###########################################

# A. Initial Vectors

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

# B. Let a dataset "data.raw", n x p matrix with n analysis and p analytes. 
data.raw <- read_csv("~/Documents/MStatistics/MA2/Thesis/Repository/data/raw/GeoPT48 -84Ra.csv")

# Elementary cleaning, removing columns with no values in the major elements

df.majors <- select(data.raw,all_of(ox.majors))
df.traces <- select(data.raw,all_of(traces))

cutoff.major <- 1/2
cutoff.trace <- 3/4

# Cut off rows with more than 50 % of major analytes are dropped this step is 
# important to avoid infinity values in the logit space
m <- which(rowSums(is.na(df.majors)) < cutoff.major*ncol(df.majors))

# Cut off rows with more than 75 % trace elements are dropped
n <- which(rowSums(is.na(df.traces)) < cutoff.trace *ncol(df.traces))

index <- intersect(m,n)

data <- data.raw[index,]



# Unit Conformity 
data[traces] <- data[traces] *10E-4

# Create columns X (sum of major elements), OX (sum of oxygen), Tr (sum of traces)

data <- data %>% mutate(X = rowSums(data %>% select(ox.majors)*frac_el,na.rm = T),
                OX = rowSums(data %>% select(ox.majors)*frac_ox,na.rm = T),
                Tr = rowSums(data %>% select(traces),na.rm = T))

# Feasibility matrix

f.mat <- data %>% select(X,OX,Tr)
# Scaling this matrix to unit sum
f.mat <- as_tibble(clo(f.mat))

f.mat <- f.mat %>% mutate(blrX = log(X/(1-X)),
                 blrTr = log(Tr/(1-Tr)),
                 blrOX = log(OX/(1-OX)),
                 R = 100-X-OX-Tr
                )
blrX.m <- mean(f.mat$blrX)
blrTr.m <- mean(f.mat$blrTr)
blrOx.m <- mean(f.mat$blrOX)

Inv.X <- (exp(blrX.m) )/ ( 1 + exp(blrX.m))
Inv.Tr <- (exp(blrTr.m) )/ ( 1 + exp(blrTr.m))
Inv.Ox <- (exp(blrOx.m) )/ ( 1 + exp(blrOx.m))
# R is :
Inv.X + Inv.Tr + Inv.Ox


# Next we impute missing values by their means after logit transform. 
# Convergence criterion : Euclidean distance of successive estimates X' and T'

# Should I proceed column (element-wise) or row-wise (analyte after analyte) ?

# Replace iteratively NAs values columnwise

df <- 1/100*data %>% select(ox.majors,traces)

#for i in 1:ncol(df)


# computing mean of all columns using apply()

logit <- function(x){log((x)/(1-x))}

blr.mean <- function(x){
  logit.t <- sapply(x,logit)
  m <- mean(logit.t,na.rm=T)
  result <- (exp(m))/( 1 + exp(m))
  return(result)
}

blr.mean.v <- c()
for (i in 1:ncol(df)){
  blr.mean.v[i] <- blr.mean(df[i])
}
names(blr.mean.v) <- names(df)
# imputing NA with the mean calculated 
replace_na(as_tibble(df),as.list(blr.mean.v))

# Iteratively ? 