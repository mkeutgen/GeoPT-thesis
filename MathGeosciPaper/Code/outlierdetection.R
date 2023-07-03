#########################
### OUTLIER DETECTION ###
#########################

setwd("MathGeosciPaper/Code/")

library(MASS)

library(tidyverse)
library(compositions)
library(xtable)
library(rootSolve)
library(EnvStats)
library(ggpubr)
library(quadprog)
library(missMethods)

#### Import True Datasets with GEOPT Missing Values ####

file.list <- list.files("sim/")
file.list
mean.compo <- readRDS("~/Documents/GeoPTManuscript/MathGeosciPaper/Code/meancompo.RData")
list.sim <- lapply(paste("sim/",file.list,sep = ""),readRDS)
names(list.sim) <- file.list

source("importgeoptdata.R")




mean.compo$`GeoPT25 .csv`
# Contaminate it with 50 \% of missing values :
list.sim$sim_list_GeoPT25.RData[[1]]

fun.mval <- function(tibble,propmissing){
  output <- gen.mval.fun(tibble[-ncol(tibble)],prop = propmissing)
  output$U <- output$sum.vec + tibble$U
  output$sum.vec <- NULL
  return(output)
}



# And a function to impute those missing values 
imputed.df.fun <- function(dataset){
  # Imputed df in logit space
  dataset <- dataset[!rowSums(dataset>1,na.rm = T),]
  imputed.df.logit <- logit(dataset)  %>% impute_EM() 
  imputed.df <- logitInv(imputed.df.logit)
  sum.imp <- c()
  for (i in 1:nrow(dataset)){
    # For each row, sum the missing elements which are imputed and substract them from U
    sum.imp[i] <- sum(imputed.df[i,which(is.na(dataset[i,]))] )
  }
  # Yields a set of realistic compositions, from compositions with missing values 
  imputed.df$U <- abs(imputed.df$U - sum.imp)
  return(imputed.df)
}

gen.mval.fun <- function(df,prop){
  # inner function, select random elements and take their sum
  random_row_sum <- function(row) {
    # select 2 val
    selected_indices <- sample(length(row), floor(prop*length(row)))
    selected_values <- row[selected_indices]
    sum(selected_values)
    return(list(sum(selected_values),selected_indices))
  }
  
  list.sim$sim_list_GeoPT2.RData[[4]] %>% names() 
  
  
  # create a sum vector 
  sum.vec <- c()
  for (i in 1:nrow(df)){
    random.row.o <- random_row_sum(df[i,])
    sum.vec[i] <- random.row.o[[1]]
    df[i,c(random.row.o[[2]])] <- NA
  }
  output <- cbind(df,sum.vec) %>% as_tibble()
  return(output)
}

optimizing.blr.mean <- function(dataset, min = -10, max = 10) {
  # mu_star is the initial mean estimate from the blr transformed data 
  logit.df <- logit(dataset)
  
  mu_star <- logit.df %>% colMeans()
  
  # sigma 
  sigma <- logit.df %>% apply(MARGIN = 2, sd)
  
  # Optimization: sum(logitInv(mu_star + alpha * sigma)) = 1
  f2 <- function(alpha) 1-sum(logitInv(mu_star + alpha * sigma)) 
  
  alpha <- tryCatch({
    uniroot(f2, c(min, max))$root
  }, error = function(e) {
    if (e$message == "f() values at end points not of opposite sign") {
      mid <- (min + max) / 2
      if (f2(mid) > 0) {
        max <- mid
      } else {
        min <- mid
      }
      uniroot(f2, c(min, max))$root
    } else {
      stop(e)
    }
  })
  
  mu_optimized <- logitInv(mu_star + alpha * sigma)
  
  return(mu_optimized)
}

#### Only select 20 first elements

# Outlier Detection Method 



df <- select(list.df$`GeoPT1 -83R.csv`,all_of(names.full.wO))

df


df <- df %>%
  filter(rowSums(is.na(.)) / ncol(.) <= 0.5)  %>%
  select(where(~ sum(is.na(.)) / length(.) <= 0.5))

df %>% length()

df[1:10] <- 1E-2 * df[1:10]
df[11:17] <- 1E-6 * df[11:17]
df$U <- ifelse(1-rowSums(df,na.rm = T)>0,1-rowSums(df,na.rm=T),NA)

imputed.df <- df %>% logit() %>% impute_EM() %>% logitInv() %>% as_tibble()

estimated.mean <- optimizing.blr.mean(imputed.df)
format(estimated.mean,digits = 4)

# Mahalanobis distance
maha.vec <- mahalanobis(logit(imputed.df),logit(estimated.mean),cov(logit(imputed.df)))
# Chi2 quantile
chisq <- qchisq(p=0.975, nrow(imputed.df) )

cbind(maha.vec,chisq) %>% as_tibble() %>% ggplot(aes(x=maha.vec))+
  geom_histogram()+theme_bw()+
  geom_vline(xintercept = chisq,color="red")+labs(x="Mahalanobis distances between the estimated mean and the chemical analyses")





order(maha.vec)
imputed.df[27,]
estimated.mean




fun.mval(list.sim$sim_list_GeoPT25.RData[[1]],0.5)
mval.df <- fun.mval(list.sim$sim_list_GeoPT25.RData[[1]],0.5)
imputed.df <- mval.df %>% imputed.df.fun()
pca <- prcomp(imputed.df %>% logit())  
pca$center %>% logitInv()

pca7 <- prcomp(logit(imputed.df),rank. = 7)


pca_result <- prcomp(logit(imputed.df), scale = TRUE)

# Step 2: Extract the first k principal components
k <- 7  # Number of components to keep
pcs <- pca_result$x[, 1:k]

# Step 3: Project the data onto the selected principal components
projected_data <- t(as.matrix(logit(imputed.df))) %*% pcs



X <- logit(imputed.df)
Xpca = prcomp(X,scale. = T)
mu = colMeans(X)
col.sd <- apply(X,2,sd)
nComp <- 7
Xhat = Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
Xhat <- Xhat %*% diag(col.sd)
Xhat = scale(Xhat, center = -mu, scale = FALSE)



distances <- mahalanobis(Xhat,center = logit(estimated.mean),cov = var(Xhat))



estimated.mean <- imputed.df %>% optimizing.blr.mean()
mean.compo$`GeoPT25 .csv`

distances <- mahalanobis(logit(imputed.df),center = logit(estimated.mean),cov = var(logit(imputed.df)))
order(distances,decreasing = T)




