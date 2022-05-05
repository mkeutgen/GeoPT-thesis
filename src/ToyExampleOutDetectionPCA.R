#### Toy Example Outlier Detection ####

# Toy example, reconstruction of data matrix based on PCA on correlation matrix (scaled)
X = iris[,1:4]


RMSEfun <- function(X,k=1){
  Xpca = prcomp(X,scale. = T)
  mu = colMeans(X)
  col.sd <- apply(X,2,sd)
  
  nComp <- k
  
  Xhat = Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
  
  Xhat <- Xhat %*% diag(col.sd)
  
  Xhat = scale(Xhat, center = -mu, scale = FALSE)
  # Euclidean distance between each observations X and Xhat : 
  ED <- sqrt(rowSums((X-Xhat)^2))
  RMSE <- sqrt(mean(ED^2))
  output <- list(RMSE,ED)
  names(output) <- c("RMSE","L2 distance between X and Xhat")
  return(output)
}


RMSEfun(X,k = 3)[1]
