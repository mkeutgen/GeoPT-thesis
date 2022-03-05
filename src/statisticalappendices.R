# Statistical Appendices 

# Simulation for the length of an IQR with 99 % of coveraging probabilities :

# What should be the multiple of IQR such that 99 % of observations are included in the asymptotic limiting normal distribution ?
# Average length of the interval which contains 99 % of observations :


# IQR for a standard normal
n.sim <- 10E3 # 1000 simulations
n <- 10E3 # 1000 observations in each simulation
set.seed(17)
X <- matrix(rnorm(n*n.sim), n)


IQRnorm <- mean(apply(X,MARGIN = 2,IQR))
IQRnorm
# What should be the multiple of IQR such that 99 % of observations are included in the asymptotic limiting normal distribution ?
# Average length of the interval which contains 99 % of observations :
X <- NULL
X <- matrix(rnorm(n*n.sim), n)


length99percent <- function(x) {abs(quantile(x,probs=c(0.005,0.995))[1])+quantile(x,probs=c(0.005,0.995))[2]}
vec <- c()
for (i in 1:ncol(X)){
  vec[i] <- length99percent(X[,i])  
}

coef.99 <- mean(vec)/IQRnorm
coef.99 # multiplying coeff found 3.81387


## Arguing against the KS test, artificially diminshes p-values. 
n.sim <- 1e4
n <- 50
set.seed(17)
X <- matrix(rnorm(n*n.sim), n)

f <- function(x) ks.test(x, "pnorm")$p.value
ks.1 <- apply(X, 2, f)
ks.2 <- apply(scale(X), 2, f)


# Figure 1: Histograms
par(mfrow=c(1,2))
b <- seq(0, 1, by=0.05)
hist(ks.1, breaks=b, freq=FALSE, col=gray(.9), main="Non-standardized", xlab="p-value")
abline(h=1, lwd=2, col=hsv(0,1,3/4))
hist(ks.2, breaks=b, freq=FALSE, col=gray(.9), main="Standardized", xlab="p-value")
abline(h=1, lwd=2, col=hsv(0,1,3/4))
par(mfrow=c(1,1))

# Figure 2: Scatterplot
plot(ks.1, ks.2, pch=21, bg=gray(0, alpha=.05), col=gray(0, alpha=.2), cex=.5,
     xlab="Non-standardized p-value", ylab="Standardized p-value", asp=1)

