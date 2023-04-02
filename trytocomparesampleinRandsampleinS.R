# Define the lognormal-on-I density
library(compositions)
lognormal.on.I <- function(x, sigma, mu) {
  (1/(sqrt(2*pi))) * exp(-((log(x/(1-x))-mu)^2)/(2*sigma^2))
}

logit <- function(x){log((x)/(1-x))}
logitInv <- function(x){1/(1+exp(-x))}



# Define the proposal distribution (uniform between 0 and 1)
proposal <- function(n) {
  runif(n)
}

# Set the parameters
sigma <- 0.5
mu <- -2

# Define the acceptance ratio
acceptance <- function(x, sigma, mu) {
  lognormal.on.I(x, sigma, mu) / dunif(x, min = 0, max = 1)
}

# Generate samples using rejection sampling
n <- 1E6 # number of samples
samples <- numeric(n)
i <- 1
while (i <= n) {
  x <- proposal(1)
  u <- runif(1)
  if (u <= acceptance(x, sigma, mu)) {
    samples[i] <- x
    i <- i + 1
  }
}

# Plot the resulting samples
hist(samples, breaks = 30, freq = FALSE, main = "Rejection Sampling")


hist(logitInv(rnorm(1000,1,0.5)))


logit.in.R <- logitInv(rnorm(1E6,-2,0.5))
logit.in.I <- samples

logit.in.I %>% mean()
logit.in.R %>% mean()
logit.in.I %>% median()

mean.I <- geometricmean((logit.in.I)/(1-logit.in.I))/(1+geometricmean(logit.in.I/(1-logit.in.I)))
mean.R <- geometricmean((logit.in.R)/(1-logit.in.R))/(1+geometricmean(logit.in.R/(1-logit.in.R)))
logitInv(-2)


tibble(logit.in.I,logit.in.R) %>% pivot_longer(cols=everything()) %>%
  ggplot(aes(x=value,color=name))+geom_density()

max(logit.in.I)
max(logit.in.R)

my_mode <- function(x) {                     # Create mode function 
  unique_x <- unique(x)
  tabulate_x <- tabulate(match(x, unique_x))
  unique_x[tabulate_x == max(tabulate_x)]
}

my_mode(logit.in.I)[1]
my_mode(logit.in.R)[1]

mean(logit.in.I)
mean(logit.in.R)

5/3

exp(log(5)-log(3))
5-3
exp(log(5/3))

curve(lognormal.on.I(x, sigma, mu), add = TRUE, col = "red", lwd = 2)



logit(rnorm(1000,0.5,1)) %>% hist()


