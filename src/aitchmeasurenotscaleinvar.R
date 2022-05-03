# Aitchison measure is not scale invariant.


compo1 <- c(0.6,0.4)

compo2 <- c(0.65,0.35)
# Translate it

compo3 <- c(0.05,0.95)
compo4 <- c(0.1,0.9)

# Lebesgue measure between compo 1 and compo 2 is translation invariant

abs(compo1 - compo2) %>% prod()
abs(compo3 - compo4) %>% prod()

# Aitchison measure is not translation invariant
abs(ilr(compo1)-ilr(compo2)) %>% prod() 
abs(ilr(compo3)-ilr(compo4)) %>% prod() 

# BLR MEAN Versus ALR MEAN
c(0.99,0.02)
df <- cbind(rnorm(10,0.99,sd = 0.01),rnorm(10,0.02,sd = 0.01)) %>% as_tibble()
alr(df) %>% mean() %>% alrInv() 
apply(df,2,blr.mean)
