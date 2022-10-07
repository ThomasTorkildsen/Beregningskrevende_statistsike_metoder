getwd()
setwd("/Users/thoma/Documents/skole/beregningskrevende/comp3")
source("additionalFiles/probAhelp.R")
source("additionalFiles/probAdata.R")

plot(data3A$x)

x = data3A$x
betas = ARp.beta.est(x, p = 2)
betaLS = betas$LS
betaLA = betas$LA
betas
epsLA = ARp.resid(x, betaLA)
epsLS = ARp.resid(x, betaLS)

plot(eLA)
plot(eLS)
x
betaLA = as.double(betaLA)
filter(x,c(1,-betaLA),sides=1)[-(1:length(beta))]
filter(x, filter = c(1), sides = 1)

#install.packages("matrixStats")
library(matrixStats)

Te = 100
B = 1500

boot_resLA <- boot_resLS <- matrix(nrow = Te, ncol = B)
boot_AR_LA <- boot_AR_LS <- matrix(nrow = Te, ncol = B)
boot_betaLA <- boot_betaLS <- matrix(nrow = 2, ncol = B)

for (i in 1:B) {
  boot_resLA[,i]= sample(epsLA, size=Te, replace=TRUE)
  boot_resLS[,i] = sample(epsLS, size=Te, replace=TRUE)
  
  ind = sample(1:(Te-1), size = 2, replace = TRUE)
  boot_AR_LA[,i] = ARp.filter(x[ind[1]:(ind[1]+1)],betaLA, boot_resLA[,i])[3:(Te+2)]
  boot_AR_LS[,i] = ARp.filter(x[ind[2]:(ind[2]+1)],betaLS, boot_resLS[,i])[3:(Te+2)]
  
  boot_betaLA[,i] = ARp.beta.est(boot_AR_LA[,i], p = 2)$LA
  boot_betaLS[,i] = ARp.beta.est(boot_AR_LS[,i], p = 2)$LS
  
}

betas = data.frame(LA = betaLA, mean_LA_boot = rowMeans(boot_betaLA), 
                     var_LA_boot = rowVars(boot_betaLA), 
                     mean_LS_boot =rowMeans(boot_betaLS), 
                     var_LA_boot = rowVars(boot_betaLS))

row.names(betas) = c("beta_1", "beta_2")
hist(boot_betaLA[1,], n = 50, freq = T)
abline(v = betasLA$est_betaLA_bt[1], col = "red")
abline(v = betasLA$betaLA[1], col = "green")
hist(boot_betaLS[1,], n = 50, freq = T)
abline(v = betasLS$est_betaLS_bt[1], col = "red")
abline(v = betasLS$betaLS[1], col = "green")



?mean()

## Copper-nickel alloy, page 291 in Givens & Hoeting book

y <- c(127.6, 124.0, 110.8, 103.9, 101.5, 130.1, 
       122.0, 92.3, 113.1, 83.7, 128.0, 91.4, 86.2)
x <- c(0.01, 0.48, 0.71, 0.95, 1.19, 0.01, 0.48, 
       1.44, 0.71, 1.96, 0.01, 1.44, 1.96)

ndata <- length(y)


## this function compute the ratio of interest beta_1/beta_0
mylm <- function(y, x){
  
  fit <- lm(y ~ x)
  coefs <- coef(fit)
  theta <- as.numeric(coefs[2]/coefs[1])
  return(theta)
}

originalModel <- lm(y ~ x)
plot(fitted(originalModel), y)
abline(0,1, col=2)
plot(residuals(originalModel))
abline(h=0,lty=2)
plot(fitted(originalModel),residuals(originalModel))
abline(h=0,lty=2)


originalFit =  mylm(y,x)


#### boottrap  residuals
B <- 1000

res = originalModel$residuals
fit = originalModel$fitted.values

idx = matrix(sample(1:ndata, ndata*B, replace=T),
             ncol=B, byrow=F)

idx[,1]

bootest1 = sapply(1:B, function(i){
  mylm(y = fit[idx[,i]]+res[idx[,i]], x = x[idx[,i]])})

?sapply()

## bootstrap pairs
# decide which observation to include 
#(at random and with replacement)
obs.idx <- matrix(sample(1:ndata, ndata*B, replace=T), ncol=B, byrow=F)
# apply the regression model to each data subset
ob
bootest2 <- sapply(1:B, function(i){
  mylm(y[obs.idx[,i]], x[obs.idx[,i]])
}
)



library(MASS)
par(mfrow=c(1,2))
truehist(bootest1, prob=TRUE, 
         ylab="Density", xlab=expression(theta), 
         col="lightblue", xlim=range(c(bootest1, bootest2)),
         main="resample residuals")
abline(v=originalFit, col=2, lwd=3)

truehist(bootest2, prob=TRUE, 
         ylab="Density", xlab=expression(theta), 
         col="lightblue", xlim=range(c(bootest1, bootest2)),
         main="resample pairs")
abline(v=originalFit, col=2, lwd=3)

sd(bootest1)
sd(bootest2)


hist(bootest1)
abline(v = originalFit, col=2, lwd=3)
abline(v = mean(bootest1), col=3, lwd=3)


# determine bias
mean(bootest1) - originalFit
mean(bootest2) - originalFit

# bias corrected estimate
originalFit - mean(bootest1 - originalFit)
originalFit - mean(bootest2 - originalFit)


# compute 95% CI based on percentile method.
par(mfrow=c(1,1))
plot(sort(bootest1),seq(1:B)/B)

abline(h = c(0.025, 0.975))

abline(v = quantile(bootest1, prob=c(0.025, 0.975)))

