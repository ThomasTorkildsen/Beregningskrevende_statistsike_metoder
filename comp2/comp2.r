
library(ggplot2)
require(smooth)

#Loading file

getwd()
setwd("/Users/thoma/Documents/skole/beregningskrevende/comp2/")
load(file = "rain.rda") #loading data


summary(rain)

r = rain #redefining rain as r


r["prob_rain"] = r$n.rain/r$n.years #calculating the estimated probability of rain


# fix initial values, since log(0) gives error
init_fix <- function(p){
  p = ifelse(p == 0, 0.001, p)
  return(p)
}

r$prob_rain = init_fix(r$prob_rain)

plot(r$day, r$prob_rain)
r["tau"] = log(r$prob_rain) - log(1-r$prob_rain)



#log likelihood function for a joint binomial
log_likelihood <- function(n,y,p) {
  lbin = sum(lchoose(n,y) + y*log(p) + (n-y)*log(1-p))
  return(lbin)
}

#test
log_likelihood(r$n.years, r$n.rain, r$prob_rain)

#defining Q, without sigma
Q_def <- function(n=366){
  Q <- diag(2, n)
  Q[1,1] <- Q[n,n] <- 1
  Q[abs(row(Q)-col(Q)) == 1] = -1
  return(Q)
}

Q <- Q_def(6) #test function

#install.packages("invgamma")
library(invgamma) #package for sampling inverse gamma

# Task e, the MCMC function
mcmc_RW = function(ntimes = 1000)
{
  
  te = 366          #number of days
  Q = Q_def(te)     #366x366 precition matrix, excluding sigma_u
  
  # initial values
  init_pi = r$prob_rain #init probability of rain
  init_n = r$n.years #number of years for each day
  init_y = r$n.rain #number of days with rain more than 1mm
  init_tau = r$tau #initial tau value
  tau = matrix(nrow = ntimes, ncol = te) #creating a n x T matrix
  
  tau[1,] = init_tau #insering the initial tau into the matrix
  sigmas = vector(length = ntimes-1)
  
  alpha = 2 #given alpha value
  beta = 0.05 #given beta value
  count = 0 #counting the number of accepted values
  
  
  for(i in 2:ntimes) #starting the for loop
  {
    
    sigmas[i-1] = rinvgamma(1, shape = alpha + (te-1)/2, rate = 1/2*sum(diff(tau[i-1,])^2)+beta) #finding sigma, using gibbs
    
    #loop for estimating each tau
    for(j in 1:te){
      #print(tau[i-1,-j])
      if(j == 1) {
        mu = tau[i-1,2]
        var = sigmas[i-1]
      }
      else if(j == te) {
        mu = tau[i-1,te-1]
        var = sigmas[i-1]
      }
      else {
        mu = -1/2*(Q[j,j-1]*tau[i-1,j-1] + Q[j,j+1]*tau[i-1,j+1])
        var = sigmas[i-1]/2
      }
      #mu = -solve(Q[j,j]) %*% Q[j,-j] %*% tau[i-1,-j] #expected value of the conditional normal dist.
      #var = 1/sigmas[i-1]*Q[j,j] #variance of the conditional normal dist
      
      #print(var)
      temp = rnorm(1, mean = mu, sd = sqrt(var))
      #print(temp)
      temp_pi = 1/(1+exp(-temp))
      ratio = log_likelihood(init_n[j], init_y[j], temp_pi) - log_likelihood(init_n[j], init_y[j], 1/(1+exp(-tau[i-1,j])))
      log_alpha = min(0, ratio)
      #print("log_alpha ",log_alpha)
      #print("log_first: ", log_likelihood(init_n[j], init_y[j], temp_pi))
      #print("log_sec:", log_likelihood(init_n[j], init_y[j], 1/(1+exp(-tau[i-1,j]))))
      if(log(runif(1))< log_alpha){
        tau[i,j] = temp
        count = count + 1 #counting total number of times accepted
      }
      else
        tau[i,j] = tau[i-1,j]
    }
  }
  acceptrate = count/(ntimes*te) #general acceptance rate 
  ret = list(tau = tau, rate = acceptrate, sigmas = sigmas, pi = 1/(1+exp(-tau)))
  return(ret)
}


st = proc.time()[3]
test = mcmc_RW(n = 1000)
sto = proc.time()[3]
sto-st
test$rate
hist(test$pi[,1], freq = TRUE, n = 100)





mcmc_RW2 = function(ntimes = 1000, M = 3)
{
  if (M == 1){
    return(mcmc_RW(ntimes))
  }
  
  M_need = M        #saving M for later
  te = 366          #number of days
  Q = Q_def(te)     #366x366 precition matrix, excluding sigma_u
  
  # initial values
  init_pi = r$prob_rain #init probability of rain
  init_n = r$n.years #number of years for each day
  init_y = r$n.rain #number of days with rain more than 1mm
  init_tau = r$tau #initial tau value
  tau = matrix(nrow = ntimes, ncol = te) #creating a n x T matrix
  
  tau[1,] = init_tau #insering the initial tau into the matrix
  sigmas = vector(length = ntimes-1) #initialising a vector of sigmas
  
  alpha = 2 #given alpha value
  beta = 0.05 #given beta value
  count = 0 #counting the number of accepted values
  
  loopseq = seq(1,te,M) #creating loop for updating 
  
  #precomputing
  Qinv_j1 = solve(Q[1:(1+M-1),1:(1+M-1)]) #Q inverse for j = 1
  Qinv_jm = solve(Q[2:(2+M-1),2:(2+M-1)]) #Q inverse for j = m
  #QinvQ_j1 = Qinv_j1 %*% Q[1:(1+M-1),-(1:(1+M-1))] #QAA mult QAB for j = 1
  QinvQ_j1 = matrix(0, nrow = M, ncol = te-M)
  QinvQ_j1[,1]=-1
  #QinvQ_jm = Qinv_jm %*% Q[2:(2+M-1),-(2:(2+M-1))] #QAA mult QAB for j = m
  l = loopseq[length(loopseq)]
  if((l+M-1) > te) { #In case te mod M != 0
    M = te-l+1
  }
  #Q inverse for j = te
  Qinv_jM = solve(Q[l:(l+M-1),l:(l+M-1)]) #Q inverse for j = M
  
  QinvQ_jM = Qinv_jM %*% Q[l:(l+M-1),-(l:(l+M-1))] #QAA mult. QAB for j = M
  
  M = M_need #resetting M to original value
  
  
  
  for(i in 2:ntimes) #starting the for loop
  {
    
    sigmas[i-1] = rinvgamma(1, shape = alpha + (te-1)/2, rate = 1/2*sum(diff(tau[i-1,])^2)+beta) #finding sigma, using gibbs
    
    #loop for estimating each tau
    for(j in loopseq){
      if(j == 1){ 
        mu = -QinvQ_j1%*%tau[i-1,-(j:(j+M-1))]
        var = sigmas[i-1]*Qinv_j1 
      }
      else if ((j+M-1) < te){
        #mu = -Qinv_jm %*% Q[j:(j+M-1),-(j:(j+M-1))] %*% tau[i-1,-(j:(j+M-1))]
        mult = rep(0, M)
        mult[1] = -1*tau[i-1,j-1] 
        mult[M] = -1*tau[i-1,j+M]
        mu = -Qinv_jm %*% mult
        var = sigmas[i-1]*Qinv_jm
      }
      else { #In case te mod M != 0
        if((l+M-1) >= te) { #In case te mod M != 0
          M = te-l+1
        }
        mu = -QinvQ_jM %*% tau[i-1, -(j:(j+M-1))]
        var = sigmas[i-1]*Qinv_jM
      }
      #print(tau[i-1,-j])
      #Q_inv = solve(Q[j:(j+M-1),j:(j+M-1)])
      #mu = -Q_inv %*% Q[j:(j+M-1),-(j:(j+M-1))] %*% tau[i-1,-(j:(j+M-1))] #expected value of the conditional normal dist.
 #covariance matrix
      
      #print(var)
      temp = rmvnorm(1, mean = mu, sigma = var) #sample from mulitvariate normall
      #print(temp)
      temp_pi = 1/(1+exp(-temp))
      ratio = log_likelihood(init_n[j:(j+M-1)], init_y[j:(j+M-1)], temp_pi) - log_likelihood(init_n[j:(j+M-1)], init_y[j:(j+M-1)], 1/(1+exp(-tau[i-1,j:(j+M-1)])))
      log_alpha = min(0, ratio)
      
      #print("log_alpha ",log_alpha)
      #print("log_first: ", log_likelihood(init_n[j], init_y[j], temp_pi))
      #print("log_sec:", log_likelihood(init_n[j], init_y[j], 1/(1+exp(-tau[i-1,j]))))
      if(log(runif(1))< log_alpha){
        tau[i,j:(j+M-1)] = temp
        count = count + M #counting total number of times accepted
      }
      else
        tau[i,j:(j+M-1)] = tau[i-1,j:(j+M-1)]
    }
    
    M = M_need #In case te mod M != 0
  }
  acceptrate = count/(ntimes*(ceiling(te))) #general acceptance rate 
  ret = list(tau = tau, rate = acceptrate, sigmas = sigmas, pi = 1/(1+exp(-tau)))
  return(ret)
}

start = proc.time()[3]
x = mcmc_RW2(n=1000, M=3)
stop = proc.time()[3]
time = stop-start
x$rate
1/as.numeric(time)
hist(x$pi[,7], n = 10)


#plots...
plot(ts(x$tau[,1]))
hist(x$tau[,1], freq = TRUE, n = 20)
hist(x$sigmas, freq = TRUE, n = 20)
solve(Q[1,1]) %*% Q[1,-1] %*% T
## Task 2
a
init_tau = r$tau
tau = matrix(nrow = 1000, ncol = 366)
tau[1,] = init_tau  
tau  

#install.packages("INLA",repos=c(getOption("repos"),
#                                INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)

#fit the same model as in the previous problem with
control.inla = list(strategy="simplified.laplace", int.strategy="ccd")
hyper = list(prec = list(prior = "loggamma", param = c(2,0.05)))
mod <- inla(n.rain ~ -1 + f(day, model="rw1", constr=FALSE, hyper = hyper),
            data=r, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)


mod3 <- inla(n.rain ~ f(day, model="rw1", constr=TRUE, hyper = hyper),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)

mod3$summary.fixed
mod3$summary.fitted.values
INLA_fit = mod$summary.fitted.values #fitted values

#mean and 95% CI for pi1
pi1_inla = data.frame(mean = INLA_fit$mean[1], CI_lower = INLA_fit$`0.025quant`[1], CI_upper = INLA_fit$`0.975quant`[1])

#mean and 95% CI for pi201
pi201_inla = data.frame(mean = INLA_fit$mean[201], CI_lower = INLA_fit$`0.025quant`[201], CI_upper = INLA_fit$`0.975quant`[201])

#mean and 95% CI for pi366
pi366_inla = data.frame(mean = INLA_fit$mean[366], CI_lower = INLA_fit$`0.025quant`[366], CI_upper = INLA_fit$`0.975quant`[366])

t1 = rbind(pi1_inla, pi201_inla, pi366_inla)

plot(mod3$summary.fixed)
mod3$model.matrix
mod3$summary.fitted.values
plot(mod3$summary.linear.predictor$mean, type="l")
lines(mod$summary.fitted.values$mean)

INLA_fit = mod3$summary.fitted.values #fitted values

#mean and 95% CI for pi1
pi1_inla3 = data.frame(mean = INLA_fit$mean[1], CI_lower = INLA_fit$`0.025quant`[1], CI_upper = INLA_fit$`0.975quant`[1])

#mean and 95% CI for pi201
pi201_inla3 = data.frame(mean = INLA_fit$mean[201], CI_lower = INLA_fit$`0.025quant`[201], CI_upper = INLA_fit$`0.975quant`[201])

#mean and 95% CI for pi366
pi366_inla3 = data.frame(mean = INLA_fit$mean[366], CI_lower = INLA_fit$`0.025quant`[366], CI_upper = INLA_fit$`0.975quant`[366])

t2 = (rbind(pi1_inla3, pi201_inla3, pi366_inla3))

mod2 <- inla(n)(n.rain ~ -1 + f(day, model="rw1", constr=FALSE, hyper = hyper),
                data=r, Ntrials=n.years, control.compute=list(config = TRUE),
                family="binomial", verbose=TRUE)

mod3 <- inla(n.rain ~ f(day, model="rw1", constr=TRUE, hyper = hyper),
             data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
             family="binomial", verbose=TRUE, control.inla=control.inla)


summary(mod3)
mod3$summary.fitted.values$mean[1]
mod$summary.fitted.values$mean[1]
mod$marginals.hyperpar
?control.fixed

2+365/2

mean(x$pi[,1])
print(mod$summary.fitted.values - mean(x$pi[,200]))
?f
inla.doc("posterior")
?control.inla
?control.fixed
?inla.models
#available likelihood functions: names ( inla . models ()$ likelihood )
inla.doc("rw1")
#I Use -1 if you don???t want an automatic intercept

inla.list.models("prior")
# Given that the precision matrix is sparse we can use INLA

mod3$summary.fixed
mod
?inla()

inla.doc("f")

hei = mcmc_RW(50000)
p = mod3$summary.linear.predictor

summary.linear.predictor
m = mean(hei$pi[,1])
l = m - 1.96*sqrt(1/50000*var(hei$pi[,1]))
u = m + 1.96*sqrt(1/50000*var(hei$pi[,1]))

p$kld
mod3
t = exp(p)/(1+exp(p))

plot(t$mean)
?control.inla
#bruke conditional prior