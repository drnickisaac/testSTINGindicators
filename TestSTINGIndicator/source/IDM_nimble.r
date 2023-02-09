
require(nimble)

modelcode <- nimbleCode({ # BUGS code 
  
  ######################### State model
  for (i in 1:nsite){ 
    z[i] ~ dbern(psi[i]) 
    log(lambda[i]) <- cloglog(psi[i])
    logit(psi[i]) <- eta[i] 
  } 
  
  # random effect for site
  for (i in 1:nsite) {
    eta[i] ~ dnorm(0, tau2)       
  } 
  
  # hyperpriors on site effect
  tau2 <- 1/(sigma2 * sigma2) 
  sigma2 ~ T(dt(0, 1, 1), 0, Inf)  # Half Cauchy
  
  ######################### Obs model - pan traps
  
  ### Observation Model (occupancy detection model)
  for(i in 1:nsite) {
    for(j in 1:nvisit) {
      y1[i,j] ~ dbinom(Py[i,j], nT[i,j])
      Py[i,j]<- z[i] * p1[i,j]
      logit(p1[i,j]) <-  alpha.p + log(lambda[i]) + beta3 * f_JD[JulDate[i,j]]
    }}
  
  # Observation model priors 
  alpha.p ~ dnorm(mu.lp1, tau.lp1) # constant       
  
  # hyperpriors
  mu.lp1 ~ dnorm(-2, 0.01)
  tau.lp1 <- 1 / (sd.lp1 * sd.lp1)                 
  sd.lp1 ~ T(dt(0, 1, 1), 0, Inf)  # Half Cauchy
  
  ######################### Obs model - transects
  
  ### Observation Model
  for(i in 1:nsite) {
    for(j in 1:nvisit) {
      y2[i,j] ~ dpois(LamThin[i,j])
      LamThin[i,j] <- z[i] * lambda[i] * p2[i,j]
      logit(p2[i,j]) <-  alpha.p2 + beta3 * f_JD[JulDate[i,j]]
    } }
  
  # Observation model priors 
  alpha.p2 ~ dnorm(mu.lp2, tau.lp2) # constant       
  
  # hyperpriors
  mu.lp2 ~ dnorm(-2, 0.01)
  tau.lp2 <- 1 / (sd.lp2 * sd.lp2)                 
  sd.lp2 ~ T(dt(0, 1, 1), 0, Inf) # Half Cauchy
  
  
  ######################### Seasonality shared effect
  
  beta1 ~ dunif(1, 366)
  beta2 ~ T(dt(0, 1, 1), 0, Inf) # Half Cauchy
  beta3 ~ dnorm(0, 0.0001)
  for (d in 1:366){
    f_JD[d] <- 1/((2*3.14159265359)^0.5 * beta2) * exp(-((d - beta1)^2 / (2* beta2^2)))
  }
  
  ######################### 
  # Derived parameters: occupancy
  psi.fs <- mean(z[1:nsite])
  
})
