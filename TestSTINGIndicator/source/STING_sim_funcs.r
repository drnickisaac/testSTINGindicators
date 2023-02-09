#


# create some data
createSimData <- function(nSite = 75,
                          psi.tru = 0.5, 
                          meanLambda = 4,
                          nVis = 3, # visits per site per year 
                          nTraps = 5, # pan traps per site 
                          pDetPan = 0.5
                          ){

  # define true occupancy - random process
  occupancy <- rbinom(n = nSite, size = 1, prob = psi.tru)
  
  # define the intensity of the point process. 
  # assume it's homogenous, i.e. all drawn from same distribution
  # but only occupied sites have intensity > 0
  intensity <- rpois(n = nSite, lambda = meanLambda) * occupancy
  
  # when do the visit occur?
  # assume evenly spaced between day 100 (10 April) & 269 (25 September)
  # so the middle visit is on day 181 (30 June)
  # we further assume that each site is sampled on the same day (not realistic!)
  when <- round(seq(from = 100, to = 269, length.out = nVis))

  # to add: Julian date modifier on abundance
  # simple function such that detection decays expoentially away from midsummer
  phenoMod <- exp(-abs(when - 180)/90)
  
  # pan traps
  # I want to create a data table with nSite * nVis 
  # in which the numbers are the number of times the species of interest was observed
  # first calculate a matrix of detection probabilities
  
  # the detection probability is given by 
  #  - intensity (transformed to a scaled probability), because detection is higher where populations are abundant
  #  - time of year (phenoMod)

  # start on the unbounded scale?
  # define the shape of the output
  
  # site level detection is determined by the intensity
  siteDet <- pDetPan * boot::inv.logit(log(intensity))^3
  
  # this is modified by the phenoMod to create visDet: the detection probability per visit
  visDet <- t(sapply(siteDet, function(x) x * phenoMod))
    
  # now there should be nTraps Bernoulli trials with probablities defined by these probabilities
  panTrapData <- t(apply(visDet, 1, function(p) {
                            rbinom(n = length(p), size = nTraps, prob = p)
    }))
  
  transDet <- t(sapply(intensity, function(x) x * phenoMod))
  transectData <- t(apply(transDet, 1, function(lam) {
      rpois(n = length(lam), lambda = lam)
    }))

  return(list(sites = data.frame(id = 1:nSite,
                                  occupancy = occupancy,
                                  intensity = intensity),
              dates =  t(replicate(nSite, when)),
              panTrapData = panTrapData,
              transectData = transectData,
              nTraps = nTraps
              )
         )
}

##############################################

summariseSimData <- function(data){
  occSites <- with(data, apply(cbind(panTrapData, transectData), 1, max)>0)
  c(
    naiveOcc = mean(occSites),
    reportingRate = mean(data$panTrapData), # per pan trap day, not per trap
    meanCount = mean(data$transectData),
    reportingRate_1 = mean(data$panTrapData[occSites,]), # per pan trap day, not per trap
    meanCount_1 = mean(data$transectData[occSites,])
  )
  # maybe add other summary stats (e.g. quantiles) on the counts
}


##############################################

getBUGSdata <- function(data){
  with(data, list(y1 = panTrapData,
                  y2 = transectData,
                  JulDate = dates,
                  nsite = dim(panTrapData)[1],
                  nvisit = dim(panTrapData)[2],
                  nT = matrix(nTraps, nr=dim(panTrapData)[1], nc=dim(panTrapData)[2])
  ))
}

##############################################

# fit the model
# for now use JAGS.
runModel <- function(data, 
                     nit = 1000,
                     nb = 200, 
                     nchain = 3, 
                     nthin = 5,
                     params = NULL){
  
  # coerce simulated data into a BUGS data object
  bugs_data <- getBUGSdata(data)

  # initial values
  initiate <- function(z, nsite) {
    init <- list (z = as.numeric(z),
                  alpha.p = runif(1, -2, 2),
                  alpha.p2 = runif(1, -2, 2)
                  #eta = rep(runif(1, -2, 2), nsite)
                  )
  }
  
  zst <- rowSums(cbind(bugs_data$y1, bugs_data$y2)) > 0
  init.vals <- replicate(nchain, initiate(z = zst, bugs_data$nsite),
                         simplify = F)
  
  if(is.null(params)) params <- c(
              "psi.fs", "alpha.p", "alpha.p2", 
              "beta1", "beta2", "beta3",
              "psi.fs", "lambda")
  
  # main event
  out <- R2jags::jags(data = bugs_data, 
               parameters.to.save = params,
               model.file = "source/IDM.txt",
               inits = init.vals,
               n.iter = nit,
               n.burnin = nb,
               n.chains = nchain,
               n.thin = nthin
  )
  
  return(out)
}


