require(nimble)

### Set of functions needed in the simulation study to compute connectivity from occupancy data ###

inv.logit <- function(x){
  y <- exp(x)/(1+exp(x))
  return(y)
}

Simulate_occ_data <- function(M, J, T, psi1, alpha2, sigma, gamma0, phi, p, method){
  # Simulate dynamic site occupancy data from least cost path or circuit occupancy model
  
  # M = Number of sites
  # J = Number of secondary sampling occasion
  # T = Number primary sampling occasions = seasons
  # psi1 = initial occupancy probability
  # alpha2 = resistance parameter
  # sigma = scale parameter
  # gamma0 = baseline colonization probability 
  # phi = persistence probability = 1 - extinction probability
  # p = detection probability
  # method = method used to compute distance between site: 
  # - "LCP" = Least Cost Path (Howell et al. 2018), 
  # - "Circuit" = Circuit theory commute time (this study) 
  
  # Initialization 
  muZ <- gamma <- z <- array(NA, dim = c(M, T)) # initialize muZ, gamma : colonization probability, z : latent state of occurrence
  y <- array(NA, dim = c(M, J, T)) # detection histories
  
  # Compute non-Euclidean distance
  if (method == "LCP"){
    distSq <- leastcostpath(alpha2) # least cost path
  } else {
    if (method == "Circuit"){
      distSq <- CircuitDistance(alpha2) # circuit distance
    }
  }
  distSqReduced <- array(NA, dim = c(M, M)) # reduced to help convergence
  distSqReduced[1:M,1:M] <- distSq[1:M,1:M]/sd(distSq[1:M,1:M])
  
  gammaDistPairs <- array(NA, dim = c(M, M, T-1)) # pairwise colonization probability
  
  # Generate latent states of occurrence
  ## First year
  #z[1:M,1] <- rbinom(M, 1, psi1) # initial occupancy state
  z[1:M,1] <- 0 # initial occupancy state
  z[c(37,48, 49, 60, 70),1] <- 1
  #z[c(140,141,142,294,295,250,251,252, 253, 227,230,186,207,162, 420),1] <- 1 # 500m grid
  
  
  ## Following seasons
  for(t in 2:T){                                  # loop over seasons
    for(i in 1:M){                                # loop over sites
      for(n in 1:M){                              # loop over sites
        ### Pairwise colonization probability
        gammaDistPairs[i,n,t-1] <- 1 - gamma0 * exp(-(distSqReduced[i,n]) / (2*sigma^2)) * z[n,t-1]
      }
      ### Colonization probability
      gamma[i,t-1] <- 1 - prod(gammaDistPairs[i,1:M,t-1])
      muZ[i,t-1] <- gamma[i,t-1]*(1-z[i,t-1]) + phi* z[i,t-1]
      z[i,t] <- rbinom(1, 1, muZ[i,t-1])
    }
  }
  
  # Generate detection/non-detection data
  for(i in 1:M){
    for(t in 1:T){
      muP <- z[i,t] * p
      for(j in 1:J){
        y[i,j,t] <- rbinom(1, 1, muP)
      }
    }
  }
  return(list(z = z, y = y))
}
# 1. Mackenzie Model --------------------------------------------------------------
mackenzie.model <- nimbleCode({
  
  # Priors: ecological parameters
  psi1 ~ dunif(0, 1)  # initial occupancy probability
  gamma ~ dunif(0, 1) # colonization probability
  phi ~ dunif(0, 1)   # extinction probability
  p ~ dunif(0, 1)     # detection probability
  
  # ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    
    for (t in 2:nseasons){
      muZ[i,t]<- z[i,t-1] * phi + (1 - z[i,t-1]) * gamma
      z[i,t] ~ dbern(muZ[i,t])
    }
  }
  
  # observation model
  for (i in 1:nsites){
    for (t in 1:nseasons){
      y[i,t] ~ dbinom(size = nsurveys, prob = z[i,t] * p)
    }
  }
})

# 2. Euclidean Model --------------------------------------------------------------
euclidean.model <- nimble::nimbleCode({
  
  # Priors: ecological parameter
  psi1 ~ dunif(0, 1)   # first session occupancy probability
  gamma0 ~ dunif(0, 1) # baseline colonization probability
  sigma ~ dgamma(shape = 1, rate = 1) # scale parameter
  phi ~ dunif(0, 1)   # extinction probability
  p ~ dunif(0, 1)     # detection probability
  
  # ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    
    for(t in 2:nseasons) {
      for(n in 1:nsites) {
        # pairwise colonization probability
        gammaDistPairs[i,n,t-1] <- 1 - gamma0 * exp(-distSq[i,n] / (2 * sigma^2)) * z[n,t-1]
      }
      # colonization probability
      gamma[i,t-1] <- 1 - prod(gammaDistPairs[i,1:nsites,t-1])
      ## Pr(z=1|z=1) = 1 - Pr(extinction)*Pr(not colonized) = 1 - ((1-phi)*(1-gamma))
      muz[i,t-1] <- gamma[i,t-1] * (1 - z[i,t-1]) + phi * z[i,t-1]
      z[i,t] ~ dbern(muz[i,t-1])
    }
  }
  # observation model
  for (i in 1:nsites){
    for (t in 1:nseasons){
      y[i,t] ~ dbinom(size = nsurveys, prob = z[i,t] * p)
    }
  }
})

# 3. Lcp Model --------------------------------------------------------------

leastcostpath <- function(alpha){
  # Compute least cost path between a set of points given a resistance parameter alpha and a resistance surface rcov
  # alpha is a scalar
  # rcov is a RasterLayer
  # xy is a matrix of dims nsites x 2
  
  # obtain resistance surface
  cost <- exp(alpha * rcov)
  # calculate conductances among neighbours 
  # probability of transition from one cell to another
  tr <- gdistance::transition(x = cost, 
                              transitionFunction = function(x) 1/mean(x),
                              directions = 8) # class TransitionLayer
  # adjust diag.conductances
  trCorrC <- gdistance::geoCorrection(x = tr, 
                                      type = "c", 
                                      multpl = FALSE, 
                                      scl = FALSE) # class TransitionLayer
  # compute ecological distance
  d <- gdistance::costDistance(x = trCorrC,
                               fromCoords = xy, 
                               toCoords = xy) # class Matrix
  d^2
}

# Create the nimble funtion to compute distance between sites 
LCP <- nimbleRcall(function(alpha = double(0)){},
                   Rfun = 'leastcostpath',
                   returnType = double(2))
## 3.1 Normal ----------------------------------------------------------------
lcp.model <- nimble::nimbleCode({
  
  # Priors: ecological parameter
  psi1 ~ dunif(0, 1)   # first session occupancy probability
  gamma0 ~ dunif(0, 1) # baseline colonization probability
  alpha2 ~ dunif(-5, 5) # resistance parameter
  sigma ~ dgamma(shape = 1, rate = 1) # scale parameter
  phi ~ dunif(0, 1)   # extinction probability
  p ~ dunif(0, 1)     # detection probability
  
  # computation of distances for a given resistance
  distSq[1:nsites,1:nsites] <- LCP(alpha2) 
  distSqreduced[1:nsites,1:nsites] <- distSq[1:nsites,1:nsites] / sd(distSq[1:nsites,1:nsites]) # reduced to help convergence
  
  # ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    
    for(t in 2:nseasons) {
      for(n in 1:nsites) {
        # pairwise colonization probability
        gammaDistPairs[i,n,t-1] <- 1 - gamma0 * exp(-distSqreduced[i,n] / (2 * sigma^2)) * z[n,t-1]
      }
      # colonization probability
      gamma[i,t-1] <- 1 - prod(gammaDistPairs[i,1:nsites,t-1])
      muz[i,t-1] <- gamma[i,t-1] * (1 - z[i,t-1]) + phi * z[i,t-1]
      z[i,t] ~ dbern(muz[i,t-1])
    }
  }
  # observation model
  for (i in 1:nsites){
    for (t in 1:nseasons){
      y[i,t] ~ dbinom(size = nsurveys, prob = z[i,t] * p)
    }
  }
})

## 3.2 Cached ----------------------------------------------------------------
lcp.model.cached <- nimble::nimbleCode({
  
  # Priors: ecological parameter
  psi1 ~ dunif(0, 1)   # first session occupancy probability
  gamma0 ~ dunif(0, 1) # baseline colonization probability
  alpha2 ~ dunif(-5, 5) # resistance parameter
  sigma ~ dgamma(shape = 1, rate = 1) # scale parameter
  one ~ dconstraint(alpha2 >= -5 & alpha2 <= 5) # be sure that alpha2 is not outside cached array
  phi ~ dunif(0, 1)   # extinction probability
  p ~ dunif(0, 1)     # detection probability
  
  # computation of distances for a given resistance
  AlphaID <- trunc(alpha2*10+51) # formula (round(alpha2,1)*10+51)  to find were the distances values are recorded for these alpha 
  distSq[1:nsites,1:nsites] <- cached[1:nsites, 1:nsites,AlphaID] / sd(cached[1:nsites, 1:nsites,AlphaID])# retrive distances from the cached matrix 
  
  # ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    
    for(t in 2:nseasons) {
      for(n in 1:nsites) {
        # pairwise colonization probability
        gammaDistPairs[i,n,t-1] <- 1 - gamma0 * exp(-distSq[i,n] / (2 * sigma^2)) * z[n,t-1]
      }
      # colonization probability
      gamma[i,t-1] <- 1 - prod(gammaDistPairs[i,1:nsites,t-1])
      muz[i,t-1] <- gamma[i,t-1] * (1 - z[i,t-1]) + phi * z[i,t-1]
      z[i,t] ~ dbern(muz[i,t-1])
    }
  }
  # observation model
  for (i in 1:nsites){
    for (t in 1:nseasons){
      y[i,t] ~ dbinom(size = nsurveys, prob = z[i,t] * p)
    }
  }
})

# 4. Circuit model -----------------------------
circuitdistance <- function(alpha){
  # alpha is a scalar
  # rcov is a RasterLayer
  # xy is a matrix of dims nsites x 2
  # obtain resistance surface
  cost <- exp(alpha * rcov)    
  # calculate conductances among neighbours 
  # probability of transition from one cell to another
  tr <- gdistance::transition(x = cost, 
                              transitionFunction = function(x) 1/mean(x),
                              directions = 8) # class TransitionLayer
  # adjust diag.conductances
  trCorrC <- gdistance::geoCorrection(x = tr, 
                                      type = "r", # for randomwalk
                                      scl = TRUE) # class TransitionLayer 
  
  #Circuit theory = random walk with theta = 0
  circuitDist <- gdistance::commuteDistance(x = trCorrC, coords = xy) %>%  
    as.matrix() 
  dist <- circuitDist
  return(dist^2)
}

CircuitDistance <- nimbleRcall(function(alpha = double(0)){},
                               Rfun = 'circuitdistance',
                               returnType = double(2))

## 4.1 Normal ------------------------------------------------------------------
circuit.model <- nimble::nimbleCode({
  
  # Priors: ecological parameter
  psi1 ~ dunif(0, 1)   # first session occupancy probability
  gamma0 ~ dunif(0, 1) # baseline colonization probability
  alpha2 ~ dunif(-5, 5) # resistance parameter
  sigma ~ dgamma(shape = 1, rate = 1) # scale parameter
  phi ~ dunif(0, 1)   # extinction probability
  p ~ dunif(0, 1)     # detection probability
  
  # computation of distances for a given resistance
  distSq[1:nsites,1:nsites] <- CircuitDistance(alpha2)
  distSqreduced[1:nsites,1:nsites] <- distSq[1:nsites,1:nsites] / sd(distSq[1:nsites,1:nsites])
  
  # ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    
    for(t in 2:nseasons) {
      for(n in 1:nsites) {
        # pairwise colonization probability
        gammaDistPairs[i,n,t-1] <- 1 - gamma0 * exp(-distSqreduced[i,n] / (2 * sigma^2)) * z[n,t-1]
      }
      # colonization probability
      gamma[i,t-1] <- 1 - prod(gammaDistPairs[i,1:nsites,t-1])
      muz[i,t-1] <- gamma[i,t-1] * (1 - z[i,t-1]) + phi * z[i,t-1]
      z[i,t] ~ dbern(muz[i,t-1])
    }
  }
  # observation model
  for (i in 1:nsites){
    for (t in 1:nseasons){
      y[i,t] ~ dbinom(size = nsurveys, prob = z[i,t] * p)
    }
  }
})

## 4.2 Cached ------------------------------------------------------------------

circuit.model.cached <- nimble::nimbleCode({
  # Priors: ecological parameter
  psi1 ~ dunif(0, 1)   # first session occupancy probability
  gamma0 ~ dunif(0, 1) # baseline colonization probability
  alpha2 ~ dunif(-5, 5) # resistance parameter
  one ~ dconstraint(alpha2 >= -5 & alpha2 <= 5) # be sure that alpha2 is not outside cached array
  sigma ~ dgamma(shape = 1, rate = 1) # scale parameter
  phi ~ dunif(0, 1)   # extinction probability
  p ~ dunif(0, 1)     # detection probability
  
  # computation of distances for a given resistance
  AlphaID <- trunc(alpha2*10+51) # formula (round(alpha2,1)*10+51)  to find were the distances values are recorded for these alpha 
  distSq[1:nsites,1:nsites] <- cached[1:nsites, 1:nsites,AlphaID]/ sd(cached[1:nsites, 1:nsites,AlphaID]) # retrieve distances from the cached matrix 
  
  # ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    
    for(t in 2:nseasons) {
      for(n in 1:nsites) {
        # pairwise colonization probability
        gammaDistPairs[i,n,t-1] <- 1 - gamma0 * exp(-distSq[i,n] / (2 * sigma^2)) * z[n,t-1]
      }
      # colonization probability
      gamma[i,t-1] <- 1 - prod(gammaDistPairs[i,1:nsites,t-1])
      muz[i,t-1] <- gamma[i,t-1] * (1 - z[i,t-1]) + phi * z[i,t-1]
      z[i,t] ~ dbern(muz[i,t-1])
    }
  }
  # observation model
  for (i in 1:nsites){
    for (t in 1:nseasons){
      y[i,t] ~ dbinom(size = nsurveys, prob = z[i,t] * p)
    }
  }
})
