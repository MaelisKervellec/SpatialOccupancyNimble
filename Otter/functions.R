require(nimble)

inv.logit <- function(x){
  y <- exp(x)/(1+exp(x))
  return(y)
}

# 2. Circuit model -----------------------------
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


model.otter <- nimble::nimbleCode({
  # Priors: ecological parameter
  alpha_psi1 ~ dnorm(0, sd = 10) # initial occupancy probability
  beta_psi1 ~ dnorm(0, sd = 10) # initial occupancy probability
  for(l in 1:nseasons){
    p[l] ~ dunif(0, 1)} # detection probability per session
  gamma0 ~ dunif(0, 1) # baseline colonization probability
  alpha2 ~ dunif(1, 5) # resistance parameter ASSUME TO BE POSITIVE and STRONG
  one ~ dconstraint(alpha2 >= 0 & alpha2 < 5) # be sure that alpha2 is not outside cached array
  sigma ~ dgamma(shape = 1, rate = 1) # scale parameter
  phi ~ dunif(0, 1) # 1 - extinction probability
  
  # computation of distances for a given resistance
  #AlphaID <- trunc(alpha2*10+51) # formula (round(alpha2,1)*10+51)  to find were the distances values are recorded for these alpha 
  AlphaID <- trunc(alpha2*10)+1
  distSq[1:nsites,1:nsites] <- cached[1:nsites, 1:nsites,AlphaID] / sd(cached[1:nsites, 1:nsites,AlphaID]) # retrieve distances from the cached matrix 
  
  
  # ecological submodel
  for (i in 1:nsites){
    logit(psi1[i]) <- beta_psi1*elev[i] + alpha_psi1
    z[i,1] ~ dbern(psi1[i])
    
    for(t in 2:nseasons) {
      for(n in 1:nsites) {
        # pairwise colonization probability
        gammaDistPairs[i,n,t-1] <- 1 - gamma0 * exp(-distSq[i,n] / (2 * sigma^2)) * z[n,t-1]
      }
      # colonization probability
      gamma[i,t-1] <- 1 - prod(gammaDistPairs[i,1:nsites,t-1])
      ## Pr(z=1|z=1) = 1 - Pr(extinction)*Pr(not colonized) = 1 - ((1-phi)*(1-gamma))
      #muz[i,t-1] <- gamma[i,t-1] * (1 - z[i,t-1]) + (1 - (1 - phi) * (1 - gamma[i,t-1])) * z[i,t-1] # with rescue effect
      muz[i,t-1] <- gamma[i,t-1] * (1 - z[i,t-1]) + phi * z[i,t-1] # without rescue effect
      z[i,t] ~ dbern(muz[i,t-1])
    }
  }
  # observation model
  for (i in 1:nsites){
    for (j in 1:nsurveys){
      for (t in 1:nseasons){
        y[i,j,t] ~ dbern(z[i,t] * p[t])
      }
    }
  }
})

