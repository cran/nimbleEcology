## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  ,eval = TRUE ## uncomment this to build quickly without running code.
)

## ---- results='hide', messages=FALSE,warnings=FALSE---------------------------
library(nimble)
library(nimbleEcology)

## -----------------------------------------------------------------------------
occupancy_code <- nimbleCode({
  psi ~ dunif(0,1)
  p ~ dunif (0,1)
  for(i in 1:nSites) {
    z[i] ~ dbern(psi)
    for(j in 1:nVisits) {
      y[i, j] ~ dbern(z[i] * p)
    }
  }
})

## -----------------------------------------------------------------------------
occupancy_code_new <- nimbleCode({
  psi ~ dunif(0,1)
  p ~ dunif (0,1)
  for(i in 1:nSites) {
    y[i, 1:nVisits] ~ dOcc_s(probOcc = psi, probDetect = p, len = nVisits)
  }
})

## -----------------------------------------------------------------------------
occupancy_model <- nimbleModel(occupancy_code,
                               constants = list(nSites = 50, nVisits = 5))

## -----------------------------------------------------------------------------
occupancy_model$psi <- 0.7
occupancy_model$p <- 0.15
simNodes <- occupancy_model$getDependencies(c("psi", "p"), self = FALSE)
occupancy_model$simulate(simNodes)
occupancy_model$z
head(occupancy_model$y, 10) ## first 10 rows
occupancy_model$setData('y') ## set "y" as data

## -----------------------------------------------------------------------------
MCMCconf <- configureMCMC(occupancy_model)
MCMC <- buildMCMC(occupancy_model)

## -----------------------------------------------------------------------------
## These can be done in one step, but many people
## find it convenient to do it in two steps.
Coccupancy_model <- compileNimble(occupancy_model)
CMCMC <- compileNimble(MCMC, project = occupancy_model)

## ---- results='hide', messages=FALSE,warnings=FALSE---------------------------
samples <- runMCMC(CMCMC, niter = 10000, nburnin = 500, thin = 10)

## ---- results='hide', messages=FALSE,warnings=FALSE---------------------------
occupancy_model_new <- nimbleModel(occupancy_code_new,
                                   constants = list(nSites = 50, nVisits = 5),
                                   data = list(y = occupancy_model$y),
                                   inits = list(psi = 0.7, p = 0.15))
MCMC_new <- buildMCMC(occupancy_model_new) ## This will use default call to configureMCMC.
Coccupancy_model_new <- compileNimble(occupancy_model_new)
CMCMC_new <- compileNimble(MCMC_new, project = occupancy_model_new)
samples_new <- runMCMC(CMCMC_new, niter = 10000, nburnin = 500, thin = 10)

## ---- echo = FALSE, results='hide', messages=FALSE,warnings=FALSE-------------

{
plot(density(as.data.frame(samples)$psi), col = "red", main = "psi")
points(density(as.data.frame(samples_new)$psi), col = "blue", type = "l")
}

{
plot(density(as.data.frame(samples)$p), col = "red", main = "p")
points(density(as.data.frame(samples_new)$p), col = "blue", type = "l")
}


## -----------------------------------------------------------------------------
CalcLogLik <- nimbleFunction(
  setup = function(model, nodes)
    calcNodes <- model$getDependencies(nodes, self = FALSE),
  run = function(v = double(1)) {
    values(model, nodes) <<- v
    return(model$calculate(calcNodes))
    returnType(double(0))
  }
)
OccLogLik <- CalcLogLik(occupancy_model_new, c("psi", "p"))
COccLogLik <- compileNimble(OccLogLik, project = occupancy_model_new)
optim(c(0.5, 0.5), COccLogLik$run, control = list(fnscale = -1))$par

