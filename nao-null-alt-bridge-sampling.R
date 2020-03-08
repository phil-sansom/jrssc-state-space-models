#####################
## Bridge sampling ##
#####################

## Load libraries
library(compiler)
library(MASS)
library(parallel)
library(triangle)
library(bridgesampling)
library(coda)

## Load trajectories
load("nao-null-alt.RData")
options(mc.cores = mc.cores)

## Determine when adaptation ended
psrfs = sapply(parameters, function(x) max(gelman.diag(x, autoburnin = FALSE)$psrf[,1]))
alpha = min(which(psrfs < 1.10)) + 1
gamma = length(parameters)
blocks = gamma - alpha + 1

## Extract parameters
samples = list()
for (i in 1:n.chains) {
  samples[[i]] = parameters[[alpha]][[i]]
  for (j in (alpha + 1):gamma) {
    samples[[i]] = rbind(samples[[i]],parameters[[j]][[i]])
  } ## j
  samples[[i]] = samples[[i]][seq(n.chains*blocks,nrow(samples[[i]]),n.chains*blocks),]
} ## i
for (i in 1:n.chains)
  samples[[i]] = mcmc(samples[[i]], start = n.chains*blocks,
                      thin = n.chains*blocks)
samples = as.mcmc.list(samples)

lb = c(-32,-32,-64,-2.5,-1.0,-3.0,-64)
ub = c(  0,- 2,-18,+3.5,+5.0,+3.0,- 8)
names(lb) = c("v","wmu","wbeta","wxa","wxb","wxc","wphi")
names(ub) = c("v","wmu","wbeta","wxa","wxb","wxc","wphi")

log.posterior = function(x, data) {

    psi = as.list(x)
    names(psi) = c("v","wmu","wbeta","wxa","wxb","wxc","wphi")

    model = fit.model(psi, data, n, k, p, m0, C0, omega)

    zz = sum(logLik.ekf(model)) + log.prior(psi) + log.jacobian(psi)

    return(zz)

}

## Compute marginal likelihood
marginal.likelihood =
    bridge_sampler(samples = samples, log_posterior = log.posterior,
                   data = nao, lb = lb, ub = ub, cores = mc.cores)

## Save results
save(marginal.likelihood, file = "nao-null-alt-marginal-likelihood.RData")
