## Load libraries
library(compiler)
library(coda)
library(parallel)
library(triangle)
mc.cores = 4
options(mc.cores = mc.cores)

## MCMC setup
chunk.size = 1000
n.chains   = 4
n.samples  = 1000
file.name  = "nao-null-alt.RData"

## Load data
data = read.csv("era-interim.csv")
nao  = data$nao
time = as.Date(data$time)
rm(data)

y = nao
days = match(format(time, "%m-%d"), format(time[1:365], "%m-%d"))
days[is.na(days)] = 59.5
omega = 2*pi/365.25 ## Observation frequency

## Load source
source("make-model.R")
source("ekf.R")


#####################
## Model structure ##
#####################

n  = 2      ## Order of polynomial growth model for mean level
k  = c(1,2) ## Harmonics to include in seasonal component
p  = 5      ## Order of TVAR residual component

## Dimensions
tt = length(y)
kk = length(k)

## State prior mean and variance
m0 = c(6, 0, 3.6, 1.0, 1.3, 0.7, rep(0,p), c( 1.8,-1.3, 0.7,-0.3, 0.1))
C0 = diag(c(1, 2e-3, 1, 1.5, 0.9, 1.3, rep(10,p), rep(0.2,p))^2)
dim(m0) = c(length(m0),1)

## TVAR innovation variance function
make.Wx = function(x,omega)
    Wx = function(t)
        x[1] + x[2]*sin(omega*t) + x[3]*cos(omega*t)

## Parameter names
parameter.names = c("v","wmu","wbeta","wxa","wxb","wxc","wphi")

## Initial proposal distribution
Sigma0 = diag(c(1,1,1,1e-3,1e-3,1e-3,1))

## Evaluate log prior
log.prior = function (psi) {

    v      = dnorm(psi$v     , 0,  8, log = TRUE)
    wmu    = dnorm(psi$wmu   , 0,  9, log = TRUE)
    wbeta  = dnorm(psi$wbeta , 0, 17, log = TRUE)
    wphi   = dnorm(psi$wphi  , 0, 12, log = TRUE)

    return(sum(v,wmu,wbeta,wphi))

}
log.prior = cmpfun (log.prior)

## Sample initial conditions from prior
r.prior = function() {

    v      = rnorm(1, -10, 3)
    wmu    = rnorm(1, -12, 3)
    wbeta  = rnorm(1, -28, 3)
    wxa    = rnorm(1, 0.5, 1)
    wxb    = rnorm(1, 2.0, 1)
    wxc    = rnorm(1, 0.0, 1)
    wphi   = rnorm(1, -18, 3)

    return(list(v = v, wmu = wmu, wbeta = wbeta, wxa = wxa, wxb = wxb,
                wxc = wxc, wphi = wphi))

}

## Fit model
fit.model = function (psi, y, n, k, p, m0, C0, omega) {

    ## Extract hyper-parameters
    V      = exp(psi$v     )
    wmu    = exp(psi$wmu   )
    wbeta  = exp(psi$wbeta )
    wxa    = psi$wxa
    wxb    = psi$wxb
    wxc    = exp(psi$wxc)
    wphi   = exp(psi$wphi  )

    ## Weather variance
    wx = sqrt(wxa^2 + wxb^2) + wxc
    Wx = make.Wx(c(wx,wxa,wxb), omega)

    ## Innovation covariance
    W  = diag(c(wmu,wbeta,rep(wmu,2*kk),1,rep(0,p-1),rep(wphi,p)))

    model = make.model(y, n, k, p, V, W, Wx, omega)
    model = filter.ekf(y, model$Fn, model$Gn, model$Vn, model$Wn, m0, C0)

    return(model)

}
fit.model = cmpfun (fit.model)

## Evaluate log jacobian
log.jacobian = function (psi) {

    0

}
log.jacobian = cmpfun (log.jacobian)

## Sample from proposal distribution
proposal = function (psi, sigma) {

    psi. = mvrnorm(1, unlist(psi), sigma)

    v      = psi.[1]
    wmu    = psi.[2]
    wbeta  = psi.[3]
    wxa    = psi.[4]
    wxb    = psi.[5]
    wxc    = psi.[6]
    wphi   = psi.[7]

    return(list(v = v, wmu = wmu, wbeta = wbeta, wxa = wxa, wxb = wxb,
                wxc = wxc, wphi = wphi))

}
proposal = cmpfun (proposal)

## Metropolis-Hastings
source("metropolis-hastings.R")
