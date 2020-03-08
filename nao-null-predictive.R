##################################
## Posterior predictive samples ##
##################################

## Load libraries
library(coda)
library(compiler)
library(MASS)
library(parallel)
mc.cores = 4
options(mc.cores = mc.cores)

## Load trajectories
load("nao-null.RData")
load("nao-null-trajectories.RData")

## Initialisation and length of posterior predictive simulations
t0s = which(time == "1987-12-31")
ks  = length(y) - t0s

## Extract results
theta0s = trajectories[[1]][,t0s]
for (i in 2:length(trajectories))
    theta0s = rbind(theta0s, trajectories[[i]][,t0s])

## Garbage collection
rm(trajectories)
gc()

sample.model = function(j) {

    ## Initial state
    theta0 = theta0s[j,]

    ## Extract hyper-parameters
    V      = exp(parameters[j,1])
    wmu    = exp(parameters[j,2])
    wbeta  = exp(parameters[j,3])
    wxa    =     parameters[j,4]
    wxb    =     parameters[j,5]
    wxc    = exp(parameters[j,6])
    wphi   = exp(parameters[j,7])

    ## Weather variance
    wx = sqrt(wxa^2 + wxb^2) + wxc
    Wx = make.Wx(c(wx,wxa,wxb), omega)

    ## Innovation covariance
    W  = diag(c(wmu,wbeta,rep(wmu,2*kk),1,rep(0,p-1),rep(wphi,p)))

    ## Make model
    model = make.model(y, n, k, p, V, W, Wx, omega)

    ## Initialise storage
    yy = numeric(ks)

    ## Posterior predictive sampling
    at = theta0
    for (i in 1:ks) {

        ## Sample new state
        Gt = model$Gn(at, 0, t0s+i)
        Wt = model$Wn(at, 0, t0s+i)
        at = mvrnorm(1, Gt$a, tcrossprod(Wt$Wj %*% Wt$W, Wt$Wj))

        ## Sample new observation
        Ft = model$Fn(at, 0, t0s+i)
        Vt = model$Vn(at, 0, t0s+i)
        ft = rnorm(1, Ft$f, Vt$Vj * sqrt(Vt$V))

        ## Store samples
        yy[i] = ft

    } ## i

    ## Return samples
    return(yy)

}

wrapper = function(j) {

    ## Set random seed
    set.seed(j)

    jj = 0
    loop = TRUE
    while(loop & jj < 10) {
        jj = jj + 1
        zz = try(sample.model(j), silent = TRUE)
        if (is.numeric(zz)) {
            loop = FALSE
        } else {
            zz = rep(NA,ks)
        }
    }

    return(zz)

}

predictive = mclapply(1:nrow(parameters), wrapper)
predictive = simplify2array(predictive)

save(predictive, file = "nao-null-predictive.RData")
