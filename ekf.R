############################
## Extended Kalman Filter ##
############################

## Load libraries
library(compiler)
library(MASS)

## Forward filtering
step.ekf = function(y, Fn, Gn, Vn, Wn, m0, C0, t) {

    ## Mean vector
    zero = rep(0, nrow(C0))

    ## Single update
    Gt = Gn(m0  , zero, t)
    Wt = Wn(m0  , zero, t)
    Ft = Fn(Gt$a, 0   , t)
    Vt = Vn(Gt$a, 0   , t)

    R = tcrossprod(Gt$Gj %*% C0  , Gt$Gj) +
        tcrossprod(Wt$Wj %*% Wt$W, Wt$Wj)
    Q = tcrossprod(Ft$Fj %*% R, Ft$Fj) + Vt$Vj^2 * Vt$V

    A = tcrossprod(R, Ft$Fj) / drop(Q)
    e = y - Ft$f

    m = Gt$a + A * drop(e)
    C = R - tcrossprod(A) * drop(Q)

    logLik = dnorm(y, Ft$f, sqrt(Q), TRUE)

    ## Return results
    results = list(a = Gt$a, R = R, f = Ft$f, Q = Q, m = m, C = C,
                   logLik = logLik)
    return(results)

}
step.ekf = cmpfun(step.ekf)


## Forward filtering
filter.ekf = function(y, Fn, Gn, Vn, Wn, m0, C0, t0 = 0, t1 = length(y),
                      keep.pred = TRUE, keep.jac = TRUE) {

    ## Dimensions
    tt = t1 - t0   ## Length of observation vector
    nn = nrow(m0)  ## Length of state vector

    ## Create arrays to store results
    a  = array(NA, dim = c(nn, 1,tt))  ## state prior mean
    R  = array(NA, dim = c(nn,nn,tt))  ## state prior scale

    if (keep.pred) {
        f  = array(NA, dim = c( 1, 1,tt))  ## forecast mean
        Q  = array(NA, dim = c( 1, 1,tt))  ## forecast scale
        m  = array(NA, dim = c(nn, 1,tt))  ## state posterior mean
        C  = array(NA, dim = c(nn,nn,tt))  ## state posterior scale
    }

    if (keep.jac) {
        Gj = array(NA, dim = c(nn,nn,tt))  ##
        Wj = array(NA, dim = c(nn,nn,tt))  ##
        Fj = array(NA, dim = c( 1,nn,tt))  ##
        Vj = array(NA, dim = c( 1, 1,tt))  ##
    }

    ## Mean vector
    zero = rep(0, nn)

    ## Initialisation
    mt = m0
    Ct = C0

    ## Loop over time
    for (t in 1:tt) {

        ## Prediction step
        Gt = Gn(mt  , zero, t0 + t)
        Wt = Wn(mt  , zero, t0 + t)
        Ft = Fn(Gt$a, 0   , t0 + t)
        Vt = Vn(Gt$a, 0   , t0 + t)

        Rt = tcrossprod(Gt$Gj %*% Ct  , Gt$Gj) +
             tcrossprod(Wt$Wj %*% Wt$W, Wt$Wj)
        Qt = tcrossprod(Ft$Fj %*% Rt, Ft$Fj) + Vt$Vj^2 * Vt$V

        ## Update step
        if (is.na(y[t])) {

            mt = Gt$a
            Ct = Rt

        } else {

            A = tcrossprod(Rt, Ft$Fj) / Qt[1,1]
            e = y[t0 + t] - Ft$f

            mt = Gt$a + A * e[1,1]
            Ct = Rt - tcrossprod(A) * Qt[1,1]

        }

        ## Store results
        m[,,t] = mt
        C[,,t] = Ct

        if (keep.pred) {
            a[,,t] = Gt$a
            R[,,t] = Rt
            f[,,t] = Ft$f
            Q[,,t] = Qt
        }

        if (keep.jac) {
            Gj[,,t] = Gt$Gj
            Wj[,,t] = Wt$Wj
            Fj[,,t] = Ft$Fj
            Vj[,,t] = Vt$Vj
        }

    } ## t

    ## Return results
    results = list(m = m, C = C, y = y, Fn = Fn, Gn = Gn, Vn = Vn, Wn = Wn,
                   m0 = m0, C0 = C0, t0 = t0, t1 = t1)
    if (keep.pred)
        results = c(results, list(a = a, R = R, f = f, Q = Q))
    if (keep.jac)
        results = c(results, list(Gj = Gj, Wj = Wj, Fj = Fj, Vj = Vj))

    return(results)

}
filter.ekf = cmpfun(filter.ekf)


## Compute log-likelihood
logLik.ekf = function(object) {

    ## Expects the following inputs:
    ## object - an object generated by filter.ekf()

    results = dnorm(object$y[(object$t0+1):object$t1],
                    mean = object$f[1,1,], sd = sqrt(object$Q[1,1,]),
                    log = TRUE)

    ## Return result
    return(results)

}
logLik.udlm = cmpfun(logLik.ekf)


## Backward smoothing
smooth.ekf = function(object, from = NULL, to = NULL) {

    ## Expects the following inputs:
    ## object - an object generated by filter.ekf()

    ## Dimensions
    t0 = object$t0
    t1 = object$t1
    tt = t1 - t0
    nn = nrow(object$m0)  ## Length of state vector

    from = min(from-t0,tt)
#    to   = max(to,t0+nn) - t0
    to   = max(1,to)

    ## Create arrays to store results
    at = array(NA, dim = c(nn, 1,tt))  ## state prior mean
    Rt = array(NA, dim = c(nn,nn,tt))  ## state prior scale matrix

    ## Mean vector
    zero = rep(0, nn)

    ## Final conditions
    at[,,from] = object$m[,,from]
    Rt[,,from] = object$C[,,from]

    ## Loop over remaining times
    for (t in (from-1):to) {

        Bt = t(solve(object$R[,,t+1],
                     tcrossprod(object$Gj[,,t+1], object$C[,,t])))

        at[,,t] = object$m[,,t] + Bt %*% (at[,,t+1] - object$a[,,t+1])
        Rt[,,t] = object$C[,,t] + Bt %*% (Rt[,,t+1] - object$R[,,t+1]) %*% t(Bt)

    }

    ## Return results
    results = list(at = at, Rt = Rt)
    return(results)

}
smooth.ekf = cmpfun(smooth.ekf)


## Backward sampling
sample.ekf = function(object, from = NULL, to = 1, n.samples = 1e3) {

    ## Dimensions
    t0 = object$t0
    tt = length(object$y)  ## Length of observation vector
    nn = nrow(object$m0)   ## Length of state vector

    from = min(from,tt)
    to   = max(to,object$t0+1)

    ## Mean vector
    zero = rep(0, nn)

    ## Create arrays to store results
    theta = array(NA, dim = c(n.samples,nn,tt))

    ## Sample final state
    theta[,,from] = mvrnorm(n.samples, object$m[,,from], object$C[,,from])

    ## Loop over remaining times
    for (t in (from-1):to) {

        B  = t(solve(object$R[,,t+1],
                     tcrossprod(object$Gj[,,t+1], object$C[,,t])))
        h  = t(object$m[,,t] + B %*% (t(theta[,,t+1]) - object$a[,,t+1]))
        H  = object$C[,,t] - tcrossprod(B %*% object$R[,,t+1], B)

        theta[,,t] = h + mvrnorm(n.samples, zero, H)

    } ## t

    ## Return results
    return(theta)

}
sample.ekf = cmpfun(sample.ekf)


## Backward sampling
sample0.ekf = function(object, from = NULL, to = 1) {

    ## Dimensions
    t0 = object$t0
    tt = length(object$y)  ## Length of observation vector
    nn = nrow(object$m0)   ## Length of state vector

    from = min(from,tt)
#    to   = max(to,object$t0+nn)
    to   = max(to,object$t0+1)

    ## Create arrays to store results
    theta = array(NA, dim = c(nn,tt))

    ## Sample final state
    theta[,from] = mvrnorm(1, object$m[,,from], object$C[,,from])

    ## Loop over remaining times
    for (t in (from-1):to) {

        B  = t(solve(object$R[,,t+1],
                     tcrossprod(object$Gj[,,t+1], object$C[,,t])))
        h  = object$m[,,t] + B %*% (theta[,t+1] - object$a[,,t+1])
        H  = object$C[,,t] - tcrossprod(B %*% object$R[,,t+1], B)

        buffer = try(mvrnorm(1, h, H), TRUE)
        if (class(buffer) == "try-error") {
            break
        } else {
            theta[,t] = buffer
        }

    } ## t

    ## Return results
    return(theta)

}
sample0.ekf = cmpfun(sample0.ekf)


## Prediction
predict.ekf = function(object, t, k = 1) {

    ## t - reference time
    ## k - lead time

    ## Dimensions
    nn = nrow(object$m0)

    ## Create arrays to hold forecasts
    a = array(NA, c(nn,k))
    R = array(NA, c(nn,nn,k,k))

    f = numeric(k)
    Q = array(NA, c(k,k))

    Gj = array(NA, dim = c(nn,nn,k))
    Fj = array(NA, dim = c( 1,nn,k))

    ## Mean vector
    zero = rep(0, nn)

    ## Initialise forecasts
    i = 1

    ## Forecast new state
    Gt = object$Gn(object$m[,,t], zero, t+i)
    Wt = object$Wn(object$m[,,t], zero, t+i)
    Ct = object$C[,,t]
    at = Gt$a
    Rt = tcrossprod(Gt$Gj %*% Ct, Gt$Gj) + tcrossprod(Wt$Wj %*% Wt$W, Wt$Wj)

    ## Forecast new observation
    Ft = object$Fn(at, 0, t+i)
    Vt = object$Vn(at, 0, t+i)
    f[i]   = Ft$f
    Q[i,i] = tcrossprod(Ft$Fj %*% Rt, Ft$Fj) + Vt$Vj^2 * Vt$V

    ## Store results
    a[,i]    = at
    R[,,i,i] = Rt

    ## Loop over remaining forecasts
    if (k > 1) {

        ## Fill in forecast covariances
        Gj[,,i] = Gt$Gj; Fj[,,i] = Ft$Fj

        for (i in 2:k) {

            ## Forecast new state
            Gt = object$Gn(a[,i-1], zero, t+i)
            Wt = object$Wn(a[,i-1], zero, t+i)
            at = Gt$a
            Rt = tcrossprod(Gt$Gj %*% Rt  , Gt$Gj) +
                 tcrossprod(Wt$Wj %*% Wt$W, Wt$Wj)

            ## Forecast new observation
            Ft = object$Fn(at, 0, t+i)
            Vt = object$Vn(at, 0, t+i)
            f[i]   = Ft$f
            Q[i,i] = tcrossprod(Ft$Fj %*% Rt, Ft$Fj) + Vt$Vj^2 * Vt$V

            ## Store results
            a[,i]    = at
            R[,,i,i] = Rt

            ## Fill in forecast covariances
            Gj[,,i] = Gt$Gj; Fj[,,i] = Ft$Fj
            for (j in 1:(i-1)) {

                ## State covariances
                R[,,i,j] = Gj[,,i] %*% R[,,i-1,j]
                R[,,j,i] = t(R[,,i,j])

                ## Observation covariances
                Q[i,j] = t(Fj[,,i]) %*% R[,,i,j] %*% Fj[,,j]
                Q[j,i] = t(Q[i,j])

            } ## j

        } ## i

    } ## k > 1

    ## Return results
    results = list(f = f, Q = Q, a = a, R = R, t = t, k = k)
    return(results)

}
predict.ekf = cmpfun(predict.ekf)


## Simulation
sim.ekf = function(object, t, k = 1, n.samples = 1e3) {

    ## t - reference time
    ## k - lead time

    ## Dimensions
    nn = nrow(object$m0)  ## Length of state vector

    ## Create arrays to hold forecasts
    a = array(0, dim = c(nn,n.samples,k))
    f = array(0, dim = c( 1,n.samples,k))

    ## Mean vector
    zero = rep(0, nn)

    for (j in 1:n.samples) {

        ## Sample initial state
        at = mvrnorm(1, object$m[,,t], object$C[,,t])

        for (i in 1:k) {

            ## Sample new state
            Gt = object$Gn(at, zero, t+i)
            Wt = object$Wn(at, zero, t+i)
            at = mvrnorm(1, Gt$a, tcrossprod(Wt$Wj %*% Wt$W, Wt$Wj))

            ## Sample new observation
            Ft = object$Fn(at, 0, t+i)
            Vt = object$Vn(at, 0, t+i)
            ft = rnorm(1, Ft$f, Vt$Vj * sqrt(Vt$V))

            ## Store samples
            a[,j,i] = at
            f[,j,i] = ft

            atm1 = at

        } ## i

    } ## j

    ## Return results
    results = list(a = a, f = f, t = t, k = k)
    return(results)

}
sim.ekf = cmpfun(sim.ekf)
