## Metropolis-Hastings sampler
metropolis.hastings = function(n.samples, initial.values, random.seed, Sigma,
                               adaptive = FALSE, jj = 1, mu = NULL,
                               parameter.names = NULL) {

    ## Set seed
    assign(".Random.seed", random.seed, .GlobalEnv)

    ## Initialise adaptation
    if (adaptive) {
        lambda = 2.38^2/nrow(Sigma)
        M2     = jj*Sigma
    } else {
        lambda = 1
    }

    ## Initialise storage
    acceptance = numeric(n.samples)
    prior      = numeric(n.samples)
    jacobian   = numeric(n.samples)
    likelihood = numeric(n.samples)
    parameters = array(NA, dim = c(n.samples,nrow(Sigma)))
    if (length(parameter.names) == nrow(Sigma))
        colnames(parameters) = parameter.names

    ## Initial values
    psi   = initial.values
    model = fit.model(psi, y, n, k, p, m0, C0, omega)
    lls   = logLik.ekf(model)
    ll    = sum(lls, na.rm = TRUE)
    lp    = log.prior(psi)
    lj    = log.jacobian(psi)

    ## Sampling loop
    for (j in 1:n.samples) {

        ## Metropolis-Hastings
        psi.   = proposal(psi, lambda*Sigma)
        model. = fit.model(psi., y, n, k, p, m0, C0, omega)
        lls.   = logLik.ekf(model.)
        ll.    = sum(lls., na.rm = TRUE)
        lp.    = log.prior(psi.)
        lj.    = log.jacobian(psi.)
        a      = min(1, exp(ll. + lp. + lj. - ll - lp - lj))
        u      = runif(1)
        if (u < a) {

            model = model.
            psi   = psi.
            lls   = lls.
            ll    = ll.
            lp    = lp.
            lj    = lj.

        }

        ## Store results
        acceptance[j ] = a
        likelihood[j ] = ll
        prior     [j ] = lp
        jacobian  [j ] = lj
        parameters[j,] = unlist(psi)

        ## Update proposal
        if (adaptive) {
            jj     = jj + 1
            delta  = unlist(psi) - mu
            mu     = mu + delta/jj
            delta2 = unlist(psi) - mu
            M2     = M2 + tcrossprod(delta,delta2)
            Sigma  = M2/jj
        }

    } ## j

    ## Return results
    return(list(acceptance = acceptance, likelihood = likelihood,
                prior = prior, jacobian = jacobian,
                parameters = parameters, random.seed = .Random.seed,
                mu = mu, Sigma = Sigma))

} ## metropolis.hastings

## Wrapper for sampling
sampler = function(i, n.samples, initial.values, random.seeds, Sigma, ...)
    metropolis.hastings(n.samples, initial.values[[i]], random.seeds[[i]],
                        Sigma, ...)

## Initialise storage
parameters = list()
acceptance = list()
jacobian   = list()
likelihood = list()
prior      = list()

## Initialise chains
initial.values = list()
random.seeds   = list()
for (i in 1:n.chains) {
    set.seed(i)
    initial.values[[i]] = r.prior()
    random.seeds  [[i]] = .Random.seed
} ## i

## Continue sampling until approximate convergence
rho = 0
not.converged = TRUE
while (not.converged) {

    ## Increment counter
    rho = rho + 1
    print(rho)

    ## Sampling
    buffer = mclapply(1:n.chains, sampler, n.samples = chunk.size,
                      initial.values = initial.values,
                      random.seeds = random.seeds, Sigma = Sigma0,
                      parameter.names = parameter.names)

    ## Extract results
    parameters.buffer = lapply(buffer, function(x) mcmc(x$parameters))
    parameters.buffer = as.mcmc.list(parameters.buffer)
    acceptance.buffer = lapply(buffer, function(x) x$acceptance)
    likelihood.buffer = lapply(buffer, function(x) x$likelihood)
    prior.buffer      = lapply(buffer, function(x) x$prior)
    jacobian.buffer   = lapply(buffer, function(x) x$jacobian)
    random.seeds      = lapply(buffer, function(x) x$random.seed)

    ## Store results
    parameters[[rho]] = parameters.buffer
    acceptance[[rho]] = acceptance.buffer
    likelihood[[rho]] = likelihood.buffer
    prior     [[rho]] = prior.buffer
    jacobian  [[rho]] = jacobian.buffer

    ## Update initial values
    initial.values = lapply(parameters.buffer,
                            function(x) as.list(x[chunk.size,]))

    ## Check approximate convergence
    gelman = gelman.diag(parameters.buffer, autoburnin = FALSE)
    not.converged = any(gelman$psrf[,1] > 2.00)

} ## not.converged


## Adaptive Sampling
not.converged = TRUE
while (not.converged) {

    ## Extract proposal distribution
    mu    = apply(simplify2array(parameters.buffer), 2, mean)
    Sigma = as.matrix(parameters.buffer[[1]])
    for (i in 2:n.chains)
        Sigma = rbind(Sigma,as.matrix(parameters.buffer[[i]]))
    mu    = apply(Sigma, 2, mean)
    Sigma = cov(Sigma)
    jj    = as.integer(min(effectiveSize(parameters.buffer)))

    ## Increment counter
    rho = rho + 1
    print(rho)

    ## Sample
    buffer = mclapply(1:n.chains, sampler, n.samples = chunk.size,
                      initial.values = initial.values,
                      random.seeds = random.seeds, Sigma = Sigma,
                      adaptive = TRUE, jj = jj, mu = mu,
                      parameter.names = parameter.names)

    ## Extract results
    parameters.buffer = lapply(buffer, function(x) mcmc(x$parameters))
    parameters.buffer = as.mcmc.list(parameters.buffer)
    acceptance.buffer = lapply(buffer, function(x) x$acceptance)
    likelihood.buffer = lapply(buffer, function(x) x$likelihood)
    prior.buffer      = lapply(buffer, function(x) x$prior)
    jacobian.buffer   = lapply(buffer, function(x) x$jacobian)
    random.seeds      = lapply(buffer, function(x) x$random.seed)

    ## Store results
    parameters[[rho]] = parameters.buffer
    acceptance[[rho]] = acceptance.buffer
    likelihood[[rho]] = likelihood.buffer
    prior     [[rho]] = prior.buffer
    jacobian  [[rho]] = jacobian.buffer

    ## Update initial values
    initial.values = lapply(parameters.buffer,
                            function(x) as.list(x[chunk.size,]))

    ## Check convergence
    gelman = gelman.diag(parameters.buffer, autoburnin = FALSE)
    not.converged = any(gelman$psrf[,1] > 1.10)

} ## not.converged

## Extract proposal distribution
Sigma = as.matrix(parameters.buffer[[1]])
for (i in 2:n.chains)
    Sigma = rbind(Sigma,as.matrix(parameters.buffer[[i]]))
Sigma = cov(Sigma)*2.38^2/ncol(Sigma)

## Final sampling
sample.size = rep(0,nrow(Sigma))
blocks = 0
while (any(sample.size < n.samples)) {

    ## Increment counter
    rho = rho + 1
    blocks = blocks + 1
    print(rho)

    ## Sample
    buffer = mclapply(1:n.chains, sampler, n.samples = chunk.size,
                      initial.values = initial.values,
                      random.seeds = random.seeds, Sigma = Sigma,
                      parameter.names = parameter.names)

    ## Extract results
    parameters.buffer = lapply(buffer, function(x) mcmc(x$parameters))
    parameters.buffer = as.mcmc.list(parameters.buffer)
    acceptance.buffer = lapply(buffer, function(x) x$acceptance)
    likelihood.buffer = lapply(buffer, function(x) x$likelihood)
    prior.buffer      = lapply(buffer, function(x) x$prior)
    jacobian.buffer   = lapply(buffer, function(x) x$jacobian)
    random.seeds      = lapply(buffer, function(x) x$random.seed)

    ## Store results
    parameters[[rho]] = parameters.buffer
    acceptance[[rho]] = acceptance.buffer
    likelihood[[rho]] = likelihood.buffer
    prior     [[rho]] = prior.buffer
    jacobian  [[rho]] = jacobian.buffer

    ## Update initial values
    initial.values = lapply(parameters.buffer,
                            function(x) as.list(x[chunk.size,]))

    ## Check sample size
    sample.size = sample.size + effectiveSize(parameters.buffer)

} ## sample.size

## Save output
save.image(file = file.name)
