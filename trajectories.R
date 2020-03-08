## Load samples
load(input.file)

## Determine when adaptation ended
psrfs = sapply(parameters, function(x) max(gelman.diag(x, autoburnin = FALSE)$psrf[,1]))
alpha = min(which(psrfs < 1.10)) + 1
gamma = length(parameters)
blocks = gamma - alpha + 1

## Extract parameters
buffers = list()
for (i in 1:n.chains) {
    buffers[[i]] = parameters[[alpha]][[i]]
    for (j in (alpha + 1):gamma) {
        buffers[[i]] = rbind(buffers[[i]],parameters[[j]][[i]])
    } ## j
    buffers[[i]] = buffers[[i]][seq(n.chains*blocks,nrow(buffers[[i]]),n.chains*blocks),]
} ## i
parameters = buffers[[1]]
for (i in 2:n.chains)
    parameters = rbind(parameters,buffers[[i]])

# ## Extract parameters
# buffer = matrix(NA, 0, ncol(Sigma))
# for (i in alpha:gamma)
#     for (j in 1:n.chains)
#         buffer = rbind(buffer, parameters[[i]][[j]])

## Sampling function
fun = function(i) {

    set.seed(i)
    model = fit.model(as.list(parameters[i,]), y, n, k, p, m0, C0, omega)
    traj  = sample0.ekf(model)
    # mask  = 1:nrow(traj)
    # if (nrow(traj) > n+2*kk+2*p+1)
    #     mask = mask[-(n+2*kk+2*p+1)]
    # mask  = mask[-(n+2*kk+2:p)]
    # traj = traj[mask,]
    return(traj)

}

## Parallel sampling
trajectories = mclapply(1:nrow(parameters), fun)

## Save trajectories
save(parameters, trajectories, file = output.file)
