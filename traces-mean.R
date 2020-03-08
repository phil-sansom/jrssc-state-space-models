## Load libraries
library(coda)

## Load data
load("nao-mean.RData")

## PSRFs
psrfs = sapply(parameters, function(x) max(gelman.diag(x, autoburnin = FALSE)$psrf[,1]))
burnin     = min(which(psrfs < 2.0))*chunk.size
adaptation = min(which(psrfs < 1.1))*chunk.size

## Extract parameters
samples = array(NA, c(chunk.size*length(parameters),nrow(Sigma0),n.chains))
for (i in 1:length(parameters))
    for(j in 1:n.chains)
        samples[(i-1)*1000 + 1:1000,,j] = parameters[[i]][[j]]
samples[, 9,] = 1/(1+exp(-samples[, 9,]))
samples[,12,] = 1/(1+exp(-samples[,12,]))

ylabs = list("log(V)",expression(log(W[mu])),expression(log(W[beta])),
             "a","b","log(c)",
             expression(log(W[phi])), expression(log(W[delta])),
             expression(phi[delta]),
             expression(paste("Start date ", alpha)),
             expression(paste("Length ", gamma, " (days)")),
             expression(paste("Tapering ", rho)))

yats = apply(samples, 2, range)
yats = yats[1,] + 0.05*(yats[2,] - yats[1,])

width  = 33*12/72.27
height = 4*width/3*9/16

postscript(file = "traces-mean.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)


## Plotting parameters
par(ann = FALSE, las = 1, mar = c(1.5,2.5,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3, xaxs = "i")

layout(matrix(1:dim(samples)[2], 4, 3, byrow = TRUE))

for (i in 1:dim(samples)[2]) {
    plot(samples[,i,1], type = "l",
         xlim = c(0,nrow(samples)), ylim = range(samples[,i,]))
    title(ylab = ylabs[[i]], line = 1.5)
    for (j in 2:dim(samples)[3]) {
        lines(samples[,i,j], col = j)
    } ## j
    abline(v = burnin, lty = "dashed")
    abline(v = adaptation, lty = "dashed")
    text(0.5*burnin, yats[i], "Burn-in", cex = 0.5)
    text(0.5*(burnin + adaptation), yats[i], "Adaptation", cex = 0.5)
    text(0.5*(adaptation + nrow(samples)), yats[i], "Sampling", cex = 0.5)
} ## i

dev.off()

