############################
## Posterior Traces: Mean ##
############################

## Load libraries
library(coda)
library(triangle)

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
samples = samples[(adaptation+1):nrow(samples),,]

## Device width
width  = 33*12/72.27
height = 4*width/3*9/16

priors = list(
    v      = rnorm(1e6, -10, 3),
    wmu    = rnorm(1e6, -12, 3),
    wbeta  = rnorm(1e6, -28, 3),
    wxa    = rnorm(1e6, 0.5, 1),
    wxb    = rnorm(1e6,   2, 1),
    wxc    = rnorm(1e6,   0, 1),
    wphi   = rnorm(1e6, -18, 3),
    wdelta = rnorm(1e6,  -8, 4),
    pdelta = rbeta(1e6,   4, 1),
    cdate  = floor(rtriangle(1e6, 120, 365 + 120, 305)%%365) + 1,
    clen   = floor(rtriangle(1e6, 0, 365, 180)),
    cprop  = rbeta(1e6,   4, 6)
)


xlabs = list("log(V)",expression(log(W[mu])),expression(log(W[beta])),
             "a","b","log(c)",
             expression(log(W[phi])), expression(log(W[delta])),
             expression(phi[delta]),
             expression(paste("Start date ", alpha)),
             expression(paste("Length ", gamma, " (days)")),
             expression(paste("Tapering ", rho)))
xlims = list(c(-19,+1), c(-21,-3), c(-37,-19),c(-1.5,+2.5),c(0,4),c(-2,+2),
             c(-27,-9),c(-16,0),c(0.9,1.0),c(0,365),c(0,365), c(0,1))
xbreaks = list(seq(floor(min(samples[,1,])),ceiling(max(samples[,1,])),+1),
               seq(floor(min(samples[,2,])),ceiling(max(samples[,2,])),+1),
               seq(floor(min(samples[,3,])),ceiling(max(samples[,3,])),+1),
               seq(-1.5,+2.5,+0.05), seq( 0.0,+4.0,+0.05), seq(-2.0,+2.0,+0.05),
               seq(floor(min(samples[,7,])),ceiling(max(samples[,7,])),+1.0),
               seq(floor(min(samples[,8,])),ceiling(max(samples[,8,])),+0.1),
               seq(0.9,+1.0,+0.001), seq(0,365,10), seq(0,365,10),
               seq(0.0,1.0,0.05))

## Open device
postscript(file = "posterior-mean.eps", width = height, height = width,
           horizontal = FALSE, pointsize = 8, paper = "special")

par(las = 1, mar = c(2.5,3,0.5,0.5), mgp = c(1.5,0.5,0.0),
    ps = 8, tcl= -1/3, xaxs = "i")

layout(matrix(1:12, 4, 3, byrow = TRUE))

for (i in 1:12) {

  hist(samples[,i,], breaks = xbreaks[[i]], freq = FALSE, col = "darkgrey",
       main = NULL, xlim = xlims[[i]], xlab = xlabs[[i]], ylab = NULL)
  lines(density(priors[[i]]), lwd = 1)
    title(ylab = "Density", line = 2)

}

dev.off()

