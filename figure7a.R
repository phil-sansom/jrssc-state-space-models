#########################################
## Figure 7a: Correlation of forecasts ##
#########################################

## Load libraries
library(coda)
library(compiler)
library(MASS)
library(parallel)

## Load posterior samples
load("nao-mean.RData")

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

## Earliest forecast start times
tas = which(format(time, "%m-%d") == "10-31" &
                "1987-01-01" < time & time < "2017-01-01")

## Forecast function
fun = function(j) {

    ## Set seed
    set.seed(j)

    ## Fit model
    psi = as.list(parameters[j,])
    names(psi) = c("v","wmu","wbeta","wxa","wxb","wxc","wphi",
                   "wdelta","pdelta","cdate","clen","cprop")
    model = fit.model(psi, y, n, k, p, m0, C0, omega)

    ## Forecast
    forecasts = list()
    for (i in 0:92) {
        buffer = matrix(NA, 120-i, length(tas))
        for (j in 1:length(tas)) {
            buffer[,j] = sim.ekf(model, tas[j]+i, 120-i, 1)$f
        } ## i
        forecasts[[i+1]] = buffer
    } ## j

    return(forecasts)

}

## Parallel forecasting
correlation.forecasts = mclapply(1:1000, fun, mc.cores = 4)

correlation.forecasts = lapply(correlation.forecasts,
                               function(x) sapply(x, function(z) apply(z, 2, mean)))
correlation.forecasts = simplify2array(correlation.forecasts)
correlation.forecasts = apply(correlation.forecasts, c(1,2), mean)

yy = matrix(NA, 30, 93)
tas = which(format(time, "%m-%d") == "10-31" &
                "1987-01-01" < time & time < "2017-01-01")
for (i in 1:30)
    for (j in 0:92)
        yy[i,j+1] = mean(nao[tas[i] + j:120])


cors = numeric(93)
for (j in 0:92)
    cors[j+1] = cor(yy[,j+1],correlation.forecasts[,j+1])

## Plotting parameters
width  = 33*12/72.27/2
height = width*9/16

xat = which(format(time[tas[1]+0:92+1], "%d") == "01")

## Open eps device
postscript(file = "figure7a.eps",
           width = width, height = height, horizontal = FALSE,
           pointsize = 8, paper = "special")

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2,2.5,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

plot(cors, type = "p", pch = 16, xlim = c(1,93), ylim = c(0.0,0.6),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
axis(side = 1, at = xat,  labels = c("Nov","Dec","Jan","Feb"))
axis(side = 2, at = seq(0.0,0.6,0.1), hadj = 1.2)
title(xlab = "Date", line = 1)
title(ylab = "Correlation", line = 1.5)

## Close device
dev.off()
