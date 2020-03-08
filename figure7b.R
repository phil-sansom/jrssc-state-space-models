###################################
## Figure 7b: Seasonal forecasts ##
###################################

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

## Forecast times
t0s = which(format(time, "%m-%d") == "11-30" & time >= "1987-11-30")
t1s = which(format(time, "%m-%d") == "03-01" & time >= "1987-11-30") - 1
t0s = t0s[-length(t0s)]
ks  = t1s - t0s

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
    forecasts = array(NA, c(max(ks),length(t0s)))
    for (i in 1:length(t0s))
        forecasts[1:ks[i],i] = sim.ekf(model, t0s[i], ks[i], 1)$f

    ## Return forecasts
    return(forecasts)

}

## Parallel forecasting
seasonal.forecasts = mclapply(1:1000, fun, mc.cores = 4)
seasonal.forecasts = simplify2array(seasonal.forecasts)
dimnames(seasonal.forecasts) = list(NULL, year = format(time[t0s],"%Y"), NULL)

means = apply(seasonal.forecasts, 2, median, na.rm = TRUE)
lwr   = apply(seasonal.forecasts, 2, quantile, probs = 0.05, na.rm = TRUE)
upr   = apply(seasonal.forecasts, 2, quantile, probs = 0.95, na.rm = TRUE)

t0s = which(format(time, "%m-%d") == "11-30" & time >= "1987-11-30")
t1s = which(format(time, "%m-%d") == "03-01" & time >= "1987-11-30") - 1
t0s = t0s[-length(t0s)]
ks  = t1s - t0s

data.means = sapply(1:30, function(x) mean(nao[t0s[x]+1:ks[x]]))

xx = as.numeric(names(means))

## Plotting parameters
width  = 33*12/72.27/2
height = width*9/16

## Open eps device
postscript(file = "figure7b.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2,2.5,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

plot(xx, data.means, type = "b", xlim = c(1986,2017), ylim = c(-5,+20),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", pch = 4)
lines(xx, means, type = "b", pch = 16, lwd = 2)
axis(side = 1, at = seq(1990,2015,5))
axis(side = 2, at = seq(-5,+20,+5), hadj = 1.4)
title(xlab = "Year", line = 1)
title(ylab = "NAO Index (hPa)", line = 1.5)

## Close device
dev.off()
