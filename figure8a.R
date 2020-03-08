#####################
## What happenned? ##
#####################

## Load trajectories
load("nao-mean.RData")
load("nao-mean-trajectories.RData")

## Extract parameters and trajectories
trajectories = simplify2array(trajectories)
gc()

eta = trajectories[1,,] + trajectories[3,,] + trajectories[5,,]
eta = apply(eta, 1, median, na.rm = TRUE)

fit = apply(trajectories[17,,], 1, median, na.rm = TRUE)
lwr = apply(trajectories[17,,], 1, quantile, probs = 0.05, na.rm = TRUE)
upr = apply(trajectories[17,,], 1, quantile, probs = 0.95, na.rm = TRUE)

## Plotting parameters
width  = 33*12/72.27/2
height = width*9/16

## Open eps device
postscript(file = "figure8a.eps",
           width = width, height = height, horizontal = FALSE,
           pointsize = 8, paper = "special")

xlim = as.Date(c("2009-11-01","2010-04-01"))
xat  = xlim[1] + which(format(seq(xlim[1],xlim[2],by="day"), "%d") == "01") - 1

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

plot(time, y - eta, type = "n", xlim = xlim, ylim = c(-25,+15),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", lwd = 2)
for (i in sample.int(1e3,1e1))
    lines(time, trajectories[17,,i], col = "lightgrey")
abline(h = 0, lty = "dashed")
lines(time, fit, lwd = 2, col = "darkgrey")
lines(time, lwr, col = "darkgrey", lty = "dashed")
lines(time, upr, col = "darkgrey", lty = "dashed")
lines(time, y - eta, lwd = 2)
axis(side = 1, at = xat, labels = c("Nov","Dec","Jan","Feb","Mar","Apr"))
axis(side = 2, at = seq(-25,+15,+5), hadj = 1, mgp = c(2.5,1/2,0))
title(xlab = "Date", line = 1)
title(ylab = "NAO Index (hPa)", line = 2.0)

## Close device
dev.off()
