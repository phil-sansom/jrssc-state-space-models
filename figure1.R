##############################################
## Figure 1: The North Atlantic Oscillation ##
##############################################

## Load data
data = read.csv("era-interim.csv")
nao  = data$nao
time = as.Date(data$time)
rm(data)

## Plotting parameters
width  = 2*33*12/72.27
height = width*9/16

## Set filtering bandwidth
h = 90

## x axis
xlim = as.Date(c("1979-01-01","2018-01-01"))
xat  = as.Date(c("1980-01-01","1985-01-01","1990-01-01",
                 "1995-01-01","2000-01-01","2005-01-01",
                 "2010-01-01","2015-01-01"))
xlab = seq(1980,2015,5)

## Inter-annual SD
months = format(time, "%b")
years  = format(time, "%Y")
monthly.means = tapply(nao, list(months, years), mean, na.rm = TRUE)
monthly.means = monthly.means[month.abb,]
monthly.sds = apply(monthly.means, 1, sd, na.rm = TRUE)

## ACF
source("acf.R")
tt = length(nao)
omega = 2*pi/365.25
c1 = cos(  omega*1:tt)
s1 = sin(  omega*1:tt)
c2 = cos(2*omega*1:tt)
s2 = sin(2*omega*1:tt)
naos = residuals(lm(nao ~ c1 + s1 + c2 + s2))
month = format(time, "%b")
season = format(time, "%m")
season[! season %in% c("12","01","02","03")] = "jja"
season[  season %in% c("12","01","02","03")] = "djf"
acfm = myacf(naos, 30, month)
acfs = myacf(naos, 30, season)


## Open eps device
postscript(file = "figure1.eps", width = width, height = height,
           horizontal = FALSE, pointsize = 16, paper = "special")

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2,2.5,0.25,0.5), mgp = c(2.5,1/3,0),
    ps = 16, tcl = -1/3)

layout(matrix(c(1,2,1,3),2,2), heights = c(0.6,0.4))

## Time series
plot (time, nao, type = "l", xlim = xlim, xaxt ="n", ylim = c(-20,+35),
      col = "darkgrey", xaxs = "i", yaxs = "i")
lines(time, filter(nao, rep(1/h,h), sides = 2), col = "black", lwd = 1)
axis (side = 1, at = xat, labels = xlab)
title(xlab = "Date", line = 1)
title(ylab = "NAO Index (hPa)", line = 1.5)
# mtext("(a)", side = 3, adj = -0.08, padj = 1, cex = 0.8)

## Inter-annual SD
plot(1:12, monthly.sds, type = "b", lwd = 1, col = "black", pch = 19,
     xaxt = "n", ylim = c(0,7), yaxs = "i", yaxt = "n")
axis(side = 1, at = seq(1,12,2),  labels = month.abb[seq(1,12,2)])
axis(side = 1, at = seq(2,12,2),  labels = month.abb[seq(2,12,2)])
axis(side = 2, at = seq(0,7,1))
title(xlab = "Month", line = 1.)
title(ylab = "Inter-annual SD (hPa)", line = 1.5)
# mtext("(b)", side = 3, adj = -0.12, padj = 1, cex = 0.8)

## ACF
plot(0:30, acfm$acf[,1], type = "n", xaxs = "i", yaxs = "i",
     ylim = c(-0.2,+1.0), las = 1, xaxt = "n", cex = 1.5)
axis(1, seq(0,30,5))
title(xlab = "Lag (days)", line = 1.)
title(ylab =  "Autocorrelation function", line = 1.5)
abline(h = 0, lty = "dashed")
for (i in 1:12)
    lines(0:30, acfm$acf[,i], col = "darkgrey")
lines(0:30, acfs$acf[,1], type = "b", pch = 16)
lines(0:30, acfs$acf[,2], type = "b", pch = 17)
legend("topright", legend = c("Dec-Mar","Apr-Nov"), col = c(1,1),
       lty = c(1,1), pch = c(16,17), bty = "n")
# mtext("(c)", side = 3, adj = -0.12, line = -0.5)

dev.off()
