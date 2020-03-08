## Load source
source("durbin-levinson.R")

## Set seed
set.seed(0)

## Simulation length
n = 30*365

## Observation uncertainty
v = 1e-1

## Mean and trend
wmu   = 1e-1
wbeta = 1e-4
beta  = rnorm(1, 0, 0.002) + cumsum(rnorm(n, 0, wbeta))
mu    = rnorm(1, 0, 5) + cumsum(beta) + cumsum(rnorm(n, 0, wmu))

## Seasonal component
omega = 2*pi/365
wpsi  = 1e-1
psi   = matrix(0, 4, n)
G     = matrix(0, 4, 4)
G[1:2,1:2] = c( cos(1*omega),-sin(1*omega), sin(1*omega), cos(1*omega))
G[3:4,3:4] = c( cos(2*omega),-sin(2*omega), sin(2*omega), cos(2*omega))
psi[,1] = rnorm(4, 0, 5)
for (i in 2:n)
    psi[,i] = G %*% psi[,i-1] + rnorm(4, 0, wpsi)

## TVAR coefficients
p    = 5
wrho = 1e-2
rho  = matrix(0, p, n)
phi  = matrix(0, p, n)
rho[,1] = runif(p, -1, +1)
phi[,1] = pacf2ar(rho[,1])
for (i in 2:n) {
    rho[,i] = rho[,i-1] + rnorm(p, 0, wrho)
    for (j in 1:p)
        rho[j,i] = sign(rho[j,i])*min(abs(rho[j,i]),1)
    phi[,i] = pacf2ar(rho[,i])
}
wphi = median(sapply(1:5, function(x) sd(diff(phi[x,]))))

## Irregular component
wx = 5
x  = rnorm(n+p, 0, wx)
for (i in 1:n)
    x[p+i] = phi[,i] %*% x[p+i-1:p] + x[p+i]

## Create influence function
days = rep(1:365,30)
cprop = 0.2
clen  = 90
cdate = 180
f  = numeric(length(days))
if (clen > 0) {
  c0 = floor(0.5*clen*cprop)
  if (c0 > 0) {
    fs = c(seq( 1/(c0+1),c0/(c0+1), 1/(c0+1)),rep(1,clen-2*c0),
           seq(c0/(c0+1), 1/(c0+1),-1/(c0+1)))
  } else {
    fs = rep(1,clen)
  }
  starts = which(days == cdate)
  for (i in starts)
    f[i-1+1:clen] = fs
  if(cdate - 1 + clen > 365)
    f[1:(cdate - 1 + clen - 365)] = fs[(365 - cdate + 2):clen]
}
make.Bp = function(f)
  function(t) f[t]
Bp = make.Bp(f)

## TVAR interventions
wdelta = 1e-3
pdelta = 1
rhod = matrix(0, p, n)
phid = matrix(0, p, n)
i = 1
rhod[,i] = rnorm(rho[,i], 0, 0.2)
for (j in 1:p)
  rhod[j,i] = sign(rhod[j,i])*min(abs(rhod[j,i]),1)
phid[,1] = pacf2ar(rhod[,1])
for (i in 2:n) {
  rhod[,i] = rhod[,i-1] + rnorm(p, 0, wdelta)
  for (j in 1:p)
    rhod[j,i] = sign(rhod[j,i])*min(abs(rhod[j,i]),1)
  phid[,i] = pacf2ar(rhod[,i])
}
phid = phid - phi
wphid = median(apply(phid, 1, function(x) sd(diff(x))))

## TVAR intervention effect
delta = numeric(n)
for (i in 1:n)
    delta[i] = phid[,i] %*% x[p+i-1:p]

x = x[p+1:n]

## Observations
y = mu + psi[1,] + psi[3,] + x + f*delta + rnorm(n, 0, v)

## Load source
source("make-model.R")
source("ekf.R")

## State prior mean and variance
m0 = c(0,0,rep(0,4),rep(0,p),0.5,rep(0,p-1),0,rep(0,p))
C0 = diag(c(5,wbeta,rep(5,4),rep(1,p),rep(1,p),10,rep(1,p))^2)
dim(m0) = c(length(m0),1)

V = matrix(v^2, 1, 1)
W = diag(c(wmu,wbeta,rep(wpsi,4),1,rep(0,p-1),rep(wphi,p),0,rep(wphid,p))^2)
Wx = function(t) wx^2

model = make.model(y, 2, c(1,2), p, V, W, Wx, omega, dp = pdelta, Bp = Bp)
model = filter.ekf(y, model$Fn, model$Gn, model$Vn, model$Wn, m0, C0)
models = smooth.ekf(model)

## Plotting parameters
width  = 33*12/72.27/2
height = width*9/16

# ## Open eps device
# postscript(file = "sim-tvar-delta.eps", width = width, height = height,
#            pointsize = 8, paper = "special", horizontal = FALSE)
#
# ## Plotting parameters
# par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
#     ps = 8, tcl = -1/3)
#
# i   = 17
# fit = models$at[i,1,]
# lwr = fit + qnorm(0.025)*sqrt(models$Rt[i,i,])
# upr = fit + qnorm(0.975)*sqrt(models$Rt[i,i,])
# plot(delta, type = "l", col = "red",  xaxt = "n",
#      xlim = c(0,n), ylim = range(delta,fit,lwr,upr), xaxs = "i", yaxs = "i")
# lines(fit, lwd = 2)
# lines(lwr, lty = "dashed")
# lines(upr, lty = "dashed")
# title(xlab = "Year", line = 1.0)
# axis(side = 1, at = seq(0,n,365), labels = 0:30)
# title(ylab = expression(paste("Effect ", delta, " (hPa)")), line = 2.0)
#
# dev.off()


## Open eps device
postscript(file = "sim-tvar-mu.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

i   = 1
fit = models$at[i,1,]
lwr = fit + qnorm(0.025)*sqrt(models$Rt[i,i,])
upr = fit + qnorm(0.975)*sqrt(models$Rt[i,i,])
plot(mu, type = "l", col = "red",  xaxt = "n",
     xlim = c(0,n), ylim = range(fit,lwr,upr), xaxs = "i", yaxs = "i")
lines(fit, lwd = 2)
lines(lwr, lty = "dashed")
lines(upr, lty = "dashed")
title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,n,365), labels = 0:30)
title(ylab = expression(paste("Mean ", mu, " (hPa)")), line = 2.0)

dev.off()

## Open eps device
postscript(file = "sim-tvar-beta.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

i   = 2
fit = 365*10*models$at[i,1,]
lwr = fit + 365*10*qnorm(0.025)*sqrt(models$Rt[i,i,])
upr = fit + 365*10*qnorm(0.975)*sqrt(models$Rt[i,i,])
plot(365*10*beta, type = "l", col = "red",  xaxt = "n",
     xlim = c(0,n), ylim = range(fit,lwr,upr), xaxs = "i", yaxs = "i")
lines(fit, lwd = 2)
lines(lwr, lty = "dashed")
lines(upr, lty = "dashed")
title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,n,365), labels = 0:30)
title(ylab = expression(paste("Trend ", beta, " (hPa/10yr)")), line = 2.0)

dev.off()

## Open eps device
postscript(file = "sim-tvar-psi1.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

i   = 3
fit = models$at[i,1,]
lwr = fit + qnorm(0.025)*sqrt(models$Rt[i,i,])
upr = fit + qnorm(0.975)*sqrt(models$Rt[i,i,])
plot(psi[1,], type = "l", col = "red",  xaxt = "n",
     xlim = c(0,n), ylim = range(fit,lwr,upr), xaxs = "i", yaxs = "i")
lines(fit, lwd = 2)
lines(lwr, lty = "dashed")
lines(upr, lty = "dashed")
title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,n,365), labels = 0:30)
title(ylab = expression(paste("Annual cycle ", psi[1], " (hPa)")), line = 2.0)

dev.off()


## Open eps device
postscript(file = "sim-tvar-psi2.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

i   = 5
fit = models$at[i,1,]
lwr = fit + qnorm(0.025)*sqrt(models$Rt[i,i,])
upr = fit + qnorm(0.975)*sqrt(models$Rt[i,i,])
plot(psi[3,], type = "l", col = "red",  xaxt = "n",
     xlim = c(0,n), ylim = range(fit,lwr,upr), xaxs = "i", yaxs = "i")
lines(fit, lwd = 2)
lines(lwr, lty = "dashed")
lines(upr, lty = "dashed")
title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,n,365), labels = 0:30)
title(ylab = expression(paste("Semi-annual cycle ", psi[2], " (hPa)")), line = 2.0)

dev.off()

## Coefficients
## Open eps device
postscript(file = "sim-tvar-phi.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

j   = 1
i   = 11+j
fit = models$at[i,1,]
lwr = fit + qnorm(0.025)*sqrt(models$Rt[i,i,])
upr = fit + qnorm(0.975)*sqrt(models$Rt[i,i,])
plot(phi[j,], type = "l", col = "red",  xaxt = "n",
     xlim = c(0,n), ylim = range(phi), xaxs = "i", yaxs = "i")
lines(fit, lwd = 2)
lines(lwr, lty = "dashed")
lines(upr, lty = "dashed")
for (j in 2:p) {
  i   = 11+j
  fit = models$at[i,1,]
  lwr = fit + qnorm(0.025)*sqrt(models$Rt[i,i,])
  upr = fit + qnorm(0.975)*sqrt(models$Rt[i,i,])
  lines(phi[j,], type = "l", col = j+1)
  lines(fit, lwd = 2)
  lines(lwr, lty = "dashed")
  lines(upr, lty = "dashed")
}

title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,n,365), labels = 0:30)
title(ylab = expression(paste("TVAR coefficients ", phi[p])), line = 2.0)

dev.off()

## Effects
## Open eps device
postscript(file = "sim-tvar-delta.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

j   = 1
i   = 17+j
fit = models$at[i,1,]
lwr = fit + qnorm(0.025)*sqrt(models$Rt[i,i,])
upr = fit + qnorm(0.975)*sqrt(models$Rt[i,i,])
plot(phid[j,], type = "l", col = "red",  xaxt = "n",
     xlim = c(0,n), ylim = range(phid), xaxs = "i", yaxs = "i")
lines(fit, lwd = 2)
lines(lwr, lty = "dashed")
lines(upr, lty = "dashed")
for (j in 2:p) {
  i   = 17+j
  fit = models$at[i,1,]
  lwr = fit + qnorm(0.025)*sqrt(models$Rt[i,i,])
  upr = fit + qnorm(0.975)*sqrt(models$Rt[i,i,])
  lines(phid[j,], type = "l", col = j+1)
  lines(fit, lwd = 2)
  lines(lwr, lty = "dashed")
  lines(upr, lty = "dashed")
}

title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,n,365), labels = 0:30)
title(ylab = expression(paste("TVAR effects ", delta[p])), line = 2.0)

dev.off()




library(parallel)
buffer = mclapply(1:4, function(x) {set.seed(x); sample.ekf(model, n.samples = 250)},
                  mc.cores = 4)
sims = array(NA, c(1000,22,n))
for (i in 1:4)
  sims[(i-1)*250 + 1:250,,] = buffer[[i]]
rm(buffer); gc()

## Intervention effect
deltaf      = f*delta
dim(deltaf) = c(365,30)
deltaf = apply(deltaf, 2, sum)/clen
means = matrix(NA, 1000, 30)
for (i in 1:1000) {
    buffer = f*sims[i,17,]
    dim(buffer) = c(365,30)
    means[i,] = apply(buffer, 2, sum)/clen
}
mylist = list()
mylist$stats = apply(means, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975))
mylist$n     = apply(means, 2, length)
mylist$names = as.character(1:ncol(means))

## Open eps device
postscript(file = "sim-tvar-delta-mean.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

bxp(mylist, outline = FALSE, xlim = c(1,30), xaxs = "i", xaxt = "n",
    ylim = c(-12.5,+12.5), yaxs = "i")
points(deltaf, col = "red", pch = 16)
title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,30,1))
title(ylab = expression(paste("Mean intervention ", delta, " (hPa)")), line = 2.0)

dev.off()


## Open eps device
postscript(file = "sim-tvar-amplitude.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

A1 = sqrt(sims[,3,]^2 + sims[,4,]^2)
A2 = sqrt(sims[,5,]^2 + sims[,6,]^2)
fit1 = apply(A1, 2, mean)
lwr1 = apply(A1, 2, quantile, probs = 0.025)
upr1 = apply(A1, 2, quantile, probs = 0.975)
fit2 = apply(A2, 2, mean)
lwr2 = apply(A2, 2, quantile, probs = 0.025)
upr2 = apply(A2, 2, quantile, probs = 0.975)
a1   = sqrt(psi[1,]^2 + psi[2,]^2)
a2   = sqrt(psi[3,]^2 + psi[4,]^2)

plot(a1, type = "l", col = "red",  xaxt = "n",
     xlim = c(0,n), ylim = range(fit1,lwr1,upr1,fit2,lwr2,upr2),
     xaxs = "i", yaxs = "i")
lines(fit1, lwd = 2)
lines(lwr1, lty = "dashed")
lines(upr1, lty = "dashed")
lines(a2, col = "blue")
lines(fit2, lwd = 2)
lines(lwr2, lty = "dashed")
lines(upr2, lty = "dashed")
title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,n,365), labels = 0:30)
title(ylab = "Amplitude (hPa)", line = 2.0)

dev.off()


phase1 = atan(sims[,3,]/sims[,4,])
for (i in 1:ncol(phase1))
  phase1[,i] = (phase1[,i] - (i*omega)%%pi)%%pi

rho1 = atan(psi[1,]/psi[2,])
rho1 = (rho1 - ((1:n)*omega)%%pi)%%pi

fit1 = apply(phase1, 2, mean, na.rm = TRUE)
lwr1 = apply(phase1, 2, quantile, probs = 0.025, na.rm = TRUE)
upr1 = apply(phase1, 2, quantile, probs = 0.975, na.rm = TRUE)

# fit1 = sapply(1:length(fit1),
#               function(i) which.max(sin(omega*1:(365*24)/24 + fit1[i]))/24)
# lwr1 = sapply(1:length(lwr1),
#               function(i) which.max(sin(omega*1:(365*24)/24 + lwr1[i]))/24)
# upr1 = sapply(1:length(upr1),
#               function(i) which.max(sin(omega*1:(365*24)/24 + upr1[i]))/24)
#
# fit1 = fit1 - 365
# lwr1 = lwr1 - 365
# upr1 = upr1 - 365

phase2 = atan(sims[,5,]/sims[,6,])
for (i in 1:ncol(phase2))
  phase2[,i] = (phase2[,i] - (i*2*omega)%%pi)%%pi

rho2 = atan(psi[3,]/psi[4,])
rho2 = (rho2 - ((1:n)*2*omega)%%pi)%%pi

fit2 = apply(phase2, 2, mean, na.rm = TRUE)
lwr2 = apply(phase2, 2, quantile, probs = 0.025, na.rm = TRUE)
upr2 = apply(phase2, 2, quantile, probs = 0.975, na.rm = TRUE)

# fit2 = sapply(1:length(fit2),
#               function(i) which.max(sin(2*omega*1:(182.5*24)/24 + fit2[i]))/24)
# lwr2 = sapply(1:length(lwr2),
#               function(i) which.max(sin(2*omega*1:(182.5*24)/24 + lwr2[i]))/24)
# upr2 = sapply(1:length(upr2),
#               function(i) which.max(sin(2*omega*1:(182.5*24)/24 + upr2[i]))/24)

## Open eps device
postscript(file = "sim-tvar-phase.eps", width = width, height = height,
           pointsize = 8, paper = "special", horizontal = FALSE)

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3)

plot(rho1, type = "l", col = "red",  xaxt = "n",
     xlim = c(0,n), ylim = range(fit1,lwr1,upr1,fit2,lwr2,upr2),
     xaxs = "i", yaxs = "i")
lines(fit1, lwd = 2)
lines(lwr1, lty = "dashed")
lines(upr1, lty = "dashed")
lines(rho2, col = "blue")
lines(fit2, lwd = 2)
lines(lwr2, lty = "dashed")
lines(upr2, lty = "dashed")
title(xlab = "Year", line = 1.0)
axis(side = 1, at = seq(0,n,365), labels = 0:30)
title(ylab = "Phase", line = 2.0)

dev.off()

