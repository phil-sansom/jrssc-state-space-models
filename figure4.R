###########################################
## Figure 4: Posterior predictive checks ##
###########################################

## Loop over intervention types
for (type in c("mean","tvar","mean-alt","tvar-alt")) {

  ## Load samples
  load(paste0("nao-",type,".RData"))
  load(paste0("nao-",type,"-predictive.RData"))

  ## Plotting parameters
  width  = 33*12/72.27/2
  height = width*9/16


  ##################
  ## Monthly Mean ##
  ##################

  tt = 1:length(nao)
  data.mean = tapply(lm(nao ~ tt)$residuals, format(time, "%m"), mean)
  sim.means = apply(predictive, 2,
                    function(x) tapply(lm(x ~ tt)$residuals, format(time, "%m"), mean))

  fit = apply(sim.means, 1, median, na.rm = TRUE)
  lwr = apply(sim.means, 1, quantile, probs = 0.05, na.rm = TRUE)
  upr = apply(sim.means, 1, quantile, probs = 0.95, na.rm = TRUE)

  polyx = c(1:12,12:1)
  polyy = c(upr,rev(lwr))

  ## Open eps device
  postscript(file = paste0("figure4a-",type,".eps"),
             width = width, height = height, horizontal = FALSE,
             pointsize = 8, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.0,2.5,0.5,0.5), mgp = c(2.5,1/3,0),
      ps = 8, tcl = -1/3)

  plot(1:12, data.mean, type = "n", ylim = c(-8,+8), xlim = c(1,12),
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
  polygon(polyx, polyy, border = NA, col = "lightgrey")
  lines(1:12, fit, lwd = 2, col = "darkgrey")
  lines(1:12, data.mean, type = "b", pch = 16)
  box()

  axis(side = 1, at = seq(1,11,2), labels = month.abb[seq(1,11,2)])
  axis(side = 1, at = seq(2,12,2), labels = month.abb[seq(2,12,2)])
  axis(side = 2, at = seq(-8,+8,2), hadj = 1.4)
  title(xlab = "Month", line = 1)
  title(ylab = "Mean NAO Index (hPa)", line = 1.5)

  dev.off()


  #########################
  ## Day-on-day Variance ##
  #########################

  data.var = tapply(diff(lm(nao ~ tt)$residuals), format(time[-1], "%m"), sd)
  sim.vars = apply(predictive, 2,
                   function(x) tapply(diff(lm(x ~ tt)$residuals), format(time[-1], "%m"), sd))

  fit = apply(sim.vars, 1, median, na.rm = TRUE)
  lwr = apply(sim.vars, 1, quantile, probs = 0.05, na.rm = TRUE)
  upr = apply(sim.vars, 1, quantile, probs = 0.95, na.rm = TRUE)

  polyx = c(1:12,12:1)
  polyy = c(upr,rev(lwr))

  ## Open eps device
  postscript(file = paste0("figure4b-",type,".eps"),
             width = width, height = height, horizontal = FALSE,
             pointsize = 8, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
      ps = 8, tcl = -1/3)

  plot(1:12, data.var, type = "n", ylim = c(0,4), xlim = c(1,12),
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
  polygon(polyx, polyy, border = NA, col = "lightgrey")
  lines(1:12, fit, lwd = 2, col = "darkgrey")
  lines(1:12, data.var, type = "b", pch = 16)
  box()

  axis(side = 1, at = seq(1,11,2), labels = month.abb[seq(1,11,2)])
  axis(side = 1, at = seq(2,12,2), labels = month.abb[seq(2,12,2)])
  axis(side = 2, at = seq(0,4,0.5), hadj = 1.2)
  title(xlab = "Month", line = 1)
  title(ylab = "Day-on-day SD (hPa)", line = 2.0)

  dev.off()


  ###########################
  ## Inter-annual Variance ##
  ###########################

  year  = format(time, "%Y")
  month = format(time, "%m")

  inter.var = function(x) {

    x = tapply(lm(x ~ tt)$residuals, list(year = year, month = month), mean)
    return(apply(x, 2, sd))

  }

  data.var = inter.var(nao)
  sim.vars = apply(predictive, 2, inter.var)

  fit = apply(sim.vars, 1, median, na.rm = TRUE)
  lwr = apply(sim.vars, 1, quantile, probs = 0.05, na.rm = TRUE)
  upr = apply(sim.vars, 1, quantile, probs = 0.95, na.rm = TRUE)

  polyx = c(1:12,12:1)
  polyy = c(upr,rev(lwr))

  ## Open eps device
  postscript(file = paste0("figure4c-",type,".eps"),
             width = width, height = height, horizontal = FALSE,
             pointsize = 8, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.0,2.5,0.5,0.5), mgp = c(2.5,1/3,0),
      ps = 8, tcl = -1/3)

  plot(1:12, data.var, type = "n", ylim = c(0,8), xlim = c(1,12),
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
  polygon(polyx, polyy, border = NA, col = "lightgrey")
  lines(1:12, fit, lwd = 2, col = "darkgrey")
  lines(1:12, data.var, type = "b", pch = 16)
  box()

  axis(side = 1, at = seq(1,11,2), labels = month.abb[seq(1,11,2)])
  axis(side = 1, at = seq(2,12,2), labels = month.abb[seq(2,12,2)])
  axis(side = 2, at = seq(0,8,1), hadj = 1.2)
  title(xlab = "Month", line = 1)
  title(ylab = "Inter-annual SD (hPa)", line = 1.5)

  dev.off()


  ###############################
  ## Auto-correlation function ##
  ###############################

  source("acf.R")

  tt = length(nao)
  omega = 2*pi/365.25
  ts = 1:tt
  c1 = cos(  omega*ts)
  s1 = sin(  omega*ts)
  c2 = cos(2*omega*ts)
  s2 = sin(2*omega*ts)

  season = format(time, "%m")
  season[! season %in% c("12","01","02","03")] = "jja"
  season[  season %in% c("12","01","02","03")] = "djf"

  acf.seas = function(x) {

    if (any(is.na(x))) {
      zz = matrix(NA, 31, 2)
    } else {
      xx = lm(x ~ ts + c1 + s1 + c2 + s2)$residuals
      zz = myacf(xx, 30, season)$acf
    }
    return(zz)

  }

  data.acf = acf.seas(nao)
  sim.acfs = lapply(1:1e3, function(j) acf.seas(predictive[,j]))
  sim.acfs = simplify2array(sim.acfs)

  fit = apply(sim.acfs, c(1,2), median, na.rm = TRUE)
  lwr = apply(sim.acfs, c(1,2), quantile, probs = 0.05, na.rm = TRUE)
  upr = apply(sim.acfs, c(1,2), quantile, probs = 0.95, na.rm = TRUE)

  polyx.djf = c(0:30,30:0)
  # polyy.djf = c(upr[,1],rev(lwr[,1]))
  polyy.djf = c(upr[,1],rev(upr[,2]))

  polyx.mid = c(0:30,30:0)
  polyy.mid = c(lwr[,1],rev(upr[,2]))

  polyx.jja = c(0:30,30:0)
  polyy.jja = c(lwr[,1],rev(lwr[,2]))


  ## Open eps device
  postscript(file = paste0("figure4d-",type,".eps"),
             width = width, height = height, horizontal = FALSE,
             pointsize = 8, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.0,3.0,0.5,0.5), mgp = c(2.5,1/3,0),
      ps = 8, tcl = -1/3)

  lightgrey = rgb(211,211,211,211,"lightgrey",255)
  grey      = rgb(190,190,190,190,"grey",255)

  plot(0:30, data.acf[,1], type = "n", ylim = c(-0.1,1), xlim = c(0,30),
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
  polygon(polyx.djf, polyy.djf, border = NA, col = "lightgrey")
  polygon(polyx.mid, polyy.mid, border = NA, col = "grey")
  polygon(polyx.jja, polyy.jja, border = NA, col = "lightgrey")
  lines(0:30, fit[,1], lwd = 2, col = "darkgrey")
  lines(0:30, fit[,2], lwd = 2, col = "darkgrey")
  lines(0:30, data.acf[,1], type = "b", pch = 16)
  lines(0:30, data.acf[,2], type = "b", pch = 17)
  abline(h = 0, lty = 3)
  box()

  axis(side = 1, at = seq(0,30,5))
  axis(side = 2, at = seq(-0.1,1,0.1), hadj = 1.2)
  title(xlab = "Lag", line = 1)
  title(ylab = "Autocorrelation", line = 2.0)

  legend("topright", legend = c("Dec-Mar","Apr-Nov"), col = c(1,1),
         lty = c(1,1), pch = c(16,17), bty = "n")

  dev.off()

  ## Garbage collection
  rm(list = ls()[ls() != "type"])
  gc()

} ## type
