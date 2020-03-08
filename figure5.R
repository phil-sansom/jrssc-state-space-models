################################
## Figure 5: Long-term trends ##
################################

## Loop over intervention types
for (type in c("mean","tvar","mean-alt","tvar-alt")) {

  ## Load trajectories
  load(paste0("nao-",type,".RData"))
  load(paste0("nao-",type,"-trajectories.RData"))

  trajectories = simplify2array(trajectories)
  gc()

  ## Plotting parameters
  width  = 33*12/72.27
  height = width*9/16


  ################
  ## Mean level ##
  ################

  fit = apply(trajectories[1,,], 1, median, na.rm = TRUE)
  lwr = apply(trajectories[1,,], 1, quantile, probs = 0.025, na.rm = TRUE)
  upr = apply(trajectories[1,,], 1, quantile, probs = 0.975, na.rm = TRUE)

  xat = as.Date(paste(seq(1980,2015,5),"01","01",sep="-"))
  xlab = format(xat, "%Y")

  mask = which(format(time,"%m-%d") == "01-01")
  masks = which(format(time,"%d") == "01")

  ## Open PDF device
  postscript(file = paste0("figure5a-",type,".eps"),
             width = width, height = width/2, horizontal = FALSE,
             pointsize = 16, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.5,3.0,0.5,0.5)+0.1, mgp = c(2.0,1/2,0),
      ps = 16, tcl = -1/3, xaxs = "i", yaxs = "i")

  plot (time, fit, type = "n", xaxt = "n", yaxt = "n", lwd = 2,
        xlim = as.Date(c("1979-01-01","2018-01-01")),
        # ylim = c(0,10), bty = "n")
        ylim = c(4.0,7.0), bty = "n")
  for (i in sample.int(1e3,1e1))
    lines(time[masks], trajectories[1,masks,i], col = "lightgrey")
  lines(time[mask], fit[mask], lwd = 2)
  lines(time[mask], lwr[mask], lty = 2)
  lines(time[mask], upr[mask], lty = 2)
  box()
  axis(side = 1, at = xat, labels = xlab)
  axis(side = 2, at = seq(3.0,7.0,0.25),
       labels = c("3.0","","3.5","","4.0","","4.5","","5.0","","5.5","",
                  "6.0","","6.5","","7.0"),
       mgp = c(2.5,1/2,0))
  title(ylab = "NAO Index (hPa)", line = 2.0)
  title(xlab = "Date", line = 1.5)

  #mtext("(a)", side = 3, line = -0.5, at = -10200)

  dev.off()


  ###########
  ## Trend ##
  ###########

  fit = apply(trajectories[2,,]*365.25, 1, median,
              na.rm = TRUE)
  lwr = apply(trajectories[2,,]*365.25, 1, quantile,
              probs = 0.025, na.rm = TRUE)
  upr = apply(trajectories[2,,]*365.25, 1, quantile,
              probs = 0.975, na.rm = TRUE)

  xat = as.Date(paste(seq(1980,2015,5),"01","01",sep="-"))
  xlab = format(xat, "%Y")

  mask = which(format(time,"%m-%d") == "01-01")
  masks = which(format(time,"%d") == "01")

  ## Open PDF device
  postscript(file = paste0("figure5b-",type,".eps"),
             width = width, height = width/2, horizontal = FALSE,
             pointsize = 16, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.5,4.0,0.5,0.5)+0.1, mgp = c(2.0,1/2,0),
      ps = 16, tcl = -1/3, xaxs = "i", yaxs = "i")

  plot (time, fit, type = "n", xaxt = "n", yaxt = "n", lwd = 2,
        xlim = as.Date(c("1979-01-01","2018-01-01")),
        # ylim = c(0,10), bty = "n")
        ylim = c(-0.1,+0.1), bty = "n")
  for (i in sample.int(1e3,1e1))
    lines(time[masks], trajectories[2,masks,i]*365.25, col = "lightgrey")
  lines(time[mask], fit[mask], lwd = 2)
  lines(time[mask], lwr[mask], lty = 2)
  lines(time[mask], upr[mask], lty = 2)
  abline(h = 0, lty = 3)
  box()
  axis(side = 1, at = xat, labels = xlab)
  axis(side = 2, at = seq(-0.1,+0.1,0.02),
       #     labels = c("3.0","","3.5","","4.0","","4.5","","5.0","","5.5","",
       #                "6.0","","6.5","","7.0"),
       mgp = c(2.5,1/2,0))
  title(ylab = "Trend (hPa/yr)", line = 3.0)
  title(xlab = "Date", line = 1.5)

  #mtext("(a)", side = 3, line = -0.5, at = -10200)

  dev.off()


  ##################
  ## Annual cycle ##
  ##################

  xat = as.Date(paste(seq(1980,2015,5),"01","01",sep="-"))
  xlab = format(xat, "%Y")

  samples1 = sqrt(trajectories[3,,]^2 + trajectories[4,,]^2)
  fit1 = apply(samples1, 1, median, na.rm = TRUE)
  lwr1 = apply(samples1, 1, quantile, probs = 0.025, na.rm = TRUE)
  upr1 = apply(samples1, 1, quantile, probs = 0.975, na.rm = TRUE)

  samples2 = sqrt(trajectories[5,,]^2 + trajectories[6,,]^2)
  fit2 = apply(samples2, 1, median, na.rm = TRUE)
  lwr2 = apply(samples2, 1, quantile, probs = 0.025, na.rm = TRUE)
  upr2 = apply(samples2, 1, quantile, probs = 0.975, na.rm = TRUE)

  mask = which(format(time,"%m-%d") == "01-01")
  masks = which(format(time,"%d") == "01")

  ## Open PDF device
  postscript(file = paste0("figure5c-",type,".eps"),
             width = width, height = width/2, horizontal = FALSE,
             pointsize = 16, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.5,3.0,0.5,0.5)+0.1, mgp = c(2.0,1/2,0),
      ps = 16, tcl = -1/3, xaxs = "i", yaxs = "i")

  plot (time[mask], fit1[mask], type = "n", xaxt = "n", yaxt = "n", lwd = 2,
        xlim = as.Date(c("1979-01-01","2018-01-01")), ylim = c(0,6), bty = "n")
  for (i in sample.int(1e3,1e1)) {
    lines(time[masks], samples1[masks,i], col = "lightgrey")
    lines(time[masks], samples2[masks,i], col = "lightgrey")
  }
  lines(time[mask], fit1[mask], lwd = 2)
  lines(time[mask], lwr1[mask], lty = 2)
  lines(time[mask], upr1[mask], lty = 2)
  lines(time[mask], fit2[mask], lwd = 2)
  lines(time[mask], lwr2[mask], lty = 2)
  lines(time[mask], upr2[mask], lty = 2)


  text(1e4, fit1[1], expression(psi[1]))
  text(1e4, fit2[1], expression(psi[2]))

  box()
  axis(side = 1, at = xat, labels = xlab)
  axis(side = 2, at = seq(0,6,0.5),
       mgp = c(2.5,1/2,0))
  title(ylab = "Amplitude (hPa)", line = 2.0)
  title(xlab = "Date", line = 1.5)

  #mtext("(c)", side = 3, line = -0.5, at = -10200)

  dev.off()


  ###########
  ## Phase ##
  ###########

  xat = as.Date(paste(seq(1980,2015,5),"01","01",sep="-"))
  xlab = format(xat, "%Y")

  samples1 = atan(trajectories[3,,]/trajectories[4,,])
  for (i in 1:nrow(samples1))
    samples1[i,] = (samples1[i,] - (i*omega)%%pi)%%pi

  fit1 = apply(samples1, 1, median, na.rm = TRUE)
  lwr1 = apply(samples1, 1, quantile, probs = 0.025, na.rm = TRUE)
  upr1 = apply(samples1, 1, quantile, probs = 0.975, na.rm = TRUE)

  fit1 = sapply(1:length(fit1),
                function(i) which.max(sin(omega*1:(365*24)/24 + fit1[i]))/24)
  lwr1 = sapply(1:length(lwr1),
                function(i) which.max(sin(omega*1:(365*24)/24 + lwr1[i]))/24)
  upr1 = sapply(1:length(upr1),
                function(i) which.max(sin(omega*1:(365*24)/24 + upr1[i]))/24)

  fit1 = fit1 - 365
  lwr1 = lwr1 - 365
  upr1 = upr1 - 365

  samples2 = atan(trajectories[5,,]/trajectories[6,,])
  for (i in 1:nrow(samples2))
    samples2[i,] = (samples2[i,] - (i*2*omega)%%pi)%%pi

  fit2 = apply(samples2, 1, mean, na.rm = TRUE)
  lwr2 = apply(samples2, 1, quantile, probs = 0.025, na.rm = TRUE)
  upr2 = apply(samples2, 1, quantile, probs = 0.975, na.rm = TRUE)

  fit2 = sapply(1:length(fit2),
                function(i) which.max(sin(2*omega*1:(182.5*24)/24 +
                                            fit2[i]))/24)
  lwr2 = sapply(1:length(lwr2),
                function(i) which.max(sin(2*omega*1:(182.5*24)/24 +
                                            lwr2[i]))/24)
  upr2 = sapply(1:length(upr2),
                function(i) which.max(sin(2*omega*1:(182.5*24)/24 +
                                            upr2[i]))/24)

  mask = which(format(time,"%m-%d") == "01-01")
  masks = which(format(time,"%d") == "01")

  p01 = which.max(sin(  omega*1:(365  *24)/24))/24
  p02 = which.max(sin(2*omega*1:(182.5*24)/24))/24

  ## Open PDF device
  postscript(file = paste0("figure5d-",type,".eps"),
             width = width, height = width/2, horizontal = FALSE,
             pointsize = 16, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.5,3.5,0.5,0.5), mgp = c(2.0,1/2,0),
      ps = 16, tcl = -1/3, xaxs = "i", yaxs = "i")

  plot (time[mask], fit1[mask], type = "n", xaxt = "n", yaxt = "n", lwd = 2,
        xlim = as.Date(c("1979-01-01","2018-01-01")), ylim = c(-31,+31),
        bty = "n")
  # lines(time, lwr, lty = "solid", lwd = 2, col = "darkgrey")
  # lines(time, upr, lty = "solid", lwd = 2, col = "darkgrey")
  for (i in sample.int(1e3,1e1))
    lines(time[masks], p01 - 182.5*samples1[masks,i]/pi, col = "lightgrey")
  for (i in sample.int(1e3,1e1))
    lines(time[masks], p02 - 365/4*samples2[masks,i]/pi, col = "lightgrey")
  lines(time[mask], fit1[mask], lwd = 2)
  lines(time[mask], lwr1[mask], lty = 2)
  lines(time[mask], upr1[mask], lty = 2)
  lines(time[mask], fit2[mask], lwd = 2)
  lines(time[mask], lwr2[mask], lty = 2)
  lines(time[mask], upr2[mask], lty = 2)


  text(1e4, fit1[1], expression(psi[1]))
  text(1e4, fit2[1], expression(psi[2]))

  box()
  axis(side = 1, at = xat, labels = xlab)
  axis(side = 2, at = c(-31,0,+31), labels = c("Dec","Jan","Feb"),
       mgp = c(2.5,1/2,0))
  title(ylab = "Phase", line = 2.5)
  title(xlab = "Date", line = 1.5)

  #mtext("(d)", side = 3, line = -0.5, at = -10200)

  dev.off()


  ###################
  ## AR Parameters ##
  ###################

  fit = apply(trajectories[12:16,,], c(1,2), median,
              na.rm = TRUE)
  lwr = apply(trajectories[12:16,,], c(1,2), quantile,
              probs = 0.025, na.rm = TRUE)
  upr = apply(trajectories[12:16,,], c(1,2), quantile,
              probs = 0.975, na.rm = TRUE)

  xat = as.Date(paste(seq(1980,2015,5),"01","01",sep="-"))
  xlab = format(xat, "%Y")

  mask = which(format(time,"%m-%d") == "01-01")
  masks = which(format(time,"%d") == "01")

  ## Open PDF device
  postscript(file = paste0("figure5e-",type,".eps"),
             width = width, height = width/2, horizontal = FALSE,
             pointsize = 16, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2.5,3.0,0.5,0.5)+0.1, mgp = c(2.0,1/2,0),
      ps = 16, tcl = -1/3, xaxs = "i", yaxs = "i")

  plot (time[mask], fit[1,mask], type = "n", xaxt = "n", yaxt = "n", lwd = 2,
        xlim = as.Date(c("1979-01-01","2018-01-01")), ylim = c(-2,+2),
        bty = "n")
  # lines(time, lwr, lty = "solid", lwd = 2, col = "darkgrey")
  # lines(time, upr, lty = "solid", lwd = 2, col = "darkgrey")
  for (i in sample.int(1e3,1e1))
    for (j in 1:p)
      lines(time[masks], trajectories[11+j,masks,i], col = "lightgrey")
  for (j in 1:p) {
    lines(time[mask], fit[j,mask], lwd = 2)
    lines(time[mask], lwr[j,mask], lty = 2)
    lines(time[mask], upr[j,mask], lty = 2)
  }
  box()
  axis(side = 1, at = xat, labels = xlab)
  axis(side = 2, at = seq(-2,+2,+0.5))
  title(ylab = "Autoregressive coefficient", line = 2.0)
  title(xlab = "Date", line = 1.5)
  for (j in 1:p)
    text(1e4, fit[j,length(y)], substitute(phi[j], list(j = j)))

  #mtext("(b)", side = 3, line = -0.5, at = -10200)

  dev.off()

  ## Garbage collection
  rm(list = ls()[ls() != "type"])
  gc()

} ## type
