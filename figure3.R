##############################
## Figure 3: Coupled period ##
##############################

## Load library
library(coda)

## Loop over intervention types
for (type in c("mean","tvar","mean-alt","tvar-alt")) {

  ## Load posterior samples
  load(paste0("nao-",type,".RData"))

  ## Determine when adaptation ended
  psrfs  = sapply(parameters,
                  function(x) max(gelman.diag(x, autoburnin = FALSE)$psrf[,1]))
  alpha  = min(which(psrfs < 1.10)) + 1
  gamma  = length(parameters)
  blocks = gamma - alpha + 1

  ## Extract parameters
  buffers = list()
  for (i in 1:n.chains) {
    buffers[[i]] = parameters[[alpha]][[i]]
    for (j in (alpha + 1):gamma) {
      buffers[[i]] = rbind(buffers[[i]],parameters[[j]][[i]])
    } ## j
  } ## i
  parameters = buffers[[1]]
  for (i in 2:n.chains)
    parameters = rbind(parameters,buffers[[i]])
  parameters = parameters[,10:12]

  make.f = function(x) {

    cdate  = floor(x[1]) + 1
    clen   = floor(x[2])
    cprop  = 1/(1+exp(-x[3]))

    f  = numeric(365)
    c0 = floor(0.5*clen*cprop)
    if (c0 > 0) {
      fs = c(seq( 1/(c0+1),c0/(c0+1), 1/(c0+1)),rep(1,clen-2*c0),
             seq(c0/(c0+1), 1/(c0+1),-1/(c0+1)))
    } else {
      fs = rep(1,clen)
    }
    f[cdate-1+1:clen] = fs
    f = f[1:365]
    if(cdate - 1 + clen > 365)
      f[1:(cdate - 1 + clen - 365)] = fs[(365 - cdate + 2):clen]

    return(f)

  }

  fs = array(NA, c(nrow(parameters),365))
  for (i in 1:nrow(parameters))
    fs[i,] = make.f(parameters[i,])
  pm = apply(fs, 2, function(x) mean(x, na.rm = TRUE))

  ## Plotting parameters
  width  = 33*12/72.27/2
  height = width*9/16

  xat = which(format(time[1:365], "%d") == "01")

  ## Open eps device
  postscript(file = paste0("figure3-",type,".eps"),
             width = width, height = height, horizontal = FALSE,
             pointsize = 8, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2,2.5,0.5,0.5), mgp = c(2.5,1/2,0),
      ps = 8, tcl = -1/3)

  plot (pm, type = "n", xlim = c(1,365), xaxt ="n", ylim = c(0,1),
        xaxs = "i", frame.plot = FALSE)
  for (i in sample.int(nrow(fs),1e1))
    lines(fs[i,], col = "darkgrey")
  lines(pm, lwd = 2)
  axis (side = 1, at = xat[seq(1,11,2)], labels = month.abb[seq(1,11,2)],
        mgp = c(2.5,1/3,0))
  axis (side = 1, at = xat[seq(2,12,2)], labels = month.abb[seq(2,12,2)],
        mgp = c(2.5,1/3,0))
  title(xlab = "Date", line = 1)
  title(ylab = expression(paste("Intervention ",lambda[t])), line = 1.5)
  box()

  dev.off()

  ## Garbage collection
  rm(list = ls()[ls() != "type"])
  gc()

} ## type
