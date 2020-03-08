###########################
## Figure 6: Attribution ##
###########################

## Loop over intervention types
for (type in c("mean","tvar","mean-alt","tvar-alt")) {

  ## Load trajectories
  load(paste0("nao-",type,".RData"))
  load(paste0("nao-",type,"-trajectories.RData"))

  ## Extract parameters and trajectories
  parameters = t(parameters)
  trajectories = simplify2array(trajectories)
  gc()

  ## Define seasons
  month = format(time, "%b")
  seas = c("DJF","DJF","MAM","MAM","MAM","JJA",
           "JJA","JJA","SON","SON","SON","DJF")
  names(seas) = month.abb
  season = seas[month]
  names(season) = NULL

  year = as.numeric(format(time, "%Y"))
  year[month %in% c("Jan","Feb")] = year[month %in% c("Jan","Feb")] - 1

  parameters = parameters[10:12,]

  ## Extract driver weights
  make.f = function(x) {

    cdate  = floor(x[1]) + 1
    clen   = floor(x[2])
    cprop  = 1/(1+exp(-x[3]))

    f  = numeric(length(y))
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
    f = f[1:length(y)]

    return(f)

  }

  fs = array(NA, c(length(y),1000))
  for (i in 1:1000)
    fs[,i] = make.f(parameters[,i])

  ## Extract forecast components
  mu      = trajectories[ 1,,]
  s1      = trajectories[ 3,,]
  s2      = trajectories[ 5,,]
  x       = trajectories[ 7,,]
  delta   = trajectories[17,,]
  eta     = mu + s1 + s2
  delta   = fs * delta
  yy      = array(y, c(length(y),1e3))
  epsilon = yy - eta - x - delta
  rm(trajectories)
  gc()

  ## Compute seasonal means
  eta     = tapply(apply(eta    , 1, mean, na.rm = TRUE),
                   list(year,season), mean, na.rm = TRUE)
  x       = tapply(apply(x      , 1, mean, na.rm = TRUE),
                   list(year,season), mean, na.rm = TRUE)
  delta   = tapply(apply(delta  , 1, mean, na.rm = TRUE),
                   list(year,season), mean, na.rm = TRUE)
  epsilon = tapply(apply(epsilon, 1, mean, na.rm = TRUE),
                   list(year,season), mean, na.rm = TRUE)

  seasonal.samples = list(eta = eta, delta = delta, x = x, epsilon = epsilon)
  seasonal.samples = simplify2array(seasonal.samples)
  seasonal.samples = seasonal.samples[-1,,]
  seasonal.samples["2017","DJF",] = NA
  for (i in 1:4)
    seasonal.samples[,i,"eta"] = seasonal.samples[,i,"eta"] -
      mean(seasonal.samples[,i,"eta"], na.rm = TRUE)

  rm(mu,s1,s2,x,delta,eta,epsilon)
  gc()

  ## Bar plotting function
  myrect = function(x, theta, col) {

    col   = c(rev(col[theta < 0]),col[theta > 0])
    theta = c(rev(cumsum(theta[theta < 0])),0,cumsum(theta[theta > 0]))

    rect(x - 1/3, theta[1], x + 1/3, theta[2], col = col[1])
    rect(x - 1/3, theta[2], x + 1/3, theta[3], col = col[2])
    rect(x - 1/3, theta[3], x + 1/3, theta[4], col = col[3])
    rect(x - 1/3, theta[4], x + 1/3, theta[5], col = col[4])

  }

  col   = c("white","lightgrey","darkgrey","black")

  nn = nrow(seasonal.samples)

  ## Plotting parameters
  width  = 33*12/72.27
  height = width*9/16

  ## Open PDF device
  postscript(file = paste0("figure6-",type,".eps"),
             width = width, height = height, horizontal = FALSE,
             pointsize = 8, paper = "special")

  ## Plotting parameters
  par(ann = FALSE, las = 1, mar = c(2,2.5,0.5,0.5), mgp = c(2.5,1/3,0),
      ps = 8, tcl = -1/3, xaxs = "i", yaxs = "i")

  xlab = seq(1980, 2015, 5)

  plot(1:nn, 1:nn, type = "n", xlim = c(1/3,nn+2/3), xaxt = "n",
       ylim = c(-12,+10), yaxt = "n")
  for (i in 1:nn)
    myrect(i, seasonal.samples[i,"DJF",], col)
  axis (side = 1, at = seq(2,nn,5), labels = xlab)
  axis (side = 2, at = seq(-12,+10,2), hadj = 1.2)
  title(xlab = "Date", line = 1)
  title(ylab = paste("Contribution to","DJF","mean (hPa)"), line = 1.5)

  legend("bottom", legend = c("Mean","Coupling","Irregular","Obs error"),
         fill = col, bty = "n", horiz = TRUE)

  dev.off()

  ## Garbage collection
  rm(list = ls()[ls() != "type"])
  gc()

} ## type
