###########
## ANOVA ##
###########

## Loop over intervention types
for (type in c("mean","mean-alt")) {

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

  ## Extract driver weights
  parameters = parameters[10:12,]
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

  ## Compute seasonal means for each posterior sample
  yy      = tapply(y, list(year, season), mean)
  eta     = lapply(1:ncol(eta    ),
                   function(j) tapply(eta    [,j], list(year, season), mean,
                                      na.rm = TRUE))
  x       = lapply(1:ncol(x      ),
                   function(j) tapply(x      [,j], list(year, season), mean,
                                      na.rm = TRUE))
  delta   = lapply(1:ncol(delta  ),
                   function(j) tapply(delta  [,j], list(year, season), mean,
                                      na.rm = TRUE))
  epsilon = lapply(1:ncol(epsilon),
                   function(j) tapply(epsilon[,j], list(year, season), mean,
                                      na.rm = TRUE))
  yy      = array(yy, c(dim(yy),length(eta)))
  eta     = simplify2array(eta    )
  x       = simplify2array(x      )
  delta   = simplify2array(delta  )
  epsilon = simplify2array(epsilon)
  dimnames(yy) = dimnames(eta)

  yy      = yy     [rownames(yy)      != "1978",,]
  eta     = eta    [rownames(eta)     != "1978",,]
  x       = x      [rownames(x)       != "1978",,]
  delta   = delta  [rownames(delta)   != "1978",,]
  epsilon = epsilon[rownames(epsilon) != "1978",,]

  yy     ["2017","DJF",] = NA
  eta    ["2017","DJF",] = NA
  x      ["2017","DJF",] = NA
  delta  ["2017","DJF",] = NA
  epsilon["2017","DJF",] = NA

  years   = as.numeric(rownames(eta))
  seasons = colnames(eta)

  ## ANOVA on posterior samples
  ss = array(NA, c(length(seasons),1000,4))
  for (s in 1:length(seasons))
    for (j in 1:1000) {
      buffer =
        anova(lm(yy[,s,j] ~ eta[,s,j] + delta[,s,j] + x[,s,j]))$`Sum Sq`
      if (length(buffer) == 4) {
        ss[s,j,] = buffer / sum(buffer, na.rm = TRUE)
      } else {
        ss[s,j,] = c(buffer[1],0,buffer[2:3]) / sum(buffer, na.rm = TRUE)
      }
    }

  atab = apply(ss, c(1,3), quantile, probs = c(0.05,0.50,0.95), na.rm = TRUE)
  dimnames(atab)[[3]] = c("eta","delta","x","epsilon")
  dimnames(atab)[[2]] = seasons

  mam = round(t(atab[c(2,1,3),"MAM",]),2)
  jja = round(t(atab[c(2,1,3),"JJA",]),2)
  son = round(t(atab[c(2,1,3),"SON",]),2)
  djf = round(t(atab[c(2,1,3),"DJF",]),2)

  paste(type)
  paste(t(mam), collapse = " ")
  paste(t(jja), collapse = " ")
  paste(t(son), collapse = " ")
  paste(t(djf), collapse = " ")

} ## type
