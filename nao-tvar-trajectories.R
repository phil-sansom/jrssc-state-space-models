## Load libraries
library(coda)
library(compiler)
library(MASS)
library(parallel)
mc.cores = 4
options(mc.cores = mc.cores)

 input.file = "nao-tvar.RData"
output.file = "nao-tvar-trajectories.RData"

source("trajectories.R")
