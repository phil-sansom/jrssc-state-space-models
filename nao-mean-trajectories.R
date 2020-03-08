## Load libraries
library(coda)
library(compiler)
library(MASS)
library(parallel)
mc.cores = 4
options(mc.cores = mc.cores)

 input.file = "nao-mean.RData"
output.file = "nao-mean-trajectories.RData"

source("trajectories.R")
