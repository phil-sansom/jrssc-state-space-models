## Load libraries
library(coda)
library(compiler)
library(MASS)
library(parallel)
mc.cores = 4
options(mc.cores = mc.cores)

 input.file = "nao-mean-alt.RData"
output.file = "nao-mean-alt-trajectories.RData"

source("trajectories.R")
