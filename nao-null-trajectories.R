## Load libraries
library(coda)
library(compiler)
library(MASS)
library(parallel)
mc.cores = 4
options(mc.cores = mc.cores)

 input.file = "nao-null.RData"
output.file = "nao-null-trajectories.RData"

source("trajectories.R")
