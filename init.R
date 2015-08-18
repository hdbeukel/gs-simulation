# load required libraries
library(BGLR)     # version used: 1.0.4
library(coda)     # version used: 0.17-1
library(rrBLUP)   # version used: 4.3
library(synbreed) # version used: 0.10-5
library(hypred)   # version used: 0.5
library(gdata)    # version used: 2.17.0
library(Hmisc)    # version used: 3.16-0
library(rJava)    # version used: 0.10-5
library(setRNG)   # version used: 2013.9-1

# set path to BEAGLE used by synbreed
beaglePath <- "beaglejar"

# load scripts
source("data/dataprocessing.R")
source("sim/hypred.R")
source("sim/prediction.R")
source("sim/simulation.R")
source("sim/selection.R")
source("sim/java.R")