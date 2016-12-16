# load required libraries
library(BGLR)        # version used: 1.0.4
library(coda)        # version used: 0.17-1
library(rrBLUP)      # version used: 4.3
library(hypred)      # version used: 0.5
library(gdata)       # version used: 2.17.0
library(Hmisc)       # version used: 3.17-2
library(rJava)       # version used: 0.9-8
library(setRNG)      # version used: 2013.9-1
library(animation)   # version used: 2.4
library(synbreed)    # version used: 0.10-5

# set path to BEAGLE used by synbreed
beaglePath <- "beaglejar"

# load scripts
source("R/dataprocessing.R")
source("R/hypred.R")
source("R/prediction.R")
source("R/simulation.R")
source("R/plots.R")
source("R/selection.R")
source("R/java.R")