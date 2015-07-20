##############################################################
# Run one iteration of a simulation & output seasons to file #
##############################################################

# command line arguments:
#  1: simulation type (PS, GS, WGS)
#  2: heritability
#  3: iteration number
# output file written to: out/<1>/<2>/<3>.Rdata

# load scripts
suppressMessages(source("scripts.R"))
# load founder dataset
load("data/ProcessedData.RData")

# get command line arguments
args = commandArgs(trailingOnly = TRUE)
sim.function.name = args[1]
sim.function = get(sim.function.name)
heritability = as.numeric(args[2])
it = as.numeric(args[3])

# run simulation
seasons = sim.function(founders, heritability)

# write output
out.dir = paste("out", sim.function.name, heritability, sep="/")
file.name = paste(it, "RDS", sep=".")
dir.create(out.dir, showWarning = FALSE, recursive = TRUE)
full.path = paste(out.dir, file.name, sep="/")
saveRDS(seasons, file = full.path)
