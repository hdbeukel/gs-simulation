##############################################################
# Run one iteration of a simulation & output seasons to file #
##############################################################

# command line arguments:
#  1: simulation type (PS, GS, WGS)
#  2: number of seasons (even)
#  3: heritability
#  4: effect sampling method (normal, jannink)
#  5: iteration number
# output file written to: out/<1>/<2>-seasons/h2-<3>/<4>-effects/<5>.RDS

# load scripts
suppressMessages(source("scripts.R"))
# load founder dataset
load("data/ProcessedData.RData")

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
sim.function.name <- args[1]
sim.function = get(sim.function.name)
num.seasons <- args[2]
heritability <- as.numeric(args[3])
QTL.effects <- args[4]
it <- as.numeric(args[5])

# run simulation
seasons <- sim.function(founders, heritability,
                        num.seasons = num.seasons,
                        QTL.effects = QTL.effects)

# write output
out.dir <- sprintf("out/%s/%d-seasons/h2-%.1f/%s-effects",
                   sim.function.name, num.seasons,
                   heritability, QTL.effects)
file.name = paste(it, "RDS", sep=".")
dir.create(out.dir, showWarning = FALSE, recursive = TRUE)
full.path = paste(out.dir, file.name, sep="/")
saveRDS(seasons, file = full.path)
