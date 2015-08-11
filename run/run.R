##############################################################
# Run one iteration of a simulation & output seasons to file #
##############################################################

# command line arguments:
#  1: simulation type ("PS", "GS", "WGS", "CGS")
#  2: number of seasons
#  3: heritability
#  4: additional TP size, ignored for PS
#  5: effect sampling method ("normal", "jannink")
#  6: GP method ("RR", "BRR"), ignored for PS
#  7: diversity measure ("MR", "HE"; ignored for all types except CGS) 
#  8: diversity weight (ignored for all types except CGS) 
#  9: CGS type ("index", "split"; ignored for all types except CGS) 
#  10: random generator seed (empty string for random seed)
#  11: iteration number
# output file written to:
#  - GS, WGS: out/<1>/<2>-seasons/h2-<3>/addTP-<4>/<5>-effects/<6>/<11>.RDS
#  - CGS:     out/<1>/<2>-seasons/h2-<3>/addTP-<4>/<5>-effects/<6>/<7>-<8>/<9>/<11>.RDS
#  - PS:      out/<1>/<2>-seasons/h2-<3>/<5>-effects/<11>.RDS

# load scripts
suppressMessages(source("scripts.R"))
# set Java options
options(java.parameters = "-Xmx1g")

# load founder dataset
load("data/ProcessedData.RData")

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
sim.function.name <- args[1]
sim.function = get(sim.function.name)
num.seasons <- as.numeric(args[2])
heritability <- as.numeric(args[3])
add.TP <- as.numeric(args[4])
QTL.effects <- args[5]
gp.method <- args[6]
if(sim.function.name == "PS"){
  gp.method <- ""
}
div.measure <- args[7]
div.weight <- as.numeric(args[8])
CGS.type <- args[9]
seed <- as.numeric(args[10])
it <- as.numeric(args[11])

# set seed
if(is.na(seed)){
  seed <- ceiling(runif(1, 0, 2^31-1))
}
set.seed(seed)
message(sprintf("Seed: %d", seed))

# run simulation
seasons <- sim.function(founders, heritability,
                        num.seasons = num.seasons,
                        QTL.effects = QTL.effects,
                        gp.method = gp.method,
                        add.TP = add.TP,
                        div.measure = div.measure,
                        div.weight = div.weight,
                        type = CGS.type)

# print warnings
warnings()

# write output
if(sim.function.name == "PS"){
  out.dir <- sprintf("out/%s/%d-seasons/h2-%.1f/%s-effects",
                     sim.function.name, num.seasons,
                     heritability, QTL.effects)
} else if(sim.function.name == "CGS") {
  out.dir <- sprintf("out/%s/%d-seasons/h2-%.1f/addTP-%d/%s-effects/%s/%s-%.2f/%s",
                     sim.function.name, num.seasons, heritability,
                     add.TP, QTL.effects, gp.method,
                     div.measure, div.weight, CGS.type)
} else {
  out.dir <- sprintf("out/%s/%d-seasons/h2-%.1f/addTP-%d/%s-effects/%s",
                     sim.function.name, num.seasons, heritability,
                     add.TP, QTL.effects, gp.method)
}
file.name <- paste(it, "RDS", sep=".")
dir.create(out.dir, showWarning = FALSE, recursive = TRUE)
full.path <- paste(out.dir, file.name, sep="/")
saveRDS(seasons, file = full.path)








