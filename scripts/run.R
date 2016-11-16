
##############################################################
# Run one iteration of a simulation & output seasons to file #
##############################################################

# should be run from project root directory, with 'Rscript scripts/run.R <args>'

# command line arguments:
#  1:  simulation type ("PS", "GS", "WGS", "WGS2", "CGS", "OC1", "OC2")
#  2:  number of seasons
#  3:  heritability
#  4:  additional TP size, ignored for PS
#  5:  effect sampling method ("normal", "jannink")
#  6:  GP method ("RR", "BRR"), ignored for PS
#  7:  diversity measure ("LOGall", "OC", "HEall", "HEfav", "LOGfav"; CGS only) 
#  8:  diversity weight (CGS only) 
#  9:  CGS type ("index", "split"; CGS only)
#  10: desired rate of inbreeding delta F (OC only)
#  11: iteration number
#  12: store full marker data of all intermediate populations (logical)
#      when missing or empty string, defaults to FALSE
#  13: random generator seed (numerical)
#      when missing or empty string, a random seed is set
#      if set, the seed is included in the output file name
#  14: fixed base population (numerical)
#      when missing or empty string, a base population is generated when starting the simulation
#      if set, a fixed bp is loaded from 'data/fixed-bp' based on heritability and the given bp number,
#      and the bp number is included in the output file name
#
# output file written to:
#  - GS, WGS:  out/<1>/<2>-seasons/h2-<3>/addTP-<4>/<5>-effects/<6>/*<11>.RDS
#  - CGS:      out/<1>/<2>-seasons/h2-<3>/addTP-<4>/<5>-effects/<6>/<7>-<8>/<9>/*<11>.RDS
#  - OC1, OC2: out/<1>/<2>-seasons/h2-<3>/addTP-<4>/<5>-effects/<6>/dF-<10>/*<11>.RDS
#  - PS:       out/<1>/<2>-seasons/h2-<3>/<5>-effects/*<11>.RDS
# if seed or BP index are specified, they are prepended to the file name

# load scripts
suppressMessages(source("init.R"))
# set Java options
options(java.parameters = "-Xmx1g")

# load founder dataset
founders <- readRDS("data/founders.RDS")

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
delta.F <- as.numeric(args[10])
it <- as.numeric(args[11])

# store full marker data?
store.all.pops <- as.logical(args[12])
if(is.na(store.all.pops)){
  store.all.pops <- FALSE
}

# set seed
custom.seed <- as.numeric(args[13])
if(is.na(custom.seed)){
  seed <- ceiling(runif(1, 0, 2^31-1))
} else {
  seed <- custom.seed
}
set.seed(seed, kind = "default", normal.kind = "default")
message(sprintf("RNG Seed: %d", seed))
message(sprintf("RNG Kind: %s", getRNG()$kind))
message(sprintf("RNG Normal kind: %s", getRNG()$normal.kind))

# load fixed base population if requested
fixed.bp <- as.numeric(args[14])
if(!is.na(fixed.bp)){
  base.pop <- readRDS(sprintf("data/fixed-bp/bp-%d.RDS", fixed.bp))
} else {
  base.pop <- NULL
}

# run simulation
seasons <- sim.function(founders, heritability, base.pop,
                        num.seasons = num.seasons,
                        QTL.effects = QTL.effects,
                        gp.method = gp.method,
                        add.TP = add.TP,
                        div.measure = div.measure,
                        div.weight = div.weight,
                        type = CGS.type,
                        store.all.pops = store.all.pops,
                        delta.F = delta.F)

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
} else if(sim.function.name == "OC") {
  out.dir <- sprintf("out/%s/%d-seasons/h2-%.1f/addTP-%d/%s-effects/%s/dF-%.5f",
                     sim.function.name, num.seasons, heritability,
                     add.TP, QTL.effects, gp.method, delta.F)
} else {
  # WGS or WGS2
  out.dir <- sprintf("out/%s/%d-seasons/h2-%.1f/addTP-%d/%s-effects/%s",
                     sim.function.name, num.seasons, heritability,
                     add.TP, QTL.effects, gp.method)
}
file.name <- paste(it, "RDS", sep=".")
if(!is.na(custom.seed)){
  file.name <- sprintf("seed-%d-%s", custom.seed, file.name)
}
if(!is.na(fixed.bp)){
  file.name <- sprintf("bp-%d-%s", fixed.bp, file.name)
}
dir.create(out.dir, showWarning = FALSE, recursive = TRUE)
full.path <- paste(out.dir, file.name, sep="/")
saveRDS(seasons, file = full.path)








