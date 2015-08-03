library(Hmisc)

#################################
# SELECTION STRATEGY SIMULATION #
#################################

# simulate phenotypic selection for a number of seasons:
#  - preprocessing: assign QTL, infer genetic values
#  - season 0: randomly mate founders & inbreed (DH) to create base population,
#              fix initial heritability
#  - odd seasons >= 1: evaluate (phenotypes) and select
#  - even seasons >= 2: randomly mate selected individuals & inbreed (DH)
# default selection criterion = pure phenotypic mass selection (highest phenotype value)
PS = function(founders, heritability,
              num.QTL=100, QTL.effects = c("normal", "jannink"),
              F1.size=200, num.select=20, num.seasons=26,
              selection.criterion=select.highest.score,
              ...){
  
  # check input
  if(missing(founders)){
    stop("founder population is required")
  }
  if(missing(heritability)
     || is.null(heritability)
     || !is.double(heritability)
     || heritability < 0.0
     || heritability > 1.0){
    stop("heritability is required (real value in [0,1])")
  }
  if(num.seasons < 2 || num.seasons %% 2 != 0){
    stop("number of seasons should be an even number >= 2")
  }
  
  # initialize results list (one entry per season 0-n)
  seasons = as.list(rep(NA, num.seasons+1))
  
  # PREPROCESSING
  
  # assign QTL (if not yet assigned)
  if(is.null(founders$hypred$realQTL)){
    founders = assign.qtl(founders, num.QTL, method = QTL.effects)
  } else {
    message("QTL already assigned in founder population, using existing effects")
  }
  # (re-)infer genetic values
  founders = infer.genetic.values(founders)
  
  # SEASON 0
  
  message("Season 0: Cross & inbreed founders")
  
  # mate founders to create base population
  base.pop = mate.founders(founders, F1.size, "bp")
  # fix heritability
  base.pop = set.heritability(base.pop, heritability)
  
  # store input/output populations of crossing & inbreeding in season 0
  season.0 = list(cross.inbreed = list(pop.in = founders, pop.out = base.pop))
  seasons[[1]] = season.0
      
  # simulate subsequent seasons (1-N)
  for(s in 1:num.seasons){
    start.time = Sys.time()
    if(s %% 2 == 0){
      # even season: cross & inbreed
      message(paste("Season ", s, ": Cross & inbreed selected population", sep=""))
      parents = seasons[[s]]$select$pop.out
      offspring = mate.dh(parents, F1.size, paste("s", s, sep=""))
      # store
      new.season = list(cross.inbreed = list(pop.in = parents, pop.out = offspring))
      seasons[[s+1]] = new.season
    } else {
      # odd season: evaluate & select
      message(paste("Season ", s, ": Evaluate & select", sep=""))
      pop = seasons[[s]]$cross.inbreed$pop.out
      # evaluate (phenotypes)
      evaluated.pop = infer.phenotypes(pop)
      # select
      selected.names = selection.criterion(evaluated.pop$pheno, num.select)
      selected.pop = restrict.population(evaluated.pop, selected.names)
      # store
      new.season = list(evaluate = list(pop = evaluated.pop),
                        select = list(pop.in = evaluated.pop, pop.out = selected.pop))
      seasons[[s+1]] = new.season
    }
    stop.time = Sys.time()
    time = as.numeric(stop.time - start.time, units = "secs")
    message("|- ", time, " seconds elapsed")
  }
  
  # return simulated seasons metadata
  return(extract.metadata(seasons))
  
}

# simulate genomic selection (possibly weighted by favourable allele frequencies)
#  - season 0: cross & inbreed founders
#  - season 1: evaluate offspring, train GP & select (on predicted values)
#  - season 2: cross, inbreed & select (on predicted values, no model update)
#  - season >= 3: (1) evaluate previous offspring, update GP model
#               + (2) cross, inbreed & select (on predicted values)
WGS = function(founders, heritability,
               num.QTL=100, QTL.effects = c("normal", "jannink"),
               F1.size=200, add.TP=0, num.select=20, num.seasons=26,
               selection.criterion=select.highest.score,
               gp.method = c("RR", "BRR")){
  return(GS(founders, heritability, num.QTL, QTL.effects, F1.size,
            add.TP, num.select, num.seasons, selection.criterion,
            gp.method, weighted = TRUE))
}
GS = function(founders, heritability,
              num.QTL=100, QTL.effects = c("normal", "jannink"),
              F1.size=200, add.TP=0, num.select=20, num.seasons=26,
              selection.criterion=select.highest.score,
              gp.method = c("RR", "BRR"), weighted = FALSE){
  
  # check input
  if(missing(founders)){
    stop("founder population is required")
  }
  if(missing(heritability)
     || is.null(heritability)
     || !is.double(heritability)
     || heritability < 0.0
     || heritability > 1.0){
    stop("heritability is required (real value in [0,1])")
  }
  if(num.seasons < 3){
    stop("number of seasons should be >= 3")
  }
  if(add.TP < 0){
    stop("additional training population size (add.TP) should be >= 0")
  }
  
  # set weigthed/unweighted prediction
  if(weighted){
    # weighted GS
    predict.values = function(model, pop) gp.predict(model, gp.design.matrix(pop), infer.weights(model, pop))
  } else {
    # unweighted GS
    predict.values = function(model, pop) gp.predict(model, gp.design.matrix(pop))
  }
  
  # initialize results list (one entry per season 0-n)
  seasons = as.list(rep(NA, num.seasons+1))
  
  # PREPROCESSING
  
  # assign QTL (if not yet assigned)
  if(is.null(founders$hypred$realQTL)){
    founders = assign.qtl(founders, num.QTL, method = QTL.effects)
  } else {
    message("QTL already assigned in founder population, using existing effects")
  }
  # (re-)infer genetic values
  founders = infer.genetic.values(founders)
  
  # SEASON 0: mate founders to create base population (+ additional training population, if requested)
  
  message("Season 0: Cross & inbreed founders")
  
  # mate founders to create base population
  base.pop = mate.founders(founders, F1.size, "bp")
  # fix heritability (error variation is inferred)
  base.pop = set.heritability(base.pop, heritability)
  # generate additional TP if requested
  if(add.TP > 0){
    add.TP.pop = mate.founders(founders, add.TP, "add-tp")
    # set same error variance as base population
    add.TP.pop = set.error.variance(add.TP.pop, base.pop$errorVar)
  } else {
    add.TP.pop = NULL
  }
  
  # store season
  season.0 = list(cross.inbreed = list(pop.in = founders, pop.out = base.pop))
  seasons[[1]] = season.0
  
  # SEASON 1: evaluate base population (and additional TP, if any)
  
  message("Season 1: Evaluate & select (on predicted values)")
  
  # evaluate base population
  evaluated.base.pop = infer.phenotypes(base.pop)
  # evaluate additional TP, if any
  if(!is.null(add.TP.pop)){
    evaluated.add.TP.pop = infer.phenotypes(add.TP.pop)
    # combine all training data
    tp = merge.populations(evaluated.base.pop, evaluated.add.TP.pop)
  } else {
    evaluated.add.TP.pop = NULL
    # only base population itself serves as TP
    tp = evaluated.base.pop
  }
  
  # train GP
  gp.trained.model = gp.train(pheno = tp$pheno, Z = gp.design.matrix(tp), method = gp.method)
  
  # select based on estimated values
  est.values = predict.values(gp.trained.model, evaluated.base.pop)
  selected.names = selection.criterion(est.values, num.select)
  selected.pop = restrict.population(evaluated.base.pop, selected.names)
  
  # store season
  season.1 = list(evaluate = list(pop = evaluated.base.pop, add.tp = evaluated.add.TP.pop),
                  select = list(pop.in = evaluated.base.pop, pop.out = selected.pop),
                  tp = list(pop = tp))
  seasons[[2]] = season.1
  
  # SEASON 2: cross & inbreed selected population, select (based on predictions, no model update)
  
  message("Season 2: Cross, inbreed & select (on predicted values, no model update)")
  
  # cross & inbreed
  offspring = mate.dh(selected.pop, F1.size, "s2")
  # select based on estimated values
  est.values = predict.values(gp.trained.model, offspring)
  selected.names = selection.criterion(est.values, num.select)
  selected.offspring = restrict.population(offspring, selected.names)
  
  # store season
  season.2 = list(cross.inbreed = list(pop.in = selected.pop, pop.out = offspring),
                  select = list(pop.in = offspring, pop.out = selected.offspring),
                  tp = list(pop = tp)) # TP did not change in season 2
  seasons[[3]] = season.2
  
  # iterate over subsequent seasons (3-N)
  for(s in 3:num.seasons){
    # record start time
    start.time = Sys.time()
    message("Season ", s, ": Evaluate (previous) + Cross, inbreed & select")
    # evaluate offspring from previous season
    prev.offspring = seasons[[s]]$cross.inbreed$pop.out
    evaluated.prev.offspring = infer.phenotypes(prev.offspring)
    # add evaluated population to TP
    tp = merge.populations(tp, evaluated.prev.offspring)
    # cross & inbreed selection from previous offspring
    parents = seasons[[s]]$select$pop.out
    offspring = mate.dh(parents, F1.size, paste("s", s, sep=""))
    # update GP model
    gp.trained.model = gp.train(pheno = tp$pheno, Z = gp.design.matrix(tp), method = gp.method)
    # select from offspring based on estimated values using updated GP model
    est.values = predict.values(gp.trained.model, offspring)
    selected.names = selection.criterion(est.values, num.select)
    selected.offspring = restrict.population(offspring, selected.names)
    # store season
    new.season = list(evaluate = list(add.tp = evaluated.prev.offspring),
                      cross.inbreed = list(pop.in = parents, pop.out = offspring),
                      select = list(pop.in = offspring, pop.out = selected.offspring),
                      tp = list(pop = tp)) # enlarged TP
    seasons[[s+1]] = new.season
    # record stop time
    stop.time = Sys.time()
    # print time spent on this generation
    time = as.numeric(stop.time - start.time, units = "secs")
    message("|- ", time, " seconds elapsed")
  }
  
  # return simulated seasons metadata
  return(extract.metadata(seasons))
}

infer.weights = function(gp.trained.model, pop){
  Z = gp.design.matrix(pop)
  # weigh according to favourable allele frequencies (inversely proportional)
  weights = 1/sqrt(get.favourable.allele.frequencies(gp.trained.model, Z))
  # set Inf's to zero (fixed alleles no longer influence GP anyway)
  weights[is.infinite(weights)] = 0
  return(weights)
}

# replicate simulation of any selection strategy
replicate.simulation = function(num.rep = 100, simulate){
  # initialize output list
  replicates = as.list(rep(NA, num.rep))
  # run simulations
  for(r in 1:num.rep){
    message("---------------\nReplication ", r, "\n---------------")
    replicates[[r]] = simulate()
  }
  return(replicates)
}

#########
# PLOTS #
#########

# plot genetic gain, with one of these scales:
# 1) in terms of number of standard deviations of genetic value in founder population
# 2) normalized wrt maximal genotypic value possible, as in the paper by Jannink
#  --> input:   list of replicates (each replicate is a list of simulated seasons)
#  --> plotted: average genetic gain with confidence intervals (by default 95%)
#               calculated from normal distribution
plot.genetic.gain <- function(replicates,
                              type = c("generations", "seasons"),
                              scale = c("jannink", "sd"),
                              ci=0.95, add=FALSE, pch=23,
                              bg="black", lty=2){
  # get selected options
  type <- match.arg(type)
  scale <- match.arg(scale)
  # infer number of replicates and (maximum) number of seasons
  num.rep <- length(replicates)
  num.seasons <- max(sapply(replicates, length))-1
  # initialize matrix (rep x seasons) to store gains
  gains <- matrix(NA, num.rep, num.seasons+1)
  # go through replicates and fill gains matrix
  for(i in 1:num.rep){
    
    # extract seasons of current replicate
    seasons <- replicates[[i]]
    
    # extract founder data
    founders <- seasons[[1]]
    
    # extract mean genetic value of founder population
    gains[i,1] <- mean(founders$geneticValues)
    # extract mean genetic value of selected populations during simulation
    for(s in 1:num.seasons){
      season <- seasons[[s+1]]
      # compute mean genetic value in selected population
      if(!is.null(season$geneticValues)){
        gains[i,s+1] <- mean(season$geneticValues)
      }
    }
    
    if (scale == "sd"){
      # subtract values from initial value and divide by intial sd
      gains[i,] <- gains[i,] - gains[i,1]
      gains[i,] <- gains[i,] / sd(founders$geneticValues)
    } else if (scale == "jannink") {
      # normalize genetic values to [-1,1] based on minimum and maximum possible value
      gains[i,] <- normalize.genetic.values(gains[i,], founders$qtl.effects)
      # subtract values from initial value
      gains[i,] <- gains[i,] - gains[i,1]
    } else {
      stop(sprintf("Unknown scale option %s (should not happen)", scale))
    }
    
    # when plotting generations, remove NA's introduced because of season timing
    if(type == "generations"){
      gen.gains <- gains[i,!is.na(gains[i,])]
      gains[i,] <- NA
      gains[i,1:length(gen.gains)] <- gen.gains
    }
    
  }
  # compute averages and standard error + CI across replicates
  gain.avg <- colMeans(gains)
  gain.std.err <- apply(gains, 2, sd)/sqrt(num.rep)
  gain.ci.halfwidth <- qnorm(ci+(1-ci)/2)*gain.std.err
  # only plot CI when large enough to be visible ???
  # gain.ci.halfwidth[gain.ci.halfwidth < 0.15] = NA
  gain.ci.top <- gain.avg + gain.ci.halfwidth
  gain.ci.bottom <- gain.avg - gain.ci.halfwidth
  # retain only non NA values
  non.na <- !is.na(gain.avg)
  gain.avg <- gain.avg[non.na]
  gain.ci.top <- gain.ci.top[non.na]
  gain.ci.bottom <- gain.ci.bottom[non.na]
  x <- (0:num.seasons)[non.na]
  # set x-axis label
  if(type == "generations"){
    xlab <- "Generation"
  } else {
    xlab <- "Season"
  }
  # first plot CI bars (points and lines are plotted on top)
  final.gain <- ceil(gain.avg[length(gain.avg)])
  errbar(x, gain.avg, gain.ci.top, gain.ci.bottom, type="n",
         xlab=xlab, ylab="Mean genetic gain from selection",
         xaxp=c(0,num.seasons,num.seasons/2),
         add=add)
  points(x, gain.avg, type="o", pch=pch, bg=bg, lty=lty)
}

################
# LOAD RESULTS #
################

# create list containing simulated breeding cycles read from all RDS files found in the given directory
load.simulation.results <- function(dir) {
  # get file paths
  files <- Sys.glob(sprintf("%s/*.RDS", dir))
  # load list of simulated breeding cycles
  breeding.cycles <- lapply(files, readRDS)
  # return results
  return(breeding.cycles)
}

#####################################################
# EXTRACT SIMULATION METADATA FROM FULL OUTPUT DATA #
#####################################################

extract.metadata <- function(seasons){
  
  # initialize result list
  num.seasons <- length(seasons)-1
  metadata <- lapply(1:(num.seasons+1), function(i) {list()} )
  
  # store genetic values of founders and QTL effects
  founders <- seasons[[1]]$cross.inbreed$pop.in
  metadata[[1]]$qtl.effects <- get.qtl.effects(founders)
  metadata[[1]]$geneticValues <- founders$geneticValues

  # go through all seasons
  for(s in 1:num.seasons){
    
    # get season
    season <- seasons[[s+1]]
    
    # extract selected parents (if seasons involves crossing & inbreeding)
    selection <- NULL
    if(!is.null(season$cross.inbreed)){
      selection <- season$cross.inbreed$pop.in
    }
    
    # store genetic values of selected parents (if applicable)
    if(!is.null(selection)){
      metadata[[s+1]]$geneticValues <- selection$geneticValues
    }
    
  }
  
  return(metadata)
  
}

##################################################
# CONVERT FULL OUTPUT DATA TO EXTRACTED METADATA #
##################################################

# converts all .RDS files in the given directory (overwrites original files)
convert.to.metadata <- function(dir){
  
  # get file paths
  files <- Sys.glob(sprintf("%s/*.RDS", dir))
  
  for(file in files){
    message(sprintf("Converting file %s ...", file))
    # read file
    seasons <- readRDS(file)
    # extract metadata
    metadata <- extract.metadata(seasons)
    # clear seasons
    rm(seasons)
    # overwrite file
    saveRDS(metadata, file = file)
  }
  
}














