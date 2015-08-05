library(Hmisc)

#################################
# SELECTION STRATEGY SIMULATION #
#################################

# simulate phenotypic selection for a number of seasons:
#  - season 0: randomly mate founders & inbreed (DH) to create base population,
#              assign QTL, infer genetic values, fix initial heritability
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
  
  # SEASON 0
  
  message("Season 0: Cross & inbreed founders")
  
  # mate founders to create base population
  base.pop = mate.founders(founders, F1.size, "bp")
  # assign QTL (if not yet assigned)
  if(is.null(base.pop$hypred$realQTL)){
    base.pop = assign.qtl(base.pop, num.QTL, method = QTL.effects)
  } else {
    message("QTL already assigned in founder population, using existing effects")
  }
  # infer genetic values
  base.pop = infer.genetic.values(base.pop)
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
#  - season 0: cross & inbreed founders, assign QTL, infer genetic values and fix heritability
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
    predict.values <- function(model, pop) {
      est.values <- gp.predict(model, gp.design.matrix(pop), infer.weights(model, pop))
      # store in population
      pop$estGeneticValues <- est.values
      # return modified population
      return(pop)
    }
  } else {
    # unweighted GS
    predict.values <- function(model, pop) {
      est.values <- gp.predict(model, gp.design.matrix(pop))
      # store in population
      pop$estGeneticValues <- est.values
      # return modified population
      return(pop)
    }
  }
  
  # wrap GP training to include season number in warning when all markers are fixed 
  train <- function(pheno, Z, method, s){
    tryCatch(gp.train(pheno, Z, method),
             warning = function(w){
               if(grepl("all markers fixed", w)){
                 warning(sprintf("all markers fixed in season %d", s))
               } else {
                 # throw back other warnings
                 warning(w)
               }
             })
  }
  
  # initialize results list (one entry per season 0-n)
  seasons = as.list(rep(NA, num.seasons+1))
  
  # SEASON 0: mate founders to create base population (+ additional training population, if requested)
  
  message("Season 0: Cross & inbreed founders")
  
  # mate founders to create base population
  base.pop = mate.founders(founders, F1.size, "bp")
  # assign QTL (if not yet assigned)
  if(is.null(base.pop$hypred$realQTL)){
    base.pop = assign.qtl(base.pop, num.QTL, method = QTL.effects)
  } else {
    message("QTL already assigned in founder population, using existing effects")
  }
  # infer genetic values
  base.pop = infer.genetic.values(base.pop)
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
  gp.trained.model = train(pheno = tp$pheno, Z = gp.design.matrix(tp), method = gp.method, 1)
  
  # select based on estimated values
  evaluated.base.pop = predict.values(gp.trained.model, evaluated.base.pop)
  selected.names = selection.criterion(evaluated.base.pop$estGeneticValues, num.select)
  selected.pop = restrict.population(evaluated.base.pop, selected.names)
  
  # store season
  season.1 = list(evaluate = list(pop = evaluated.base.pop, add.tp = evaluated.add.TP.pop),
                  select = list(pop.in = evaluated.base.pop, pop.out = selected.pop),
                  gp = list(model = gp.trained.model, tp = tp))
  seasons[[2]] = season.1
  
  # SEASON 2: cross & inbreed selected population, select (based on predictions, no model update)
  
  message("Season 2: Cross, inbreed & select (on predicted values, no model update)")
  
  # cross & inbreed
  offspring = mate.dh(selected.pop, F1.size, "s2")
  
  # select based on estimated values
  offspring = predict.values(gp.trained.model, offspring)
  selected.names = selection.criterion(offspring$estGeneticValues, num.select)
  selected.offspring = restrict.population(offspring, selected.names)
  
  # store season
  season.2 = list(cross.inbreed = list(pop.in = selected.pop, pop.out = offspring),
                  select = list(pop.in = offspring, pop.out = selected.offspring),
                  gp = list(model = gp.trained.model, tp = tp)) # GP model was not updated in season 2
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
    gp.trained.model = train(pheno = tp$pheno, Z = gp.design.matrix(tp), method = gp.method, s)
    # select from offspring based on estimated values using updated GP model
    offspring = predict.values(gp.trained.model, offspring)
    selected.names = selection.criterion(offspring$estGeneticValues, num.select)
    selected.offspring = restrict.population(offspring, selected.names)
    # store season
    new.season = list(evaluate = list(pop = evaluated.prev.offspring),
                      cross.inbreed = list(pop.in = parents, pop.out = offspring),
                      select = list(pop.in = offspring, pop.out = selected.offspring),
                      gp = list(model = gp.trained.model, tp = tp)) # updated GP model from enlarged TP
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

# plot a certain variable extracted from the simulations:
#  --> input:   list of replicates (each replicate is a list of simulated seasons)
#  --> plotted: values of variable extracted with the given function 'extract.values',
#               averaged over all replicates, with confidence intervals (by default 95%)
#               calculated from a normal distribution
# the function 'extract.values' should take a single argument, i.e. the list of seasons
# produced from a single simulation run, and return a vector with one extracted value per season
# (for seasons where the value is not reported, NA should be returned)
plot.simulation.variable <- function(replicates,
                                     extract.values,
                                     ylab, type = c("generations", "seasons"),
                                     ci=0.95, add=FALSE, pch=23,
                                     bg="black", lty=2,
                                     ...){
 
  # check input
  if(!is.function(extract.values)){
    stop("'extract.values' should be a function")
  }
  # get selected type
  type <- match.arg(type)
  
  # infer number of replicates and (maximum) number of seasons
  num.rep <- length(replicates)
  num.seasons <- max(sapply(replicates, length))-1
  # initialize matrix (rep x seasons) to store values
  values <- matrix(NA, num.rep, num.seasons+1)
  
  # go through replicates and fill value matrix
  for(i in 1:num.rep){
    
    # extract seasons of current replicate
    rep.seasons <- replicates[[i]]
    
    # extract variable values from season
    rep.values <- extract.values(rep.seasons)
    
    # convert to generations if requested, by removing NA
    # values, i.e. seasons where no value was reported
    if(type == "generations"){
      rep.values <- rep.values[!is.na(rep.values)]
    }
    
    # store
    values[i, 1:length(rep.values)] <- rep.values
    
  }
  
  # compute averages and standard error + CI across replicates
  value.avg <- colMeans(values)
  value.std.err <- apply(values, 2, sd)/sqrt(num.rep)
  value.ci.halfwidth <- qnorm(ci+(1-ci)/2)*value.std.err
  # only plot CI when large enough to be visible ???
  # value.ci.halfwidth[value.ci.halfwidth < 0.15] = NA
  value.ci.top <- value.avg + value.ci.halfwidth
  value.ci.bottom <- value.avg - value.ci.halfwidth
  # plot non NA values only
  non.na <- !is.na(value.avg)
  value.avg <- value.avg[non.na]
  value.ci.top <- value.ci.top[non.na]
  value.ci.bottom <- value.ci.bottom[non.na]
  x <- (0:num.seasons)[non.na]
  # set x-axis label
  xlab <- ifelse(type == "generations", "Generation", "Season")
  # first plot CI bars (points and lines are plotted on top)
  final.value <- ceil(value.avg[length(value.avg)])
  errbar(x, value.avg, value.ci.top, value.ci.bottom, type="n",
         xlab=xlab, ylab=ylab,
         xaxp=c(0,num.seasons,num.seasons/2),
         add=add,
         ...)
  points(x, value.avg, type="o", pch=pch, bg=bg, lty=lty)
   
}

# plot genetic gain, with one of these scales:
#  1) in terms of number of standard deviations of genetic value in founder population
#  2) normalized wrt maximal genotypic value possible, as in the paper by Jannink
plot.genetic.gain <- function(replicates,
                              scale = c("jannink", "sd"),
                              ylab = "Genetic gain from selection",
                              ...){
  
  # get selected scale
  scale <- match.arg(scale)
  
  # set function to extract genetic gains from a simulated list of seasons
  extract.gain <- function(seasons){
    # extract number of seasons
    num.seasons <- length(seasons)-1
    # initialize gain vector
    gains <- rep(NA, length(seasons))
    # extract general variables
    general <- seasons[[1]]$general
    # extract base pop variables
    base.pop <- seasons[[1]]$candidates
    
    # set mean genetic value of base population
    gains[1] <- mean(base.pop$geneticValues)
    # set mean genetic value of selected populations during simulation
    for(s in 1:num.seasons){
      season <- seasons[[s+1]]
      # compute mean genetic value in selected population
      if(!is.null(season$selection$geneticValues)){
        gains[s+1] <- mean(season$selection$geneticValues)
      }
    }
    
    if (scale == "sd"){
      # subtract values from initial value and divide by intial sd
      gains <- gains - gains[1]
      gains <- gains / sd(base.pop$geneticValues)
    } else {
      # scale according to Jannink: normalize genetic values to [-1,1]
      # based on minimum and maximum possible value
      gains <- normalize.genetic.values(gains, general$qtl.effects)
      # subtract values from initial value
      gains <- gains - gains[1]
    }
    
    # return extracted gains
    return(gains)
    
  }
 
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values = extract.gain, ylab = ylab, ...)
  
}

# plot genetic standard deviation among selection candidates
plot.genetic.standard.deviation <- function(replicates,
                                            ylab = "Genetic standard deviation",
                                            ...){
  
  # set function to extract genetic standard deviation
  extract.genetic.sd <- function(seasons){
    # initialize result vector
    genetic.sd <- rep(NA, length(seasons))
    # extract general variables
    general <- seasons[[1]]$general
    # extract sd for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether selection candidates have been produced in this season
      if(!is.null(season$candidates)){
        # normalize genetic values to [-1,1] based on minimum and maximum possible value
        genetic.values <- normalize.genetic.values(season$candidates$geneticValues, general$qtl.effects)
        # extract and store genetic sd
        genetic.sd[s] <- sd(genetic.values)
      }
    }
    return(genetic.sd)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.genetic.sd, ylab = ylab, ...)
  
}

# plot average inbreeding coefficient among selection candidates
plot.mean.inbreeding <- function(replicates,
                                 ylab = "Mean inbreeding coefficient",
                                 ...){
  
  # set function to extract mean inbreeding coefficient
  extract.mean.inbr <- function(seasons){
    # initialize result vector
    mean.inbr <- rep(NA, length(seasons))
    # extract mean inbreeding for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether selection candidates have been produced in this season
      if(!is.null(season$candidates)){
        # extract and store mean inbreeding coefficient
        mean.inbr[s] <- mean(season$candidates$inbreeding)
      }
    }
    return(mean.inbr)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.mean.inbr, ylab = ylab, ...)
  
}

# plot QTL - marker LD in selection candidates, averaged over all polymorphic QTL
plot.mean.QTL.marker.LD <- function(replicates,
                                    ylab = "Mean QTL - marker LD",
                                    ...){
  
  # set function to extract mean LD
  extract.mean.LD <- function(seasons){
    # initialize result vector
    mean.LD <- rep(NA, length(seasons))
    # extract mean LD for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether selection candidates have been produced in this season
      if(!is.null(season$candidates)){
        # extract and store mean QTL - marker LD
        mean.LD[s] <- mean(season$candidates$QTL.marker.LD$LD, na.rm = TRUE)
      }
    }
    return(mean.LD)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.mean.LD, ylab = ylab, ...)
  
}

# plot marker effect estimation accuracy
plot.effect.estimation.accuracy <- function(replicates,
                                            corrected = TRUE,
                                            ylab = "Marker effect accuracy",
                                            ...){
  
  # set function to extract accuracy
  extract.accuracy <- function(seasons){
    # initialize result vector
    accuracies <- rep(NA, length(seasons))
    # extract accuracy for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether a GP model has been trained in this season
      if(!is.null(season$gp)){
        # extract and store effect estimation accuracy
        acc <- season$gp$effect.estimation.accuracy
        accuracies[s] <- ifelse(corrected, acc$corrected, acc$plain)
      }
    }
    return(accuracies)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.accuracy, ylab = ylab, ...)
  
}

# plot number of favourable QTL lost
plot.num.fav.QTL.lost <- function(replicates,
                                  ylab = "Number favourable QTL lost",
                                  ...){

  # set function to extract number of favourable QTL lost
  extract.num.fav.QTL.lost <- function(seasons){
    # initialize result vector
    num.lost <- rep(NA, length(seasons))
    # extract number lost for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether selection candidates have been produced in this season
      if(!is.null(season$candidates)){
        # extract and store number
        num.lost[s] <- season$candidates$num.fav.QTL.lost
      }
    }
    return(num.lost)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.num.fav.QTL.lost, ylab = ylab, ...)
  
}

# plot mean QTL favourable allele frequency in selection candidates, averaged over all QTL
plot.mean.QTL.fav.allele.freq <- function(replicates,
                                         ylab = "Mean QTL favourable allele frequency",
                                         ...){
  
  # set function to extract mean QTL favourable allele frequency
  extract.mean.QTL.fav.allele.freq <- function(seasons){
    # initialize result vector
    mean.QTL.fav.allele.freq <- rep(NA, length(seasons))
    # extract
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether selection candidates have been produced in this season
      if(!is.null(season$candidates)){
        mean.QTL.fav.allele.freq[s] <- mean(season$candidates$fav.QTL.allele.freqs)
      }
    }
    return(mean.QTL.fav.allele.freq)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.mean.QTL.fav.allele.freq, ylab = ylab, ...)
  
}

# plot mean SNP favourable allele frequency in selection candidates, averaged over all SNP
plot.mean.marker.fav.allele.freq <- function(replicates,
                                             ylab = "Mean marker favourable allele frequency",
                                             ...){
  
  # set function to extract mean marker favourable allele frequency
  extract.mean.marker.fav.allele.freq <- function(seasons){
    # initialize result vector
    mean.marker.fav.allele.freq <- rep(NA, length(seasons))
    # extract
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether season involves GP
      if(!is.null(season$gp)){
        mean.marker.fav.allele.freq[s] <- mean(season$gp$fav.marker.allele.freqs)
      }
    }
    return(mean.marker.fav.allele.freq)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.mean.marker.fav.allele.freq, ylab = ylab, ...)
  
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
  
  # store general variables
  base.pop <- seasons[[1]]$cross.inbreed$pop.out
  # 1) QTL effects
  metadata[[1]]$general$qtl.effects <- get.qtl.effects(base.pop)

  # go through all seasons
  for(s in 0:num.seasons){
    
    # get current and next season (if any)
    season <- seasons[[s+1]]
    next.season <- NULL
    if(s+2 <= num.seasons){
      next.season <- seasons[[s+2]]
    }
    
    ################################################################################
    # TRACK VARIABLES FOR SELECTION CANDIDATES PRODUCED BY CROSSING AND INBREEDING #
    ################################################################################
    
    # check whether season involves crossing and inbreeding
    if(!is.null(season$cross.inbreed)){
      
      # extract produced selection candidates
      candidates <- season$cross.inbreed$pop.out
      # extract evaluated version (from next season, if any)
      evaluated.candidates <- NULL
      if(!is.null(next.season)){
        evaluated.candidates <- next.season$evaluate$pop
      }
      
      # store variables
      
      # 1) genetic values
      metadata[[s+1]]$candidates$geneticValues <- candidates$geneticValues
      # 2) simulated phenotypes (if available)
      if(!is.null(evaluated.candidates)){
        metadata[[s+1]]$candidates$pheno <- evaluated.candidates$pheno
      }
      # 3) estimated genetic values (if any)
      if(!is.null(candidates$estGeneticValues)){
        metadata[[s+1]]$candidates$estGeneticValues <- candidates$estGeneticValues
      }
      # 4) inbreeding coefficients
      metadata[[s+1]]$candidates$inbreeding <- inbreeding.coefficients(candidates)
      # 5) QTL favourable allele frequencies
      metadata[[s+1]]$candidates$fav.QTL.allele.freqs <- get.favourable.qtl.allele.frequencies(candidates)
      # 6) QTL - marker LD
      metadata[[s+1]]$candidates$QTL.marker.LD <- QTL.marker.highest.LD(candidates)
      # 7) number of favourable QTL lost
      num.lost <- sum(get.favourable.qtl.allele.frequencies(candidates) == 0)
      metadata[[s+1]]$candidates$num.fav.QTL.lost <- num.lost
      
    }
    
    #############################################################
    # TRACK VARIABLES FOR SELECTED PARENTS FOR FUTURE CROSSINGS #
    #############################################################
    
    # check whether season involves selection
    if(!is.null(season$select)){
      
      # extract selection
      selection <- season$select$pop.out
      # extract names of selected individuals
      selection.names <- names(selection$geneticValues)
      # extract previously inferred metadata of selection candidates from which selection was made
      if(!is.null(metadata[[s+1]]$candidates)){
        # selection was made from candidates produced in same season (GS/WGS)
        candidates <- metadata[[s+1]]$candidates
      } else {
        # selection was made from candidates produced in previous generation (PS + first season of GS/WGS) 
        candidates <- metadata[[s]]$candidates
      }
      
      # store variables
      
      # 1) genetic values
      metadata[[s+1]]$selection$geneticValues <- selection$geneticValues
      # 2) simulated phenotypes (extract from previously inferred metadata of  selection candidates)
      if(!is.null(candidates$pheno)){
        metadata[[s+1]]$selection$pheno <- candidates$pheno[selection.names]
      }
      # 3) estimated genetic values (extract from previously inferred metadata of  selection candidates)
      if(!is.null(candidates$estGeneticValues)){
        metadata[[s+1]]$selection$estGeneticValues <- candidates$estGeneticValues[selection.names]
      }

    }
    
    ################################
    # TRACK VARIABLES FOR GP MODEL #
    ################################

    # check whether season involves GP    
    if(!is.null(season$gp)){
      
      # extract GP model and TP
      gp.model <- season$gp$model
      gp.tp <- season$gp$tp
      # extract selection candidates for which model is used
      candidates <- season$select$pop.in
      # selection candidates design matrix
      Z <- gp.design.matrix(candidates)
                  
      # store variables
      
      # 1) estimated marker effects
      metadata[[s+1]]$gp$effects <- gp.get.effects(gp.model)
      # 2) marker favourable allele frequencies in selection candidates
      metadata[[s+1]]$gp$fav.marker.allele.freqs <- get.favourable.allele.frequencies(gp.model, Z)
      # 3) marker effect estimation accuracy: computed as correlation between polymorphic QTL effects
      #    and estimated effects of SNP in highest LD, corrected by dividing by average actual LD (in TP)
      QTL.marker.LD <- QTL.marker.highest.LD(gp.tp)
      QTL.marker.LD <- QTL.marker.LD[!is.na(QTL.marker.LD$LD), ]
      # compute mean QTL - marker LD in TP
      mean.LD <- mean(QTL.marker.LD$LD)
      # retrieve polymorphic QTL and corresponding marker effects
      # (!! based on marker *names*, not indices as the latter
      #  include dummies for which no effect was estimated)
      all.names <- rownames(gp.tp$map)
      all.marker.effects <- gp.get.effects(gp.model)
      all.qtl.effects <- get.qtl.effects(gp.tp)
      marker.names <- all.names[QTL.marker.LD$marker.index]
      marker.effects <- all.marker.effects[marker.names]
      qtl.names <- all.names[QTL.marker.LD$QTL.index]
      qtl.effects <- all.qtl.effects[qtl.names]
      # compute and store accuracy (both plain and corrected)
      plain.acc <- cor(marker.effects, qtl.effects)
      metadata[[s+1]]$gp$effect.estimation.accuracy$plain <- plain.acc
      metadata[[s+1]]$gp$effect.estimation.accuracy$corrected <- plain.acc / mean.LD
      
    }
    
  }
  
  return(metadata)
  
}

#############################################
# UTILITY FUNCTIONS FOR METADATA EXTRACTION #
#############################################

# make genomic relationship matrix G from marker matrix Z (0/1 DHs)
genomic.relationship.matrix <- function (Z){
  nSNP <- ncol(Z)
  pfreq <- colMeans(Z)
  Zt <- t(apply(Z, 1, function(x, pfreq) { x-pfreq }, pfreq))
  G <- (Zt%*%t(Zt))/(sum(pfreq*(1-pfreq)))
  return(G)
}

# get coefficients of inbreeding for each individual
inbreeding.coefficients <- function(pop){
  Z <- gp.design.matrix(pop)
  G <- genomic.relationship.matrix(Z)
  coeff <- diag(G)-1
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














