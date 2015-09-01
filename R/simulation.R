
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
              F1.size=200, num.select=20, num.seasons=30,
              selection.criterion=select.highest.score,
              extract.metadata = TRUE,
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
      selected.names = selection.criterion(n = num.select,
                                           values = evaluated.pop$pheno)
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
  
  if(extract.metadata){
    # return simulated seasons metadata
    return(extract.metadata(seasons))
  } else {
    # return full data
    return(seasons)
  }
  
}

WGS <- function(founders, heritability,
               num.QTL=100, QTL.effects = c("normal", "jannink"),
               F1.size=200, add.TP=0, num.select=20, num.seasons=30,
               selection.criterion=select.highest.score,
               gp.method = c("BRR", "RR"), extract.metadata = TRUE,
               ...){
  return(GS(founders, heritability, num.QTL, QTL.effects, F1.size,
            add.TP, num.select, num.seasons, selection.criterion,
            gp.method, extract.metadata, weighted = TRUE, ...))
}
CGS <- function(founders, heritability,
               num.QTL=100, QTL.effects = c("normal", "jannink"),
               F1.size=200, add.TP=0, num.select=20, num.seasons=30,
               gp.method = c("BRR", "RR"), extract.metadata = TRUE,
               div.weight, div.measure = c("MR", "HE", "HEadj"),
               type = c("index", "split"), ...){
  
  # get selected type
  type <- match.arg(type)
  
  if(type == "index"){
    # maximize weighted index
    sel.crit <- function(n, values, markers, fav.alleles){
      select.weighted.index(n, values, markers, div.weight, div.measure, fav.alleles)
    }
  } else {
    # split & combine: highest quality + most diverse individuals
    stop("deprecated")
  }
  
  # run GS with combinatorial selection strategy
  return(GS(founders, heritability, num.QTL, QTL.effects, F1.size,
            add.TP, num.select, num.seasons, selection.criterion = sel.crit,
            gp.method, extract.metadata, weighted = FALSE, ...))
  
}
# simulate genomic selection (possibly weighted by favourable allele frequencies)
#  - season 0: cross & inbreed founders, assign QTL, infer genetic values and fix heritability
#  - season 1: evaluate offspring, train GP & select (on predicted values)
#  - season 2: cross, inbreed & select (on predicted values, no model update)
#  - season >= 3: (1) evaluate previous offspring, update GP model
#               + (2) cross, inbreed & select (on predicted values)
GS <- function(founders, heritability,
              num.QTL=100, QTL.effects = c("normal", "jannink"),
              F1.size=200, add.TP=0, num.select=20, num.seasons=30,
              selection.criterion=select.highest.score,
              gp.method = c("BRR", "RR"), extract.metadata = TRUE,
              weighted = FALSE, ...){
  
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
  
  # initialize results list (one entry per season 0-n)
  seasons = as.list(rep(NA, num.seasons+1))
  
  # SEASON 0: mate founders to create base population (+ additional training population, if requested)
  
  message("Season 0: Cross & inbreed founders")
  
  # mate founders to create base population
  message("|- Cross & inbreed founders: generate base population")
  base.pop = mate.founders(founders, F1.size, "bp")
  # assign QTL (if not yet assigned)
  if(is.null(base.pop$hypred$realQTL)){
    message("|- Assign QTL")
    base.pop = assign.qtl(base.pop, num.QTL, method = QTL.effects)
  } else {
    message("|- QTL already assigned in founder population, using existing effects")
  }
  message("|- Fix heritability to ", heritability)
  # infer genetic values
  base.pop = infer.genetic.values(base.pop)
  # fix heritability (error variation is inferred)
  base.pop = set.heritability(base.pop, heritability)
  # generate additional TP if requested
  if(add.TP > 0){
    message("|- Cross & inbreed founders: generate additional TP")
    # cross founders
    add.TP.pop = mate.founders(founders, add.TP, "add-tp")
    # set same QTL effects as in base population
    add.TP.pop$hypred = base.pop$hypred
    # infer genetic values
    add.TP.pop = infer.genetic.values(add.TP.pop)
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
  
  message("|- Evaluate base population")
  # evaluate base population
  evaluated.base.pop = infer.phenotypes(base.pop)
  # evaluate additional TP, if any
  if(!is.null(add.TP.pop)){
    message("|- Evaluate additional TP")
    evaluated.add.TP.pop = infer.phenotypes(add.TP.pop)
    # combine all training data
    tp = enlarge.tp(evaluated.base.pop, evaluated.add.TP.pop)
  } else {
    evaluated.add.TP.pop = NULL
    # only base population itself serves as TP
    tp = evaluated.base.pop
  }
  
  # train GP
  message("|- Train GP model")
  gp.trained.model = gp.train(pheno = tp$pheno, Z = gp.design.matrix(tp), method = gp.method)
  
  # select based on estimated values
  message("|- Select")
  evaluated.base.pop = predict.values(gp.trained.model, evaluated.base.pop)
  selected.names = selection.criterion(n = num.select,
                                       values = evaluated.base.pop$estGeneticValues,
                                       markers = gp.design.matrix(evaluated.base.pop),
                                       fav.alleles = get.favourable.alleles(gp.get.effects(gp.trained.model)))
  selected.pop = restrict.population(evaluated.base.pop, selected.names)

  # store season
  season.1 = list(evaluate = list(pop = evaluated.base.pop, add.tp = evaluated.add.TP.pop),
                  select = list(pop.in = evaluated.base.pop, pop.out = selected.pop),
                  gp = list(model = gp.trained.model, tp = tp))
  seasons[[2]] = season.1
  
  # SEASON 2: cross & inbreed selected population, select (based on predictions, no model update)
  
  message("Season 2: Cross, inbreed & select (on predicted values, no model update)")
  
  # cross & inbreed
  message("|- Cross & inbreed")
  offspring = mate.dh(selected.pop, F1.size, "s2")

  # select based on estimated values
  message("|- Select (no model update)")
  offspring = predict.values(gp.trained.model, offspring)
  selected.names = selection.criterion(n = num.select,
                                       values = offspring$estGeneticValues,
                                       markers = gp.design.matrix(offspring),
                                       fav.alleles = get.favourable.alleles(gp.get.effects(gp.trained.model)))
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
    message("|- Evaluate selection candidates from previous generation")
    prev.offspring = seasons[[s]]$cross.inbreed$pop.out
    evaluated.prev.offspring = infer.phenotypes(prev.offspring)
    # add evaluated population to TP
    message("|- Enlarge TP")
    tp = enlarge.tp(tp, evaluated.prev.offspring)
    # update GP model
    message("|- Update GP model")
    gp.trained.model = gp.train(pheno = tp$pheno, Z = gp.design.matrix(tp), method = gp.method)
    # cross & inbreed selection from previous offspring
    message("|- Cross & inbreed parents selected in previous generation")
    parents = seasons[[s]]$select$pop.out
    offspring = mate.dh(parents, F1.size, paste("s", s, sep=""))
    # select from offspring based on estimated values using updated GP model
    message("|- Select")
    offspring = predict.values(gp.trained.model, offspring)
    selected.names = selection.criterion(n = num.select,
                                         values = offspring$estGeneticValues,
                                         markers = gp.design.matrix(offspring),
                                         fav.alleles = get.favourable.alleles(gp.get.effects(gp.trained.model)))
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
  
  if(extract.metadata){
    # return simulated seasons metadata
    return(extract.metadata(seasons))
  } else {
    # return full data
    return(seasons)
  }
}

infer.weights <- function(gp.trained.model, pop){
  Z <- gp.design.matrix(pop)
  # weigh according to favourable allele frequencies (inversely proportional)
  weights <- 1/sqrt(get.favourable.allele.frequencies(gp.get.effects(gp.trained.model), Z))
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

#####################################################
# EXTRACT SIMULATION METADATA FROM FULL OUTPUT DATA #
#####################################################

extract.metadata <- function(seasons){
  
  message("Extracting metadata ...")
  
  # initialize result list
  num.seasons <- length(seasons)-1
  metadata <- lapply(1:(num.seasons+1), function(i) {list()})
  
  # initialize complete pedigree
  pedigree <- list()
  
  # store general variables
  base.pop <- seasons[[1]]$cross.inbreed$pop.out
  # 1) QTL effects
  metadata[[1]]$general$qtl.effects <- get.qtl.effects(base.pop)
  # 2) complete pedigree --> after processing each season (see below)
  
  # go through all seasons
  for(s in 0:num.seasons){
    
    # get current and next season (if any)
    season <- seasons[[s+1]]
    next.season <- NULL
    if(s+2 <= (num.seasons+1)){
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
      } else if(!is.null(evaluated.candidates$estGeneticValues)){
        metadata[[s+1]]$candidates$estGeneticValues <- evaluated.candidates$estGeneticValues
      }
      
      # 4) inbreeding coefficients
      
      # 4A) genomic
  
      metadata[[s+1]]$candidates$inbreeding$geno <- genomic.inbreeding.coefficients(candidates)
      
      # 4B) pedigree
      
      # update pedigree data
      pedigree$ID <- c(pedigree$ID, geno.names(candidates))
      pedigree$par1 <- c(pedigree$par1, candidates$pedigree$par1)
      pedigree$par2 <- c(pedigree$par2, candidates$pedigree$par2)
      # convert to synbreed data
      candidates.ped <- create.pedigree(pedigree$ID, pedigree$par1, pedigree$par2, add.ancestors = TRUE)
      candidates.ped <- create.gpData(pedigree = candidates.ped)
      # compute A-matrix
      A <- kin(candidates.ped, ret = "add")
      # infer all inbreeding coefficients
      all.inbreeding.coefs <- diag(A) - 1
      # erase relationship matrix class
      class(all.inbreeding.coefs) <- NULL
      # store inbreeding coefficients of current selection candidates only
      cur.gener <- max(candidates.ped$pedigree$gener)
      metadata[[s+1]]$candidates$inbreeding$ped <- all.inbreeding.coefs[candidates.ped$pedigree$gener == cur.gener]
      
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
      metadata[[s+1]]$gp$fav.marker.allele.freqs <- get.favourable.allele.frequencies(gp.get.effects(gp.model), Z)
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
  
  # store complete pedigree
  metadata[[1]]$general$pedigree <- pedigree
  
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

# get genomic coefficients of inbreeding for each individual
genomic.inbreeding.coefficients <- function(pop){
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














