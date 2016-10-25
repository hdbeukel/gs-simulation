
#################################
# SELECTION STRATEGY SIMULATION #
#################################

# simulate phenotypic selection for a number of seasons:
#  - season 0: randomly mate founders & inbreed (DH) to create base population,
#              assign QTL, infer genetic values, fix initial heritability
#  - odd seasons >= 1: evaluate (phenotypes) and select
#  - even seasons >= 2: randomly mate selected individuals & inbreed (DH)
# default selection criterion = pure phenotypic mass selection (highest phenotype value)
PS <- function(founders, heritability, base.pop = NULL,
               num.QTL=1000, QTL.effects = c("normal", "jannink"),
               F1.size=200, num.select=20, num.seasons=30,
               selection.criterion=select.highest.score,
               extract.metadata = TRUE,
               store.all.pops = FALSE,
               ...){
  
  # check input
  if(missing(founders)){
    stop("founder population is required")
  }
  if(!missing(heritability) && !is.null(heritability) && !is.null(base.pop)){
    message("both 'base.pop' and 'heritability' specified, ignoring 'heritability'")
  }
  if(is.null(base.pop) && (missing(heritability)
                           || is.null(heritability)
                           || !is.double(heritability)
                           || heritability < 0.0
                           || heritability > 1.0)){
    stop("either 'base.pop' or 'heritability' (real value in [0,1]) is required")
  }
  if(num.seasons < 2 || num.seasons %% 2 != 0){
    stop("number of seasons should be an even number >= 2")
  }
  
  # initialize results list (one entry per season 0-n)
  seasons <- as.list(rep(NA, num.seasons+1))
  
  # SEASON 0
  
  message("Season 0: Cross & inbreed founders")
  
  if(is.null(base.pop)){
    # create base population by crossing founders + inbreeding
    base.pop <- create.base.population(founders, F1.size, num.QTL, QTL.effects)
  }
  
  # store input/output populations of crossing & inbreeding in season 0
  season.0 <- list(cross.inbreed = list(pop.in = founders, pop.out = base.pop))
  seasons[[1]] <- season.0
      
  # simulate subsequent seasons (1-N)
  for(s in 1:num.seasons){
    start.time <- Sys.time()
    if(s %% 2 == 0){
      # even season: cross & inbreed
      message(paste("Season ", s, ": Cross & inbreed selected population", sep=""))
      parents <- seasons[[s]]$select$pop.out
      offspring <- mate.dh(parents, F1.size, paste("s", s, sep=""))
      # store
      new.season <- list(cross.inbreed = list(pop.in = parents, pop.out = offspring))
      seasons[[s+1]] <- new.season
    } else {
      # odd season: evaluate & select
      message(paste("Season ", s, ": Evaluate & select", sep=""))
      pop <- seasons[[s]]$cross.inbreed$pop.out
      # evaluate (phenotypes)
      evaluated.pop <- infer.phenotypes(pop)
      # select
      selected.names <- selection.criterion(n = num.select,
                                           values = evaluated.pop$pheno)
      selected.pop <- restrict.population(evaluated.pop, selected.names)
      # store
      new.season <- list(evaluate = list(pop = evaluated.pop),
                        select = list(pop.in = evaluated.pop, pop.out = selected.pop))
      seasons[[s+1]] <- new.season
    }
    stop.time <- Sys.time()
    time <- as.numeric(stop.time - start.time, units = "secs")
    message("|- ", time, " seconds elapsed")
  }
  
  if(extract.metadata){
    # return simulated seasons metadata
    return(extract.metadata(seasons, store.all.pops))
  } else {
    # return full data
    return(seasons)
  }
  
}

# weighted genomic selection (cfr. Jannink)
WGS <- function(founders, heritability, base.pop = NULL,
               num.QTL=1000, QTL.effects = c("normal", "jannink"),
               F1.size=200, add.TP=0, num.select=20, num.seasons=30,
               gp.method = c("BRR", "RR"), extract.metadata = TRUE,
               store.all.pops = FALSE, ...){
  return(GS(founders, heritability, base.pop, num.QTL, QTL.effects, F1.size,
            add.TP, num.select, num.seasons, selection.criterion = select.highest.score,
            gp.method, extract.metadata, store.all.pops, weights = weights.jannink, ...))
}
# weighted genomic selection (cfr. Liu and Woolliams)
WGS2 <- function(founders, heritability, base.pop = NULL,
                num.QTL=1000, QTL.effects = c("normal", "jannink"),
                F1.size=200, add.TP=0, num.select=20, num.seasons=30,
                gp.method = c("BRR", "RR"), extract.metadata = TRUE,
                store.all.pops = FALSE, ...){
  return(GS(founders, heritability, base.pop, num.QTL, QTL.effects, F1.size,
            add.TP, num.select, num.seasons, selection.criterion = select.highest.score,
            gp.method, extract.metadata, store.all.pops, weights = weights.liu.woolliams, ...))
}

# combinatorial genomic selection
CGS <- function(founders, heritability, base.pop = NULL,
               num.QTL=1000, QTL.effects = c("normal", "jannink"),
               F1.size=200, add.TP=0, num.select=20, num.seasons=30,
               gp.method = c("BRR", "RR"), extract.metadata = TRUE,
               store.all.pops = FALSE, div.weight,
               div.measure = c("HEall", "HEfav", "LOGall", "LOGfav"),
               type = c("index", "split"), ...){
  
  # get selected type
  type <- match.arg(type)
  
  if(type == "index"){
    # maximize weighted index
    sel.crit <- function(n, values, markers, fav.alleles, ...){
      select.weighted.index(n = n,
                            values = values,
                            markers = markers,
                            div.weight = div.weight,
                            div.measure = div.measure,
                            fav.alleles = fav.alleles)
    }
  } else {
    # split & combine: highest quality + most diverse individuals
    stop("deprecated")
  }
  
  # run GS with combinatorial selection strategy
  return(GS(founders, heritability, base.pop, num.QTL, QTL.effects, F1.size,
            add.TP, num.select, num.seasons, selection.criterion = sel.crit,
            gp.method, extract.metadata, store.all.pops, ...))
  
}

# optimal contributions simulation
OC <- function(founders, heritability, base.pop = NULL,
               num.QTL=1000, QTL.effects = c("normal", "jannink"),
               F1.size=200, add.TP=0, num.select=20, num.seasons=30,
               gp.method = c("BRR", "RR"), extract.metadata = TRUE,
               store.all.pops = FALSE, delta.F, verbose = FALSE, ...){

  sel.crit <- function(n, values, markers, generation, ...){
    select.fixed.size.oc(n = n,
                        values = values,
                        markers = markers,
                        generation = generation,
                        delta.F = delta.F,
                        verbose = verbose)
  }
  
  return(GS(founders, heritability, base.pop, num.QTL, QTL.effects, F1.size,
            add.TP, num.select, num.seasons, selection.criterion = sel.crit,
            gp.method, extract.metadata, store.all.pops, ...))
  
}

# simulate genomic selection (possibly weighted by favourable allele frequencies)
#  - season 0: cross & inbreed founders, assign QTL, infer genetic values and fix heritability
#  - season 1: evaluate offspring, train GP & select (on predicted values)
#  - season 2: cross, inbreed & select (on predicted values, no model update)
#  - season >= 3: (1) evaluate previous offspring, update GP model
#               + (2) cross, inbreed & select (on predicted values)
GS <- function(founders, heritability, base.pop = NULL,
              num.QTL=1000, QTL.effects = c("normal", "jannink"),
              F1.size=200, add.TP=0, num.select=20, num.seasons=30,
              selection.criterion = select.highest.score,
              gp.method = c("BRR", "RR"), extract.metadata = TRUE,
              store.all.pops = FALSE, weights = weights.none,
              mating.probability = equal.contributions, ...){
  
  # check input
  if(missing(founders)){
    stop("founder population is required")
  }
  if(missing(heritability)
       || is.null(heritability)
       || !is.double(heritability)
       || heritability < 0.0
       || heritability > 1.0){
    stop("'heritability' is required (real value in [0,1])")
  }
  if(num.seasons < 3){
    stop("number of seasons should be >= 3")
  }
  if(add.TP < 0){
    stop("additional training population size (add.TP) should be >= 0")
  }
  
  # wrap prediction to store predicted values
  predict.values <- function(model, pop, generation) {
    est.values <- gp.predict(model, gp.design.matrix(pop), weights(model, pop, t = generation, N = num.seasons))
    # store in population
    pop$estGeneticValues <- est.values
    # return modified population
    return(pop)
  }
  
  # initialize results list (one entry per season 0-n)
  seasons <- as.list(rep(NA, num.seasons+1))
  
  # SEASON 0: mate founders to create base population (+ additional training population, if requested)
  
  message("Season 0: Cross & inbreed founders")
  
  if(is.null(base.pop)){
    # create base population by crossing founders + inbreeding
    base.pop <- create.base.population(founders, F1.size, num.QTL, QTL.effects)
  }
  
  # generate additional TP if requested
  if(add.TP > 0){
    message("|- Cross & inbreed founders: generate additional TP")
    # cross founders
    add.TP.pop <- mate.founders(founders, add.TP, "add-tp")
    # set same QTL effects as in base population
    add.TP.pop$hypred <- base.pop$hypred
    # infer genetic values
    add.TP.pop <- infer.genetic.values(add.TP.pop)
  } else {
    add.TP.pop <- NULL
  }
  
  # store season
  season.0 <- list(cross.inbreed = list(pop.in = founders, pop.out = base.pop))
  seasons[[1]] <- season.0
  
  # SEASON 1: evaluate base population (and additional TP, if any)
  
  message("Season 1: Evaluate & select (on predicted values)")
  
  message("|- Evaluate base population")
  # evaluate base population
  message("|- Fix heritability to ", heritability)
  # fix heritability (error variation is inferred)
  base.pop <- set.heritability(base.pop, heritability)
  # simulate phenotypes
  evaluated.base.pop <- infer.phenotypes(base.pop)
  # init TP from base population
  tp <- init.tp(evaluated.base.pop)
  # evaluate and add additional TP, if any
  if(!is.null(add.TP.pop)){
    message("|- Evaluate additional TP")
    # set same error variance as base population
    add.TP.pop <- set.error.variance(add.TP.pop, base.pop$errorVar)
    # simulate phenotypes
    evaluated.add.TP.pop <- infer.phenotypes(add.TP.pop)
    # combine all training data
    tp <- enlarge.tp(tp, evaluated.add.TP.pop)
  } else {
    evaluated.add.TP.pop <- NULL
  }
  
  # train GP
  message("|- Train GP model")
  gp.trained.model <- gp.train(pheno = tp$pheno, Z = gp.design.matrix(tp), method = gp.method)
  
  # select based on estimated values
  message("|- Select")
  evaluated.base.pop <- predict.values(gp.trained.model, evaluated.base.pop, generation = 1)
  selected.names <- selection.criterion(n = num.select,
                                        values = evaluated.base.pop$estGeneticValues,
                                        markers = gp.design.matrix(evaluated.base.pop),
                                        fav.alleles = get.favourable.alleles(gp.get.effects(gp.trained.model)),
                                        generation = 1)
  selected.pop <- restrict.population(evaluated.base.pop, selected.names)

  # store season
  season.1 <- list(evaluate = list(pop = evaluated.base.pop, add.tp = evaluated.add.TP.pop),
                  select = list(pop.in = evaluated.base.pop, pop.out = selected.pop),
                  gp = list(model = gp.trained.model, tp = tp))
  seasons[[2]] <- season.1
  
  # SEASON 2: cross & inbreed selected population, select (based on predictions, no model update)
  
  message("Season 2: Cross, inbreed & select (on predicted values, no model update)")
  
  # cross & inbreed
  message("|- Cross & inbreed")
  offspring <- mate.dh(selected.pop, F1.size, "s2",
                       probs = mating.probability(
                         values = selected.pop$estGeneticValues,
                         markers = gp.design.matrix(selected.pop),
                         generation = 1
                       ))

  # select based on estimated values
  message("|- Select (no model update)")
  offspring <- predict.values(gp.trained.model, offspring, generation = 2)
  selected.names <- selection.criterion(n = num.select,
                                       values = offspring$estGeneticValues,
                                       markers = gp.design.matrix(offspring),
                                       fav.alleles = get.favourable.alleles(gp.get.effects(gp.trained.model)),
                                       generation = 2)
  selected.offspring <- restrict.population(offspring, selected.names)
  
  # store season
  season.2 <- list(cross.inbreed = list(pop.in = selected.pop, pop.out = offspring),
                  select = list(pop.in = offspring, pop.out = selected.offspring),
                  gp = list(model = gp.trained.model, tp = tp)) # GP model was not updated in season 2
  seasons[[3]] <- season.2
  
  # iterate over subsequent seasons (3-N)
  for(s in 3:num.seasons){
    # record start time
    start.time <- Sys.time()
    message("Season ", s, ": Evaluate (previous) + Cross, inbreed & select")
    # evaluate offspring from previous season
    message("|- Evaluate selection candidates from previous generation")
    prev.offspring <- seasons[[s]]$cross.inbreed$pop.out
    evaluated.prev.offspring <- infer.phenotypes(prev.offspring)
    # add evaluated population to TP
    message("|- Enlarge TP")
    tp <- enlarge.tp(tp, evaluated.prev.offspring)
    # update GP model
    message("|- Update GP model")
    gp.trained.model <- gp.train(pheno = tp$pheno, Z = gp.design.matrix(tp), method = gp.method)
    # cross & inbreed selection from previous offspring
    message("|- Cross & inbreed parents selected in previous generation")
    parents <- seasons[[s]]$select$pop.out
    offspring <- mate.dh(parents, F1.size, paste("s", s, sep=""),
                         probs = mating.probability(
                           values = parents$estGeneticValues,
                           markers = gp.design.matrix(parents),
                           generation = s-1
                         ))
    # select from offspring based on estimated values using updated GP model
    message("|- Select")
    offspring <- predict.values(gp.trained.model, offspring, generation = s)
    selected.names <- selection.criterion(n = num.select,
                                         values = offspring$estGeneticValues,
                                         markers = gp.design.matrix(offspring),
                                         fav.alleles = get.favourable.alleles(gp.get.effects(gp.trained.model)),
                                         generation = s)
    selected.offspring <- restrict.population(offspring, selected.names)
    # store season
    new.season <- list(evaluate = list(pop = evaluated.prev.offspring),
                      cross.inbreed = list(pop.in = parents, pop.out = offspring),
                      select = list(pop.in = offspring, pop.out = selected.offspring),
                      gp = list(model = gp.trained.model, tp = tp)) # updated GP model from enlarged TP
    seasons[[s+1]] <- new.season
    # record stop time
    stop.time <- Sys.time()
    # print time spent on this generation
    time <- as.numeric(stop.time - start.time, units = "secs")
    message("|- ", time, " seconds elapsed")
  }
  
  if(extract.metadata){
    # return simulated seasons metadata
    return(extract.metadata(seasons, store.all.pops))
  } else {
    # return full data
    return(seasons)
  }
}

create.base.population <- function(founders, num.ind, num.QTL, QTL.effects){
  # mate founders to create base population
  message("|- Cross & inbreed founders: generate base population")
  base.pop <- mate.founders(founders, num.ind, "bp")
  # assign QTL (if not yet assigned)
  if(is.null(base.pop$hypred$realQTL)){
    message("|- Assign QTL")
    base.pop <- assign.qtl(base.pop, num.QTL, method = QTL.effects)
  } else {
    message("|- QTL already assigned in founder population, using existing effects")
  }
  # infer genetic values
  base.pop <- infer.genetic.values(base.pop)
  return(base.pop)
}

weights.none <- function(gp.trained.model, ...){
  effects <- gp.get.effects(gp.trained.model)
  weights <- rep(1, times = length(effects))
  return(weights)
}

weights.jannink <- function(gp.trained.model, pop, ...){
  Z <- gp.design.matrix(pop)
  # weigh according to Jannink (2010)
  p <- get.favourable.allele.frequencies(gp.get.effects(gp.trained.model), Z)
  weights <- 1/sqrt(p)
  # set Inf's to zero (fixed markers no longer influence GP anyway)
  weights[is.infinite(weights)] = 0
  return(weights)
}

weights.liu.woolliams <- function(gp.trained.model, pop, t, N){
  Z <- gp.design.matrix(pop)
  # weigh according to formulas by Liu and Woolliams (2010)
  p <- get.favourable.allele.frequencies(gp.get.effects(gp.trained.model), Z)
  B <- asin(1 - 2*p)
  A <- - 1/(N-t+1) * (pi/2 + B)
  p.next <- 1/2 * (1 - sin(A + B))
  weights <- (p.next - p)/sqrt(1/2 * p * (1-p))
  # set Inf's and NaN's to zero (fixed markers no longer influence GP anyway)
  weights[is.infinite(weights) | is.nan(weights)] = 0
  return(weights)
}

# replicate simulation of any selection strategy
replicate.simulation <- function(num.rep = 100, simulate){
  # initialize output list
  replicates <- as.list(rep(NA, num.rep))
  # run simulations
  for(r in 1:num.rep){
    message("---------------\nReplication ", r, "\n---------------")
    replicates[[r]] <- simulate()
  }
  return(replicates)
}

# returns vector of equal probabilities [1/n, ..., 1/n]'
equal.contributions <- function(values, ...){
  n <- length(values)
  probs <- rep(1/n, n)
  return(probs)
}

# returns vector with optimal contributions, maximizing expected response
# with a predefined rate of inbreeding according to Meuwissen 1997
# (adjusted for plant breeding with constraint 1'c_t = 1 instead of Q'c_t = 1/2)
optimal.contributions <- function(values, markers, C, cmin = 0, cmax = 1, verbose = FALSE){
  
  message("[OC] Maximum expected coancestry = ", C)
  
  tol <- 1e-10
  
  all.values <- values
  all.markers <- markers

  n.all <- length(values)
  
  # indices of all individuals
  i.all <- 1:n.all
  # indices of individuals whose contribution is optimized (initially all)
  i.opt <- i.all
  # indices of individuals whose contribution has been
  # fixed to the allowed maximum or zero (initially none)
  i.fixed.max <- numeric(0)
  i.fixed.zero <- numeric(0)
  
  # init vector to store computed optimal contributions
  c <- rep(-1, n.all)

  # compute G
  G <- genomic.relationship.matrix(all.markers)
  
  # continue until cmin and cmax are satisfied
  while(length(i.opt) > 0 && (any(c > cmax + tol) || any(c < 0) || any(c[c > 0] < cmin - tol))){
    
    # progress message
    if(verbose){
      message(sprintf(
        "Remaining: %d - Fixed to zero: %d - Fixed to cmax: %d",
        length(i.opt), length(i.fixed.zero), length(i.fixed.max)
      ))
    }
    
    # set previously fixed contributions
    i.fixed <- c(i.fixed.max, i.fixed.zero)
    c.fixed <- c(rep(cmax, length(i.fixed.max)), rep(0, length(i.fixed.zero)))
    
    if(length(i.opt) == 1){
      
      # one individual left: assign remaining contribution
      c.opt <- 1 - sum(c.fixed)
      
    } else {
    
      # retrieve values and marker data of individuals
      # whose contribution remains to be optimized
      values <- all.values[i.opt]
      markers <- all.markers[i.opt,]
      n <- length(values)
      m <- ncol(markers)
      
      if(isTRUE(all.equal(apply(markers, 2, sd), rep(0, m), check.names = FALSE))){
        
        # candidates fixed at all markers: equally divide remaining contribution
        c.remaining <- 1 - sum(c.fixed)
        c.opt <- rep(c.remaining/n, n)
        
      } else {
        
        # subset genomic relationship matrix
        Goo <- G[i.opt, i.opt, drop = F]
        Gof <- G[i.opt, i.fixed, drop = F]
        Gfo <- G[i.fixed, i.opt, drop = F]
        Gff <- G[i.fixed, i.fixed, drop = F]
        # invert Goo
        Goo.inv <- solve(make.positive.definite(Goo))
        
        # precompute some values
        K <- as.numeric(2*C - t(c.fixed) %*% Gff %*% c.fixed)
        s <- 1 - sum(c.fixed)
        Goo.inv.sum <- sum(Goo.inv)
        ones <- rep(1, n)
        P <- Goo.inv - (Goo.inv %*% ones %*% t(ones) %*% Goo.inv) / Goo.inv.sum
        
        # compute lambda 0
        a <- as.numeric(t(values) %*% P %*% values)
        b <- as.numeric(
          K + t(c.fixed) %*% Gfo %*% P %*% Gof %*% c.fixed
          - s^2/Goo.inv.sum - (2*s * t(ones) %*% Goo.inv %*% Gof %*% c.fixed)/Goo.inv.sum
        )
        lambda.0.sq <- 1/4 * a/b
        if(lambda.0.sq < 0){
          stop("Constraints can not be satisfied")
        }
        lambda.0 <- sqrt(lambda.0.sq)
        
        # compute lambda
        a <- as.numeric(t(ones) %*% Goo.inv %*% (values - 2*lambda.0 * Gof %*% c.fixed) - 2*lambda.0*s)
        lambda <- a / Goo.inv.sum
        
        # compute optimal contributions
        c.opt <- Goo.inv %*% (values - 2*lambda.0 * Gof %*% c.fixed - lambda) / (2*lambda.0)
        
      }
      
    }
    
    # merge fixed and optimized contributions
    c[i.opt] <- c.opt
    c[i.fixed] <- c.fixed
    
    # eliminate most negative contribution, if any
    if(any(c < 0)){
      elim <- which.min(c)
      i.fixed.zero <- c(i.fixed.zero, elim)
      i.opt <- setdiff(i.opt, elim)
    } else {
      # truncate highest contribution if above maximum allowed
      if(any(c > cmax + tol)){
        trunc <- which.max(c)
        i.fixed.max <- c(i.fixed.max, trunc)
        # reset those fixed to zero after truncating too large contribution
        # (so that others are reconsidered for inclusion to decrease inbreeding)
        i.fixed.zero <- numeric(0)
        i.opt <- setdiff(i.all, i.fixed.max)
      } else {
        # eliminate smallest positive contribution below minimum, if any
        if(any(c[c > 0] < cmin - tol)){
          c.pos <- rep(NA, length(c))
          c.pos[c > 0] <- c[c > 0]
          elim <- which.min(c.pos)
          i.fixed.zero <- c(i.fixed.zero, elim)
          i.opt <- setdiff(i.opt, elim)
        }
      }
    }
          
  }
  
  # print some messages
  if(length(i.fixed.zero) > 0){
    message(sprintf("[OC] Set contribution of %d/%d individuals to zero", length(i.fixed.zero), n.all))
  }
  if(length(i.fixed.max) > 0){
    message(sprintf(
      "[OC] Truncated contribution of %d/%d individuals to cmax = %.2f",
      length(i.fixed.max), n.all, cmax
    ))
  }
  
  # return final contributions
  names(c) <- names(all.values)
  return(c)
  
}

# make genomic relationship matrix G from marker matrix M (0/1 DHs)
genomic.relationship.matrix <- function (M){
  # convert to 0/2 format
  M <- 2*M
  # center
  pfreq <- colMeans(M)/2
  Z <- t(apply(M, 1, function(row) { row - 2*pfreq }))
  # compute G
  # G <- Z %*% t(Z) / (2*sum(pfreq*(1-pfreq))) # Van Raden 2008
  G <- Z %*% t(Z) / ncol(M) # Sonesson 2012
  return(G)
}

make.positive.definite <- function(X, tol = 1e-6) {
  eig <- eigen(X, symmetric = TRUE)
  rtol <- tol * eig$values[1]
  if(eig$values[length(eig$values)] < rtol) {
    vals <- eig$values
    vals[vals < rtol] <- rtol
    Srev <- eig$vectors %*% (vals * t(eig$vectors))
    dimnames(Srev) <- dimnames(X)
    return(Srev)
  } else {
    return(X)
  }
}

# compute inbreeding rate delta F based on markers (cfr. Sonesson)
inbreeding.rate.sonesson <- function(pop, prev.pop){
  # get minor allele frequencies of of both populations
  Z.cur <- gp.design.matrix(pop)
  Z.prev <- gp.design.matrix(prev.pop)
  freqs.cur <- maf(Z.cur, encoding = "dh")
  freqs.prev <- maf(Z.prev, encoding = "dh")
  # compute average homozygosity
  hom <- function(f){
    f^2 + (1-f)^2
  }
  hom.cur <- mean(sapply(freqs.cur, hom))
  hom.prev <- mean(sapply(freqs.prev, hom))
  # delta F is relative increase in average homozygosity
  delta.F <- (hom.cur - hom.prev) / (1.0 - hom.prev)
  return(delta.F)
}

# compute inbreeding rate delta F based on markers (cfr. Jannink)
inbreeding.rate.jannink <- function(pop, prev.pop){
  cur.fixed <- proportion.fixed.markers(pop)
  prev.fixed <- proportion.fixed.markers(prev.pop)
  delta.F <- (cur.fixed - prev.fixed) / (1.0 - prev.fixed)
  return(delta.F)
}

# select inbreeding rate formula
inbreeding.rate <- inbreeding.rate.sonesson

proportion.fixed.markers <- function(pop){
  M <- gp.design.matrix(pop)
  polymorphic <- sum(apply(M, 2, sd) > 0)/ncol(M)
  fixed <- 1.0 - polymorphic
  return(fixed)
}

#####################################################
# EXTRACT SIMULATION METADATA FROM FULL OUTPUT DATA #
#####################################################

extract.metadata <- function(seasons, store.all.pops = FALSE){
  
  message("Extracting metadata ...")
  
  # initialize result list
  num.seasons <- length(seasons)-1
  metadata <- lapply(1:(num.seasons+1), function(i) {list()})
  
  # store general variables
  base.pop <- seasons[[1]]$cross.inbreed$pop.out
  # 1) QTL effects
  metadata[[1]]$general$qtl.effects <- get.qtl.effects(base.pop)

  # go through all seasons
  prev.candidates <- NULL
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
      
      # 4) inbreeding rate (marker based, cfr. Jannink)
      if(!is.null(prev.candidates)){
        metadata[[s+1]]$candidates$inbreeding <- inbreeding.rate(candidates, prev.candidates)
      }
      
      # 5) QTL favourable allele frequencies
      metadata[[s+1]]$candidates$fav.QTL.allele.freqs <- get.favourable.qtl.allele.frequencies(candidates)
      # 6) QTL - marker LD
      metadata[[s+1]]$candidates$QTL.marker.LD <- QTL.marker.highest.LD(candidates)
      # 7) number of favourable QTL lost
      num.lost <- sum(get.favourable.qtl.allele.frequencies(candidates) == 0)
      metadata[[s+1]]$candidates$num.fav.QTL.lost <- num.lost
      
      # 8) store full marker data (only if requested to store all populations)
      if(store.all.pops){
        metadata[[s+1]]$candidates$markers <- gp.design.matrix(candidates)
        metadata[[s+1]]$candidates$qtl <- get.qtl.allele.matrix(candidates)
      }
      
      # store for inbreeding rate calculation in next generation
      prev.candidates <- candidates
      
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
      
      # 4) store full marker data final selection (or of all selections if requested)
      if(store.all.pops || s == num.seasons){
        metadata[[s+1]]$selection$markers <- gp.design.matrix(selection)
        metadata[[s+1]]$selection$qtl <- get.qtl.allele.matrix(selection)
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
      # (!! based on marker *names*, *not* indices, as the latter
      #  include QTL and dummies which were dropped for GP analysis)
      all.names <- rownames(gp.tp$map)
      all.marker.effects <- gp.get.effects(gp.model)
      all.qtl.effects <- get.qtl.effects(gp.tp)
      marker.names <- all.names[QTL.marker.LD$marker.index] # convert indices to names
      marker.effects <- all.marker.effects[marker.names]
      qtl.names <- all.names[QTL.marker.LD$QTL.index] # convert indices to names
      qtl.effects <- all.qtl.effects[qtl.names]
      # extend QTL-marker LD table with effects and names
      QTL.marker.LD$QTL.name <- qtl.names
      QTL.marker.LD$marker.name <- marker.names
      QTL.marker.LD$QTL.effect <- qtl.effects
      QTL.marker.LD$marker.effect <- marker.effects
      # store QTL-marker LD table
      metadata[[s+1]]$gp$qtl.marker.ld <- QTL.marker.LD
      # compute and store accuracy (both plain and corrected)
      plain.acc <- cor(marker.effects, qtl.effects)
      metadata[[s+1]]$gp$effect.estimation.accuracy$plain <- plain.acc
      metadata[[s+1]]$gp$effect.estimation.accuracy$corrected <- plain.acc / mean.LD
      
      # 4) proportion of effect sign mismatches
      sign.diff <- sign(marker.effects) - sign(qtl.effects)
      sign.mismatches <- sum(sign.diff != 0) / length(marker.effects)
      metadata[[s+1]]$gp$sign.mismatches <- sign.mismatches
      
      # 5) TP size (after possible filtering)
      metadata[[s+1]]$gp$tp.size <- gp.tp$numGenotypes
      
    }
    
  }
  
  return(metadata)
  
}

###########################
# OUTPUT DATA CONVERSIONS #
###########################

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

# strip full marker data from output file to reduce memory footprint (except for final selection)
# NOTE: information needed for detailed snapshot MDS plots and movies is lost
strip.marker.data <- function(file){
 
  message("Processing file ", file)
  
  # load data
  seasons <- readRDS(file)
  
  # discard full marker data
  for(s in 1:length(seasons)){
    seasons[[s]]$candidates$markers <- NULL
    seasons[[s]]$candidates$qtl <- NULL
    # retain full data for final selection
    if(s < length(seasons)){
      seasons[[s]]$selection$markers <- NULL
      seasons[[s]]$selection$qtl <- NULL
    }
  }
  
  # overwrite data
  saveRDS(seasons, file = file)
   
}














