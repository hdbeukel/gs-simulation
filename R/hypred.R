
###############################
# ADD AND STRIP DUMMY MARKERS #
###############################

add.dummies = function(data){
  # check if dummies have not already been added
  if(data$hypred$dummiesAdded){
    stop("dummmies already added")
  }
  # initialize extended genotype matrices and map
  num.genotypes = data$numGenotypes
  num.markers = data$numChroms * data$hypred$markersPerChrom
  if(!is.null(data$geno)){
    geno = matrix(NA, num.genotypes, num.markers)
    rownames(geno) = rownames(data$geno)
    colnames(geno) = seq(1, num.markers)
  }
  if(!is.null(data$phasedGeno)){
    phased.hap1 = matrix(NA, num.genotypes, num.markers)
    phased.hap2 = matrix(NA, num.genotypes, num.markers)
    rownames(phased.hap1) = rownames(data$phasedGeno$hap1)
    colnames(phased.hap1) = seq(1, num.markers)
    rownames(phased.hap2) = rownames(data$phasedGeno$hap2)
    colnames(phased.hap2) = seq(1, num.markers)
  }
  if(!is.null(data$dh)){
    dh = matrix(NA, num.genotypes, num.markers)
    rownames(dh) = rownames(data$dh)
    colnames(dh) = seq(1, num.markers)
  }
  map = matrix(NA, num.markers, 2)
  rownames(map) = seq(1, num.markers)
  colnames(map) = colnames(data$map)
  # fill extended genotype matrices and map (per chromosome)
  for(i in 1:data$numChroms){
    # compute marker index bounds (real part + dummy part)
    start.real = (i-1)*data$hypred$markersPerChrom + 1
    end.real = start.real + data$chrNumMarkers[i] - 1
    start.dummy = end.real + 1
    end.dummy = i*data$hypred$markersPerChrom
    orig.start = data$hypred$chromBounds[i,1]
    orig.end = data$hypred$chromBounds[i,2]
    message("Extending chromosome ", i, ":", " start.real = ", start.real,
                                             ", end.real = ", end.real,
                                             ", start.dummy = ", start.dummy,
                                             ", end.dummy = ", end.dummy,
                                             ", orig.start = ", orig.start,
                                             ", orig.end = ", orig.end)
    # copy real genotype data
    if(exists("geno")){
      geno[,start.real:end.real] = data$geno[, orig.start:orig.end]
      colnames(geno)[start.real:end.real] = colnames(data$geno)[orig.start:orig.end]
    }
    if(exists("phased.hap1")){
      phased.hap1[,start.real:end.real] = data$phasedGeno$hap1[, orig.start:orig.end]
      phased.hap2[,start.real:end.real] = data$phasedGeno$hap2[, orig.start:orig.end]
      colnames(phased.hap1)[start.real:end.real] = colnames(data$phasedGeno$hap1)[orig.start:orig.end]
      colnames(phased.hap2)[start.real:end.real] = colnames(data$phasedGeno$hap2)[orig.start:orig.end]
    }
    if(exists("dh")){
      dh[,start.real:end.real] = data$dh[, orig.start:orig.end]
      colnames(dh)[start.real:end.real] = colnames(data$dh)[orig.start:orig.end]
    }
    # set dummy genotype data
    if(data$hypred$chrNumDummies[i] > 0){
      dummy.names = paste("dummy", i, seq(1, data$hypred$chrNumDummies[i]), sep="-")
      if(exists("geno")){
        geno[,start.dummy:end.dummy] = 0
        colnames(geno)[start.dummy:end.dummy] = dummy.names
      }
      if(exists("phased.hap1")){
        phased.hap1[,start.dummy:end.dummy] = 0
        phased.hap2[,start.dummy:end.dummy] = 0
        colnames(phased.hap1)[start.dummy:end.dummy] = dummy.names
        colnames(phased.hap2)[start.dummy:end.dummy] = dummy.names
      }
      if(exists("dh")){
        dh[,start.dummy:end.dummy] = 0
        colnames(dh)[start.dummy:end.dummy] = dummy.names
      }
    }
    # update map (dummies at 1 cM distance)
    map[start.real:end.real,] = as.matrix(data$map[orig.start:orig.end,])
    rownames(map)[start.real:end.real] = rownames(data$map)[orig.start:orig.end]
    if(data$hypred$chrNumDummies[i] > 0){
      # chromosome number
      map[start.dummy:end.dummy,1] = i
      # position
      map[start.dummy:end.dummy,2] = data$chrLengths[i] + seq(1,data$hypred$chrNumDummies[i])
      # names
      rownames(map)[start.dummy:end.dummy] = dummy.names
    }
  }
  # set extended genotype matrices and map
  if(exists("geno")){
    data$geno = geno
  }
  if(exists("phased.hap1")){
    data$phasedGeno = list(hap1 = phased.hap1, hap2 = phased.hap2)
  }
  if(exists("dh")){
    data$dh = dh
  }
  data$map = as.data.frame(map)
  
  # update metadata
  data$numMarkers = num.markers
  data$chrLengths = data$chrLengths + data$hypred$chrNumDummies
  data$chrNumMarkers[] = data$hypred$markersPerChrom
  data$hypred$chromBounds[,1] = data$hypred$chromBounds[,1] + c(0, cumsum(data$hypred$chrNumDummies[1:(data$numChroms-1)]))
  data$hypred$chromBounds[,2] = data$hypred$chromBounds[,2] + cumsum(data$hypred$chrNumDummies)
  
  # attach hypred genome definition if not yet present
  if(is.null(data$hypred$genome)){
    # IMPORTANT: positions and lengths in Morgan (not cM)
    hypred.genome = hypredGenome(num.chr = data$numChroms,
                                 len.chr = data$chrLengths/100.0,
                                 num.snp.chr = data$hypred$markersPerChrom)
    hypred.genome = hypredNewMap(hypred.genome, data$map$pos/100.0)
    data$hypred$genome = hypred.genome
  }
  
  # change dummy flag
  data$hypred$dummiesAdded = TRUE
  
  return(data)
}

strip.dummies = function(data){
  # check if dummies have been added
  if(!data$hypred$dummiesAdded){
    stop("dummmies currently not added")
  }
  # initialize stripped genotype matrices and map
  num.genotypes = data$numGenotypes
  num.markers = data$numMarkers - sum(data$hypred$chrNumDummies)
  if(!is.null(data$geno)){
    geno = matrix(NA, num.genotypes, num.markers)
    rownames(geno) = rownames(data$geno)
    colnames(geno) = seq(1, num.markers)
  }
  if(!is.null(data$phasedGeno)){
    phased.hap1 = matrix(NA, num.genotypes, num.markers)
    phased.hap2 = matrix(NA, num.genotypes, num.markers)
    rownames(phased.hap1) = rownames(data$phasedGeno$hap1)
    colnames(phased.hap1) = seq(1, num.markers)
    rownames(phased.hap2) = rownames(data$phasedGeno$hap2)
    colnames(phased.hap2) = seq(1, num.markers)
  }
  if(!is.null(data$dh)){
    dh = matrix(NA, num.genotypes, num.markers)
    rownames(dh) = rownames(data$dh)
    colnames(dh) = seq(1, num.markers)
  }
  map = matrix(NA, num.markers, 3)
  rownames(map) = seq(1, num.markers)
  colnames(map) = colnames(data$map)
  # strip genotype matrices and map (per chromosome)
  strippedNumMarkers = data$chrNumMarkers - data$hypred$chrNumDummies
  for(i in 1:data$numChroms){
    # compute marker index bounds
    start.stripped = sum(strippedNumMarkers[1:i]) - strippedNumMarkers[i] + 1
    end.stripped = sum(strippedNumMarkers[1:i])
    start.extended = data$hypred$chromBounds[i,1]
    end.extended = data$hypred$chromBounds[i,2] - data$hypred$chrNumDummies[i]
    message("Stripping chromosome ", i, ":", " start.stripped = ", start.stripped,
                                             ", end.stripped = ", end.stripped,
                                             ", start.extended = ", start.extended,
                                             ", end.extended = ", end.extended)
    # strip genotype data
    if(exists("geno")){
      geno[,start.stripped:end.stripped] = data$geno[,start.extended:end.extended]
      colnames(geno)[start.stripped:end.stripped] = colnames(data$geno)[start.extended:end.extended]
    }
    if(exists("phased.hap1")){
      phased.hap1[,start.stripped:end.stripped] = data$phasedGeno$hap1[,start.extended:end.extended]
      phased.hap2[,start.stripped:end.stripped] = data$phasedGeno$hap2[,start.extended:end.extended]
      colnames(phased.hap1)[start.stripped:end.stripped] = colnames(data$phasedGeno$hap1)[start.extended:end.extended]
      colnames(phased.hap2)[start.stripped:end.stripped] = colnames(data$phasedGeno$hap2)[start.extended:end.extended]
    }
    if(exists("dh")){
      dh[,start.stripped:end.stripped] = data$dh[,start.extended:end.extended]
      colnames(dh)[start.stripped:end.stripped] = colnames(data$dh)[start.extended:end.extended]
    }
    # strip map
    map[start.stripped:end.stripped,] = as.matrix(data$map[start.extended:end.extended,])
    rownames(map)[start.stripped:end.stripped] = rownames(data$map)[start.extended:end.extended]
  }
  # set stripped genotype matrices and map
  if(exists("geno")){
    data$geno = geno
  }
  if(exists("phased.hap1")){
    data$phasedGeno = list(hap1 = phased.hap1, hap2 = phased.hap2)
  }
  if(exists("dh")){
    data$dh = dh
  }
  data$map = as.data.frame(map)
  
  # update metadata
  data$numMarkers = num.markers
  data$chrLengths = data$chrLengths - data$hypred$chrNumDummies
  data$chrNumMarkers = data$chrNumMarkers - data$hypred$chrNumDummies
  data$hypred$chromBounds[,1] = data$hypred$chromBounds[,1] - c(0, cumsum(data$hypred$chrNumDummies[1:(data$numChroms-1)]))
  data$hypred$chromBounds[,2] = data$hypred$chromBounds[,2] - cumsum(data$hypred$chrNumDummies)
  
  # change dummy flag
  data$hypred$dummiesAdded = FALSE
  
  return(data)
}

get.dummy.indices <- function(pop){
  # assert: dummies added
  if(!pop$hypred$dummiesAdded){
    stop("no dummies added in pop")
  }
  dummy.ends <- pop$hypred$chromBounds[,2]
  dummy.starts <- pop$hypred$chromBounds[,2] - pop$hypred$chrNumDummies + 1
  dummy.bounds <- cbind(dummy.starts, dummy.ends)
  dummy.indices <- unlist(apply(dummy.bounds, 1, function(row){
    if(row[1] <= row[2]){
      seq(row[1], row[2])
    }
  }))
  return(dummy.indices)
}

#######################################
# ASSIGN QTL AND INFER GENETIC VALUES #
#######################################

# assign QTL:
# - each QTL is randomly assigned to a chromosome with probability proportional
#   to the chromosome length
# - only polymorphic markers are chosen as QTL
# - dummy QTLs with effect of 0 are added because hypred requires same number of
#   QTL per chromosome
# - real QTLs effects are determined using one of two provided methods:
#   1) "normal":  effects are sampled from a normal distribution, avoiding the tails
#                 (i.e. large effects) by sampling only from a certain inner region
#                 of the distribution (based on parameter 'random.effect.coverage')
#   2) "jannink": effects are scaled to the inverse of the standard deviation of the
#                 allelic states 0/1/2, with random sign
assign.qtl <- function(pop, num.qtl,
                       method = c("normal", "jannink"),
                       random.effect.coverage = 0.9){
  # process arguments
  method <- match.arg(method)
  # assert: dummies added
  if(!pop$hypred$dummiesAdded){
    stop("first add dummy markers using add.dummies(pop)")
  }
  # assert: no QTL added yet
  if(!is.null(pop$hypred$realQTL)){
    stop("QTL already assigned")
  }
  # check input
  if(missing(num.qtl) || is.null(num.qtl)){
    stop("specify total number of QTL")
  }
  # compute num QTL per chrom
  num.chroms <- pop$numChroms
  chrom.lengths <- pop$chrLengths
  total.length <- sum(chrom.lengths)
  chrom.lengths.rel <- chrom.lengths/total.length
  qtl.assignment <- sample(1:num.chroms, replace=TRUE, size=num.qtl, prob=chrom.lengths.rel)
  real.qtl.per.chrom.table <- table(qtl.assignment)
  real.qtl.per.chrom <- rep(0, num.chroms)
  real.qtl.per.chrom[as.numeric(names(real.qtl.per.chrom.table))] <- real.qtl.per.chrom.table
  hypred.qtl.per.chrom <- max(real.qtl.per.chrom)
  dummy.qtl.per.chrom <- hypred.qtl.per.chrom - real.qtl.per.chrom
  # get ranges of non dummy markers
  non.dummy.starts <- pop$hypred$chromBounds[,1]
  non.dummy.stops <- pop$hypred$chromBounds[,2] - pop$hypred$chrNumDummies
  non.dummy.ranges <- cbind(non.dummy.starts, non.dummy.stops)
  # get genotypes
  if(!is.null(pop$dh)){
    genotypes <- pop$dh * 2
  } else {
    genotypes <- pop$geno
  }
  # infer positions of possible fixed markers
  mafs <- maf(genotypes, encoding = "012")
  fixed <- which(mafs == 0)
  # pick real QTL indices per chromosome (NOT allowed to be dummy markers)
  real.qtl.indices <- sort(unlist(sapply(1:num.chroms, function(c){
    # QTL candidates on current chromosome
    candidates <- setdiff(seq(non.dummy.ranges[c,1], non.dummy.ranges[c,2]), fixed)
    # sample real QTL indices
    sample(candidates, real.qtl.per.chrom[c])
  })))
  # pick dummy QTL indices per chromosome (can be dummy markers, don't care)
  full.chrom.ranges <- pop$hypred$chromBounds
  dummy.qtl.indices <- sort(unlist(sapply(1:num.chroms, function(c){
    # get candidiates
    candidates <- seq(full.chrom.ranges[c,1], full.chrom.ranges[c,2])
    # remove indices which are already picked as real QTL
    candidates <- setdiff(candidates, real.qtl.indices)
    # sample dummy QTL indices
    sample(candidates, dummy.qtl.per.chrom[c])
  })))
  # store real/dummy QTL indices in population object
  pop$hypred$realQTL <- real.qtl.indices
  pop$hypred$dummyQTL <- dummy.qtl.indices
  
  # generate real QTL effects
  if(method == "normal"){
    # random effects sampled from normal distribution (without tail)
    real.qtl.effects <- rep(NA, num.qtl)
    i <- 1
    batch.size <- 10
    q <- abs(qnorm((1-random.effect.coverage)/2))
    while(i <= num.qtl){
      rand.effects <- rnorm(batch.size)
      inner.rand.effects <- rand.effects[abs(rand.effects) <= q]
      num.new.effects <- length(inner.rand.effects)
      num.used.effects <- min(num.new.effects, num.qtl - i + 1)
      real.qtl.effects[i:(i+num.used.effects-1)] <- inner.rand.effects[1:num.used.effects]
      i <- i + num.used.effects
    }
  } else if (method == "jannink"){
    # effects scaled to inverse of standard deviation of QTL allelic states 0/1/2, random sign
    real.qtl.effect.scales = 1 / apply(genotypes[,real.qtl.indices], 2, sd)
    real.qtl.effects = sample(c(-1,1), length(real.qtl.effect.scales), replace = TRUE) * real.qtl.effect.scales
  } else {
    stop(sprintf("Unknown effect sampling method %s (should not happen)", method))
  }
  
  # set dummy QTL effects to ZERO
  dummy.qtl.effects <- rep(0, length(dummy.qtl.indices))

  # combine and sort real and dummy QTL IDs (sort effects accordingly)
  real.and.dummy.qtl.indices <- c(real.qtl.indices, dummy.qtl.indices)
  real.and.dummy.qtl.effects <- c(real.qtl.effects, dummy.qtl.effects)
  sorted.index.order <- order(real.and.dummy.qtl.indices)
  real.and.dummy.qtl.indices <- real.and.dummy.qtl.indices[sorted.index.order]
  real.and.dummy.qtl.effects <- real.and.dummy.qtl.effects[sorted.index.order]
  
  # update hypred genome
  genome.with.qtl <- hypredNewQTL(pop$hypred$genome,
                                 new.id.add = real.and.dummy.qtl.indices,
                                 new.eff.add = real.and.dummy.qtl.effects)
  pop$hypred$genome <- genome.with.qtl
  
  return(pop)
}

# get QTL effects (*real* QTL only)
get.qtl.effects <- function(pop){
  # assert: QTLs assigned
  if(attr(pop$hypred$genome, 'num.add.qtl.chr') == 0){
    stop("first assign QTL using assign.qtl(pop, ...)")
  }
  # retrieve positions of all QTL (including dummies)
  all.QTL.pos <- pop$hypred$genome@pos.add.qtl$ID
  # retrieve positions of real QTL
  real.QTL.pos <- pop$hypred$realQTL
  # retrieve all effects (including dummies)
  all.QTL.effects <- pop$hypred$genome@add.and.dom.eff$add
  # select only real effects
  real.QTL.effects <- all.QTL.effects[all.QTL.pos %in% real.QTL.pos]
  # retrieve and set names of real QTL in effect vector
  real.QTL.names <- rownames(pop$map)[real.QTL.pos]
  names(real.QTL.effects) <- real.QTL.names
  # return named real QTL effects
  return(real.QTL.effects)
}

get.qtl.allele.matrix <- function(pop){
  # assert: DH population
  if(is.null(pop$dh)){
    stop("pop should be a DH population")
  }
  Q <- pop$dh[, pop$hypred$realQTL]
  return(Q)
}

# get favourable QTL alleles (*real* QTL only)
get.favourable.qtl.alleles <- function(pop){
  # get effects
  qtl.eff <- get.qtl.effects(pop)
  # initialize output
  desired.qtl.alleles <- rep(0, length(qtl.eff))
  # set to 1 for positive effects
  desired.qtl.alleles[qtl.eff > 0] <- 1
  # set names
  names(desired.qtl.alleles) <- names(qtl.eff)
  
  return(desired.qtl.alleles)
}

# get favourable QTL allele frequencies (*real* QTL only)
get.favourable.qtl.allele.frequencies <- function(pop){
  # assert: DH population
  if(is.null(pop$dh)){
    stop("pop should be a DH population")
  }
  # get favourable alleles
  fav.alleles <- get.favourable.qtl.alleles(pop)
  # get QTL allele matrix
  Q <- get.qtl.allele.matrix(pop)
  # compute frequencies
  n <- nrow(Q)
  m <- ncol(Q)
  freqs <- rep(NA, m)
  for(i in 1:m){
    freqs[i] <- sum(Q[,i] == fav.alleles[i]) / n
  }
  # set names
  names(freqs) <- names(fav.alleles)
  
  return(freqs)
}

infer.genetic.values = function(pop){
  # assert: QTLs assigned
  if(attr(pop$hypred$genome, 'num.add.qtl.chr') == 0){
    stop("first assign QTL using assign.qtl(pop, ...)")
  }
  # infer genetic values
  if(!is.null(pop$dh)){
    genetic.values = as.vector(hypredTruePerformance(pop$hypred$genome, pop$dh, DH = TRUE))
    names(genetic.values) = rownames(pop$dh)
  } else {
    genotypes = interleave(pop$phasedGeno$hap1, pop$phasedGeno$hap2)
    genetic.values = as.vector(hypredTruePerformance(pop$hypred$genome, genotypes, DH = FALSE))
    names(genetic.values) = rownames(pop$geno)
  }
  pop$geneticValues = genetic.values
  # compute genetic variation
  genetic.var = var(genetic.values)
  pop$geneticVar = genetic.var
  
  return(pop)
}

# normalize genetic values to [min, max] based on the minimum
# and maximum possible genetic value (inferred from the QTL effects)
normalize.genetic.values <- function(values, effects, min = -1, max = 1){
  # compute minimum and maximum possible genetic value
  min.gen.value <- 2 * sum(effects[effects < 0])
  max.gen.value <- 2 * sum(effects[effects > 0])
  # normalize genetic values to [min, max]
  values <- min + (max-min) * (values - min.gen.value)/(max.gen.value - min.gen.value)
  return(values)
}

###########################
# CROSSING AND INBREEDING #
###########################

# create initial population by random mating of founders (non DH) + DH creation
mate.founders <- function(founders, n, name.prefix = "i"){
  # assert: dummies added
  if(!founders$hypred$dummiesAdded){
    stop("first add dummy markers using add.dummies(founders)")
  }
  # assert: no DH population!
  if(is.null(founders$geno)){
    stop("only use this function for NON-DH populations (initial founders)")
  }
  # initialize output
  offspring <- founders
  offspring$numGenotypes <- n                          # number of DHs
  offspring$dh <- matrix(NA, n, ncol(founders$geno))   # DH matrix
  colnames(offspring$dh) <- colnames(founders$geno)
  rownames(offspring$dh) <- paste(name.prefix, seq(1:n), sep="-")
  offspring$geno <- NULL                               # discard founders genotype matrix
  offspring$phasedGeno <- NULL                         # discard founders phased genotypes
  offspring$pedigree$par1 <- rep(NA, n)
  offspring$pedigree$par2 <- rep(NA, n)
  names(offspring$pedigree$par1) <- names(offspring$pedigree$par2) <- rownames(offspring$dh)
  # generate DHs
  for(i in 1:n){
    # select two distinct parents
    parents <- sample(rownames(founders$geno), 2)
    # create gametes from each parent (recombination)
    gametes <- lapply(parents, function(p) {
      recombine(founders$hypred$genome, founders$phasedGeno$hap1[p,], founders$phasedGeno$hap2[p,])
    })
    # create DH
    dh <- recombine(founders$hypred$genome, gametes[[1]], gametes[[2]])
    # store
    offspring$dh[i,] <- dh
    offspring$pedigree$par1[i] <- parents[1]
    offspring$pedigree$par2[i] <- parents[2]
  }
  
  # update genetic values if QTL assigned
  if(!is.null(offspring$hypred$realQTL)){
    offspring <- infer.genetic.values(offspring)
  }
  
  # update heritability (if set)
  if(!is.null(offspring$heritability)){
    # update heritability
    offspring$heritability <- offspring$geneticVar / (offspring$geneticVar + offspring$errorVar)
  }
  
  # mark as DH
  offspring$isDH <- TRUE
  
  # erase phenotypes (if any)
  offspring$pheno <- NULL
  
  return(offspring)
}

# create new population through random mating of individuals in a given DH population
# random crossings create n F1 individuals + inbreeding --> one DH per F1 genotype
mate.dh <- function(pop, n, name.prefix){
  # assert: DH population
  if(is.null(pop$dh)){
    stop("only use this function for DH populations")
  }
  # assert: dummies added
  if(!pop$hypred$dummiesAdded){
    stop("first add dummy markers using add.dummies(pop)")
  }
  # initialize output
  offspring <- pop
  offspring$numGenotypes <- n                          # number of DHs
  offspring$dh <- matrix(NA, n, pop$numMarkers)        # DH matrix
  colnames(offspring$dh) <- colnames(pop$dh)
  rownames(offspring$dh) <- paste(name.prefix, seq(1:n), sep="-")
  offspring$pedigree$par1 <- rep(NA, n)
  offspring$pedigree$par2 <- rep(NA, n)
  names(offspring$pedigree$par1) <- names(offspring$pedigree$par2) <- rownames(offspring$dh)
  # generate DHs
  for(i in 1:n){
    # select two distinct DH parents
    parents <- sample(rownames(pop$dh), 2, replace = T)
    gametes <- pop$dh[parents,]
    # recombine to create new DH
    dh <- recombine(pop$hypred$genome, gametes[1,], gametes[2,])
    # store
    offspring$dh[i,] <- dh
    offspring$pedigree$par1[i] <- parents[1]
    offspring$pedigree$par2[i] <- parents[2]
  }
  
  # update genetic values if QTL assigned
  if(!is.null(offspring$hypred$realQTL)){
    offspring <- infer.genetic.values(offspring)
  }
  
  # update heritability (if set)
  if(!is.null(offspring$heritability)){
    # update heritability
    offspring$heritability <- offspring$geneticVar / (offspring$geneticVar + offspring$errorVar)
  }
  
  # erase phenotypes (if any)
  offspring$pheno <- NULL
  
  return(offspring)
}

#################################################
# DETERMINE ERROR VARIANCE AND INFER PHENOTYPES #
#################################################

# set a fixed heritability (infer corresponding error variance)
set.heritability = function(pop, heritability){
  # assert: genetic values available
  if(is.null(pop$geneticValues)){
    stop("first infer genetic values using infer.genetic.values(pop)")
  }
  # determine error variation to obtain the given heritability
  error.var = (1-heritability)/heritability * pop$geneticVar
  # attach to population object
  pop$errorVar = error.var
  pop$heritability = heritability
  
  return(pop)
}

# set a fixed error variance (infer corresponding heritability)
set.error.variance = function(pop, error.var){
  # assert: genetic values available
  if(is.null(pop$geneticValues)){
    stop("first infer genetic values using infer.genetic.values(pop)")
  }
  # infer heritability
  h = pop$geneticVar / (pop$geneticVar + error.var)
  # attach to population object
  pop$errorVar = error.var
  pop$heritability = h
  
  return(pop)
}

infer.phenotypes = function(pop){
  # assert: error variance determined
  if(is.null(pop$errorVar)){
    stop("first determine error variance using set.heritability(pop, ...) or set.error.variance(pop, ...)")
  }
  # sample normally distributed errors with predetermined variance
  errors = rnorm(pop$numGenotypes, sd = sqrt(pop$errorVar))
  # sum genetic values and error to obtain phenotypes
  pheno = pop$geneticValues + errors
  # attach to population
  pop$pheno = pheno
  
  return(pop)
}

####################################################################
# RESTRICT POPULATIONS (SELECTION) AND COMBINE POPULATIONS (MERGE) #
####################################################################

# restrict population by selecting individuals with given names
restrict.population = function(pop, names){
  # check input
  if((!is.null(pop$dh) && !all(names %in% rownames(pop$dh)))
     || (!is.null(pop$geno) && !all(names %in% rownames(pop$geno)))){
    stop ("some individuals to select are not contained in the population")
  }
  # restrict DH matrix (if any)
  if(!is.null(pop$dh)){
    pop$dh = pop$dh[names,]
  }
  # restrict (phased) genotype matrices (if any)
  if(!is.null(pop$geno)){
    pop$geno = pop$geno[names,]
    pop$phasedGeno$hap1 = pop$phasedGeno$hap1[names,]
    pop$phasedGeno$hap2 = pop$phasedGeno$hap2[names,]
  }
  # update number of genotypes
  pop$numGenotypes = length(names)
  # update genetic values and variance (if any)
  if(!is.null(pop$geneticValues)){
    pop$geneticValues = pop$geneticValues[names]
    pop$geneticVar = var(pop$geneticValues)
  }
  # update heritability (if set)
  if(!is.null(pop$heritability)){
    pop$heritability = pop$geneticVar / (pop$geneticVar + pop$errorVar)
  }
  # update phenotypes (if any)
  if(!is.null(pop$pheno)){
    pop$pheno = pop$pheno[names]
  }
  # update estimated genetic values (if any)
  if(!is.null(pop$estGeneticValues)){
    pop$estGeneticValues = pop$estGeneticValues[names]
  }
  # restrict pedigree (if any)
  if(!is.null(pop$pedigree)){
    pop$pedigree$par1 = pop$pedigree$par1[names]
    pop$pedigree$par2 = pop$pedigree$par2[names]
  }
  # return restricted population
  return(pop)
}

# merge populations into one big population (only for DH populations)
# error variation is taken from first population (if specified)
# WARNING: populations should have same underlying hypred genome (QTL etc.),
#          else the produced result is meaningless
merge.populations = function(pop1, pop2){
  # check input
  if(is.null(pop1$dh) || is.null(pop2$dh)){
    stop("merging is only supported for DH populations")
  }
  if(length(intersect(rownames(pop1$dh), rownames(pop2$dh))) > 0){
    stop("overlapping genotype names in given populations")
  }
  # initialize new population as pop1, then merge pop2 into new pop
  pop = pop1
  # combine DH matrices
  pop$dh = rbind(pop1$dh, pop2$dh)
  # update number of genotypes
  pop$numGenotypes = pop1$numGenotypes + pop2$numGenotypes
  # merge genetic values and update variance (if known in both populations)
  if(!is.null(pop1$geneticValues) && !is.null(pop2$geneticValues)){
    pop$geneticValues = c(pop1$geneticValues, pop2$geneticValues)
    pop$geneticVar = var(pop$geneticValues)
  }
  # update heritability (if set in both populations)
  if(!is.null(pop1$errorVar) && !is.null(pop$geneticVar)){
    pop$errorVar = pop1$errorVar
    pop$heritability = pop$geneticVar / (pop$geneticVar + pop$errorVar)
  }
  # merge phenotypes (if known in both populations)
  if(!is.null(pop1$pheno) && !is.null(pop2$pheno)){
    pop$pheno = c(pop1$pheno, pop2$pheno)
  }
  # merge estimated genetic values (if available for both populations)
  if(!is.null(pop1$estGeneticValues) && !is.null(pop2$estGeneticValues)){
    pop$estGeneticValues = c(pop1$estGeneticValues, pop2$estGeneticValues)
  }
  # merge pedigree (if available for both populations)
  if(!is.null(pop1$pedigree) && !is.null(pop2$pedigree)){
    pop$pedigree$par1 = c(pop1$pedigree$par1, pop2$pedigree$par1)
    pop$pedigree$par2 = c(pop1$pedigree$par2, pop2$pedigree$par2)
  }
  # return merged population
  return(pop)
}

###################
# LD COMPUTATIONS #
###################

# find QTL-marker pairs in highest LD (on the same chromosome of course)
# returns data frame with three columns: QTL.index, marker.index, LD
QTL.marker.highest.LD <- function(pop){

  # assert: DH population
  if(is.null(pop$dh)){
    stop("DH population expected")
  }
  
  # get QTL indices
  QTL.indices <- pop$hypred$realQTL
  num.QTL <- length(QTL.indices)
  # get dummy marker indices
  dummy.indices <- get.dummy.indices(pop)
  # infer non QTL non dummy indices
  SNP.indices <- setdiff(1:pop$numMarkers, c(QTL.indices, dummy.indices))
  
  # initialize result vectors
  highest.LD.marker.indices <- rep(NA, num.QTL)
  QTL.marker.LD <- rep(NA, num.QTL)
  
  # find marker in highest LD with each QTL
  for(i in 1:num.QTL){
    # get QTL index
    QTL.index <- QTL.indices[i]
    # get QTL alleles in population
    QTL.alleles <- pop$dh[, QTL.index]
    # check whether QTL is polymorphic
    if(sd(QTL.alleles) > 0){
      # get QTL chromosome
      QTL.chrom <- pop$map[QTL.index, 1]
      # get SNP indices on QTL chromosome
      QTL.chrom.SNP.indices <- intersect(SNP.indices, which(pop$map$chr == QTL.chrom))
      # compute LD with each SNP on same chromosome
      marker.LD <- sapply(QTL.chrom.SNP.indices, function(marker.index){
        marker.alleles <- pop$dh[, marker.index]
        return(ifelse(sd(marker.alleles) > 0, compute.LD(QTL.alleles, marker.alleles), NA))
      })
      # check that not all markers are fixed
      if(sum(!is.na(marker.LD)) > 0){
        # find highest LD
        highest.LD <- max(marker.LD, na.rm = TRUE)
        # get indices of markers in highest LD
        highest.LD.markers <- QTL.chrom.SNP.indices[which(marker.LD == highest.LD)]
        # select highest LD SNP closest to QTL
        if(length(highest.LD.markers) == 1){
          # unique max LD
          highest.LD.marker <- highest.LD.markers[1]
        } else {
          # multiple max LD: select SNP closest to QTL
          QTL.pos <- pop$map[QTL.index, 2]
          marker.dist <- sapply(highest.LD.markers, function(marker.index){
            marker.pos <- pop$map[marker.index, 2]
            dist <- abs(marker.pos - QTL.pos)
            return(dist)
          })
          highest.LD.marker <- highest.LD.markers[which.min(marker.dist)]
        }
        # store
        highest.LD.marker.indices[i] <- highest.LD.marker
        QTL.marker.LD[i] <- highest.LD
      }
    }
  }
  
  # combine results in data frame
  results <- data.frame(QTL.index = QTL.indices, marker.index = highest.LD.marker.indices, LD = QTL.marker.LD)
  
  return(results)
  
}

compute.LD <- function(alleles.1, alleles.2){
  cor(alleles.1, alleles.2)^2
}

#############
# UTILITIES #
#############

# wrapper around hypred recombination (no mutation, no blocks)
recombine = function(genome, gamete1, gamete2){
  return(hypredRecombine(genome,
                         gamete1,
                         gamete2,
                         mutate = F,
                         block = F))
}


























