
########################
# SELECTION STRATEGIES #
########################

# no-op selection that returns the entire set of candidates
# (ignores the selection size n)
select.all <- function(n, values, ...){
  return(names(values))
}

# random selection of requested size
select.random <- function(n, values, ...){
  # check input
  if(is.null(names(values))){
    stop("values should be named by individual")
  }
  sel <- sample(names(values), size = n)
  return(sel)
}

# selection criterion that selects the n plants with the highest value
# (phenotype, estimated genetic value, weighted genetic value, ...)
select.highest.score <- function(n, values, ...){
  # check input
  if(is.null(names(values))){
    stop("values should be named by individual")
  }
  selected.names <- names(head(sort(values, decreasing = TRUE), n=n))
  return(selected.names)
}

# selects the n individuals with the highest optimal contributions according to Meuwissen 1997
select.highest.optimal.contribution <- function(n, values, markers, generation, delta.F, ...){
  # compute optimal contributions
  C <- 1 - (1 - delta.F)^generation
  c <- optimal.contributions(values, markers, C)
  # select individuals with highest contribution
  selected.names <- names(head(sort(c, decreasing = TRUE), n=n))
  return(selected.names)
}

select.fixed.size.oc <- function(n, values, markers, generation, delta.F,
                                 adaptive = FALSE, sonesson = FALSE, verbose = FALSE, ...){
  # compute fixed size OC
  # C <- 1 - (1 - delta.F)^generation
  C <- delta.F
  if(adaptive){
    freqs <- colMeans(markers)
    he <- mean(2*freqs*(1-freqs))
    C <- delta.F * 2*he
  }
  c <- optimal.contributions(values, markers, C, size = n, sonesson = sonesson, verbose = verbose)
  # retrieve and check
  selected.names <- names(c[c > 0])
  if(length(selected.names) != n){
    stop("Incorrect number of individuals selected by fixed size OC")
  }
  return(selected.names)
}

# select by maximizing weighted index of mean breeding value and diversity
select.weighted.index <- function(n, values, markers, div.weight,
                                  div.measure = c("LOGall", "OC", "HEall", "HEfav", "LOGfav"),
                                  fav.alleles = NULL, ...){
  # check input
  if(!is.numeric(n)){
    stop("n should be an integer (selection size)")
  }
  if(is.null(names(values)) || !is.numeric(values)){
    stop("values should be a named numeric vector (one entry per individual)")
  }
  if(rownames(markers) != names(values) || !is.numeric(markers)){
    stop("markers matrix should be 0/1 with rownames corresponding to individuals")
  }
  if(!is.null(fav.alleles) && !is.numeric(fav.alleles)){
    stop("fav.alleles should be a 0/1 vector")
  }
  if(!is.numeric(div.weight) || div.weight < 0 || div.weight > 1){
    stop("div.weight should be a number in [0,1]")
  }
  # compute G matrix if OC measure selected
  G <- NULL
  if(div.measure == "OC"){
    G <- genomic.relationship.matrix(markers)
  }
  # get Java object for selected diversity measure
  div.measure <- match.arg(div.measure)
  if (div.measure == "HEall"){
    div.measure <- j.HE.all()
  } else if (div.measure == "HEfav"){
    div.measure <- j.HE.fav()
  } else if (div.measure == "LOGall"){
    div.measure <- j.LOG.all()
  } else if (div.measure == "LOGfav") {
    div.measure <- j.LOG.fav()
  } else if (div.measure == "OC") {
    div.measure <- j.OC()
  } else {
    stop("Unknown diversity measure:", div.measure)
  }
  # run optimization
  selected.names <- j.max.index(n, names(values), values, markers, div.weight, div.measure, fav.alleles, G)
  return(selected.names)
}

####################################
# WEIGHTED SELECTION PARETO FRONTS #
####################################

plot.pareto.front <- function(file, title = "Pareto front", xlab = "Diversity"){
  
  # read file
  front <- read.csv(file)
  
  # extract experiment IDs
  experiments <- unique(front$ID)
  n <- length(experiments)
  
  # initialize vectors to store mean/sd of diversity and quality scores
  div.mean <- rep(NA, n)
  div.sd   <- rep(NA, n)
  val.mean <- rep(NA, n)
  val.sd   <- rep(NA, n)
  
  # extract values
  for(i in 1:length(experiments)){
    
    experiment <- front[front$ID == experiments[i], ]
    
    div.mean[i] <- mean(experiment$div)
    div.sd[i]   <- sd(experiment$div)
    val.mean[i] <- mean(experiment$value)
    val.sd[i]   <- sd(experiment$value)
    
    # store values for equal weight case separately
    if(isTRUE(all.equal(unique(experiment$divWeight), 0.5))){
      mean.div.equal.weight   <- div.mean[i]
      mean.value.equal.weight <- val.mean[i]
    }
    
  }
  
  # plot front
  plot(x = div.mean, y = val.mean,
       main = title, xlab = xlab, ylab = "Mean genetic value")
  
  # mark point at equal weights
  if(exists("mean.div.equal.weight")){
    points(x = mean.div.equal.weight, y = mean.value.equal.weight, pch = 21, bg = "red")
  }
  
}

plot.pareto.fronts <- function(dir, div.measure = c("HEall", "HEfav", "LOGall", "LOGfav")){
  
  div.measure <- match.arg(div.measure)
  dir <- sprintf("%s/%s", dir, div.measure)
  
  for(pop in 1:3){
    for(snapshot in c("early", "medium", "late")){
      file <- sprintf("snapshot-pop-%d-%s.csv", pop, snapshot)
      full.path <- sprintf("%s/%s", dir, file)
      xlab <- sprintf("Diversity (%s)", div.measure)
      title <- sprintf("Population %d (%s)", pop, snapshot)
      plot.pareto.front(full.path, xlab = xlab, title = title)
    }
  }
  
}




