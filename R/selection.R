
########################
# SELECTION STRATEGIES #
########################

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

# select by maximizing weighted index of mean breeding value and diversity
select.weighted.index <- function(n, values, markers, div.weight, div.measure = c("MR", "HE", "HEadj"), fav.alleles = NULL, ...){
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
  # get selected diversity measure
  div.measure <- match.arg(div.measure)
  if(div.measure == "MR"){
    div.measure <- j.MR.ENE()
  } else if (div.measure == "HE"){
    div.measure <- j.HE()
  } else {
    div.measure <- j.adj.HE()
  }
  # run optimization
  selected.names <- j.max.index(n, names(values), values, markers, div.weight, div.measure, fav.alleles)
  return(selected.names)
}

####################
# DIVERSITY SCORES #
####################

# for DH population (0/1)
HE <- function(Z, sel){
  
  # compute average genome of selection
  freqs <- colMeans(Z[sel, ])

  # compute HE from p
  he <- sum(sapply(freqs, function(p){
    p * (1-p)
  }))
  he <- 2*he/ncol(Z)
  
  return(he)
    
}

# for DH individuals (0/1)
MR <- function(ind1, ind2){
  if(length(ind1) != length(ind2)){
    stop("marker vector of both individuals should be of same size")
  }
  d <- sqrt(sum(abs(ind1 - ind2))/length(ind1))
  return(d)
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

plot.pareto.fronts <- function(dir, div.measure = c("MR", "HE", "HEadj", "LOG")){
  
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




