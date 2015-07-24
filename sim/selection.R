
########################
# SELECTION STRATEGIES #
########################

# selection criterion that selects the n plants with the highest score
# (phenotype, estimated genetic value, weighted genetic value, ...)
select.highest.score <- function(scores, n){
  # check input
  if(is.null(names(scores))){
    stop("scores should be named by individual")
  }
  selected.names <- names(head(sort(scores, decreasing = TRUE), n=n))
  return(selected.names)
}

####################################
# WEIGHTED SELECTION PARETO FRONTS #
####################################

plot.pareto.front <- function(file, title = "Pareto front"){
  
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
       main = title, xlab = "Diversity", ylab = "Median genetic value")
  
  # mark point at equal weights
  if(exists("mean.div.equal.weight")){
    points(x = mean.div.equal.weight, y = mean.value.equal.weight, pch = 21, bg = "red")
  }
  
}




