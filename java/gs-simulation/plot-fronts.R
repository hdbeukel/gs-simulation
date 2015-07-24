plot.front <- function(data){

  experiments <- unique(data$ID)
  mean.divs <- rep(NA, length(experiments))
  mean.values <- rep(NA, length(experiments))
  
  for(i in 1:length(experiments)){
    
    experiment <- data[data$ID == experiments[i], ]
    mean.divs[i] <- mean(experiment$div)
    mean.values[i] <- mean(experiment$value)
    
    if(isTRUE(all.equal(unique(experiment$divWeight), 0.5))){
      mean.div.equal.weight <- mean.divs[i]
      mean.value.equal.weight <- mean.values[i]
    }
    
  }
  
  plot(x = mean.divs, y = mean.values, xlab = "Diversity", ylab = "Median genetic value")
  
  if(exists("mean.div.equal.weight")){
    points(x = mean.div.equal.weight, y = mean.value.equal.weight, pch = 21, bg = "red")
  }

}