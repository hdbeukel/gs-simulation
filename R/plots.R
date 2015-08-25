################
# LOAD RESULTS #
################

# create list containing simulated breeding cycles, read from all
# files found in the given directory that match the given file pattern
load.simulation.results <- function(dir, file.pattern = "*.RDS") {
  # get file paths
  files <- Sys.glob(sprintf("%s/%s", dir, file.pattern))
  # load list of simulated breeding cycles
  breeding.cycles <- lapply(files, readRDS)
  # return results
  return(breeding.cycles)
}

#############################
# AUTOMATED PLOT GENERATION #
#############################

create.pdf <- function(file, width, height, plot.fun){
  
  pdf(file, width = width, height = height)
  plot.fun()
  invisible(dev.off())
  
}

get.plot.functions <- function(ylims){
  plot.functions <- list(
    list(f = plot.genetic.gain, name = "gain", title = "Genetic gain", legend = "bottomright"),
    list(f = plot.ratio.fixed.QTL, name = "QTL-fixed", title = "Ratio of fixed QTL", legend = "bottomright"),
    list(f = plot.mean.QTL.fav.allele.freq, name = "QTL-fav-allele-freq", title = "Mean QTL favourable allele frequency", legend = "bottomright"),
    list(f = plot.mean.QTL.marker.LD, name = "LD", title = "Mean polymorphic QTL - marker LD", legend = "bottomleft"),
    list(f = plot.mean.inbreeding, name = "inbreeding", title = "Mean inbreeding in selection candidates", legend = "bottomright"),
    list(f = plot.genetic.standard.deviation, name = "genetic-sd", title = "Genetic standard deviation", legend = "topright"),
    list(f = plot.num.fav.QTL.lost, name = "fav-QTL-lost", title = "Number of favourable QTL lost", legend = "bottomright"),
    list(f = plot.effect.estimation.accuracy, name = "eff-acc", title = "Effect estimation accuracy", legend = "bottomright"),
    list(f = function(...){ 
      plot.effect.estimation.accuracy(..., corrected = TRUE) 
    }, name = "eff-acc-corrected", title = "Corrected effect estimation accuracy", legend = "bottomright")
  )
  # set ylims
  plot.functions <- lapply(plot.functions, function(pf){
    c(pf, list(ylim = ylims[[pf$name]]))
  })
  return(plot.functions)
}

# stores PDF plots in "figures/simulation/CGS";
# organized into subfolders according to the applied
# diversity measure (HE, MR) and range of diversity weights
plot.CGS <- function(div.weights = c(0.25, 0.50, 0.75), div.measures = c("HE", "MR"), file.pattern = "*.RDS", xlim = c(0,30)){
  
  min.w <- min(div.weights)
  max.w <- max(div.weights)
  fig.dir <- "figures/simulation/CGS"
  
  # create plots for each diversity measure
  for(div in div.measures){
    
    message("Diversity measure: ", div)
    
    fig.subdir <- sprintf("%s/%s-%.2f-%.2f", fig.dir, div, min.w, max.w)
    
    if(!dir.exists(fig.subdir)){
      message(sprintf("|- Create output directory \"%s\"", fig.subdir))
      dir.create(fig.subdir, recursive = T)
    }
    
    # heritability/TP settings
    settings <- list(
      list(h = 0.2, add.tp = 0),
      list(h = 0.2, add.tp = 800),
      list(h = 1.0, add.tp = 0),
      list(h = 1.0, add.tp = 800)
    )
    
    # load data
    message("|- Load data ...")
    
    all.data <- lapply(settings, function(setting){
      
      h <- setting$h
      add.tp <- setting$add.tp
      tp <- add.tp + 200
      
      message(sprintf(" |- Heritability: %.1f, additional TP: %d", h, add.tp))
      
      # determine data directory names
      dir.template <- sprintf("out/%%s/30-seasons/h2-%.1f/addTP-%d/normal-effects/BRR", h, add.tp)
      # CGS results with specified diversity weights
      CGS.dir.template <- paste(sprintf(dir.template, "CGS"), sprintf("%s-%%.2f/index", div), sep = "/")
      CGS.dirs <- sapply(div.weights, function(w){
        sprintf(CGS.dir.template, w)
      })
      # corresponding GS/WGS results
      GS.dir <- sprintf(dir.template, "GS")
      WGS.dir <- sprintf(dir.template, "WGS")
      
      # load:
      #      1) GS
      # 2..n-1) CGS with specified series of weights
      #      n) WGS
      data <- lapply(c(GS.dir, CGS.dirs, WGS.dir), load.simulation.results, file.pattern)
      
      return(list(data = data, h = h, tp = tp))
      
    })
    
    # set graphical parameters
    params <- as.list(rep(NA, length(div.weights)+2))
    
    # GS
    params[[1]] <- list(lty = 1, bg = "black", pch = 24)
    # CGS
    colors <- gray.colors(length(div.weights), start = 0.5, end = 0.9)
    params[2:(length(params)-1)] <- lapply(colors, function(col){
      c(list(bg = col), list(lty = 1, pch = 21))
    })
    # WGS
    params[[length(params)]] <- list(lty = 2, bg = "white", pch = 25)
    
    # set curve names
    names <- as.list(rep(NA, length(params)))
    
    # GS
    names[[1]] <- bquote(GS)
    # CGS
    names[2:(length(names)-1)] <- sapply(div.weights, function(w){
      bquote(CGS ~ (alpha == .(sprintf("%.2f", w))))
    })
    # WGS
    names[[length(names)]] <- bquote(WGS)
    # convert to expressions
    names <- as.expression(names)
    
    # setup plot functions
    ylims <- list(
      'gain' = c(0, 0.75),
      'QTL-fixed' = c(0, 1.1),
      'QTL-fav-allele-freq' = c(0.48, 0.81),
      'LD' = c(0.1, 0.9),
      'inbreeding' = c(0, 0.68),
      'genetic-sd' = c(0, 0.11),
      'fav-QTL-lost' = c(0, 36),
      'eff-acc' = c(0.35, 0.72),
      'eff-acc-corrected' = c(0.45, 1.02)
    )
    plot.functions <- get.plot.functions(ylims)
    # init function to extend plot title
    make.title <- function(title, h, tp){
      h.formatted <- sprintf("%.1f", h)
      bquote(.(title) ~ (.(substitute(list(.div, h^2 == .h2, TP == .tp), list(.div = div, .h2 = h.formatted, .tp = tp)))))
    }
    
    # create plots
    message("|- Create plots ...")
    
    for(plot.fun in plot.functions){
      
      file <- sprintf("%s/%s.pdf", fig.subdir, plot.fun$name)

      create.pdf(file, width = 16, height = 12,  function(){
        
        # combine plots for different heritabilities and TP sizes
        par(mfrow = c(2,2))
        
        for(data in all.data){
          # plot
          plot.multi(data$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim)
          # extend title (include heritability and TP size)
          title(make.title(plot.fun$title, data$h, data$tp))
          add.legend(names, params, pos = plot.fun$legend) 
        }
        
      })
      
    }
    
  }
  
}

# stores PDF plots in "figures/simulation/GS-WGS"
plot.GS.vs.WGS <- function(file.pattern = "*.RDS", xlim = c(0,30)){
  
  fig.dir <- "figures/simulation/GS-WGS"
  
  message("Load data ...")
  
  # load data:
  #  1)  GS, h2 = 0.2, add TP = 0
  #  2)  GS, h2 = 0.2, add TP = 800
  #  3)  GS, h2 = 1.0, add TP = 0
  #  4)  GS, h2 = 1.0, add TP = 800
  #  5) WGS, h2 = 0.2, add TP = 0
  #  6) WGS, h2 = 0.2, add TP = 800
  #  7) WGS, h2 = 1.0, add TP = 0
  #  8) WGS, h2 = 1.0, add TP = 800
  dirs <- sort(Sys.glob("out/[GW]*S/30-seasons/h2-*/addTP-*/normal-effects/BRR"))
  data <- lapply(dirs, load.simulation.results, file.pattern)
  # group small and large TP results
  small.TP <- data[c(1,3,5,7)]
  large.TP <- data[c(2,4,6,8)]
  results <- list(
    list(data = small.TP, title.suffix = "(TP = 200)"),
    list(data = large.TP, title.suffix = "(TP = 1000)")
  )
  
  # set graphical parameters
  params <- list(
    # GS, h2 = 0.2
    list(lty = 1, bg = "black", pch = 24),
    # GS, h2 = 1.0
    list(lty = 1, bg = "black", pch = 21),
    # WGS, h2 = 0.2
    list(lty = 2, bg = "white", pch = 24),
    # WGS, h2 = 1.0
    list(lty = 2, bg = "white", pch = 21)
  )
  # set curve names
  names <- c(
    expression(GS ~ (h^2 == 0.2)),
    expression(GS ~ (h^2 == 1.0)),
    expression(WGS ~ (h^2 == 0.2)),
    expression(WGS ~ (h^2 == 1.0))
  )
  
  message("Create plots ...")
  
  # setup plot functions
  ylims <- list(
    'gain' = c(0, 0.7),
    'QTL-fixed' = c(0, 1),
    'QTL-fav-allele-freq' = c(0.48, 0.78),
    'LD' = c(0.1, 0.9),
    'inbreeding' = c(0, 0.68),
    'genetic-sd' = c(0, 0.075),
    'fav-QTL-lost' = c(0, 36),
    'eff-acc' = c(0.35, 0.72),
    'eff-acc-corrected' = c(0.45, 0.97)
  )
  plot.functions <- get.plot.functions(ylims)

  # create plots
  for(plot.fun in plot.functions){
  
    file <- sprintf("%s/%s.pdf", fig.dir, plot.fun$name)

    create.pdf(file, width = 16, height = 6, function(){
      # combine small/large TP plots
      par(mfrow = c(1,2))
      for(res in results){
        plot.multi(res$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim)
        title(sprintf("%s %s", plot.fun$title, res$title.suffix))
        add.legend(names, params, pos = plot.fun$legend)
      }
    })
      
  }
  
}

##################
# PLOT FUNCTIONS #
##################

add.legend <- function(names, param.list = list(), pos = "bottomright", inset = c(0.02, 0.02)){
  # expand parameters
  n <- length(names)
  expanded.params <- as.list(rep(NA, n))
  for(s in 1:n){
    # retrieve graphical parameters (cyclically reused)
    p <- ((s-1) %% length(param.list)) + 1
    expanded.params[[s]] <- param.list[[p]]
  }
  # extract lty, pch and bg
  lty <- sapply(expanded.params, function(params){ params$lty })
  pch <- sapply(expanded.params, function(params){ params$pch })
  bg <- sapply(expanded.params, function(params){ params$bg })
  # add legend to active plot
  legend(x = pos, inset = inset, legend = names, lty = lty, pch = pch, pt.bg = bg)
}

plot.multi <- function(simulations, plot.function, param.list = list(), xlim, ylim, same.plot = TRUE){
  for(s in 1:length(simulations)){
    # retrieve plot function parameters (cyclically reused)
    p <- ((s-1) %% length(param.list)) + 1
    params <- param.list[[p]]
    # set simulation to plot
    sim <- simulations[[s]]
    params$replicates <- sim
    # add to same plot if requested
    add <- same.plot && s > 1
    params$add <- add
    # add ylim to params if set
    if(!missing(ylim)){
      params$ylim <- ylim
    }
    # add xlim to params if set
    if(!missing(xlim)){
      params$xlim <- xlim
    }
    # plot simulation results
    do.call(plot.function, params)
  }
}

# plot a certain variable extracted from the simulations:
#  --> input:   list of replicates (each replicate is a list of simulated seasons)
#  --> plotted: values of variable extracted with the given function 'extract.values',
#               averaged over all replicates, with Wald confidence intervals (by default 95%)
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
  value.avg <- colMeans(values, na.rm = TRUE)
  value.std.err <- apply(values, 2, sd, na.rm = TRUE)/sqrt(num.rep)
  value.ci.halfwidth <- qt(ci+(1-ci)/2, df = num.rep-1) * value.std.err
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
                                 type = c("ped", "geno"),
                                 ylab = "Mean inbreeding coefficient",
                                 ...){
  
  # get selected type
  type <- match.arg(type)
  
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
        mean.inbr[s] <- mean(season$candidates$inbreeding[[type]])
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
                                            corrected = FALSE,
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

# plot ratio of fixed QTL
plot.ratio.fixed.QTL <- function(replicates,
                                 ylab = "Ratio of fixed QTL",
                                 ...){
  
  # set function to extract ratio of fixed QTL
  extract.ratio.fixed.QTL <- function(seasons){
    # initialize result vector
    ratio.fixed <- rep(NA, length(seasons))
    # extract ratio for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether selection candidates have been produced in this season
      if(!is.null(season$candidates)){
        # extract and store ratio
        fav.QTL.allele.freqs <- season$candidates$fav.QTL.allele.freqs
        mafs <- pmin(fav.QTL.allele.freqs, 1 - fav.QTL.allele.freqs)
        ratio <- sum(mafs == 0)/length(fav.QTL.allele.freqs)
        ratio.fixed[s] <- ratio
      }
    }
    return(ratio.fixed)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.ratio.fixed.QTL, ylab = ylab, ...)
  
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