library(scales)
library(animation)

################
# LOAD RESULTS #
################

# create list containing simulated breeding cycles, read from all
# files found in the given directory that match the given file pattern
load.simulation.results <- function(dir, file.pattern = "bp-*.RDS") {
  # get file paths
  files <- Sys.glob(sprintf("%s/%s", dir, file.pattern))
  # load list of simulated breeding cycles
  breeding.cycles <- lapply(files, function(file){
    message(file)
    readRDS(file)
  })
  # return results
  return(breeding.cycles)
}

##########################################
# INFER FINAL GAIN OF VARIOUS STRATEGIES #
##########################################

# averaged over all replicates and all combinationsn of h2 + additional TP
infer.avg.final.gain <- function(type = c("GS", "WGS", "CGS"),
                                 # ignored for type GS/WGS
                                 div.measure = c("HEall", "HEfav", "LOGall", "LOGfav"), div.weight){
  
  # process arguments
  type <- match.arg(type)
  div.measure <- match.arg(div.measure)
  
  # list data directories
  if(type == "CGS"){
    # CGS
    dirs <- Sys.glob(sprintf("out/CGS/30-seasons/h2-*/addTP-*/normal-effects/BRR/%s-%.2f/index", div.measure, div.weight))
  } else {
    # GS/WGS
    dirs <- Sys.glob(sprintf("out/%s/30-seasons/h2-*/addTP-*/normal-effects/BRR", type))
  }
  
  # load data and extract final gains
  final.gains <- sapply(dirs, function(dir){
    simulations <- load.simulation.results(dir)
    sim.final.gains <- sapply(simulations, function(simulation){
      general <- simulation[[1]]$general
      base.pop <- simulation[[1]]$candidates
      initial.value <- normalize.genetic.values(mean(base.pop$geneticValues), general$qtl.effects)
      final.selection <- simulation[[length(simulation)]]$selection
      final.value <- normalize.genetic.values(mean(final.selection$geneticValues), general$qtl.effects)
      final.gain <- final.value - initial.value
      return(final.gain)
    })
    sim.avg.final.gain <- mean(sim.final.gains)
    return(sim.avg.final.gain)
  })
  # average over all replicates and h2/TP settings
  avg.final.gain <- mean(final.gains)
  
  return(avg.final.gain)
  
}

infer.CGS.avg.final.gains <- function(div.measure = c("HEall", "HEfav", "LOGall", "LOGfav"),
                                      div.weights = seq(0.35, 1.0, 0.05)){
  
  # process arguments
  div.measure <- match.arg(div.measure)
  
  avg.final.gains <- sapply(div.weights, function(w){
    message("Processing weight ", w, " ...")
    return(infer.avg.final.gain(type = "CGS", div.measure = div.measure, div.weight = w))
  })
  names(avg.final.gains) <- div.weights
  
  return(avg.final.gains)
  
}

#############################
# AUTOMATED PLOT GENERATION #
#############################

create.pdf <- function(file, plot.fun, width = 14, height = 10.5){
  
  pdf(file, width = width, height = height)
  plot.fun()
  invisible(dev.off())
  
}

get.plot.functions <- function(){
  plot.functions <- list(
    list(f = plot.genetic.gain, name = "gain", title = "Genetic gain",
         legend = "bottomright", ylim =  c(0, 0.3)),
    list(f = plot.proportion.fixed.QTL, name = "QTL-fixed", title = "Proportion of fixed QTL",
         legend = "bottomright", ylim = c(0, 1.0)),
    list(f = plot.mean.QTL.fav.allele.freq, name = "QTL-fav-allele-freq", title = "Mean QTL favourable allele frequency",
         legend = "bottomright", ylim = c(0.50, 0.625)),
    list(f = plot.mean.QTL.marker.LD, name = "LD", title = "Mean polymorphic QTL - marker LD",
    #      legend = "bottomleft", ylim = c(0.3, 0.9)),
    list(f = plot.inbreeding.rate, name = "inbreeding-rate", title = "Inbreeding rate",
         legend = "topleft", ylim = c(0, 0.8)),
    list(f = plot.genetic.standard.deviation, name = "genetic-sd", title = "Genetic standard deviation",
         legend = "topright", ylim = c(0, 0.035)),
    list(f = plot.num.fav.QTL.lost, name = "fav-QTL-lost", title = "Number of favourable QTL lost",
         legend = "bottomright", ylim = c(0, 450)),
    list(f = plot.effect.estimation.accuracy, name = "eff-acc", title = "Effect estimation accuracy",
         legend = "bottomright", ylim = c(0.1, 0.45)),
    list(
      f = function(...){
        plot.effect.estimation.accuracy(..., corrected = TRUE)
      },
      name = "eff-acc-corrected", title = "Corrected effect estimation accuracy",
      legend = "bottomright", ylim = c(0.1, 0.45)
    ),
    list(f = plot.effect.sign.mismatches, name = "sign-mismatches", title = "Proportion of effect sign mismatches",
         legend = "topright", ylim = c(0.38, 0.5)),
    list(
      f = function(...){
        plot.effect.sign.mismatches(..., max.maf = 0.10)
      },
      name = "sign-mismatches-maf-0.10", title = "Proportion of effect sign mismatches for rare alleles",
      legend = "topright", ylim = c(0.38, 0.5)
    ),
    list(
      f = function(...){
        plot.effect.sign.mismatches(..., max.maf = 0.05)
      },
      name = "sign-mismatches-maf-0.05", title = "Proportion of effect sign mismatches for rare alleles",
      legend = "topright", ylim = c(0.38, 0.5)
    ),
    list(
      f = function(...){
        plot.effect.sign.mismatches(..., eff.quantile = 0.25)
      },
      name = "sign-mismatches-eff-quant-0.25", title = "Proportion of effect sign mismatches",
      legend = "topright", ylim = c(0.38, 0.5)
    )
  )
  return(plot.functions)
}

# stores PDF plot in "figures/simulation/CGS";
# organized into subfolders according to
#  1) the two included heritabilities
#  2) the given optimal strategy name
#  3) whether confidence intervals are included (ci)
plot.CGS.opt <- function(strategy.name = "OPT-high-short-term-gain",
                         HE.all.weight = 0.35, HE.fav.weight = 0.50, LOG.all.weight = 0.35, LOG.fav.weight = 0.50,
                         heritability = c(0.2, 0.5), file.pattern = "bp-*.RDS",
                         xlim = c(0,30), ci = NA, main.plots = TRUE, MDS.methods = FALSE,
                         MDS.pops = FALSE){
  
  # check: two heritabilities
  if(length(heritability) != 2){
    stop("'heritability' should be a vector of length 2")
  }
  # extract low and high heritability
  low.h <- min(heritability)
  high.h <- max(heritability)
  
  # combine weights
  div.weights <- c(HE.all.weight, HE.fav.weight, LOG.all.weight, LOG.fav.weight)
  names(div.weights) <- c("HEall", "HEfav", "LOGall", "LOGfav")
  div.measure.labels.main <- c("HE", "HE", "LOG", "LOG")
  names(div.measure.labels.main) <- names(div.weights)
  div.measure.labels.subscript <- c("all", "fav", "all", "fav")
  names(div.measure.labels.subscript) <- names(div.weights)
  
  fig.dir <- sprintf("figures/simulation/CGS/h2-%.1f-%.1f/%s", low.h, high.h, strategy.name)
  if(!is.na(ci)){
    fig.dir <- sprintf("%s-ci-%.2f", fig.dir, ci)
  } else {
    fig.dir <- paste(fig.dir, "no-ci", sep = "-")
  }
  
  if(!dir.exists(fig.dir)){
    message(sprintf("|- Create output directory \"%s\"", fig.dir))
    dir.create(fig.dir, recursive = T)
  }
  
  # heritability/TP settings
  settings <- list(
    list(h = low.h, add.tp = 0),
    list(h = low.h, add.tp = 800),
    list(h = high.h, add.tp = 0),
    list(h = high.h, add.tp = 800)
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
    # CGS results with various diversity measures and selected weights
    CGS.dir.template <- paste(sprintf(dir.template, "CGS"), "%s-%.2f/index", sep = "/")
    CGS.dirs <- sapply(names(div.weights), function(div.measure){
      sprintf(CGS.dir.template, div.measure, div.weights[div.measure])
    })
    # corresponding GS/WGS results
    GS.dir <- sprintf(dir.template, "GS")
    WGS.dir <- sprintf(dir.template, "WGS")
    
    # load:
    #    1) GS
    # 2..5) CGS with various diversity measures and selected weights
    #    6) WGS
    data <- lapply(c(GS.dir, CGS.dirs, WGS.dir), load.simulation.results, file.pattern)
    
    return(list(data = data, h = h, tp = tp))
    
  })
  
  # set graphical parameters
  params <- list(
    list(lty = 1, bg = "black", pch = 23), # GS
    list(lty = 3, bg = "black", pch = 24), # CGS: HEall
    list(lty = 3, bg = "white", pch = 24), # CGS: HEfav
    list(lty = 4, bg = "black", pch = 25), # CGS: LOGall
    list(lty = 4, bg = "white", pch = 25), # CGS: LOGfav
    list(lty = 2, bg = "white", pch = 21)  # WGS
  )
  
  # set curve names
  names <- as.list(rep(NA, length(params)))
  
  # GS
  names[[1]] <- bquote(GS)
  # CGS
  names[2:(length(names)-1)] <- sapply(names(div.weights), function(div.measure){
    bquote(CGS ~ (.(div.measure.labels.main[div.measure])[.(div.measure.labels.subscript[div.measure])]: alpha == .(sprintf("%.2f", div.weights[div.measure]))))
  })
  # WGS
  names[[length(names)]] <- bquote(WGS)
  # convert to expressions
  names <- as.expression(names)
  
  # setup plot functions
  plot.functions <- get.plot.functions()
  # init function to extend plot title
  make.title <- function(title, h, tp){
    h.formatted <- sprintf("%.1f", h)
    bquote(.(title) ~ (.(substitute(list(h^2 == .h2, TP == .tp), list(.h2 = h.formatted, .tp = tp)))))
  }
  
  # create plots
  message("|- Create plots ...")
  
  if(main.plots){
    for(plot.fun in plot.functions){
      
      message(" |- ", plot.fun$name, " ...")
      
      # plot with all four settings
      file <- sprintf("%s/%s.pdf", fig.dir, plot.fun$name)
      
      create.pdf(file, function(){
        
        # combine plots for different heritabilities and TP sizes
        par(mfrow = c(2,2))
        
        for(data in all.data){
          # plot
          plot.multi(data$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim, ci = ci)
          # extend title (include heritability and TP size)
          title(make.title(plot.fun$title, data$h, data$tp))
          add.legend(names, params, pos = plot.fun$legend) 
        }
        
      })
      
      # plot with 2 extreme settings only
      file <- sprintf("%s/%s-extremes.pdf", fig.dir, plot.fun$name)
      
      create.pdf(file, function(){
        
        par(mfrow = c(1,2))
        
        for(data in all.data[c(1,4)]){
          # plot
          plot.multi(data$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim, ci = ci)
          # extend title (include heritability and TP size)
          title(make.title(plot.fun$title, data$h, data$tp))
          add.legend(names, params, pos = plot.fun$legend) 
        }
        
      }, height = 5.25)
      
      # plot with low heritability only
      file <- sprintf("%s/%s-h2-%.1f.pdf", fig.dir, plot.fun$name, low.h)
      
      create.pdf(file, function(){
        
        par(mfrow = c(1,2))
        
        for(data in all.data[1:2]){
          # plot
          plot.multi(data$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim, ci = ci)
          # extend title (include heritability and TP size)
          title(make.title(plot.fun$title, data$h, data$tp))
          add.legend(names, params, pos = plot.fun$legend) 
        }
        
      }, height = 5.25)
      
      # plot with high heritability only
      file <- sprintf("%s/%s-h2-%.1f.pdf", fig.dir, plot.fun$name, high.h)
      
      create.pdf(file, function(){
        
        par(mfrow = c(1,2))
        
        for(data in all.data[3:4]){
          # plot
          plot.multi(data$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim, ci = ci)
          # extend title (include heritability and TP size)
          title(make.title(plot.fun$title, data$h, data$tp))
          add.legend(names, params, pos = plot.fun$legend) 
        }
        
      }, height = 5.25)
      
    }
  }
  
  # MDS plots
  
  if(MDS.methods){
    
    # comparison of different methods, averaged over several BP
    file <- sprintf("%s/%s.pdf", fig.dir, "mds-methods")
    
    create.pdf(file, function(){
      
      # shorten names
      names[2:(length(names)-1)] <- sapply(names(div.weights), function(div.measure){
        bquote(.(div.measure.labels.main[div.measure])[.(div.measure.labels.subscript[div.measure])])
      })
      # plot
      plot.MDS.methods(all.data, names, cex = 1.5)
      title("Similarity of final selections")
      
    }, width = 8, height = 6)
    
  }
    
  # detailed populations
  if(MDS.pops){
    
    settings <- list(
      #list(h2 = 0.2, addTP = 0),
      list(h2 = 0.2, addTP = 800)
      #list(h2 = 0.5, addTP = 800)
    )
    
    types <- list(
      "qtl",
      "markers"
    )
    
    for(setting in settings){
      for(type in types){
        for(bp in 1:5){

          name <- sprintf("mds-detail-h2-%.1f-addTP-%d-BP-%d-%s", setting$h2, setting$addTP, bp, type)
          
          # plots
          message("|- Detailed MDS: plots")
          
          file <- sprintf("%s/%s.pdf", fig.dir, name)
          
          create.pdf(file, function(){
            
            par(mfrow = c(3,3))
            gen <- c(1, seq(2, 30, 4))
            plot.MDS.populations.CGS.opt(type,
                                         HE.all.weight, HE.fav.weight, LOG.all.weight, LOG.fav.weight,
                                         gen, bp, setting$h2, setting$addTP)
            
          }, width = 12, height = 12)
          
          # movie
          message("|- Detailed MDS: movies")
          
          file <- sprintf("%s/%s.mp4", fig.dir, name)
          
          saveVideo({
            
            gen <- 1:30
            plot.MDS.populations.CGS.opt(type,
                                         HE.all.weight, HE.fav.weight, LOG.all.weight, LOG.fav.weight,
                                         gen, bp, setting$h2, setting$addTP, include.gain = TRUE)
            
          }, video.name = file, ani.height = 500, ani.width = 1000, interval = 0.5, autobrowse = FALSE)
          
        }
      }
    }

  }
  
}

plot.CGS.opt.all <- function(heritability = c(0.2, 0.5), file.pattern = "bp-*.RDS", xlim = c(0,30), ci = NA,
                             main.plots = TRUE, MDS.methods = FALSE, MDS.pops = FALSE){
  
  # optimal strategies
  opt.strategies <- list(
    # short-term gain at least as high as Jannink's WGS
    list(name = "OPT-high-short-term-gain", HE.all.weight = 0.35, HE.fav.weight = 0.5, LOG.all.weight = 0.35, LOG.fav.weight = 0.5),
    # moderate short-term gain (at least that of GS in previous generation)
    list(name = "OPT-moderate-short-term-gain", HE.all.weight = 0.45, HE.fav.weight = 0.65, LOG.all.weight = 0.45, LOG.fav.weight = 0.65),
    # maximum long-term gain
    list(name = "OPT-max-long-term-gain", HE.all.weight = 0.5, HE.fav.weight = 0.75, LOG.all.weight = 0.50, LOG.fav.weight = 0.75)
  )
  
  for(strat in opt.strategies){
    message("Strategy: ", strat$name)
    plot.CGS.opt(strat$name,
                 strat$HE.all.weight, strat$HE.fav.weight, strat$LOG.all.weight, strat$LOG.fav.weight,
                 heritability, file.pattern, xlim, ci, main.plots, MDS.methods, MDS.pops)
  }
  
}

# stores PDF plots in "figures/simulation/CGS";
# organized into subfolders according to
#  1) the two included heritabilities
#  2) the applied diversity measure and weight (each combination of the given measures and weights)
plot.CGS <- function(div.weights = seq(0.35, 1.0, 0.05), div.measures = c("HEall", "HEfav", "LOGall", "LOGfav"),
                     heritability = c(0.2, 0.5), file.pattern = "bp-*.RDS", xlim = c(0,30)){
  
  # check: two heritabilities
  if(length(heritability) != 2){
    stop("'heritability' should be a vector of length 2")
  }
  # extract low and high heritability
  low.h <- min(heritability)
  high.h <- max(heritability)
  
  fig.dir <- sprintf("figures/simulation/CGS/h2-%.1f-%.1f", low.h, high.h)

  # create plots for each combination of diversity measure and weight
  for(div in div.measures){
    for(w in div.weights){
     
      message("Diversity measure: ", div, " (weight: ", w, ")")
      
      fig.subdir <- sprintf("%s/%s-%.2f", fig.dir, div, w)
      
      if(!dir.exists(fig.subdir)){
        message(sprintf("|- Create output directory \"%s\"", fig.subdir))
        dir.create(fig.subdir, recursive = T)
      }
      
      # heritability/TP settings
      settings <- list(
        #list(h = low.h, add.tp = 0),
        list(h = low.h, add.tp = 800)
        #list(h = high.h, add.tp = 0),
        #list(h = high.h, add.tp = 800)
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
        # CGS results with specified diversity weight
        CGS.dir <- paste(sprintf(dir.template, "CGS"), sprintf("%s-%.2f/index", div, w), sep = "/")
        # corresponding GS/WGS results
        GS.dir <- sprintf(dir.template, "GS")
        WGS.dir <- sprintf(dir.template, "WGS")
        
        # load:
        # 1) GS
        # 2) CGS with specified weight
        # 3) WGS
        data <- lapply(c(GS.dir, CGS.dir, WGS.dir), load.simulation.results, file.pattern)
        
        return(list(data = data, h = h, tp = tp))
        
      })
      
      # set graphical parameters
      params <- as.list(rep(NA, 3))
      
      # GS
      params[[1]] <- list(lty = 1, bg = "black", pch = 24)
      # CGS
      params[[2]] <- list(lty = 1, bg = "grey", pch = 21)
      # WGS
      params[[3]] <- list(lty = 2, bg = "white", pch = 25)
      
      # set curve names
      names <- as.list(rep(NA, length(params)))
      
      # GS
      names[[1]] <- bquote(GS)
      # CGS
      names[[2]] <- bquote(CGS ~ (alpha == .(sprintf("%.2f", w))))
      # WGS
      names[[3]] <- bquote(WGS)
      # convert to expressions
      names <- as.expression(names)
      
      # setup plot functions
      plot.functions <- get.plot.functions()
      # init function to extend plot title
      make.title <- function(title, h, tp){
        h.formatted <- sprintf("%.1f", h)
        bquote(.(title) ~ (.(substitute(list(.div, h^2 == .h2, TP == .tp), list(.div = div, .h2 = h.formatted, .tp = tp)))))
      }
      
      # create plots
      message("|- Create plots ...")
      
      for(plot.fun in plot.functions){
        
        if(plot.fun$name == "gain"){
        
          file <- sprintf("%s/%s.pdf", fig.subdir, plot.fun$name)
          
          create.pdf(file, function(){
            
            # combine plots for different heritabilities and TP sizes
            par(mfrow = c(2,2))
            
            for(data in all.data){
              # plot
              plot.multi(data$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim, ci = NA)
              # extend title (include heritability and TP size)
              title(make.title(plot.fun$title, data$h, data$tp))
              add.legend(names, params, pos = plot.fun$legend) 
            }
            
          })
          
        }
        
      }
       
    }
  }
  
}

# stores PDF plots in "figures/simulation/GS-WGS",
# within a subfolder according to the two included heritabilities
plot.GS.vs.WGS <- function(heritability = c(0.2, 0.5), file.pattern = "bp-*.RDS", xlim = c(0,30), ci = NA){
  
  # check: two heritabilities
  if(length(heritability) != 2){
    stop("'heritability' should be a vector of length 2")
  }
  # extract low and high heritabilities
  low.h <- min(heritability)
  high.h <- max(heritability)
  
  fig.dir <- sprintf("figures/simulation/GS-WGS/h2-%.1f-%.1f", low.h, high.h)
  
  # create output directory
  if(!dir.exists(fig.dir)){
    message(sprintf("Create output directory \"%s\"", fig.dir))
    dir.create(fig.dir, recursive = T)
  }
  
  message("Load data ...")
  
  # load data
  dirs.low.h2 <- sort(Sys.glob(sprintf("out/[GW]*S*/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR", low.h)))
  dirs.high.h2 <- sort(Sys.glob(sprintf("out/[GW]*S*/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR", high.h)))
  data.low.h2 <- lapply(dirs.low.h2, load.simulation.results, file.pattern)
  data.high.h2 <- lapply(dirs.high.h2, load.simulation.results, file.pattern)
  # combine data:
  # -- low h2 --
  #  1)  GS, h2: low, add TP: 0
  #  2)  GS, h2: low, add TP: 800
  #  3) WGS, h2: low, add TP: 0
  #  4) WGS, h2: low, add TP: 800
  #  5) WGS2, h2: low, add TP: 0
  #  6) WGS2, h2: low, add TP: 800
  # -- high h2 --
  #  7)  GS, h2: high, add TP: 0
  #  8)  GS, h2: high, add TP: 800
  #  9) WGS, h2: high, add TP: 0
  # 10) WGS, h2: high, add TP: 800
  # 11) WGS2, h2: high, add TP: 0
  # 12) WGS2, h2: high, add TP: 800
  data <- c(data.low.h2, data.high.h2)
  # group small and large TP results
  small.TP <- data[c(1,3,5,7,9,11)]
  large.TP <- data[c(2,4,6,8,10,12)]
  results <- list(
    list(data = small.TP, title.suffix = "(TP = 200)"),
    list(data = large.TP, title.suffix = "(TP = 1000)")
  )
  
  # set graphical parameters
  params <- list(
    # GS, h2 = low
    list(lty = 1, bg = "black", pch = 24),
    # WGS, h2 = low
    list(lty = 2, bg = "white", pch = 24),
    # WGS2, h2 = low
    list(lty = 3, bg = "grey", pch = 24),
    # GS, h2 = high
    list(lty = 1, bg = "black", pch = 21),
    # WGS, h2 = high
    list(lty = 2, bg = "white", pch = 21),
    # WGS2, h2 = high
    list(lty = 3, bg = "grey", pch = 21)
  )
  # set curve names
  names <- c(
    bquote(GS ~ (h^2 == .(low.h))),
    bquote(WGS ~ (h^2 == .(low.h))),
    bquote(WGS2 ~ (h^2 == .(low.h))),
    bquote(GS ~ (h^2 == .(high.h))),
    bquote(WGS ~ (h^2 == .(high.h))),
    bquote(WGS2 ~ (h^2 == .(high.h)))
  )
  names <- sapply(names, as.expression)
  
  message("Create plots ...")
  
  # setup plot functions
  plot.functions <- get.plot.functions()

  # create plots
  for(plot.fun in plot.functions){
  
    file <- sprintf("%s/%s.pdf", fig.dir, plot.fun$name)

    create.pdf(file, function(){
      # combine small/large TP plots
      par(mfrow = c(1,2))
      for(res in results){
        plot.multi(res$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim, ci = ci)
        title(bquote(.(plot.fun$title) ~ .(res$title.suffix)))
        add.legend(names, params, pos = plot.fun$legend)
      }
    }, height = 6)
      
  }
  
}

# stores PDF plots in "figures/simulation/GS-WGS-OC-[delta.F]",
# within a subfolder according to the two included heritabilities
plot.GS.WGS.OC <- function(heritability = c(0.2, 0.5), file.pattern = "bp-*.RDS", xlim = c(0,30), ci = NA,
                           delta.F = c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005)){
  
  # check: two heritabilities
  if(length(heritability) != 2){
    stop("'heritability' should be a vector of length 2")
  }
  # extract low and high heritabilities
  low.h <- min(heritability)
  high.h <- max(heritability)
  
  for(dF in delta.F){
    
    message(sprintf("Delta F: %.5f", dF))
    fig.dir <- sprintf("figures/simulation/GS-WGS-OC-%.5f/h2-%.1f-%.1f", dF, low.h, high.h)
    
    # create output directory
    if(!dir.exists(fig.dir)){
      message(sprintf("Create output directory \"%s\"", fig.dir))
      dir.create(fig.dir, recursive = T)
    }
    
    message("Load data ...")
    
    # load data
    dirs.low.h2 <- c(
      sort(Sys.glob(sprintf("out/[GW]*S/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR", low.h))),
      sort(Sys.glob(sprintf("out/OC/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR/dF-%.5f", low.h, dF)))
    )
    dirs.high.h2 <- c(
      sort(Sys.glob(sprintf("out/[GW]*S/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR", high.h))),
      sort(Sys.glob(sprintf("out/OC/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR/dF-%.5f", high.h, dF)))
    )
    data.low.h2 <- lapply(dirs.low.h2, load.simulation.results, file.pattern)
    data.high.h2 <- lapply(dirs.high.h2, load.simulation.results, file.pattern)
    # combine data:
    # -- low h2 --
    #  1)  GS, h2: low, add TP: 0
    #  2)  GS, h2: low, add TP: 800
    #  3) WGS, h2: low, add TP: 0
    #  4) WGS, h2: low, add TP: 800
    #  5)  OC, h2: low, add TP: 0
    #  6)  OC, h2: low, add TP: 800
    # -- high h2 --
    #  7)  GS, h2: high, add TP: 0
    #  8)  GS, h2: high, add TP: 800
    #  9) WGS, h2: high, add TP: 0
    # 10) WGS, h2: high, add TP: 800
    # 11)  OC, h2: high, add TP: 0
    # 12)  OC, h2: high, add TP: 800
    data <- c(data.low.h2, data.high.h2)
    # group small and large TP results
    small.TP <- data[c(1,3,5,7,9,11)]
    large.TP <- data[c(2,4,6,8,10,12)]
    results <- list(
      list(data = small.TP, title.suffix = "(TP = 200)"),
      list(data = large.TP, title.suffix = "(TP = 1000)")
    )
    
    # set graphical parameters
    params <- list(
      # GS, h2 = low
      list(lty = 1, bg = "black", pch = 24),
      # WGS, h2 = low
      list(lty = 2, bg = "white", pch = 24),
      # OC, h2 = low
      list(lty = 3, bg = "grey", pch = 24),
      # GS, h2 = high
      list(lty = 1, bg = "black", pch = 21),
      # WGS, h2 = high
      list(lty = 2, bg = "white", pch = 21),
      # OC, h2 = high
      list(lty = 3, bg = "grey", pch = 21)
    )
    # set curve names
    names <- c(
      bquote(GS ~ (h^2 == .(low.h))),
      bquote(WGS ~ (h^2 == .(low.h))),
      bquote(OC ~ (h^2 == .(bquote(paste(.(low.h), ", ", sep = "") ~ paste(Delta, F) == .(dF))))),
      bquote(GS ~ (h^2 == .(high.h))),
      bquote(WGS ~ (h^2 == .(high.h))),
      bquote(OC ~ (h^2 == .(bquote(paste(.(high.h), ", ", sep = "") ~ paste(Delta, F) == .(dF)))))
    )
    names <- sapply(names, as.expression)
    
    message("Create plots ...")
    
    # setup plot functions
    plot.functions <- get.plot.functions()
    
    # create plots
    for(plot.fun in plot.functions){
      
      file <- sprintf("%s/%s.pdf", fig.dir, plot.fun$name)
      
      create.pdf(file, function(){
        # combine small/large TP plots
        par(mfrow = c(1,2))
        for(res in results){
          plot.multi(res$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim, ci = ci)
          title(bquote(.(plot.fun$title) ~ .(res$title.suffix)))
          add.legend(names, params, pos = plot.fun$legend)
        }
      }, height = 6)
      
    }
    
  }
  
}

# stores PDF plots in "figures/simulation/WGS-OC-[delta.F]",
# within a subfolder according to the two included heritabilities
plot.WGS.OC <- function(heritability = c(0.2, 0.5), file.pattern = "bp-*.RDS", xlim = c(0,30), ci = NA,
                           delta.F = c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005)){
  
  # check: two heritabilities
  if(length(heritability) != 2){
    stop("'heritability' should be a vector of length 2")
  }
  # extract low and high heritabilities
  low.h <- min(heritability)
  high.h <- max(heritability)
  
  for(dF in delta.F){
    
    message(sprintf("Delta F: %.5f", dF))
    fig.dir <- sprintf("figures/simulation/WGS-OC-%.5f/h2-%.1f-%.1f", dF, low.h, high.h)
    
    # create output directory
    if(!dir.exists(fig.dir)){
      message(sprintf("Create output directory \"%s\"", fig.dir))
      dir.create(fig.dir, recursive = T)
    }
    
    message("Load data ...")
    
    # load data
    dirs.low.h2 <- c(
      sort(Sys.glob(sprintf("out/WGS/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR", low.h))),
      sort(Sys.glob(sprintf("out/OC/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR/dF-%.5f", low.h, dF)))
    )
    dirs.high.h2 <- c(
      sort(Sys.glob(sprintf("out/WGS/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR", high.h))),
      sort(Sys.glob(sprintf("out/OC/30-seasons/h2-%.1f/addTP-*/normal-effects/BRR/dF-%.5f", high.h, dF)))
    )
    data.low.h2 <- lapply(dirs.low.h2, load.simulation.results, file.pattern)
    data.high.h2 <- lapply(dirs.high.h2, load.simulation.results, file.pattern)
    # combine data:
    # -- low h2 --
    # 1) WGS, h2: low, add TP: 0
    # 2) WGS, h2: low, add TP: 800
    # 3)  OC, h2: low, add TP: 0
    # 4)  OC, h2: low, add TP: 800
    # -- high h2 --
    # 5) WGS, h2: high, add TP: 0
    # 6) WGS, h2: high, add TP: 800
    # 7)  OC, h2: high, add TP: 0
    # 8)  OC, h2: high, add TP: 800
    data <- c(data.low.h2, data.high.h2)
    # group small and large TP results
    small.TP <- data[c(1,3,5,7)]
    large.TP <- data[c(2,4,6,8)]
    results <- list(
      list(data = small.TP, title.suffix = "(TP = 200)"),
      list(data = large.TP, title.suffix = "(TP = 1000)")
    )
    
    # set graphical parameters
    params <- list(
      # WGS, h2 = low
      list(lty = 1, bg = "black", pch = 24),
      # OC, h2 = low
      list(lty = 2, bg = "white", pch = 24),
      # WGS, h2 = high
      list(lty = 1, bg = "black", pch = 21),
      # OC, h2 = high
      list(lty = 2, bg = "white", pch = 21)
    )
    # set curve names
    names <- c(
      bquote(WGS ~ (h^2 == .(low.h))),
      bquote(OC ~ (h^2 == .(bquote(paste(.(low.h), ", ", sep = "") ~ paste(Delta, F) == .(dF))))),
      bquote(WGS ~ (h^2 == .(high.h))),
      bquote(OC ~ (h^2 == .(bquote(paste(.(high.h), ", ", sep = "") ~ paste(Delta, F) == .(dF)))))
    )
    names <- sapply(names, as.expression)
    
    message("Create plots ...")
    
    # setup plot functions
    plot.functions <- get.plot.functions()
    
    # create plots
    for(plot.fun in plot.functions){
      
      file <- sprintf("%s/%s.pdf", fig.dir, plot.fun$name)
      
      create.pdf(file, function(){
        # combine small/large TP plots
        par(mfrow = c(1,2))
        for(res in results){
          plot.multi(res$data, plot.fun$f, params, ylim = plot.fun$ylim, xlim = xlim, ci = ci)
          title(bquote(.(plot.fun$title) ~ .(res$title.suffix)))
          add.legend(names, params, pos = plot.fun$legend)
        }
      }, height = 6)
      
    }
    
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

plot.multi <- function(simulations, plot.function, param.list = list(), xlim, ylim, ci, same.plot = TRUE){
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
    # add ci to params if set
    if(!missing(ci)){
      params$ci <- ci
    }
    # plot simulation results
    do.call(plot.function, params)
  }
}

# plot a certain variable extracted from the simulations:
#  --> input:   list of replicates (each replicate is a list of simulated seasons)
#  --> plotted: values of variable extracted with the given function 'extract.values',
#               averaged over all replicates, with Wald confidence intervals (by default 95%)
#               (if 'ci' is set to 'NA', no confidence intervals are included)
# the function 'extract.values' should take a single argument, i.e. the list of seasons
# produced from a single simulation run, and return a vector with one extracted value per season
# (for seasons where the value is not reported, NA should be returned)
plot.simulation.variable <- function(replicates,
                                     extract.values,
                                     ylab, shift = 0,
                                     ci = 0.95,
                                     type = c("generations", "seasons"),
                                     add=FALSE, pch=23,
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
  if(!is.na(ci)){
    value.ci.halfwidth <- qt(ci+(1-ci)/2, df = num.rep-1) * value.std.err
    value.ci.top <- value.avg + value.ci.halfwidth
    value.ci.bottom <- value.avg - value.ci.halfwidth
  }
  # plot non NA values only
  non.na <- !is.na(value.avg)
  value.avg <- value.avg[non.na]
  if(!is.na(ci)){
    value.ci.top <- value.ci.top[non.na]
    value.ci.bottom <- value.ci.bottom[non.na]
  }
  x <- (0:num.seasons)[non.na]
  # shift if requested
  if(shift > 0){
    x <- c(0:(shift-1), x+shift)
    value.avg <- c(rep(NA, shift), value.avg)
    if(!is.na(ci)){
      value.ci.top <- c(rep(NA, shift), value.ci.top)
      value.ci.bottom <- c(rep(NA, shift), value.ci.bottom)
    }
  }
  # set x-axis label
  xlab <- ifelse(type == "generations", "Generation", "Season")
  # first plot CI bars (points and lines are plotted on top)
  if(!is.na(ci)){
    errbar(x, value.avg, value.ci.top, value.ci.bottom, type="n",
           xlab=xlab, ylab=ylab,
           xaxp=c(0,num.seasons,num.seasons/2),
           add=add,
           ...)
    points(x, value.avg, type="o", pch=pch, bg=bg, lty=lty, cex = 0.75)
  } else {
    if(add){
      plot.fun <- points
    } else {
      plot.fun <- plot
    }
    plot.fun(x, value.avg, type="o", xlab=xlab, ylab=ylab, xaxp=c(0,num.seasons,num.seasons/2),
             pch=pch, bg=bg, lty=lty, cex = 0.75, ...)
  }
  
}

# set function to extract genetic gains from a simulated list of seasons
extract.gain <- function(seasons, scale = c("jannink", "sd")){
  
  scale <- match.arg(scale)
  
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

# plot genetic gain, with one of these scales:
#  1) in terms of number of standard deviations of genetic value in founder population
#  2) normalized wrt maximal genotypic value possible, as in the paper by Jannink
plot.genetic.gain <- function(replicates,
                              scale = c("jannink", "sd"),
                              ylab = "Genetic gain from selection",
                              ...){
  
  # get selected scale
  scale <- match.arg(scale)
  
  # set function to extract gain with requested scaling
  extract.gain.scaled <- function(seasons){ extract.gain(seasons, scale = scale) }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values = extract.gain.scaled, ylab = ylab, ...)
  
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

# plot inbreeding rate of selection candidates
plot.inbreeding.rate <- function(replicates,
                                 ylab = "Inbreeding rate",
                                 ...){
  
  # set function to extract inbreeding rate
  extract.inbr.rate <- function(seasons){
    # initialize result vector
    inbr.rate <- rep(NA, length(seasons))
    # extract inbreeding rate for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether selection candidates have been produced in this season
      if(!is.null(season$candidates)){
        # extract and store inbreeding rate (if available)
        inbr <- season$candidates$inbreeding
        if(!is.null(inbr)){
          inbr.rate[s] <- inbr
        }
      }
    }
    return(inbr.rate)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.inbr.rate, ylab = ylab, shift = 1, ...)
  
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
  plot.simulation.variable(replicates, extract.values = extract.accuracy, ylab = ylab, shift = 1, ...)
  
}

# plot proportion of effect sign mismatches
plot.effect.sign.mismatches <- function(replicates,
                                        ylab = "Proportion of effect sign mismatches",
                                        max.maf = 0.5, eff.quantile = 0.0,
                                        ...){
  
  # set function to extract mismatches
  extract.mismatches <- function(seasons){
    # initialize result vector
    mismatches <- rep(NA, length(seasons))
    # extract mismatches for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether a GP model has been trained in this season
      if(!is.null(season$gp)){
        # extract QTL-marker LD table
        QTL.marker.LD <- season$gp$qtl.marker.ld
        # infer MAFs from favourable allele freqs
        # NOTE: restrict to considered markers only
        #       based on *names* (indices include QTL and dummies)
        fav.freqs <- season$gp$fav.marker.allele.freqs
        fav.freqs <- fav.freqs[QTL.marker.LD$marker.name]
        mafs <- pmin(fav.freqs, 1 - fav.freqs)
        # infer effect quantile (absolute values)
        marker.effects <- QTL.marker.LD$marker.effect
        abs.effects <- abs(marker.effects)
        quant <- quantile(abs.effects, eff.quantile)
        # restrict markers
        res <- abs.effects >= quant & mafs <= max.maf
        marker.effects <- marker.effects[res]
        qtl.effects <- QTL.marker.LD$QTL.effect[res]
        # compute and store proportion of sign mismatches
        sign.diff <- sign(marker.effects) - sign(qtl.effects)
        sign.mismatches <- sum(sign.diff != 0) / length(marker.effects)
        mismatches[s] <- sign.mismatches
      }
    }
    return(mismatches)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.mismatches, ylab = ylab, shift = 1, ...)
  
}

# plot TP size (after filtering)
plot.tp.size <- function(replicates,
                        ylab = "Training population size",
                        ...){
  
  # set function to extract TP size
  extract.tp.size <- function(seasons){
    # initialize result vector
    tp.sizes <- rep(NA, length(seasons))
    # extract TP size for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether a GP model has been trained in this season
      if(!is.null(season$gp)){
        # extract and store TP size
        tp.sizes[s] <- season$gp$tp.size
      }
    }
    return(tp.sizes)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.tp.size, ylab = ylab, shift = 1, ...)
  
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

# plot proportion of fixed QTL
plot.proportion.fixed.QTL <- function(replicates,
                                      ylab = "Proportion of fixed QTL",
                                      ...){
  
  # set function to extract proportion of fixed QTL
  extract.proportion.fixed.QTL <- function(seasons){
    # initialize result vector
    proportion.fixed <- rep(NA, length(seasons))
    # extract proportion for each season
    for(s in 1:length(seasons)){
      season <- seasons[[s]]
      # check whether selection candidates have been produced in this season
      if(!is.null(season$candidates)){
        # extract and store proportion
        fav.QTL.allele.freqs <- season$candidates$fav.QTL.allele.freqs
        mafs <- pmin(fav.QTL.allele.freqs, 1 - fav.QTL.allele.freqs)
        proportion <- sum(mafs == 0)/length(fav.QTL.allele.freqs)
        proportion.fixed[s] <- proportion
      }
    }
    return(proportion.fixed)
  }
  
  # call generic variable plot function
  plot.simulation.variable(replicates, extract.values =  extract.proportion.fixed.QTL, ylab = ylab, ...)
  
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
  plot.simulation.variable(replicates, extract.values =  extract.mean.marker.fav.allele.freq, ylab = ylab, shift = 1, ...)
  
}

###################################
# Multi-dimensional scaling plots #
###################################

# for this plot, respective replicates of simulations with a different
# selection method and or h2/TP setting should have been started from
# the same base population
plot.MDS.methods <- function(settings, names, ...){

  # initialize distance matrix
  num.points <- length(settings) * length(names)
  d <- matrix(0, nrow = num.points, ncol = num.points)
  settings.ind <- rep(1:length(settings), each = length(names))
  method.ind <- rep(1:length(names), times = length(settings))
  
  # infer number of base pops (number of replicates)
  num.bp <- length(settings[[1]]$data[[1]])
  
  #message("Num BP: ", num.bp)
  # fill distance matrix
  for(bp in 1:num.bp){
    
    #message("|- BP: ", bp)
    
    for(i in 1:num.points){
      for(j in 1:num.points){
        if(i < j){
          
          #message("|-- (i,j) = (", i, ",", j, ")")
          
          # extract final selections
          final.sel <- lapply(c(i,j), function(k){
            s <- settings.ind[k]
            m <- method.ind[k]
            generations <- settings[[s]]$data[[m]][[bp]]
            final.sel <- generations[[length(generations)]]$selection$markers
            return(final.sel)
          })
          
          # compute cluster distance (average pairwise distance)
          all <- rbind(final.sel[[1]], final.sel[[2]])
          el.dist <- as.matrix(dist(all))
          clust.dist <- mean(el.dist[(nrow(final.sel[[1]])+1):nrow(all), 1:nrow(final.sel[[1]])])
          # increment overall inter-cluster distance sum
          d[i, j] <- d[i, j] + clust.dist
          d[j, i] <- d[i, j]
          
        }
      }
    }
    
  }
  
  # average
  d <- d/num.bp

  # MDS
  mds <- cmdscale(d)

  # plot
  col <- settings.ind
  pch <- method.ind + 20
  par(mar = c(2.1,2.1,4.1,2.1))
  plot(mds, col = col, bg = col, pch = pch, xlab = "", ylab = "", xaxt = "n", yaxt = "n", asp = 1, ...)
  
  # legends
  legend("bottom", legend = names, pch = (1:length(names)) + 20)
  setting.names <- lapply(settings, function(setting){
    h.formatted <- sprintf("%.1f", setting$h)
    bquote(.(substitute(list(h^2 == .h2, TP == .tp), list(.h2 = h.formatted, .tp = setting$tp))))
  })
  setting.names <- as.expression(setting.names)
  legend("top", legend = setting.names, pch = 20, col = 1:length(settings))

}

# for this plot, the simulations with each method should have been
# started from the same base population, and full marker data should
# be available for the selection candidates and selection in each generation
plot.MDS.populations <- function(type = c("markers", "qtl"), simulations, generations, method.names, include.gain = FALSE){
  
  type <- match.arg(type)
  
  # check that all simulations were started from same base population
  bp <- check.same.bp(simulations)
  
  # color-blind friendly palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00")
  
  # extract marker/qtl matrices of selection candidates and selection at each generation
  selection.candidates <- unlist(lapply(1:length(generations), function(g.i){
    g <- generations[g.i]
    methods.selection.candidates <- lapply(1:length(simulations), function(m.i){
      sim <- simulations[[m.i]]
      if(type == "markers"){
        Z <- sim[[1+g]]$candidates$markers
        if(is.null(Z)){
          # special case where crossing + inbreeding does not
          # take place in the same (but previous) generation as selection
          Z <- sim[[g]]$candidates$markers
        }
      } else {
        Z <- sim[[1+g]]$candidates$qtl
        if(is.null(Z)){
          # same special case as above
          Z <- sim[[g]]$candidates$qtl
        }
      }
      cbind(g.i, m.i, Z)
    })
  }), recursive = FALSE)
  
  selection <- unlist(lapply(1:length(generations), function(g.i){
    g <- generations[g.i]
    methods.selection <- lapply(1:length(simulations), function(m.i){
      sim <- simulations[[m.i]]
      if(type == "markers"){
        Z <- sim[[1+g]]$selection$markers
      } else {
        Z <- sim[[1+g]]$selection$qtl
      }
      cbind(g.i, m.i, Z)
    })
  }), recursive = FALSE)
  
  # combine all data
  selection.candidates <- Reduce(rbind, selection.candidates)
  selection <- Reduce(rbind, selection)
  all <- rbind(selection.candidates, selection)
  
  # MDS
  d <- dist(all[, 3:ncol(all)])
  mds <- cmdscale(d)
  # get bounds
  xlim <- c(min(mds[, 1]), max(mds[, 1]))
  ylim <- c(min(mds[, 2]), max(mds[, 2]))
  # split
  mds.cand <- mds[1:nrow(selection.candidates), ]
  mds.sel <- mds[(nrow(selection.candidates)+1):nrow(mds), ]
  
  # extract gains if requested
  if(include.gain){
    # extract
    method.gains <- lapply(simulations, function(seasons){
      gains <- extract.gain(seasons)
      gains <- gains[!is.na(gains)]
      return(gains)
    })
    # set xlim
    max.gain <- max(sapply(method.gains, max))
    gain.xlim <- c(0.0, max.gain)
  }
  
  # split across plots (1 per generation)
  for(g.i in 1:length(generations)){
    
    g <- generations[g.i]

    # graphical settings
    par(mar = c(4.1, 1.0, 4.1, 1.0))
    par(xpd = TRUE)
    sel.col <- cbPalette[1:length(method.names)]
    cand.col <- sapply(sel.col, function(c){
      alpha(c, 0.2)
    })
    # 2 plots on 1 row if gain is included
    if(include.gain){
      par(mfrow = c(1,2))
    }

    # init plot
    plot(NULL, type = "n", xlim = xlim, ylim = ylim,
         xlab = "", ylab = "", xaxt = 'n', yaxt = 'n',
         asp = 1)
    
    # 1) plot candidates
    for(m.i in 1:length(method.names)){
      m.cand <- mds.cand[selection.candidates[, 1] == g.i & selection.candidates[, 2] == m.i, ]
      points(m.cand, pch = 20, col = cand.col[m.i])
    }
    
    # 2) plot selections
    for(m.i in 1:length(method.names)){
      m.sel <- mds.sel[selection[, 1] == g.i & selection[, 2] == m.i, ]
      points(m.sel, pch = 24, bg = sel.col[m.i], col = sel.col[m.i])
    }
    
    # add title
    title(sprintf("Generation %d", g))
    # add legend
    legend(
      horiz = TRUE,
      x = "bottom",
      inset = -0.12,
      legend = method.names,
      pch = 24,
      col = sel.col,
      pt.bg = sel.col
    )
    
    # gain progress bars
    if(include.gain){
      
      par(mar = c(4.1, 2.5, 4.1, 1.0))
      barplot(unlist(lapply(method.gains, function(gains){
        gains[[g]]
      })), horiz = TRUE, col = sel.col, names.arg = method.names, xlim = gain.xlim)
      title("Genetic gain")
      
    }

  }
  
}

plot.MDS.populations.CGS.opt <- function(type = c("markers", "qtl"),
                                         HE.all.weight, HE.fav.weight, LOG.all.weight, LOG.fav.weight,
                                         generations, bp, h2, addTP, include.gain = FALSE){
  
  type <- match.arg(type)
  
  # load data
  message("Load data ...")
  # GS/WGS
  GS.data <- readRDS(
    sprintf("out/GS/30-seasons/h2-%.1f/addTP-%d/normal-effects/BRR/bp-%d-1.RDS", h2, addTP, bp)
  )
  WGS.data <- readRDS(
    sprintf("out/WGS/30-seasons/h2-%.1f/addTP-%d/normal-effects/BRR/bp-%d-1.RDS", h2, addTP, bp)
  )
  # CGS
  CGS.HE.all.data <- readRDS(
    sprintf("out/CGS/30-seasons/h2-%.1f/addTP-%d/normal-effects/BRR/HEall-%.2f/index/bp-%d-1.RDS", h2, addTP, HE.all.weight, bp)
  )
  CGS.HE.fav.data <- readRDS(
    sprintf("out/CGS/30-seasons/h2-%.1f/addTP-%d/normal-effects/BRR/HEfav-%.2f/index/bp-%d-1.RDS", h2, addTP, HE.fav.weight, bp)
  )
  CGS.LOG.all.data <- readRDS(
    sprintf("out/CGS/30-seasons/h2-%.1f/addTP-%d/normal-effects/BRR/LOGall-%.2f/index/bp-%d-1.RDS", h2, addTP, LOG.all.weight, bp)
  )
  CGS.LOG.fav.data <- readRDS(
    sprintf("out/CGS/30-seasons/h2-%.1f/addTP-%d/normal-effects/BRR/LOGfav-%.2f/index/bp-%d-1.RDS", h2, addTP, LOG.fav.weight, bp)
  )
  
  # plot
  message("Create plots/movie ...")
  
  plot.MDS.populations(
    type = type,
    list(GS.data, WGS.data, CGS.HE.all.data, CGS.HE.fav.data, CGS.LOG.all.data, CGS.LOG.fav.data),
    generations = generations, method.names = c(
      expression(GS),
      expression(WGS),
      expression(HE[all]),
      expression(HE[fav]),
      expression(LOG[all]),
      expression(LOG[fav])
    ),
    include.gain = include.gain
  )
  
}

check.same.bp <- function(simulations){
  
  bp <- simulations[[1]][[1]]$candidates$markers
  # check that all simulations started from the same base population
  if(length(simulations) >= 2){
    for(s in 2:length(simulations)){
      bp2 <- simulations[[s]][[1]]$candidates$markers
      if(!isTRUE(all.equal(bp, bp2))){
        stop("base population should be equal for all simulations")
      }
    }
  }
  
  return(bp)
  
}












