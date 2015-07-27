library(BGLR)
library(coda)
library(rrBLUP)

#########################
# PREPARE DESIGN MATRIX #
#########################

# create design matrix (0/1, strip dummies, strip real QTL)
gp.design.matrix = function(pop){
  # assert: dummies present
  if(!pop$hypred$dummiesAdded){
    stop("requires that dummies are added, first call add.dummies(pop)")
  }
  # assert: QTL assigned
  if(is.null(pop$hypred$realQTL)){
    stop("requires QTL to be assigned, first call assign.qtl(pop, ...)")
  }
  # assert: DH population
  if(is.null(pop$dh)){
    stop("pop is required to be a DH population")
  }
  # get indices of dummy markers and *real* QTL (dummy QTLs can stay, no real effect there)
  real.qtl.indices = pop$hypred$realQTL
  dummy.ends = pop$hypred$chromBounds[,2]
  dummy.starts = pop$hypred$chromBounds[,2] - pop$hypred$chrNumDummies + 1
  dummy.bounds = cbind(dummy.starts, dummy.ends)
  dummy.indices = unlist(apply(dummy.bounds, 1, function(row){
    if(row[1] <= row[2]){
      seq(row[1], row[2])
    }
  }))
  qtl.or.dummy = c(real.qtl.indices, dummy.indices)
  # get DH marker data without dummies and QTL
  matrix = pop$dh[,-qtl.or.dummy]
  
  return(matrix)
}

# restrict design matrix to subset of all individuals
# names = names of chosen individuals
# G     = full design matrix (ind x markers)
gp.restrict.design.matrix = function(names, G){
  if(!all(names %in% rownames(G))){
    stop ("not all chosen individuals have genotypes")
  }
  # select individuals
  GDesign = G[names,]
  # set row and column names
  rownames(GDesign) = names
  colnames(GDesign) = colnames(G)
  return(GDesign)
}



# FROM HERE: code using rrBLUP

# ############################
# # MARKER EFFECT ESTIMATION #
# ############################
# 
# # train model by estimating marker effects (RR-BLUP method)
# gp.train = function(pheno, Z){
#   
#   # check input
#   
#   if(missing(pheno) || is.null(pheno)){
#     stop("no phenotypes supplied")
#   }
#   if(!is.vector(pheno)){
#     stop("pheno should be a vector of phenotypes (1 value per individual)")
#   }
#   
#   if(missing(Z) || is.null(Z)){
#     stop("no marker data supplied")
#   }
#   if(!is.matrix(Z)){
#     stop("marker data should be a matrix")
#   }
#   if(nrow(Z) != length(pheno)){
#     stop("length of phenotype vector should be equal to number of rows in marker matrix")
#   }
#   
#   # estimate marker effects
#   
#   trained = mixed.solve(y = pheno, Z = Z)
#   
#   # append class name
#   
#   class(trained) = append(class(trained), "rrBLUP")
#   
#   return(trained)
# }
# 
# ####################################
# # PREDICTION OF UNKNOWN PHENOTYPES #
# ####################################
# 
# # predict value of new individuals using a trained GP model
# # if desired, marker effects can be weighted
# gp.predict = function(trained.model, Z, weights=NULL){
#   if(!is(trained.model, "rrBLUP")){
#     stop("trained model should be an rrBLUP object (call gp.train to train model)")
#   }
#   if(missing(Z) || is.null(Z) || !is.matrix(Z)){
#     stop("Z should be a matrix")
#   }
#   marker.effects = trained.model$u
#   if(length(marker.effects) != ncol(Z)){
#     stop("number of columns of Z should correspond to number of markers in trained model")
#   }
#   if(is.null(weights)){
#     weights = rep(1, length(marker.effects))
#   }
#   if(length(weights) != length(marker.effects)){
#     stop("length of weight vector should correspond to number of markers")
#   }
#   values = trained.model$beta + as.vector(Z %*% marker.effects)
#   names(values) = rownames(Z)
#   return(values)
# }
# 
# ################################################
# # DETERMINE FAVOURABLE ALLELES AND FREQUENCIES #
# ################################################
# 
# # get favourable allele at each marker (1 if positive effect, 0 if negative effect)
# get.favourable.alleles = function(trained.model){
#   if(!is(trained.model, "rrBLUP")){
#     stop("input should be an rrBLUP object (call gp.train to train model)")
#   }
#   # get effects
#   marker.effects = trained.model$u
#   # initialize output
#   desired.alleles = rep(0, length(marker.effects))
#   # set to 1 for positive effects
#   desired.alleles[marker.effects > 0] = 1
#   
#   return(desired.alleles)
# }
# 
# # get favourable allele frequencies
# get.favourable.allele.frequencies = function(trained.model, Z){
#   if(!is(trained.model, "rrBLUP")){
#     stop("input should be an rrBLUP object (call gp.train to train model)")
#   }
#   if(missing(Z) || is.null(Z) || !is.matrix(Z)){
#     stop("Z should be a matrix")
#   }
#   if(length(trained.model$u) != ncol(Z)){
#     stop("number of columns of Z should correspond to number of estimated marker effects")
#   }
#   # get favourable alleles
#   fav.alleles = get.favourable.alleles(trained.model)
#   # compute frequencies
#   n = nrow(Z)
#   m = ncol(Z)
#   freqs = rep(NA, m)
#   for(i in 1:m){
#     freqs[i] = sum(Z[,i] == fav.alleles[i]) / n
#   }
#   return(freqs)
# }

# END OF: code using rrBLUP









# FROM HERE: code using BGLR

####################
# PREPARE GP MODEL #
####################

gp.make.priors = function(modelGenetic, df0=NULL, shape0=NULL, R2=NULL,
                          rate0=NULL, S0=NULL, probIn=NULL, countIn=NULL) {
  # set priors based on chosen method
  if(modelGenetic=="BRR"){
    priors = list(df0=df0, S0=S0, R2=R2)
  }
  if(modelGenetic=="BayesA") {
    priors = list(df0=df0, shape0=shape0, R2=R2, rate0=rate0)
  }
  if(modelGenetic=="BayesB") {
    priors = list(df0=df0, shape0=shape0, probIn=probIn, counts=countIn, R2=R2, rate0=rate0)
  }
  if(modelGenetic=="BayesC") {
    priors = list(df0=df0, S0=S0, probIn=1, counts=2, R2=R2)
  }
  # retain only non null priors (return NULL if no priors set)
  priors = priors[!sapply(priors, is.null)]
  if(length(priors) == 0){
    priors = NULL
  }
  return(priors)
}

# prepare chosen model/method and parameters
gp.model = function(modelGenetic, randomGenetic=NULL, fixed=NULL, randomNuisance=NULL,
                    pedigree=NULL, df0=NULL, shape0=NULL, R2=NULL, rate0=NULL,
                    S0=NULL, probIn=NULL, countIn=NULL) {
  # check that all of these are formulas or NULL
  if(!is.null(randomGenetic) && !class(randomGenetic) == "matrix")
    stop ("provided G/Z is not a matrix")
  if(!is.null(fixed) && !class(fixed) == "formula")
    stop ("provided fixed effect are not in proper formula format")
  if(!is.null(randomNuisance) && !class(randomNuisance) == "formula")
    stop ("provided random nuisance effects are not in proper formula format")
  if(!is.null(pedigree) && !class(pedigree) == "matrix")
    stop ("provided pedigree is not a matrix")
  if(missing(modelGenetic) || is.null(modelGenetic))
    stop("no model specified")
  if(!modelGenetic %in% c("BRR", "BayesA", "BayesB", "BayesC"))
    stop("specified model not available")
  # make priors
  priors <- gp.make.priors(modelGenetic=modelGenetic, df0=df0, shape0=shape0, R2=R2,
                           rate0=rate0, S0=S0, probIn=probIn, countIn=countIn)
  # make ETAs
  if(!is.null(fixed)){
    ETAfixed = list(fixed, model="FIXED")
  } else {
    ETAfixed = NULL
  }
  if(!is.null(randomNuisance)){
    ETArandomN = list(randomNuisance, model="BRR")
  } else {
    ETArandomN = NULL
  }
  if(!is.null(randomGenetic)){
    ETArandomG = list(X=randomGenetic, model=modelGenetic)
    if (!is.null(priors)){
      ETArandomG = append(ETArandomG, priors)
    }
  } else {
    ETArandomG = NULL;
  }
  if(!is.null(pedigree)){
    ETApedigree = list(K=pedigree, model="RKHS")
  } else {
    ETApedigree = NULL
  }
  # make complete list
  ETA = list(ETAfixed, ETArandomN, ETArandomG, ETApedigree)
  ETA = ETA[!sapply(ETA, is.null)]
  return(ETA)
}


############################
# MARKER EFFECT ESTIMATION #
############################

# train model by estimating marker effects
# uses Bayesian Ridge Regression by default
gp.train = function(modelGenetic="BRR", pheno, Z=NULL,
                    fixed=NULL, randomNuisance=NULL,
                    pedigree=NULL, df0=NULL, shape0=NULL,
                    R2=NULL, rate0=NULL, S0=NULL, probIn=NULL,
                    countIn=NULL, nIter=1500, burnIn=500, thin=1){
  
  # check phenotypes
  
  if(missing(pheno) || is.null(pheno)){
    stop("no phenotypes supplied")
  }
  if(!is.vector(pheno) || is.null(names(pheno))){
    stop("pheno should be a named vector of phenotypes (1 value per individual)")
  }
  
  # prepare model (checks other parameters)
  model = gp.model(modelGenetic, Z, fixed, randomNuisance,
                   pedigree, df0, shape0, R2, rate0,
                   S0, probIn, countIn)
  
  # estimate marker effects
  
  dir.create("BGLR-tmp", showWarnings = FALSE)
  trained = BGLR(
                  y=pheno,                               # phenotypes (response variable)
                  ETA=model,                             # use chosen model
                  nIter=nIter, burnIn=burnIn, thin=thin, # MCMC parameters
                  saveAt="BGLR-tmp/",
                  verbose=FALSE
                )
  
  return(trained)
}

####################################
# PREDICTION OF UNKNOWN PHENOTYPES #
####################################

# predict value of new individuals using a trained GP model
# the mean (mu) is ignored as shifting does not influence the outcome of the selection criteria
# if desired, marker effects can be weighted
gp.predict = function(trained.model, Z, weights=NULL){
  if(class(trained.model) != "BGLR"){
    stop("input should be a BGLR object (call gp.train to train model)")
  }
  if(missing(Z) || is.null(Z) || !is.matrix(Z)){
    stop("Z should be a matrix")
  }
  marker.effects = trained.model$ETA[[1]]$b
  if(is.null(weights)){
    weights = rep(1, length(marker.effects))
  }
  if(length(weights) != length(marker.effects)){
    stop("length of weight vector should correspond to number of markers")
  }
  values = apply(Z, 1, function(genotype){
    sum(genotype * marker.effects * weights)
  })
  return(values)
}

################################################
# DETERMINE FAVOURABLE ALLELES AND FREQUENCIES #
################################################

# get favourable allele at each marker (1 if positive effect, 0 if negative effect)
get.favourable.alleles = function(trained.model){
  if(class(trained.model) != "BGLR"){
    stop("input should be a BGLR object (call gp.train to train model)")
  }
  # get effects
  marker.effects = trained.model$ETA[[1]]$b
  # initialize output
  desired.alleles = rep(0, length(marker.effects))
  # set to 1 for positive effects
  desired.alleles[marker.effects > 0] = 1
  
  return(desired.alleles)
}

get.favourable.allele.frequencies = function(trained.model, Z){
  if(class(trained.model) != "BGLR"){
    stop("input should be a BGLR object (call gp.train to train model)")
  }
  if(missing(Z) || is.null(Z) || !is.matrix(Z)){
    stop("Z should be a matrix")
  }
  # get favourable alleles
  fav.alleles = get.favourable.alleles(trained.model)
  # compute frequencies
  n = nrow(Z)
  m = ncol(Z)
  freqs = rep(NA, m)
  for(i in 1:m){
    freqs[i] = sum(Z[,i] == fav.alleles[i]) / n
  }
  return(freqs)
}

# END OF: code using BGLR



###################
# PLOT ESTIMATION #
###################

# visualize GP estimation by plotting estimated against real values
gp.plot.estimation = function(real.values, estimated.values){
  if(missing(real.values) || is.null(real.values)){
    stop("real values are required")
  }
  if(missing(estimated.values) || is.null(estimated.values)){
    stop("estimated values are required")
  }
  if(length(real.values) != length(estimated.values)){
    stop("number of real and estimated values does not correspond")
  }
  # increase margin
  old.mar = par()$mar
  par(mar=c(5,6,3,3))
  # plot estimated against real values
  x = real.values
  y = estimated.values
  plot(x=x, y=y, xlab = "real values", ylab = "estimated values")
  # draw identity line
  abline(0,1, col="red")
  # fit line
  abline(lm(y~x), col="blue")
  # add correlation
  r = cor(x, y)
  mtext(line=0.5, text = bquote(rho == .(round(r, 5))))
  # reset margin
  par(mar=old.mar)
}
















