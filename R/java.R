
###############################################################
# MAXIMIZE WEIGHTED INDEX: AVERAGE BREEDING VALUE + DIVERSITY #
###############################################################

j.max.index <- function(n, names, values, markers, div.weight, div.measure, fav.alleles = NULL, sec.without.impr = 5){
  
  # convert parameters
  n <- as.integer(n)
  values <- as.numeric(values)
  markers <- .jarray(markers, dispatch = TRUE)
  div.weight <- as.numeric(div.weight)
  div.measure <- .jcast(div.measure, j.core("problems/objectives/Objective"))
  if(is.null(fav.alleles)){
    fav.alleles <- .jnull()
  } else {
    fav.alleles <- as.integer(fav.alleles)
  }
  sec.without.impr <- as.integer(sec.without.impr)

  # call java API
  selected.names <- .jcall(j.api(),
                           "[S", "selectWeighted",
                           n, names, values, markers, fav.alleles,
                           div.weight, div.measure, sec.without.impr)
  
  # sort for reproducible results
  selected.names <- sort(selected.names)
  
  return(selected.names)
  
}

######################
# DIVERSITY MEASURES #
######################

j.HE.all <- function(){
  .jnew(j.gs("obj/HEall"))
}

j.HE.fav <- function(){
  .jnew(j.gs("obj/HEfav"))
}

j.LOG.all <- function(){
  .jnew(j.gs("obj/LOGall"))
}

j.LOG.fav <- function(){
  .jnew(j.gs("obj/LOGfav"))
}

#############
# UTILITIES #
#############

j.version <- function(){
  print(.jcall("java/lang/System", "S", "getProperty", "java.runtime.version"))
}

j.load.jar <- function(){
  .jaddClassPath("java/gs-simulation/bin/gs-simulation.jar")
}

j.api <- function(){
  .jnew(j.gs("api/API"))
}

j.gs <- function(class){
  sprintf("org/jamesframework/gs/simulation/%s", class)
}

j.core <- function(class){
  sprintf("org/jamesframework/core/%s", class)
}

#########################################################
# Start JVM and load JAR after processing all functions #
#########################################################

.jinit()
j.load.jar()


