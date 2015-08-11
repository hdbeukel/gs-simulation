library(rJava)
.jinit()

###############################################################
# MAXIMIZE WEIGHTED INDEX: AVERAGE BREEDING VALUE + DIVERSITY #
###############################################################

j.max.index <- function(n, values, markers, div.weight, div.measure){
  # ...
}

######################
# DIVERSITY MEASURES #
######################

j.MR.ENE <- function(){
  j.new("obj.EntryToNearestEntryDistance")
}

j.HE <- function(){
  j.new("obj.ExpectedProportionOfHeterozygousLoci")
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

j.new <- function(class){
  .jnew(sprintf("org.jamesframework.gs.simulation.%s", class))
}