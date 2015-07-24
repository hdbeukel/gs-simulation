#####################
# UTILITY FUNCTIONS #
#####################

write.pop.data <- function(pops, file.prefix, names = 1:length(pops)){
  
  for(i in 1:length(pops)){
    pop <- pops[[i]]
    name <- names[i]
    file.name <- paste("data/", file.prefix, "-", name, sep="")
    # generate marker matrix
    Z <- gp.design.matrix(pop)
    # write marker data
    write.table(Z, file = paste(file.name, "-markers.txt", sep=""))
    # write genetic values, normalized to [-1,1]
    normalized.values <- normalize.genetic.values(pop$geneticValues,
                                                  pop$hypred$genome@add.and.dom.eff$add)
    write.table(normalized.values, file = paste(file.name, "-values.txt", sep=""), col.names = F)
  }
  
}

###########################
# Generate data for Chloë #
###########################

# params
num.qtl <- 100
h2 <- 0.5
pop.size <- 200
num.select <- 20
num.seasons <- 10
num.simulations <- 3

for(s in 1:num.simulations){
  # load processed data
  load("data/ProcessedData.RData")
  
  # generate series of data sets produced under genomic selection
  seasons <- GS(founders, heritability = h2, F1.size = pop.size,
                num.select = num.select, num.seasons = num.seasons,
                num.QTL = num.qtl)
  
  # select 3 population snapshots at different times during simulation
  selected.pops <- as.list(rep(NA, 3))
  
  selected.pops[[1]] <- seasons[[2]]$select$pop.in
  selected.pops[[2]] <- seasons[[num.seasons/2]]$select$pop.in
  selected.pops[[3]] <- seasons[[num.seasons]]$select$pop.in
  
  # write population data
  write.pop.data(selected.pops,
                 file.prefix = sprintf("snapshot-pop-%d", s),
                 names = c("early", "medium", "late"))
}












