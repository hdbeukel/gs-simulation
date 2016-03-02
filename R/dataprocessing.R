
# check marker encoding; known encodings:
#  1)  "012":  0, 1, 2 (non DH)
#  2)   "01":  0, 1    (DH)
#  3)   "dh": alias for "01"
# if the marker matrix does not correspond to the specified
# encoding, an error is generated; else, the method returns NULL
check.marker.encoding <- function(markers, encoding = c("012", "01", "dh")){
  # retrieve encoding
  encoding <- match.arg(encoding)
  if(encoding == "dh"){
    encoding = "01"
  }
  # set expected values
  expected.values <- switch(encoding,
                             "012" = c(0, 1, 2),
                              "01" = c(0, 1))
  # check for unexpected values
  values <- levels(as.factor(markers))
  unexpected.values <- setdiff(values, expected.values)
  if(length(unexpected.values) > 0){
    stop(sprintf("values in marker matrix do not correspond to encoding '%s'; unexpected value(s): %s",
                 encoding, paste(unexpected.values, collapse = " ")))
  }
}

# convert marker data to "012" format; known input encodings:
#  2)  "012":  0, 1, 2 (non DH)
#  3)   "01":  0, 1    (DH)
#  4)   "dh": alias for "01"
convert.marker.encoding <- function(markers, input.encoding = c("012", "01", "dh")){
  # retrieve input encoding
  input.encoding <- match.arg(input.encoding)
  if(input.encoding == "dh"){
    input.encoding = "01"
  }
  # check encoding
  check.marker.encoding(markers, input.encoding)
  # convert to 0,1,2
  converted <- switch(input.encoding,
                      "012" = markers,
                      "01" = markers * 2)
  return(converted)
}

# compute frequency of missing values per marker/individual
missing.value.frequencies <- function(markers, per = c("marker", "individual")){
  # process arguments
  per <- match.arg(per)
  # infer missing value freqs
  freqs <- apply(markers, ifelse(per == "individual", 1, 2), function(values){
    sum(is.na(values))/length(values)
  })
  return(freqs)
}

# groups markers with same position (names)
markers.per.position <- function(map){
  # go through all markers
  cur.marker.index <- 1
  num.markers <- nrow(map)
  # working list
  working.list <- as.list(rep(NA, num.markers))
  working.list.index <- 1
  while(cur.marker.index <= num.markers){
    # get position of current marker
    cur.pos <- map[cur.marker.index, 2]
    # increase next pos index while markers on same position
    next.pos.index <- cur.marker.index + 1
    next.pos <- map[next.pos.index, 2]
    while(!is.na(next.pos) && isTRUE(all.equal(cur.pos, next.pos))){
      next.pos.index <- next.pos.index + 1
      next.pos <- map[next.pos.index, 2]
    }
    # store vector with all marker names on same position
    same.position <- rownames(map[cur.marker.index:(next.pos.index-1), ])
    working.list[[working.list.index]] <- same.position
    working.list.index <- working.list.index + 1
    # continue with next marker
    cur.marker.index <- next.pos.index
  }
  # copy filled in part of working list to result list
  result.list <- working.list[1:(working.list.index-1)]
  return(result.list)
}

# get names of redundant markers: same position + nearly identical
# for all individuals (allowing for some small genotyping error)
redundant.markers <- function(markers, map, allowed.error){
  groups <- markers.per.position(map)
  # skip markers with unique position
  group.sizes <- sapply(groups, length)
  groups <- groups[group.sizes > 1]
  # (heuristically) filter groups of near identical markers at same position
  num.geno <- nrow(markers)
  redundant.marker.names <- unlist(lapply(groups, function(group){
    # compute pairwise redundancy matrix
    redundant <- outer(group, group, function(m1, m2){
      pairs <- rbind(m1, m2)
      red <- apply(pairs, 2, function(pair){
        id.ratio <- sum(markers[, pair[1]] == markers[, pair[2]])/num.geno
        return(id.ratio >= (1.0 - allowed.error))
      })
    })
    # select non-redundant markers to retain (greedy heuristic)
    retain <- c()
    for(m in 1:length(group)){
      # check if marker is redundant given the current selection
      red <- redundant[m, retain]
      if(!any(red)){
        retain <- c(retain, m)
      }
    }
    # return names of redundant markers
    return(group[-retain])
  }))
  return(redundant.marker.names)
}

# check wether two markers are nearly identical across all
# individuals, allowing for some small genotyping error
identical.markers <- function(markers, marker1.name, marker2.name, allowed.error){
  # compute proportion of identical values
  identical <- sum(markers[, marker1.name] == markers[, marker2.name])/nrow(markers)
  # check: nealy identical?
  return(identical >= (1.0 - allowed.error))
}

# check wether the given group of markers contains any pair of nearly
# identical markers, allowing for some small genotyping error
contains.identical.pair <- function(markers, marker.names, allowed.error){
  if(length(marker.names) < 2){
    return(FALSE)
  }
  identical <- outer(marker.names, marker.names, function(m1, m2){
    pairs <- rbind(m1, m2)
    res <- apply(pairs, 2, function(pair){
      identical.markers(markers, pair[1], pair[2], allowed.error)
    })
    return(res)
  })
  return(any(identical[lower.tri(identical)]))
}

# arbitrarily spread out markers on same position with a given distance
spread.coinciding.markers <- function(map, max.shift, min.shift = 0.01){
  adjusted.map <- map[,1:2]
  coinciding.markers <- markers.per.position(map)
  # go through groups of markers at same position
  num.groups <- length(coinciding.markers)
  for(i in 1:num.groups){
    # adjust positions
    group <- coinciding.markers[[i]]
    deltas <- runif(length(group), min = min.shift, max = max.shift)
    deltas <- sample(c(-1,1), replace = T, size = length(group)) * deltas
    adjusted.map[group,2] <- adjusted.map[group,2] + deltas
  }
  return(adjusted.map)
}

# compute minor allele frequencies
maf <- function(markers, encoding = c("012", "01", "dh")){
  # convert to 012 encoding (also checks encoding)
  markers <- convert.marker.encoding(markers, input.encoding = encoding)
  # compute MAF
  num.genotypes <- nrow(markers)
  maf <- apply(markers, 2, function(snp){
    freq.1 <- sum(snp)/(2*num.genotypes)
    freq.0 <- 1-freq.1
    return(min(freq.1, freq.0))
  })
  return(maf)
}
# get minor alleles (0/1)
minor.alleles <- function(markers, encoding = c("012", "01", "dh")){
  # convert to 012 encoding (also checks encoding)
  markers <- convert.marker.encoding(markers, input.encoding = encoding)
  # compute minor alleles (0/1)
  num.genotypes <- nrow(markers)
  ma <- apply(markers, 2, function(snp){
    freq.1 <- sum(snp)/(2*num.genotypes)
    return(ifelse(freq.1 <= 0.5, 1, 0))
  })
  return(ma)
}

# get genotype names
geno.names <- function(pop){
  if(!is.null(pop$dh)){
    geno <- pop$dh
  } else {
    geno <- pop$geno
  }
  return(rownames(geno))
}

# get marker names
marker.names <- function(pop){
  if(!is.null(pop$dh)){
    geno <- pop$dh
  } else {
    geno <- pop$geno
  }
  return(colnames(geno))
}

##############
# CLEAN DATA #
##############

# input markers encoded as AA/AB/BB ("-" indicates missing data)
# marker matrix: ind x markers (named)
# map: one row per marker (named); 2 columns: chromosome, position (unnamed)
clean.data <- function(markers, map,
                       max.missing = 0.20,
                       max.genotype.error = 0.00,
                       max.spread.dist = 0.1,
                       min.maf = 0.00){
  
  n <- nrow(markers)
  m <- ncol(markers)
  message("Clean raw data (", n, " genotypes, ", m, " markers)")
  
  # 1: filter map
  
  # set map column names
  colnames(map) <- c("chr", "pos")
  # extract names of observed markers
  marker.names <- colnames(markers)
  # discard markers with unknown position
  unknown.chromosome <- map$chr == "UNK"
  known.chromosome <- rownames(map[!unknown.chromosome,])
  known.position <- intersect(marker.names, known.chromosome)
  map <- map[known.position, ]
  markers <- markers[, known.position]
  message("Discarded ", m - length(known.position),
          " of ", m,
          " markers with unknown position (retained ", length(known.position), " markers)")
  # convert chromosome names to numeric
  map$chr <- as.numeric(map$chr)
  # order markers according to
  # (1) chromosome
  # (2) position on chromosome
  marker.order <- order(map$chr, map$pos)
  map <- map[marker.order, ]
  markers <- markers[, marker.order]
  # mark missing values as NA
  markers[markers == '-'] <- NA
  
  # 2: remove markers with too much missing values
  marker.names <- colnames(markers)
  discarded.marker.names <- marker.names[missing.value.frequencies(markers, per = "marker") > max.missing]
  retained.marker.names <- setdiff(marker.names, discarded.marker.names)
  markers <- markers[, retained.marker.names]
  map <- map[retained.marker.names, ]
  message("Discarded ", length(discarded.marker.names),
          " of ", length(marker.names),
          " markers with more than ", max.missing * 100, "%",
          " missing values (retained ", length(retained.marker.names) , " markers)")
  
  # 3: remove genotypes with too much missing values
  geno.names <- rownames(markers)
  discarded.genotype.names <- geno.names[missing.value.frequencies(markers, per = "ind") > max.missing]
  retained.genotype.names <- setdiff(geno.names, discarded.genotype.names)
  markers <- markers[retained.genotype.names, ]
  message("Discarded ", length(discarded.genotype.names),
          " of ", length(geno.names),
          " genotypes with more than ", max.missing * 100, "%",
          " missing values (retained ", length(retained.genotype.names) , " genotypes)")
  
  # 4: recode marker data to 0,1,2 and impute missing values
  message("Recode to 0,1,2 (synbreed format) and impute missing values (BEAGLE)")
  synbreed.data <- create.gpData(geno = markers, map = map)
  synbreed.data <- codeGeno(synbreed.data, label.heter = 'AB', impute = TRUE, impute.type = "beagle")
  markers <- synbreed.data$geno
  
  # 5: remove markers with too low MAF
  marker.names <- colnames(markers)
  discarded.marker.names <- marker.names[maf(markers, "012") < min.maf]
  retained.marker.names <- setdiff(marker.names, discarded.marker.names)
  markers <- markers[, retained.marker.names]
  map <- map[retained.marker.names, ]
  message("Discarded ", length(discarded.marker.names),
          " of ", length(marker.names),
          " markers with less than ", 100*min.maf, "% minor allele frequency",
          " (retained ", length(retained.marker.names) , " markers)")
  
  # 6: remove redundant, i.e. nearly identical markers at same position (with allowance of small genotype error)
  marker.names <- colnames(markers)
  discarded.marker.names <- redundant.markers(markers, map, max.genotype.error)
  retained.marker.names <- setdiff(marker.names, discarded.marker.names)
  markers <- markers[, retained.marker.names]
  map <- map[retained.marker.names, ]
  message("Discarded ", length(discarded.marker.names),
          " of ", length(marker.names),
          " markers: identical at same position with allowance for ", 100*max.genotype.error, "% genotype error",
          " (retained ", length(retained.marker.names) , " markers)")
  
  # 7: spread remaining markers at same position
  map <- spread.coinciding.markers(map, max.spread.dist)
  message("Randomly shift remaining markers mapped to same position with  <= ", max.spread.dist, " cM")
  
  # 8: sort markers after spreading
  unordered.marker.names <- colnames(markers) 
  order <- with(map, order(chr, pos))
  map <- map[order, ]
  markers <- markers[, order]
  ordered.marker.names <- colnames(markers)
  num.repositioned <- sum(unordered.marker.names != ordered.marker.names)
  message("Sorted markers after spreading (repositioned ", num.repositioned, " markers)")
  
  # 9: shift positions so that first marker of each chromosome has position >= 0
  chromosomes <- unique(map[,1])
  shift <- sapply(chromosomes, function(c){
    -min(map[map$chr == c,2])
  })
  shift[shift <= 0] <- 0
  map$pos <- apply(map, 1, function(marker){
    marker[2] + shift[marker[1]]
  })
  message("Shifted markers so that first marker of every chromosome is mapped to a position >= 0")
  
  # 10: pack data
  cleaned.data <- list(markers = markers, map = map)
    
  message("--------------------------------------------------------------------------------------------------------------\n",
          "CLEANED DATA SET: ", nrow(markers), " genotypes, ", ncol(markers), " markers\n",
          "--------------------------------------------------------------------------------------------------------------")
  
  return(cleaned.data)
  
}

##################
# PHASE (BEAGLE) #
##################

# expected marker encoding: 0,1,2
phase <- function(markers, map){
  # run phasing through beagle
  synbreed.data <- create.gpData(geno = markers, map = map)
  suppressWarnings(invisible(capture.output(codeGeno(synbreed.data, impute=TRUE, impute.type='beagle', label.heter='1'))))
  # parse phased output files (per chromosome)
  message("Parsing phased output files ...")
  num.genotypes = nrow(synbreed.data$geno)
  num.markers = ncol(synbreed.data$geno)
  hap1.markers = matrix(NA, num.genotypes, num.markers)
  hap2.markers = matrix(NA, num.genotypes, num.markers)
  colnames(hap1.markers) = colnames(synbreed.data$geno)
  colnames(hap2.markers) = colnames(synbreed.data$geno)
  rownames(hap1.markers) = rownames(synbreed.data$geno)
  rownames(hap2.markers) = rownames(synbreed.data$geno)
  chroms = unique(synbreed.data$map$"chr")
  for(chrom in chroms){
    message(paste("|- Chromosome", chrom))
    compressed.file = paste("beagle/chr", chrom, "input.bgl.phased.gz", sep="")
    # uncompress
    system(paste("gzip -d -f", compressed.file))
    uncompressed.file = paste("beagle/chr", chrom, "input.bgl.phased", sep="")
    # remove first column
    system(paste("cut -d' ' -f2-", uncompressed.file, "> tmp.file"))
    system(paste("mv tmp.file", uncompressed.file))
    # parse
    chrom.data = t(as.matrix(read.table(uncompressed.file, header=T, check.names=F, row.names=1)))
    chrom.data[chrom.data == 'A'] = 0
    chrom.data[chrom.data == 'B'] = 1
    mode(chrom.data) = "numeric"
    odd.rows = seq(1, nrow(chrom.data), 2)
    even.rows = seq(2, nrow(chrom.data), 2)
    # store parsed data
    chrom.marker.names = colnames(chrom.data)
    hap1.markers[,chrom.marker.names] = chrom.data[odd.rows,]
    hap2.markers[,chrom.marker.names] = chrom.data[even.rows,]
  }
  message("DONE")
  # combine and return results
  haplotypes = list(hap1 = hap1.markers, hap2 = hap2.markers)
  return(haplotypes)
}

################################
# PREPARE DATA FOR SIMULATIONS #
################################

# input:
#  - map: data frame with one row per marker (named) and two columns (named):
#     1) "chr": chromosome number (integer)
#     2) "pos": position in cM
#  - hap1, hap2: 0/1 haplotype matrices (ind x markers; named)
#  - dh: TRUE for doubled haploids; ignores hap2
prepare <- function(map, hap1, hap2, dh = FALSE){
  
  # 1: process map
  
  # infer number of chromosomes
  chromosomes <- unique(map$chr)
  num.chromosomes <- length(chromosomes)
  # compute chromosome lengths
  chr.lengths <- sapply(chromosomes, function(c){
    max(map[map$chr == c, 2])
  })
  # compute number of markers per chromosome
  chr.num.markers <- sapply(chromosomes, function(c){
    sum(map$c == c)
  })
  # compute chromosome bounds
  chrom.ends <- cumsum(chr.num.markers)
  chrom.starts <- c(1, chrom.ends[1:(num.chromosomes-1)] + 1)
  chrom.bounds <- cbind(chrom.starts, chrom.ends)
  
  # 2: pack data
    
  if(dh){
    prepared.data <- list(isDH = TRUE, dh = hap1)
  } else {
    geno <- hap1 + hap2
    haplotypes <- list(hap1 = hap1, hap2 = hap2)
    prepared.data <- list(isDH = FALSE, geno = geno, phasedGeno = haplotypes)
  }
  prepared.data <- c(prepared.data, list(map = map,
                                         numGenotypes = nrow(geno),
                                         numMarkers = ncol(geno),
                                         numChroms = num.chromosomes,
                                         chrLengths = chr.lengths,
                                         chrNumMarkers = chr.num.markers))
  
  message("Packed data")
  
  # 3: prepare hypred
  
  # infer number of markers per chromosome in hypred genome
  # (hypred requires fixed marker count per chromosome)
  markers.per.chrom <- max(chr.num.markers)
  # number of dummy markers to add per chromosome to work around fixed marker count
  chr.num.dummies <- markers.per.chrom - chr.num.markers
  
  # attach hypred metadata
  hypred <- list(chromBounds = chrom.bounds, markersPerChrom = markers.per.chrom,
                 chrNumDummies = chr.num.dummies, dummiesAdded = FALSE)
  prepared.data$hypred <- hypred
  
  # add dummy markers
  prepared.data <- add.dummies(prepared.data)
  
  message("Prepared hypred (added dummy markers)")
  
  return(prepared.data)
  
}

##############################
# CREATE FOUNDER POPULATION  #
##############################

create.founders <- function(){
  
  # load raw data
  message("# Load raw data ...")
  markers <- as.matrix(read.table("data/oregon-markers.txt", check.names = FALSE))
  map <- read.table("data/barley-map.txt", row.names = 1)
  
  # clean
  message("# Clean data ...")
  cleaned <- clean.data(markers, map)
  
  # phase
  message("# Phase genotypes ...")
  haplotypes <- phase(cleaned$markers, cleaned$map)
  
  # prepare for simulations
  message("# Prepare for simulations ...")
  founders <- prepare(cleaned$map, haplotypes$hap1, haplotypes$hap2)
  
  # store and return founders
  saveRDS(founders, "data/founders.RDS")
  return(founders)
  
}













