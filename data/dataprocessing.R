# load libraries
library(synbreed) # version 0.10-5
# set path to BEAGLE used by synbreed
beaglePath <- "beaglejar"
# tolerance for floating point comparison
floating.point.tol = 1e-12

# compute ratio of inferred or heterozygous alleles per genotype (across all markers)
get.ratio.of.inferred.or.heterozygous.alleles.per.genotype = function(markers){
  num.genotypes = nrow(markers)
  num.markers = ncol(markers)
  result = rep(NA, num.genotypes)
  for(i in 1:num.genotypes){
    genotype = markers[i,]
    inferred.or.heterozygous = abs(genotype - 1) > floating.point.tol & abs(genotype + 1) > floating.point.tol
    ratio = sum(inferred.or.heterozygous)/num.markers
    result[i] = ratio
  }
  return(result)
}

# compute ratio of inferred or heterozygous alleles per marker (across all genotypes)
get.ratio.of.inferred.or.heterozygous.alleles.per.marker = function(markers){
  num.genotypes = nrow(markers)
  num.markers = ncol(markers)
  result = rep(NA, num.markers)
  for(i in 1:num.markers){
    marker = markers[,i]
    inferred.or.heterozygous = abs(marker - 1) > floating.point.tol & abs(marker + 1) > floating.point.tol
    ratio = sum(inferred.or.heterozygous)/num.genotypes
    result[i] = ratio
  }
  return(result)
}

# groups markers with same position
split.markers.on.position = function(map){
  # go through all markers
  cur.marker.index = 1
  num.markers = nrow(map)
  # working list
  working.list = vector('list', num.markers)
  working.list.index = 1
  while(cur.marker.index <= num.markers){
    # get position of current marker
    cur.pos = map[cur.marker.index, 2]
    # increase next pos index while markers on same position
    next.pos.index = cur.marker.index + 1
    next.pos = map[next.pos.index, 2]
    while(!is.na(next.pos) && isTRUE(all.equal(cur.pos, next.pos))){
      next.pos.index = next.pos.index + 1
      next.pos = map[next.pos.index, 2]
    }
    # store vector with all marker names on same position
    same.position = rownames(map[cur.marker.index:(next.pos.index-1),])
    working.list[[working.list.index]] = same.position
    working.list.index = working.list.index + 1
    # continue with next marker
    cur.marker.index = next.pos.index
  }
  # copy filled in part of working list to result list
  result.list = working.list[1:(working.list.index-1)]
  return(result.list)
}

# compute IBD matrix within groups of markers at same position
compute.coinciding.markers.ibd = function(markers, map){
  coinciding.markers = split.markers.on.position(map)
  num.groups = length(coinciding.markers)
  result = vector('list', num.groups)
  for(i in 1:num.groups){
    # get coinciding marker names
    marker.names = coinciding.markers[[i]]
    num.markers = length(marker.names)
    # compute IBD  matrix
    ibd.matrix = matrix(NA, num.markers, num.markers)
    for(m1 in 1:num.markers){
      for(m2 in m1:num.markers){
        m1.name = marker.names[m1]
        m2.name = marker.names[m2]
        m1.alleles = markers[,m1.name]
        m2.alleles = markers[,m2.name]
        ibd = sum(m1.alleles == m2.alleles)/nrow(markers)
        ibd.matrix[m1, m2] = ibd
        ibd.matrix[m2, m1] = ibd
      }
    }
    result[[i]] = ibd.matrix
  }
  return(result)
}

# compute list of all possible subset subscripts (T/F) of size n
all.subsets = function(n){
  result = matrix(NA, 2^n, n)
  result[1,1] = TRUE
  result[2,1] = FALSE
  i=2
  while(i <= n){
    part1.start = 1
    part1.stop = 2^(i-1)
    part2.start = 2^(i-1) + 1
    part2.stop = 2^i
    result[part1.start:part1.stop, i] = TRUE 
    result[part2.start:part2.stop, 1:i] = result[part1.start:part1.stop, 1:i]
    result[part2.start:part2.stop, i] = FALSE
    i = i+1
  }
  result.list = vector('list', 2^n)
  for(i in 1:(2^n)){
    result.list[[i]] = result[i,]
  }
  return(result.list)
}

# get markers to filter: same position and IBD across all lines (allowing for some error)
get.discarded.coinciding.markers = function(markers, map, allowed.error){
  coinciding.markers = split.markers.on.position(map)
  discarded.marker.names = vector()
  # go through groups of markers at same position
  num.groups = length(coinciding.markers)
  for(i in 1:num.groups){
    # consider all possible subsets
    group = coinciding.markers[[i]]
    possible.subsets = all.subsets(length(group))
    largest.valid.subset = NA
    largest.valid.subset.size = 0
    for(subset in possible.subsets){
      # make selection
      retain = group[subset]
      # verify: no IBD pairs (taking into account allowed genotyping error)
      if(!contains.ibd.pair(markers, retain, allowed.error)){
        # update if larger valid subset
        if(length(retain) > largest.valid.subset.size){
          largest.valid.subset.size = length(retain)
          largest.valid.subset = subset
        }
      }
    }
    # store names of filtered markers
    discarded.marker.names = c(discarded.marker.names, group[!largest.valid.subset])
  }
  return(discarded.marker.names)
}

# check wether two markers are IBD across all lines with allowed genotype error
is.ibd = function(markers, marker1.name, marker2.name, allowed.error){
  # compute IBD rate
  ibd.rate = sum(markers[,marker1.name] == markers[,marker2.name])/nrow(markers)
  # check: discordance within allowed error?
  return((1.0 - ibd.rate) <= allowed.error)
}

# check wether the given series of markers contains any pair of IBD (allowing error)
contains.ibd.pair = function(markers, marker.names, allowed.error){
  if(length(marker.names) < 2){
    return(FALSE)
  }
  for(m1 in 1:(length(marker.names)-1)){
    for(m2 in (m1+1):length(marker.names)){
      m1.name = marker.names[m1]
      m2.name = marker.names[m2]
      if(is.ibd(markers, m1.name, m2.name, allowed.error)){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

# arbitrarily spread out markers on same position with a given distance
spread.coinciding.markers = function(map, dist){
  adjusted.map = map[,1:2]
  coinciding.markers = split.markers.on.position(map)
  # go through groups of markers at same position
  num.groups = length(coinciding.markers)
  for(i in 1:num.groups){
    group = coinciding.markers[[i]]
    group.size = length(group)
    if(group.size %% 2 == 0){
      # even number of markers
      n = group.size - 2
      deltas = seq(dist * -n/2 - dist/2, dist * n/2 + dist/2, dist)
    } else {
      # odd number of markers
      n = group.size - 1
      deltas = seq(dist * -n/2, dist * n/2, dist)
    }
    # adjust positions
    adjusted.map[group,2] = adjusted.map[group,2] + deltas
  }
  return(adjusted.map)
}

# compute minor allele frequencies
# known encodings:
#  "dh"   ->  0,1 (homozygous) - default
#  "012"  ->  0,2 (homozygous), 1 (heterozygous)
#  "-101" -> -1,1 (homozygous), 0 (heterozygous)
compute.minor.allele.frequencies <- function(markers, encoding = c("dh", "-101", "012")){
  # check: known encoding
  encoding <- match.arg(encoding)
  # convert to 012 encoding
  if(!is.na(pmatch(encoding, "-101"))){
    markers <- markers + 1
  } else if (!is.na(pmatch(encoding, "dh"))){
    markers <- markers*2
  }
  # compute MAF
  num.genotypes <- nrow(markers)
  maf <- apply(markers, 2, function(alleles){
    freq.1 <- sum(alleles)/(2*num.genotypes)
    freq.0 <- 1-freq.1
    return(min(freq.1, freq.0))
  })
  return(maf)
}
# get minor alleles (0/1)
# known encodings:
#  "dh"   ->  0,1 (homozygous) - default
#  "012"  ->  0,2 (homozygous), 1 (heterozygous)
get.minor.alleles = function(markers, encoding = c("dh", "012")){
  # check: known encoding
  encoding <- match.arg(encoding)
  # convert to 012 encoding
  if(!is.na(pmatch(encoding, "dh"))){
    markers <- markers*2
  }
  # compute minor allel (0/1)
  num.genotypes <- nrow(markers)
  ma <- apply(markers, 2, function(alleles){
    freq.1 <- sum(alleles)/(2*num.genotypes)
    if(freq.1 <= 0.5){
      return(1)
    } else {
      return(0)
    }
  })
  return(ma)
}

# get genotype names
get.geno.names <- function(pop){
  if(!is.null(pop$dh)){
    geno <- pop$dh
  } else {
    geno <- pop$geno
  }
  return(rownames(geno))
}

# get marker names
get.marker.names <- function(pop){
  if(!is.null(pop$dh)){
    geno <- pop$dh
  } else {
    geno <- pop$geno
  }
  return(colnames(geno))
}

###################
# PREPROCESS DATA #
###################

preprocess.data = function(markers, map, max.inferred.or.heterozygous = 0.1, max.genotype.error = 0.03, spread.dist = 0.1, min.minor.allele.frequency = 0.01){
  
  colnames(map) = tolower(colnames(map))
  
  # 1: remove markers with too much inferred or heterozygous alleles across all genotypes
  marker.names = colnames(markers)
  discarded.marker.names = marker.names[get.ratio.of.inferred.or.heterozygous.alleles.per.marker(markers) > max.inferred.or.heterozygous]
  retained.marker.names = setdiff(marker.names, discarded.marker.names)
  markers = markers[,retained.marker.names]
  map = map[retained.marker.names,]
  message("Discarded ", length(discarded.marker.names),
          " of ", length(marker.names),
          " markers with more than ", max.inferred.or.heterozygous * 100, "%",
          " heterozygous or inferred alleles (retained ", length(retained.marker.names) , " markers)")
  
  # 2: remove genotypes with too much inferred or heterozygous alleles across all genotypes
  genotype.names = rownames(markers)
  discarded.genotype.names = genotype.names[get.ratio.of.inferred.or.heterozygous.alleles.per.genotype(markers) > max.inferred.or.heterozygous]
  retained.genotype.names = setdiff(genotype.names, discarded.genotype.names)
  markers = markers[retained.genotype.names,]
  message("Discarded ", length(discarded.genotype.names),
          " of ", length(genotype.names),
          " genotypes with more than ", max.inferred.or.heterozygous * 100, "%",
          " heterozygous or inferred alleles (retained ", length(retained.genotype.names) , " genotypes)")
  
  # 3: round alleles to -1, 0 or 1
  markers = round(markers, digits = 0)
  message("Rounded allele values to -1, 0 or 1")
  
  # 4: remove markers with too small minor allele frequency
  marker.names = colnames(markers)
  discarded.marker.names = marker.names[compute.minor.allele.frequencies(markers, "-101") < min.minor.allele.frequency]
  retained.marker.names = setdiff(marker.names, discarded.marker.names)
  markers = markers[, retained.marker.names]
  map = map[retained.marker.names,]
  message("Discarded ", length(discarded.marker.names),
          " of ", length(marker.names),
          " markers with less than ", 100*min.minor.allele.frequency, "% minor allele frequency",
          " (retained ", length(retained.marker.names) , " markers)")
  
  # 5: remove IBD markers at same position (with allowance of some genotype error)
  marker.names = colnames(markers)
  discarded.marker.names = get.discarded.coinciding.markers(markers, map, max.genotype.error)
  retained.marker.names = setdiff(marker.names, discarded.marker.names)
  markers = markers[,retained.marker.names]
  map = map[retained.marker.names,]
  message("Discarded ", length(discarded.marker.names),
          " of ", length(marker.names),
          " markers: IBD at same position with allowance for ", 100*max.genotype.error, "% genotype error",
          " (retained ", length(retained.marker.names) , " markers)")
  
  # 6: spread remaining markers at same position
  map = spread.coinciding.markers(map, spread.dist)
  message("Spread remaining markers mapped to same position at ", spread.dist, " cM intervals in arbitrary order")
  
  # 7: sort markers after spreading
  unordered.marker.names = colnames(markers) 
  order = with(map, order(chr, pos))
  map = map[order,]
  markers = markers[,order]
  ordered.marker.names = colnames(markers)
  num.repositioned = sum(unordered.marker.names != ordered.marker.names)
  message("Sorted markers after spreading (repositioned ", num.repositioned, " markers)")
  
  # 8: shift positions so that first marker of each chromosome has positive positin
  chromosomes = unique(map[,1])
  shift = sapply(chromosomes, function(c){
    min(map[map$chr == c,2])
  })
  shift[shift >= 0] = 0
  for(c in chromosomes){
    map[map$chr == c,2] = map[map$chr == c,2] - shift[c]
  }
  message("Shifted markers so that first marker of every chromosome is mapped to a positive position")
  
  # 9: compute chromosome lengths, marker counts and cumulative marker positions
  chromosomes = unique(map[,1])
  num.chromosomes = length(chromosomes)
  chr.lengths = sapply(chromosomes, function(i){
    max(map[map[,1] == i,2])
  })
  chr.offset = cumsum(c(0, chr.lengths[1:(num.chromosomes-1)]))
  chrom.num.markers = sapply(chromosomes, function(i){
    sum(map[,1] == i)
  })
  rep.chr.offset = rep(chr.offset, chrom.num.markers)
  cumu.pos = map[,2] + rep.chr.offset
  
  map$cumuPos = cumu.pos
  
  chr.cumu.lengths = sapply(chromosomes, function(i){
    max(map[map[,1] == i,3])
  })
    
  # finally: pack filtered data
  message("--------------------------------------------------------------------------------------------------------------\n",
          "PROCESSED DATA SET: ", nrow(markers), " genotypes, ", ncol(markers), " markers\n",
          "--------------------------------------------------------------------------------------------------------------")
  
  processed.data = list(markers = markers, map = map, numGenotypes = nrow(markers),
                        numMarkers = ncol(markers), numChroms = num.chromosomes,
                        chrLengths = chr.lengths, chrCumuLengths = chr.cumu.lengths,
                        chrNumMarkers = chrom.num.markers)
  return(processed.data)
  
}

############################
# CONVERT TO SYNBREED DATA #
############################

convert.to.synbreed.data = function(data){
  # create gpData
  synbreed.data = create.gpData(geno = data$markers, map = data$map[,c("chr", "pos")])
  # recode alleles from -1,0,1 to 0,1,2 (with 1 = minor allele)
  synbreed.data$geno = synbreed.data$geno + 1
  synbreed.data = codeGeno(synbreed.data, label.heter = '1')
  # return synbreed data
  return(synbreed.data)
}

#####################################################
# PHASE GENOTYPES USING BEAGLE (CALLED BY SYNBREED) #
#####################################################

phase.genotypes = function(synbreed.data){
  # run phasing through beagle
  message("Running BEAGLE ...")
  suppressWarnings(invisible(capture.output(codeGeno(synbreed.data, impute=TRUE, impute.type='beagle', label.heter='1'))))
  message("DONE")
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
  phased.genotypes = list(hap1 = hap1.markers, hap2 = hap2.markers)
  return(phased.genotypes)
}

prepare.hypred = function(packed.data){
  # compute chromosome bounds
  chromEnds = cumsum(packed.data$chrNumMarkers)
  chromStarts = c(1, chromEnds[1:(packed.data$numChroms-1)] + 1)
  chromBounds = cbind(chromStarts, chromEnds)
  # number of markers per chromosome (required to be fixed for hypred)
  markersPerChrom = max(packed.data$chrNumMarkers)
  # number of dummy markers per chromosome to meet marker count
  chrNumDummies = markersPerChrom - packed.data$chrNumMarkers
  
  # attach hypred data
  hypred = list(chromBounds = chromBounds, markersPerChrom = markersPerChrom,
                chrNumDummies = chrNumDummies, dummiesAdded = FALSE)
  packed.data$hypred = hypred
  
  # add dummy markers
  packed.data = add.dummies(packed.data)
  
  return(packed.data)
}

#################################################
# PACK DATA (PREPROCESS, PHASE, PREPARE HYPRED) #
#################################################

pack.data = function(data, min.maf = 0.01, max.genotype.error = 0.03){
    
  # 1) preprocess data from Jannink
  data = preprocess.data(data$markers, data$map, min.minor.allele.frequency = min.maf, max.genotype.error = max.genotype.error)
  
  # 2) convert to synbreed data
  synbreed.data = convert.to.synbreed.data(data)
    
  # 3) phase genotypes
  phased.genotypes = phase.genotypes(synbreed.data)
  
  # 4) pack data
  packed.data = list(geno = synbreed.data$geno, phasedGeno = phased.genotypes, map = data$map,
                     numGenotypes = data$numGenotypes, numMarkers = data$numMarkers,
                     numChroms = data$numChroms, chrLengths = data$chrLengths,
                     chrCumuLengths = data$chrCumuLengths, chrNumMarkers = data$chrNumMarkers)
    
  # 5) prepare hypred (metadata for dummies, genome definition, ...)
  packed.data = prepare.hypred(packed.data)
    
  return(packed.data)
  
}




















