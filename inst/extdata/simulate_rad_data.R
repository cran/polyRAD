#source('D:/GitHub/polyRAD/R/classes_methods.R')

# make a simple dataset to use in documentation examples

alleles2loc <- c(1,1,2,2,2,3,3,4,4)
alNuc <- c("A", "T", "AC", "GC", "GT", "T", "C", "A", "G")
taxa <- paste("sample", formatC(1:100, format='d', width = 3, flag = '0'), sep = "")
loci <- paste("loc", 1:4, sep = "")
realFreq <- c(0.2, 0.8, 0.1, 0.3, 0.6, 0.95, 0.05, 0.3, 0.7)
locPloidies <- c(2,2,4,2)
realGen <- matrix(0L, nrow = 100, ncol = length(alleles2loc))
depth <- matrix(NA_integer_, nrow = 100, ncol = length(alleles2loc),
                dimnames = list(taxa, paste(loci[alleles2loc], alNuc, sep = "_")))
for(L in 1:max(alleles2loc)){
  thesecol <- which(alleles2loc == L)
  thispld <- locPloidies[L]
  thesefreq <- realFreq[thesecol]
  for(s in 1:length(taxa)){
    for(i in 1:thispld){
      thisAl <- sample(thesecol, size = 1, prob = thesefreq)
      realGen[s,thisAl] <- realGen[s,thisAl] + 1
    }
    totcounts <- round(rnorm(1, mean = 20, sd = 10))
    if(totcounts < 0) totcounts <- 0
    reads <- sample(thesecol, size = totcounts, replace = TRUE, prob = realGen[s, thesecol])
    for(i in thesecol){
      depth[s,i] <- sum(reads == i)
    }
  }
}

realGen[1:10,]
depth[1:10,]

exampleRAD <- RADdata(depth, as.integer(alleles2loc), 
                      data.frame(row.names = loci, Chr = c(1,4,6,9), 
                                 Pos = c(1111020, 33069, 2637920, 5549872)),
                      list(2L, 4L), 0.001, alNuc)
exampleRAD
exampleRAD$depthRatio[1:10,]

save(exampleRAD, file = "exampleRAD.RData")

# make a simple mapping dataset to use in documentation examples
alleles2loc <- rep(1:4, each = 2)
alNuc <- c("A", "G", "A", "C", "G", "T", "T", "A")
taxa <- c("parent1", "parent2",
          paste("progeny", formatC(1:100, format='d', width = 3, flag = '0'), sep = ""))
loci <- paste("loc", 1:4, sep = "")
depth <- matrix(NA_integer_, nrow = length(taxa), ncol = length(alleles2loc),
                dimnames = list(taxa, paste(loci[alleles2loc], alNuc, sep = "_")))
# simulate a diploid BC1 pop
for(L in 1:max(alleles2loc)){
  thesecol <- which(alleles2loc == L)
  for(s in 1:length(taxa)){
    totcounts <- round(rnorm(1, mean = 20, sd = 10))
    if(totcounts < 0) totcounts <- 0
    if(s %in% 1:2){ # make sure parents have data
      while(totcounts == 0){
        totcounts <- as.integer(round(rnorm(1, mean = 20, sd = 10)))
      }
    }
    if(taxa[s] == "parent1"){
      depth[s,thesecol] <- as.integer(c(totcounts, 0L))
    }
    if(taxa[s] == "parent2"){
      depth[s,thesecol] <- as.integer(c(0L, totcounts))
    }
    if(s %in% 3:length(taxa)){
      whichgen <- sample(2, 1)
      if(whichgen == 1){
        depth[s,thesecol] <- as.integer(c(0L, totcounts))
      } else {
        reads <- sample(2, size = totcounts, replace = TRUE)
        depth[s,thesecol] <- as.integer(c(sum(reads == 1), sum(reads == 2)))
      }
    }
  }
}

depth[1:20,]

exampleRAD_mapping <- RADdata(depth, alleles2loc, 
                              data.frame(row.names = loci, Chr = c(1,4,6,9), 
                                         Pos = c(1111020, 33069, 2637920, 5549872)),
                              list(2L), 0.001, alNuc)
exampleRAD_mapping

save(exampleRAD_mapping, file = "exampleRAD_mapping.RData")
