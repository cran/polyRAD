## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 6, fig.height = 6)

## -----------------------------------------------------------------------------
library(polyRAD)
maphmcfile <- system.file("extdata", "ClareMap_HapMap.hmc.txt", 
                          package = "polyRAD")
maphmcfile

mydata <- readHMC(maphmcfile,
                  possiblePloidies = list(2, c(2, 2)))
mydata

## -----------------------------------------------------------------------------
GetTaxa(mydata)[c(1:10,293:299)]

## -----------------------------------------------------------------------------
mydata <- SetDonorParent(mydata, "Kaskade-Justin")
mydata <- SetRecurrentParent(mydata, "Zebrinus-Justin")

## -----------------------------------------------------------------------------
alignfile <- system.file("extdata", "ClareMap_alignments.csv", 
                         package = "polyRAD")

aligndata <- read.csv(alignfile, row.names = 1)
head(aligndata)

mydata$locTable$Chr <- aligndata[GetLoci(mydata), 1]
mydata$locTable$Pos <- aligndata[GetLoci(mydata), 2]
head(mydata$locTable)

## ----eval = FALSE-------------------------------------------------------------
#  mydata <- AddPCA(mydata)

## -----------------------------------------------------------------------------
load(system.file("extdata", "examplePCA.RData", package = "polyRAD"))
mydata$PCA <- examplePCA

## -----------------------------------------------------------------------------
plot(mydata)

## -----------------------------------------------------------------------------
realprogeny <- GetTaxa(mydata)[mydata$PCA[,"PC1"] > -10 &
                                 mydata$PCA[,"PC1"] < 10]
# eliminate the one doubled haploid line in this group
realprogeny <- realprogeny[!realprogeny %in% c("IGR-2011-001", "p196-150A-c",
                                               "p877-348-b")]
# also retain parents
keeptaxa <- c(realprogeny, GetDonorParent(mydata), GetRecurrentParent(mydata))

mydata <- SubsetByTaxon(mydata, taxa = keeptaxa)
plot(mydata)

## ----message = FALSE----------------------------------------------------------
myhindhe <- HindHeMapping(mydata, ploidy = 2L)
hist(colMeans(myhindhe, na.rm = TRUE), col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus")

## -----------------------------------------------------------------------------
goodMarkers <- colnames(myhindhe)[which(colMeans(myhindhe, na.rm = TRUE) < 0.5)]
mydata <- SubsetByLocus(mydata, goodMarkers)

## -----------------------------------------------------------------------------
mydata2 <- PipelineMapping2Parents(mydata, 
                                   freqAllowedDeviation = 0.06,
                                   useLinkage = FALSE,
                                   minLikelihoodRatio = 2)

## ----message = FALSE, warning = FALSE-----------------------------------------
library(qqman)
overdispersionP <- TestOverdispersion(mydata2, to_test = 6:10)
qq(overdispersionP[["6"]])
qq(overdispersionP[["7"]])
qq(overdispersionP[["8"]])
qq(overdispersionP[["9"]])
qq(overdispersionP[["10"]])

## -----------------------------------------------------------------------------
mydata <- PipelineMapping2Parents(mydata, 
                                  freqAllowedDeviation = 0.06,
                                  useLinkage = TRUE, overdispersion = 9,
                                  minLikelihoodRatio = 2)

## -----------------------------------------------------------------------------
table(mydata$alleleFreq)

## -----------------------------------------------------------------------------
mydata$alleleDepth[8,19:26]
mydata$genotypeLikelihood[[1]][,8,19:26]
mydata$genotypeLikelihood[[2]][,8,19:26]

## -----------------------------------------------------------------------------
mydata$priorProb[[1]][,19:26]
mydata$priorProb[[2]][,19:26]

## -----------------------------------------------------------------------------
mydata$ploidyChiSq[,19:26]

## -----------------------------------------------------------------------------
plot(mydata$ploidyChiSq[1,], mydata$ploidyChiSq[2,], 
     xlab = "Chi-squared for diploid model",
     ylab = "Chi-squared for tetraploid model")

## -----------------------------------------------------------------------------
mydata$posteriorProb[[1]][,10,19:26]
mydata$posteriorProb[[2]][,10,19:26]

## -----------------------------------------------------------------------------
mydata <- SubsetByPloidy(mydata, ploidies = list(2))

## -----------------------------------------------------------------------------
mydata <- RemoveUngenotypedLoci(mydata)

## -----------------------------------------------------------------------------
mywm <- GetWeightedMeanGenotypes(mydata)
round(mywm[c(276, 277, 1:5), 10:13], 3)

## -----------------------------------------------------------------------------
mydata$likelyGeno_donor[,19:26]
mydata$likelyGeno_recurrent[,19:26]

## ----message=FALSE, warning=FALSE---------------------------------------------
library(VariantAnnotation)

myVCF <- system.file("extdata", "Msi01genes.vcf", package = "polyRAD")

## ----eval=FALSE---------------------------------------------------------------
#  mybg <- bgzip(myVCF)
#  indexTabix(mybg, format = "vcf")

## -----------------------------------------------------------------------------
mydata <- VCF2RADdata(myVCF, possiblePloidies = list(2, c(2,2)),
                      expectedLoci = 100, expectedAlleles = 500)
mydata

## -----------------------------------------------------------------------------
myhindhe <- HindHe(mydata)
myhindheByLoc <- colMeans(myhindhe, na.rm = TRUE)
hist(myhindheByLoc, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus")
abline(v = 0.5, col = "blue", lwd = 2)

## -----------------------------------------------------------------------------
InbreedingFromHindHe(0.35, ploidy = 2)

## -----------------------------------------------------------------------------
set.seed(803)
ExpectedHindHe(mydata, inbreeding = 0.3, ploidy = 2, reps = 10)

## -----------------------------------------------------------------------------
mean(myhindheByLoc < 0.24) # about 29% of markers would be removed
keeploci <- names(myhindheByLoc)[myhindheByLoc >= 0.24]
mydata <- SubsetByLocus(mydata, keeploci)

## -----------------------------------------------------------------------------
overdispersionP <- TestOverdispersion(mydata, to_test = 8:10)
qq(overdispersionP[["8"]])
qq(overdispersionP[["9"]])
qq(overdispersionP[["10"]])

## ----message = FALSE----------------------------------------------------------
mydataHWE <- IterateHWE(mydata, tol = 1e-3, overdispersion = 9)

## -----------------------------------------------------------------------------
hist(mydataHWE$alleleFreq, breaks = 20, col = "lightgrey")

## ----message = FALSE----------------------------------------------------------
set.seed(3908)
mydataPopStruct <- IteratePopStruct(mydata, nPcsInit = 8, tol = 5e-03,
                                    overdispersion = 9)

## -----------------------------------------------------------------------------
hist(mydataPopStruct$alleleFreq, breaks = 20, col = "lightgrey")

## -----------------------------------------------------------------------------
plot(mydataPopStruct)

## -----------------------------------------------------------------------------
myallele <- 1
freqcol <- heat.colors(101)[round(mydataPopStruct$alleleFreqByTaxa[,myallele] * 100) + 1]
plot(mydataPopStruct, pch = 21, bg = freqcol)

## -----------------------------------------------------------------------------
plot(mydataPopStruct$ploidyChiSq[1,], mydataPopStruct$ploidyChiSq[2,], 
     xlab = "Chi-squared for diploid model",
     ylab = "Chi-squared for allotetraploid model", log = "xy")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## -----------------------------------------------------------------------------
myChiSqRat <- mydataPopStruct$ploidyChiSq[1,] / mydataPopStruct$ploidyChiSq[2,]
myChiSqRat <- tapply(myChiSqRat, mydataPopStruct$alleles2loc, mean)
allelesPerLoc <- as.vector(table(mydataPopStruct$alleles2loc))

library(ggplot2)
ggplot(mapping = aes(x = myhindheByLoc[GetLoci(mydata)], y = myChiSqRat, fill = as.factor(allelesPerLoc))) +
  geom_point(shape = 21, size = 3) +
  labs(x = "Hind/He", y = "Ratio of Chi-squared values, diploid to allotetraploid",
       fill = "Alleles per locus") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0.5) +
  scale_fill_brewer(palette = "YlOrRd")

## -----------------------------------------------------------------------------
wmgenoPopStruct <- GetWeightedMeanGenotypes(mydataPopStruct)
wmgenoPopStruct[1:10,1:5]

## ----eval = FALSE-------------------------------------------------------------
#  myHindHe <- HindHe(mydata)
#  TotDepthT <- rowSums(mydata$locDepth)

## -----------------------------------------------------------------------------
print(load(system.file("extdata", "MsaHindHe0.RData", package = "polyRAD")))

## -----------------------------------------------------------------------------
myHindHeByInd <- rowMeans(myHindHe, na.rm = TRUE)

## -----------------------------------------------------------------------------
ggplot(data.frame(Depth = TotDepthT, HindHe = myHindHeByInd,
                  Ploidy = ploidies),
  mapping = aes(x = Depth, y = HindHe, color = Ploidy)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~ Ploidy) +
  geom_hline(data = data.frame(Ploidy = c("2x", "3x", "4x"),
                               ExpHindHe = c(1/2, 2/3, 3/4)),
             mapping = aes(yintercept = ExpHindHe), lty = 2) +
  labs(x = "Read Depth", y = "Hind/He", color = "Ploidy")

## -----------------------------------------------------------------------------
myHindHe2x <- myHindHe[ploidies == "2x",]
myHindHe4x <- myHindHe[ploidies == "4x",]

## -----------------------------------------------------------------------------
myHindHeByLoc2x <- colMeans(myHindHe2x, na.rm = TRUE)
hist(myHindHeByLoc2x, breaks = 50, xlab = "Hind/He",
     main = "Distribution of Hind/He among loci in diploids",
     col = "lightgrey")
abline(v = 0.5, col = "blue", lwd = 2)

myHindHeByLoc4x <- colMeans(myHindHe4x, na.rm = TRUE)
hist(myHindHeByLoc4x, breaks = 50, xlab = "Hind/He",
     main = "Distribution of Hind/He among loci in tetraploids",
     col = "lightgrey")
abline(v = 0.75, col = "blue", lwd = 2)

## -----------------------------------------------------------------------------
goodLoci <- colnames(myHindHe)[myHindHeByLoc2x < 0.5 & myHindHeByLoc4x < 0.75]
length(goodLoci) # 3233 out of 5182 markers retained
head(goodLoci)

## ----eval = FALSE-------------------------------------------------------------
#  library(polyRAD)
#  library(VariantAnnotation)
#  
#  # Two files produced by the TASSEL-GBSv2 pipeline using two different
#  # enzyme systems.
#  NsiI_file <- "170705Msi_NsiI_genotypes.vcf.bgz"
#  PstI_file <- "170608Msi_PstI_genotypes.vcf.bgz"
#  
#  # The vector allSam was defined outside of this script, and contains the
#  # names of all samples that I wanted to import.  Below I find sample names
#  # within the VCF files that match those samples.
#  NsiI_sam <- allSam[allSam %in% samples(scanVcfHeader(NsiI_file))]
#  PstI_sam <- allSam[allSam %in% samples(scanVcfHeader(PstI_file))]
#  
#  # Import two RADdata objects, assuming diploidy.  A large yield size was
#  # used due to the computer having 64 Gb RAM; on a typical laptop you
#  # would probably want to keep the default of 5000.
#  PstI_RAD <- VCF2RADdata(PstI_file, samples = PstI_sam, yieldSize = 5e4,
#                          expectedAlleles = 1e6, expectedLoci = 2e5)
#  NsiI_RAD <- VCF2RADdata(NsiI_file, samples = NsiI_sam, yieldSize = 5e4,
#                          expectedAlleles = 1e6, expectedLoci = 2e5)
#  
#  # remove any loci duplicated across the two sets
#  nLoci(PstI_RAD)    # 116757
#  nLoci(NsiI_RAD)    # 187434
#  nAlleles(PstI_RAD) # 478210
#  nAlleles(NsiI_RAD) # 952511
#  NsiI_keeploci <- which(!GetLoci(NsiI_RAD) %in% GetLoci(PstI_RAD))
#  cat(nLoci(NsiI_RAD) - length(NsiI_keeploci),
#      file = "180522Num_duplicate_loci.txt") #992 duplicate
#  NsiI_RAD <- SubsetByLocus(NsiI_RAD, NsiI_keeploci)
#  
#  # combine allele depth into one matrix
#  PstI_depth <- PstI_RAD$alleleDepth
#  NsiI_depth <- NsiI_RAD$alleleDepth
#  total_depth <- matrix(0L, nrow = length(allSam),
#                        ncol = ncol(PstI_depth) + ncol(NsiI_depth),
#                        dimnames = list(allSam,
#                                        c(colnames(PstI_depth),
#                                          colnames(NsiI_depth))))
#  total_depth[,colnames(PstI_depth)] <- PstI_depth[allSam,]
#  total_depth[rownames(NsiI_depth),colnames(NsiI_depth)] <- NsiI_depth
#  
#  # combine other slots
#  total_alleles2loc <- c(PstI_RAD$alleles2loc,
#                         NsiI_RAD$alleles2loc + nLoci(PstI_RAD))
#  total_locTable <- rbind(PstI_RAD$locTable, NsiI_RAD$locTable)
#  total_alleleNucleotides <- c(PstI_RAD$alleleNucleotides,
#                               NsiI_RAD$alleleNucleotides)
#  
#  # build new RADdata object and save
#  total_RAD <- RADdata(total_depth, total_alleles2loc, total_locTable,
#                       list(2L), 0.001, total_alleleNucleotides)
#  #save(total_RAD, file = "180524_RADdata_NsiIPstI.RData")
#  
#  # Make groups representing pairs of chromosomes, and one group for all
#  # non-assembled scaffolds.
#  splitlist <- list(c("^01$", "^02$"),
#                    c("^03$", "^04$"),
#                    c("^05$", "^06$"),
#                    c("^07$", "^08$"),
#                    c("^09$", "^10$"),
#                    c("^11$", "^12$"),
#                    c("^13$", "^14$", "^15$"),
#                    c("^16$", "^17$"),
#                    c("^18$", "^194"), "^SCAFFOLD")
#  # split by chromosome and save seperate objects
#  SplitByChromosome(total_RAD, chromlist = splitlist,
#                    chromlist.use.regex = TRUE, fileprefix = "180524splitRAD")
#  
#  # files with RADdata objects
#  splitfiles <- grep("^180524splitRAD", list.files("."), value = TRUE)
#  
#  # list to hold markers formatted for GAPIT/FarmCPU
#  GAPITlist <- list()
#  length(GAPITlist) <- length(splitfiles)
#  
#  # loop through RADdata objects
#  for(i in 1:length(splitfiles)){
#    load(splitfiles[i])
#    splitRADdata <- IteratePopStructLD(splitRADdata)
#    GAPITlist[[i]] <- ExportGAPIT(splitRADdata)
#  }
#  #save(GAPITlist, file = "180524GAPITlist.RData")
#  
#  # put together into one dataset for FarmCPU
#  GM.all <- rbind(GAPITlist[[1]]$GM, GAPITlist[[2]]$GM, GAPITlist[[3]]$GM,
#                  GAPITlist[[4]]$GM, GAPITlist[[5]]$GM, GAPITlist[[6]]$GM,
#                  GAPITlist[[7]]$GM, GAPITlist[[8]]$GM,
#                  GAPITlist[[9]]$GM, GAPITlist[[10]]$GM)
#  GD.all <- cbind(GAPITlist[[1]]$GD, GAPITlist[[2]]$GD[,-1],
#                  GAPITlist[[3]]$GD[,-1], GAPITlist[[4]]$GD[,-1],
#                  GAPITlist[[5]]$GD[,-1], GAPITlist[[6]]$GD[,-1],
#                  GAPITlist[[7]]$GD[,-1], GAPITlist[[8]]$GD[,-1],
#                  GAPITlist[[9]]$GD[,-1], GAPITlist[[10]]$GD[,-1])
#  #save(GD.all, GM.all, file = "180525GM_GD_all_polyRAD.RData") # 1076888 markers

