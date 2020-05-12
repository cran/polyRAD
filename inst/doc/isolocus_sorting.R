## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
system.file("python", "process_sam_multi.py", package = "polyRAD")

## -----------------------------------------------------------------------------
library(polyRAD)

## ----eval = FALSE-------------------------------------------------------------
#  myRADprelim <- readProcessSamMulti("Msa_split_1_align.csv")

## ----eval = FALSE-------------------------------------------------------------
#  hh <- HindHe(myRADprelim)
#  TotDepthT <- rowSums(myRADprelim$locDepth)

## ----echo = FALSE-------------------------------------------------------------
load(system.file("extdata", "MsaHindHe.RData", package = "polyRAD"))
hh <- myHindHe
rm(myHindHe)
# get all the diploids, and a few tetraploids to demo filtering
keepind <- sort(c(which(ploidies == "2x"), 
                  match(c("KMS397", "KMS444", "UI11-00032"), rownames(hh))))
hh <- hh[keepind,]
ploidies <- ploidies[keepind]
TotDepthT <- TotDepthT[keepind]

## -----------------------------------------------------------------------------
hhByInd <- rowMeans(hh, na.rm = TRUE)

plot(TotDepthT, hhByInd, xlog = TRUE,
     xlab = "Depth", ylab = "Hind/He", main = "Samples")
abline(h = 0.5, lty = 2)

## -----------------------------------------------------------------------------
hhByInd[hhByInd > 0.5]
hh <- hh[hhByInd <= 0.5,]

## ----eval = FALSE-------------------------------------------------------------
#  writeLines(rownames(hh), con = "samples.txt")

## -----------------------------------------------------------------------------
hhByLoc <- colMeans(hh, na.rm = TRUE)

hist(hhByLoc, breaks = 100, xlab = "Hind/He", main = "Loci", col = "lightgrey")

## -----------------------------------------------------------------------------
InbreedingFromHindHe(hindhe = 0.3, ploidy = 2)

## ----eval = FALSE-------------------------------------------------------------
#  myRAD <- readProcessIsoloci("Msa_split_1_sorted.csv", min.ind.with.reads = 80,
#                              min.ind.with.minor.allele = 5)

## ----eval = FALSE-------------------------------------------------------------
#  hh2 <- HindHe(myRAD)
#  hh2ByInd <- rowMeans(hh2, na.rm = TRUE)
#  hh2ByLoc <- colMeans(hh2, na.rm = TRUE)

## ----echo = FALSE-------------------------------------------------------------
load(system.file("extdata", "MsaHindHe2.RData", package = "polyRAD"))

## -----------------------------------------------------------------------------
hist(hh2ByInd, xlab = "Hind/He", main = "Samples", breaks = 20, col = "lightgrey")
hist(hh2ByLoc, xlab = "Hind/He", main = "Loci", breaks = 50, col = "lightgrey")

## -----------------------------------------------------------------------------
mean(hh2ByLoc < 0.45) # proportion of loci retained
keeploci <- names(hh2ByLoc)[hh2ByLoc < 0.45]
head(keeploci)

## ----eval = FALSE-------------------------------------------------------------
#  myRAD <- SubsetByLocus(myRAD, keeploci)

## ----eval = FALSE-------------------------------------------------------------
#  myRAD <- IteratePopStruct(myRAD)

