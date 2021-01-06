## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
system.file("python", "process_sam_multi.py", package = "polyRAD")

## -----------------------------------------------------------------------------
library(polyRAD)

## ----eval = FALSE-------------------------------------------------------------
#  myRADprelim <- readProcessSamMulti("Msa_split_1_align.csv")

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  # subset the object to have diploids and a few tetras
#  diploids <- readLines("diploids.txt")
#  myRADprelim <- SubsetByTaxon(myRADprelim, c(diploids, "KMS397", "KMS444", "UI11-00032"))

## ----eval = FALSE-------------------------------------------------------------
#  hh <- HindHe(myRADprelim)
#  TotDepthT <- rowSums(myRADprelim$locDepth)

## ----echo = FALSE-------------------------------------------------------------
load(system.file("extdata", "MsaHindHe.RData", package = "polyRAD"))

## -----------------------------------------------------------------------------
hhByInd <- rowMeans(hh, na.rm = TRUE)

plot(TotDepthT, hhByInd, xlog = TRUE,
     xlab = "Depth", ylab = "Hind/He", main = "Samples")
abline(h = 0.5, lty = 2)

## -----------------------------------------------------------------------------
threshold <- mean(hhByInd) + 3 * sd(hhByInd)
threshold

hhByInd[hhByInd > threshold]
hh <- hh[hhByInd <= threshold,]

## ----eval = FALSE-------------------------------------------------------------
#  myRADprelim <- SubsetByTaxon(myRADprelim, rownames(hh))

## ----eval = FALSE-------------------------------------------------------------
#  writeLines(rownames(hh), con = "samples.txt")

## -----------------------------------------------------------------------------
hhByLoc <- colMeans(hh, na.rm = TRUE)

hist(hhByLoc, breaks = 50, xlab = "Hind/He", main = "Loci", col = "lightgrey")

## -----------------------------------------------------------------------------
InbreedingFromHindHe(hindhe = 0.3, ploidy = 2)

## ----eval = FALSE-------------------------------------------------------------
#  ExpectedHindHe(myRADprelim, inbreeding = 0.4, ploidy = 2)

## ----echo = FALSE-------------------------------------------------------------
message("Simulating rep 1")
message("Completed 5 simulation reps")
cat(c("Mean Hind/He: 0.293",
      "Standard deviation: 0.0881",
      "95% of observations are between 0.149 and 0.469"), sep = "\n")
load(system.file("extdata", "MsaHindHe3.RData", package = "polyRAD"))
hist(testhhdist, xlab = "Hind/He", main = "Expected distribution of Hind/He",
     breaks = 30)

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
mean(hh2ByLoc <= 0.5) # proportion of loci retained
keeploci <- names(hh2ByLoc)[hh2ByLoc <= 0.5]
head(keeploci)
hist(hh2ByLoc[keeploci], xlab = "Hind/He", main = "Loci", breaks = 50, col = "lightgrey")

## ----eval = FALSE-------------------------------------------------------------
#  myRAD <- SubsetByLocus(myRAD, keeploci)

## ----eval = FALSE-------------------------------------------------------------
#  myRAD <- IteratePopStruct(myRAD)

## ----eval = FALSE-------------------------------------------------------------
#  RADdata2VCF(myRAD, file = "Msa_test.vcf")

