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

## ----eval = FALSE-------------------------------------------------------------
#  overdispersionP <- TestOverdispersion(myRADprelim, to_test = 9:14)
#  
#  sapply(overdispersionP[names(overdispersionP) != "optimal"],
#         quantile, probs = c(0.01, 0.25, 0.5, 0.75, 0.99))

## ----echo = FALSE-------------------------------------------------------------
cat("Optimal value is 12.\n", sep = "\n")

print(structure(c(0.025793146160144, 0.253340971566736, 0.552888349753057, 
0.84390807760186, 1, 0.019449092419997, 0.241258741258741, 0.541687845916098, 
0.836394612538573, 1, 0.0147776703743338, 0.226466237547336, 
0.523304087800422, 0.829407211200752, 1, 0.0114170501660437, 
0.211650502292279, 0.505828679120272, 0.82350760103819, 1, 0.00893782083692289, 
0.197336617073809, 0.490038344650588, 0.818544664007547, 1, 0.00697809321403532, 
0.185307346326836, 0.4777134652488, 0.813849077168056, 1), .Dim = 5:6, .Dimnames = list(
    c("1%", "25%", "50%", "75%", "99%"), c("9", "10", "11", "12", 
    "13", "14"))))

## ----eval = FALSE-------------------------------------------------------------
#  ExpectedHindHe(myRADprelim, inbreeding = 0.4, ploidy = 2, overdispersion = 12)

## ----echo = FALSE-------------------------------------------------------------
message("Simulating rep 1")
message("Completed 5 simulation reps")
cat(c("Mean Hind/He: 0.286",
      "Standard deviation: 0.0891",
      "95% of observations are between 0.144 and 0.481"), sep = "\n")
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
#  myRAD <- IteratePopStruct(myRAD, overdispersion = 12)

## ----eval = FALSE-------------------------------------------------------------
#  RADdata2VCF(myRAD, file = "Msa_test.vcf")

