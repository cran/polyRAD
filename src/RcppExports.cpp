// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// AdjustAlleleFreq
NumericMatrix AdjustAlleleFreq(NumericMatrix predAl, IntegerVector alleles2loc, double minfreq);
RcppExport SEXP _polyRAD_AdjustAlleleFreq(SEXP predAlSEXP, SEXP alleles2locSEXP, SEXP minfreqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type predAl(predAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< double >::type minfreq(minfreqSEXP);
    rcpp_result_gen = Rcpp::wrap(AdjustAlleleFreq(predAl, alleles2loc, minfreq));
    return rcpp_result_gen;
END_RCPP
}
// BestGenos
IntegerMatrix BestGenos(NumericVector probs, int ploidy, int ntaxa, int nalleles);
RcppExport SEXP _polyRAD_BestGenos(SEXP probsSEXP, SEXP ploidySEXP, SEXP ntaxaSEXP, SEXP nallelesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< int >::type ntaxa(ntaxaSEXP);
    Rcpp::traits::input_parameter< int >::type nalleles(nallelesSEXP);
    rcpp_result_gen = Rcpp::wrap(BestGenos(probs, ploidy, ntaxa, nalleles));
    return rcpp_result_gen;
END_RCPP
}
// CorrectGenos
IntegerMatrix CorrectGenos(IntegerMatrix bestgenos, NumericVector probs, IntegerVector alleles2loc, int ntaxa, int ploidy, int nalleles, int nloc, bool do_correct);
RcppExport SEXP _polyRAD_CorrectGenos(SEXP bestgenosSEXP, SEXP probsSEXP, SEXP alleles2locSEXP, SEXP ntaxaSEXP, SEXP ploidySEXP, SEXP nallelesSEXP, SEXP nlocSEXP, SEXP do_correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type bestgenos(bestgenosSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< int >::type ntaxa(ntaxaSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< int >::type nalleles(nallelesSEXP);
    Rcpp::traits::input_parameter< int >::type nloc(nlocSEXP);
    Rcpp::traits::input_parameter< bool >::type do_correct(do_correctSEXP);
    rcpp_result_gen = Rcpp::wrap(CorrectGenos(bestgenos, probs, alleles2loc, ntaxa, ploidy, nalleles, nloc, do_correct));
    return rcpp_result_gen;
END_RCPP
}
// BestPloidies
IntegerVector BestPloidies(NumericMatrix chisq);
RcppExport SEXP _polyRAD_BestPloidies(SEXP chisqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type chisq(chisqSEXP);
    rcpp_result_gen = Rcpp::wrap(BestPloidies(chisq));
    return rcpp_result_gen;
END_RCPP
}
// FormatStructure
IntegerMatrix FormatStructure(IntegerMatrix genotypes, IntegerVector alleles2loc, int ploidy);
RcppExport SEXP _polyRAD_FormatStructure(SEXP genotypesSEXP, SEXP alleles2locSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type genotypes(genotypesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(FormatStructure(genotypes, alleles2loc, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// GiniSimpson
double GiniSimpson(NumericVector counts);
RcppExport SEXP _polyRAD_GiniSimpson(SEXP countsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type counts(countsSEXP);
    rcpp_result_gen = Rcpp::wrap(GiniSimpson(counts));
    return rcpp_result_gen;
END_RCPP
}
// HindHeMat
NumericMatrix HindHeMat(IntegerMatrix alleleDepth, NumericMatrix depthRatio, IntegerVector alleles2loc, int nLoci, NumericVector He);
RcppExport SEXP _polyRAD_HindHeMat(SEXP alleleDepthSEXP, SEXP depthRatioSEXP, SEXP alleles2locSEXP, SEXP nLociSEXP, SEXP HeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type alleleDepth(alleleDepthSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type depthRatio(depthRatioSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< int >::type nLoci(nLociSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type He(HeSEXP);
    rcpp_result_gen = Rcpp::wrap(HindHeMat(alleleDepth, depthRatio, alleles2loc, nLoci, He));
    return rcpp_result_gen;
END_RCPP
}
// HoOneParent
NumericVector HoOneParent(IntegerVector genotypes, IntegerVector alleles2loc, IntegerVector keeploc, double ploidy);
RcppExport SEXP _polyRAD_HoOneParent(SEXP genotypesSEXP, SEXP alleles2locSEXP, SEXP keeplocSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type genotypes(genotypesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keeploc(keeplocSEXP);
    Rcpp::traits::input_parameter< double >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(HoOneParent(genotypes, alleles2loc, keeploc, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// HoTwoParents
NumericVector HoTwoParents(IntegerVector genotypes1, IntegerVector genotypes2, IntegerVector alleles2loc, IntegerVector keeploc, double ploidy);
RcppExport SEXP _polyRAD_HoTwoParents(SEXP genotypes1SEXP, SEXP genotypes2SEXP, SEXP alleles2locSEXP, SEXP keeplocSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type genotypes1(genotypes1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type genotypes2(genotypes2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keeploc(keeplocSEXP);
    Rcpp::traits::input_parameter< double >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(HoTwoParents(genotypes1, genotypes2, alleles2loc, keeploc, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// InitHapAssign
IntegerVector InitHapAssign(IntegerMatrix NMmat);
RcppExport SEXP _polyRAD_InitHapAssign(SEXP NMmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type NMmat(NMmatSEXP);
    rcpp_result_gen = Rcpp::wrap(InitHapAssign(NMmat));
    return rcpp_result_gen;
END_RCPP
}
// Hap2SNP
List Hap2SNP(StringVector haps, std::string refhap, int pos);
RcppExport SEXP _polyRAD_Hap2SNP(SEXP hapsSEXP, SEXP refhapSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type haps(hapsSEXP);
    Rcpp::traits::input_parameter< std::string >::type refhap(refhapSEXP);
    Rcpp::traits::input_parameter< int >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(Hap2SNP(haps, refhap, pos));
    return rcpp_result_gen;
END_RCPP
}
// Hap2Hap
List Hap2Hap(StringVector haps, std::string refhap, int pos);
RcppExport SEXP _polyRAD_Hap2Hap(SEXP hapsSEXP, SEXP refhapSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type haps(hapsSEXP);
    Rcpp::traits::input_parameter< std::string >::type refhap(refhapSEXP);
    Rcpp::traits::input_parameter< int >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(Hap2Hap(haps, refhap, pos));
    return rcpp_result_gen;
END_RCPP
}
// MakeGTstrings
StringVector MakeGTstrings(IntegerMatrix genotypes, int ploidy);
RcppExport SEXP _polyRAD_MakeGTstrings(SEXP genotypesSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type genotypes(genotypesSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(MakeGTstrings(genotypes, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// PrepVCFexport
List PrepVCFexport(IntegerMatrix genotypes, IntegerVector alleles2loc, IntegerMatrix alleleDepth, StringVector alleleNucleotides, DataFrame locTable, IntegerVector ploidy, bool asSNPs);
RcppExport SEXP _polyRAD_PrepVCFexport(SEXP genotypesSEXP, SEXP alleles2locSEXP, SEXP alleleDepthSEXP, SEXP alleleNucleotidesSEXP, SEXP locTableSEXP, SEXP ploidySEXP, SEXP asSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type genotypes(genotypesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type alleleDepth(alleleDepthSEXP);
    Rcpp::traits::input_parameter< StringVector >::type alleleNucleotides(alleleNucleotidesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type locTable(locTableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< bool >::type asSNPs(asSNPsSEXP);
    rcpp_result_gen = Rcpp::wrap(PrepVCFexport(genotypes, alleles2loc, alleleDepth, alleleNucleotides, locTable, ploidy, asSNPs));
    return rcpp_result_gen;
END_RCPP
}
// simGeno
NumericMatrix simGeno(NumericVector alleleFreq, IntegerVector alleles2loc, int nsam, double inbreeding, int ploidy);
RcppExport SEXP _polyRAD_simGeno(SEXP alleleFreqSEXP, SEXP alleles2locSEXP, SEXP nsamSEXP, SEXP inbreedingSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alleleFreq(alleleFreqSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< int >::type nsam(nsamSEXP);
    Rcpp::traits::input_parameter< double >::type inbreeding(inbreedingSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(simGeno(alleleFreq, alleles2loc, nsam, inbreeding, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// simGenoMapping
NumericMatrix simGenoMapping(NumericVector donorGeno, NumericVector recurGeno, NumericMatrix progGeno, NumericVector genoProbs, IntegerVector alleles2loc, int nsam, int ploidy);
RcppExport SEXP _polyRAD_simGenoMapping(SEXP donorGenoSEXP, SEXP recurGenoSEXP, SEXP progGenoSEXP, SEXP genoProbsSEXP, SEXP alleles2locSEXP, SEXP nsamSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type donorGeno(donorGenoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type recurGeno(recurGenoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type progGeno(progGenoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type genoProbs(genoProbsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< int >::type nsam(nsamSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(simGenoMapping(donorGeno, recurGeno, progGeno, genoProbs, alleles2loc, nsam, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// simAD
IntegerMatrix simAD(IntegerMatrix locDepth, NumericMatrix genotypes, IntegerVector alleles2loc, double overdispersion, double contamRate, NumericVector alleleFreq, double errorRate);
RcppExport SEXP _polyRAD_simAD(SEXP locDepthSEXP, SEXP genotypesSEXP, SEXP alleles2locSEXP, SEXP overdispersionSEXP, SEXP contamRateSEXP, SEXP alleleFreqSEXP, SEXP errorRateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type locDepth(locDepthSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type genotypes(genotypesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alleles2loc(alleles2locSEXP);
    Rcpp::traits::input_parameter< double >::type overdispersion(overdispersionSEXP);
    Rcpp::traits::input_parameter< double >::type contamRate(contamRateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alleleFreq(alleleFreqSEXP);
    Rcpp::traits::input_parameter< double >::type errorRate(errorRateSEXP);
    rcpp_result_gen = Rcpp::wrap(simAD(locDepth, genotypes, alleles2loc, overdispersion, contamRate, alleleFreq, errorRate));
    return rcpp_result_gen;
END_RCPP
}
// ThirdDimProd
NumericMatrix ThirdDimProd(NumericVector probs, int ngen, int ntaxa);
RcppExport SEXP _polyRAD_ThirdDimProd(SEXP probsSEXP, SEXP ngenSEXP, SEXP ntaxaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< int >::type ntaxa(ntaxaSEXP);
    rcpp_result_gen = Rcpp::wrap(ThirdDimProd(probs, ngen, ntaxa));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_polyRAD_AdjustAlleleFreq", (DL_FUNC) &_polyRAD_AdjustAlleleFreq, 3},
    {"_polyRAD_BestGenos", (DL_FUNC) &_polyRAD_BestGenos, 4},
    {"_polyRAD_CorrectGenos", (DL_FUNC) &_polyRAD_CorrectGenos, 8},
    {"_polyRAD_BestPloidies", (DL_FUNC) &_polyRAD_BestPloidies, 1},
    {"_polyRAD_FormatStructure", (DL_FUNC) &_polyRAD_FormatStructure, 3},
    {"_polyRAD_GiniSimpson", (DL_FUNC) &_polyRAD_GiniSimpson, 1},
    {"_polyRAD_HindHeMat", (DL_FUNC) &_polyRAD_HindHeMat, 5},
    {"_polyRAD_HoOneParent", (DL_FUNC) &_polyRAD_HoOneParent, 4},
    {"_polyRAD_HoTwoParents", (DL_FUNC) &_polyRAD_HoTwoParents, 5},
    {"_polyRAD_InitHapAssign", (DL_FUNC) &_polyRAD_InitHapAssign, 1},
    {"_polyRAD_Hap2SNP", (DL_FUNC) &_polyRAD_Hap2SNP, 3},
    {"_polyRAD_Hap2Hap", (DL_FUNC) &_polyRAD_Hap2Hap, 3},
    {"_polyRAD_MakeGTstrings", (DL_FUNC) &_polyRAD_MakeGTstrings, 2},
    {"_polyRAD_PrepVCFexport", (DL_FUNC) &_polyRAD_PrepVCFexport, 7},
    {"_polyRAD_simGeno", (DL_FUNC) &_polyRAD_simGeno, 5},
    {"_polyRAD_simGenoMapping", (DL_FUNC) &_polyRAD_simGenoMapping, 7},
    {"_polyRAD_simAD", (DL_FUNC) &_polyRAD_simAD, 7},
    {"_polyRAD_ThirdDimProd", (DL_FUNC) &_polyRAD_ThirdDimProd, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_polyRAD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
