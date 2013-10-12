// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kdtree_build_intl
SEXP kdtree_build_intl(SEXP d, SEXP nr, SEXP nc);
RcppExport SEXP hypervolume_kdtree_build_intl(SEXP dSEXP, SEXP nrSEXP, SEXP ncSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        SEXP d = Rcpp::as<SEXP >(dSEXP);
        SEXP nr = Rcpp::as<SEXP >(nrSEXP);
        SEXP nc = Rcpp::as<SEXP >(ncSEXP);
        SEXP __result = kdtree_build_intl(d, nr, nc);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// kdtree_ball_query_multiple
SEXP kdtree_ball_query_multiple(SEXP tr, SEXP ptlist, SEXP nr, SEXP nc, SEXP r, SEXP verb);
RcppExport SEXP hypervolume_kdtree_ball_query_multiple(SEXP trSEXP, SEXP ptlistSEXP, SEXP nrSEXP, SEXP ncSEXP, SEXP rSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        SEXP tr = Rcpp::as<SEXP >(trSEXP);
        SEXP ptlist = Rcpp::as<SEXP >(ptlistSEXP);
        SEXP nr = Rcpp::as<SEXP >(nrSEXP);
        SEXP nc = Rcpp::as<SEXP >(ncSEXP);
        SEXP r = Rcpp::as<SEXP >(rSEXP);
        SEXP verb = Rcpp::as<SEXP >(verbSEXP);
        SEXP __result = kdtree_ball_query_multiple(tr, ptlist, nr, nc, r, verb);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// kdtree_range_query_multiple
SEXP kdtree_range_query_multiple(SEXP tr, SEXP pminlist, SEXP pmaxlist, SEXP nr, SEXP nc, SEXP verb);
RcppExport SEXP hypervolume_kdtree_range_query_multiple(SEXP trSEXP, SEXP pminlistSEXP, SEXP pmaxlistSEXP, SEXP nrSEXP, SEXP ncSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        SEXP tr = Rcpp::as<SEXP >(trSEXP);
        SEXP pminlist = Rcpp::as<SEXP >(pminlistSEXP);
        SEXP pmaxlist = Rcpp::as<SEXP >(pmaxlistSEXP);
        SEXP nr = Rcpp::as<SEXP >(nrSEXP);
        SEXP nc = Rcpp::as<SEXP >(ncSEXP);
        SEXP verb = Rcpp::as<SEXP >(verbSEXP);
        SEXP __result = kdtree_range_query_multiple(tr, pminlist, pmaxlist, nr, nc, verb);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}