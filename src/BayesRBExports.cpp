#include <Rcpp.h>

using namespace Rcpp;

// BayesRB
List BayesRB(int seed, int MCMC_inte, int burn_intee,int thinn,NumericMatrix X, NumericVector Y, 
	NumericVector beta_initial);
RcppExport SEXP BayesRB_BayesRB(SEXP seedSEXP, SEXP MCMC_inteSEXP, SEXP burn_inteeSEXP, 
	SEXP thinnSEXP , SEXP XSEXP, SEXP YSEXP, SEXP beta_initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_inte(MCMC_inteSEXP);
    Rcpp::traits::input_parameter< int >::type burn_intee(burn_inteeSEXP);
    Rcpp::traits::input_parameter< int >::type thinn(thinnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_initial(beta_initialSEXP);


    __result = Rcpp::wrap(BayesRB(seed, MCMC_inte, burn_intee, thinn, X, Y, beta_initial));
    return __result;
END_RCPP
}
