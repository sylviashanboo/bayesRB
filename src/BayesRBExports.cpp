#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::ArrayXXd;
using Eigen::ArrayXi;
using Eigen::ArrayXd; 

// BayesRB
List BayesRB(int seed, int MCMC_inte, int burn_intee,int thinn, ArrayXXd X_unorm, ArrayXi Y, 
             ArrayXd beta_initial);
RcppExport SEXP BayesRB_BayesRB(SEXP seedSEXP, SEXP MCMC_inteSEXP, SEXP burn_inteeSEXP, 
	SEXP thinnSEXP , SEXP XSEXP, SEXP YSEXP, SEXP beta_initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_inte(MCMC_inteSEXP);
    Rcpp::traits::input_parameter< int >::type burn_intee(burn_inteeSEXP);
    Rcpp::traits::input_parameter< int >::type thinn(thinnSEXP);
    Rcpp::traits::input_parameter< ArrayXXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< ArrayXi >::type Y(YSEXP);
    Rcpp::traits::input_parameter< ArrayXd >::type beta_initial(beta_initialSEXP);


    __result = Rcpp::wrap(BayesRB(seed, MCMC_inte, burn_intee, thinn, X, Y, beta_initial));
    return __result;
END_RCPP
}

List BayesRB_OG(int seed, int MCMC_inte, int burn_intee,int thinn, NumericMatrix X_unorm, NumericVector Y,
                NumericVector beta_initial);
RcppExport SEXP BayesRB_BayesRB_OG(SEXP seedSEXP, SEXP MCMC_inteSEXP, SEXP burn_inteeSEXP,
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


  __result = Rcpp::wrap(BayesRB_OG(seed, MCMC_inte, burn_intee, thinn, X, Y, beta_initial));
  return __result;
  END_RCPP
}
