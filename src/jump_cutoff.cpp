#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP jump_cutoff(SEXP pa_in, SEXP pb_in, SEXP xi_in, SEXP alpha_in) {
  try{
    vec xi = as<arma::vec>(xi_in);
    vec pa = as<arma::vec>(pa_in);
    vec pb = as<arma::vec>(pb_in);
    double alpha = Rcpp::as<double>(alpha_in);

    int m = pa.size();

    vec p_max(m);
    for(int i = 0; i < m; i++){
      p_max(i) = pa(i) > pb(i) ? pa(i) : pb(i);
    }

    double xi00 = xi(0), xi01 = xi(1), xi10 = xi(2);

    vec pmax_sorted = sort(p_max);

    double thr = 0, fdp = 1, thr_jump = 0;

    for(int i = m - 1; i >= 0; i--){
      thr = pmax_sorted(i);
      fdp = m * (xi00 * thr * thr + (xi01 + xi10) * thr)/(i + 1);
      if(fdp <= alpha){
        thr_jump = thr;
        break;
      }
    }

    return Rcpp::List::create(Rcpp::Named("p_max") = p_max,
                              Rcpp::Named("thr_jump") = thr_jump);

  } catch( std::exception &ex ) {
    forward_exception_to_r(ex);
    return Rcpp::List::create();
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
    return Rcpp::List::create();
  }
}
