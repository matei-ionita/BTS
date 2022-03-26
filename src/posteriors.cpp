#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// compute log(exp(x) + exp(y))
// minimizing loss of precision
double log_add(double x, double y) {
	if (x>y)
		return x + log(1 + exp(y - x));

	return y + log(1 + exp(x - y));
}

// [[Rcpp::export]]
Rcpp::List compute_post_locus(List locus_configs, arma::vec& enrich0, arma::mat& locus_annot) {
	IntegerMatrix configs=locus_configs["configs"];
	NumericVector log_bayes=locus_configs["log_bayes"];

	int n_var=locus_annot.n_rows;
	unsigned long long n_config=configs.nrow();
	int d=configs.ncol();

	// compute priors from enrichment
	// prior = 1 / (1+ exp( enrich_baseline + enrich_annot * annot ))
	arma::mat enrich = arma::conv_to<arma::mat>::from(enrich0);
	arma::rowvec exponent = arma::conv_to<arma::rowvec>::from(locus_annot * enrich);
	arma::rowvec log_prior = -log(1+exp(exponent));
	arma::rowvec log_not_prior = -log(1+exp(-exponent));
	arma::rowvec log_prior_ratios = log_prior - log_not_prior;
	log_prior_ratios.resize(n_var+1);
	log_prior_ratios(n_var)=0;

	// prior of null configuration (nothing causal)
	double prior_null = sum(log_not_prior);

	unsigned long long i;
	int var, j;
	double total_post=0, log_config_post;
	arma::vec var_post = arma::vec(n_var, arma::fill::ones) * (-arma::datum::inf);

	for (i=0; i<n_config; i++) {
		// initialize log posterior for this config
		log_config_post = log_bayes(i)+prior_null;

		// add priors for variants causal in this config
		for (j=0; j<d; j++)
			log_config_post += log_prior_ratios(configs(i,j)-1);

		// keep tally of total posterior
		total_post = log_add(total_post, log_config_post);

		// increment marginal posteriors for causal variants
		for (j=0; j<d; j++) {
			var = configs(i,j)-1;
			if (var < n_var)
				var_post(var) = log_add(var_post(var), log_config_post);
		}
	}

	// divide marginal posteriors by denominator from Bayes' theorem
	var_post = exp(var_post - total_post);
	return Rcpp::List::create(Rcpp::Named("var_post")=var_post, 
							  Rcpp::Named("total_post")=total_post);
}