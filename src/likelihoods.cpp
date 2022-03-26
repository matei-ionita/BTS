#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// enumerate over all causal configurations
// using transformation from decimal to base n
void base_change(int i, int n, arma::uvec& tup) {
  int pos=0;
  while(i > 0) {
    tup(pos++)=i%n;
    i/=n;
  }
  return;
}

// check for duplicates: configs which are
// permutations of existing ones
bool is_duplicate_config(arma::uvec& tup, int d) {
  for (int i=0; i<d-1; i++) {
    if (tup(i) > tup(i+1))
      return true;
    if (tup(i)!=tup(0) && tup(i)==tup(i+1))
      return true;
  }

  return false;
}


double compute_log_bayes(arma::mat&z, arma::mat& ld, arma::uvec& tup, double her) {
	arma::uvec unq = arma::unique(tup);
	double size = unq.n_elem;

	// subset rows and columns corresponding to causal config
	arma::mat z_sub=z.rows(unq);
	arma::mat ld_sub=ld.submat(unq,unq);

	// log Bayes factor computed using matrix inversion lemma:
	// ld^{-1} - (ld + her/size * ld^2)^{-1} = (I + her/size * ld)^{-1}
	arma::mat mat_reg= arma::eye(size,size) + her / size * ld_sub;
	arma::mat exponent=arma::trans(z_sub) * arma::inv_sympd(mat_reg) * z_sub;
	double log_bayes = (-log(arma::det(mat_reg)) + her/size * exponent(0,0))/2;

	return log_bayes;
}


// compute maximum Bayes factor over configs with d=1
double initialize_max_bayes(arma::mat&z, arma::mat& ld, double her, double n) {
	double max_log_bayes = -arma::datum::inf;
	arma::uvec tup=arma::uvec(1);

	for (unsigned j=0; j<n; j++) {
		tup(0)=j;
		double log_bayes=compute_log_bayes(z,ld,tup,her);
		if(log_bayes > max_log_bayes)
			max_log_bayes=log_bayes;
	}
	return(max_log_bayes);
}


// [[Rcpp::export]]
Rcpp::List enumerate_configs(arma::vec& z0, arma::mat& ld, 
														 double her, int d, double thresh) {
	arma::mat z = arma::conv_to<arma::mat>::from(z0);

	int n = z.size();
	unsigned long long max = pow(n,d);

	double max_log_bayes=initialize_max_bayes(z,ld,her,n);

	unsigned long long i=0;
	unsigned long long cnt=0;
	arma::uvec tup=arma::uvec(d);
	double log_bayes;

	// compute log bayes factors for all configs with up to d variants
	// count those large enough to matter
	while(i<max) {
		base_change(i,n,tup);
		if (is_duplicate_config(tup, d)) {
			i++;
			continue;
		}

		log_bayes=compute_log_bayes(z, ld, tup, her);
		if (max_log_bayes - log_bayes < thresh)
			cnt++;
		i++;
	}

	// initialize arrays of appropriate dimension
	// configs: row i is set of causal variants in config i
	// log_bayes_vec: row i is log Bayes for config i
	IntegerMatrix configs = IntegerMatrix(cnt,d);
	NumericVector log_bayes_vec = NumericVector(cnt);

	i=0;
	cnt=0;
	tup=arma::uvec(d);

	// compute log bayes factors again, this time storing relevant ones
	while(i<max) {
		base_change(i,n,tup);
		if (is_duplicate_config(tup, d)) {
			i++;
			continue;
		}
		log_bayes=compute_log_bayes(z, ld, tup, her);
		if (max_log_bayes - log_bayes >=thresh) {
			i++;
			continue;
		}

		for (int j=0; j<d; j++)
			if (j>0 && tup(j)==tup(j-1))
				configs(cnt,j)=n+1;
			else
				configs(cnt,j)=tup(j)+1;

		log_bayes_vec(cnt)=log_bayes;

		cnt++;
		i++;
	}

	return Rcpp::List::create(Rcpp::Named("configs")=configs, 
														Rcpp::Named("log_bayes")=log_bayes_vec);
}

