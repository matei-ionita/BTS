#' @useDynLib BTS, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim p.adjust pchisq
#' @importFrom magrittr "%>%"
#' @import ggplot2
#' @importFrom readr read_tsv read_delim cols
#' @importFrom stringr str_remove str_replace str_split
#' @importFrom tibble tibble
#' @importFrom dplyr arrange
NULL



#' @title Enumerate configurations; compute and store likelihoods.
#' @description Enumerate over all configurations with at most d
#' causal signals. Compute log likelihoods and store them for future
#' use in posterior computations.
#' @param path A path to data.
#' @param loci A vector of locus names; each locus must have a sub-directory
#' under the given path.
#' @param d Integer, maximal number of independent signals considered.
#' @param thresh Numeric. Configurations are not stored if their log 
#' Bayes factor is smaller than that of the size 1 configuration with 
#' largest Bayes factor, minus the threshold.
#' @return A list of two elements: an integer matrix of shape (n_config,d)
#' and a numeric vector of length n_config, where n_config is the number of 
#' configurations of d variants with large enough Bayes factor, at most
#' choose(n_variants, d). The rows of the matrix are indices of variants
#' giving a causal configuration, and the corresponding row of the vector is
#' the log Bayes factor of the configuration.
#' @export
compute_configs <- function(path, loci, d=2, thresh=12) {
  if (!(as.integer(d)==d) || d <= 0)
    stop("d must be a positive integer.")
  
  if (!is.numeric(thresh) || thresh < 0)
    stop("thresh must be a positive number.")
  
  if (thresh < 10)
    warning("Using a threshold smaller than 10 may lead to inexact results.")
  
  # Read LD matrices and GWAS summary statistics from files
  ld_matrices <- lapply(loci, FUN=read_ld_matrix, path=path, reg=1e-5)
  z_scores <- lapply(loci, FUN=read_z_scores, path=path)
  check <- Map(check_GWAS_LD_input, ld_matrices, z_scores, loci)

  # prior for the variance of effect sizes.
  pr_var <- lapply(z_scores, function(z) 16)
  
  # Store a sparse array of relevant configurations: i-th row
  # contains the variants which are causal for i-th config.
  # Separate numeric vector gives log Bayes factors of configs.
  configs <- Map(enumerate_configs, z_scores, ld_matrices,
                 pr_var, rep(d, length(z_scores)), thresh=thresh)
  return(configs)
}

check_GWAS_LD_input <- function(ld, z, locus) {
  if(any( abs(diag(ld)-1) > 1e-4 ))
    stop(paste0(locus, ": LD matrix must have 1 on diagonal."))
  
  if(max(abs( t(ld)-ld )) > 1e-10)
    stop(paste0(locus, ": LD matrix must be symmetric."))
  
  if(nrow(ld) != length(z))
    stop(paste0(locus, ": different number of variants in LD matrix
                 and z-scores."))
}



