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
#' Bayes factor is smaller than that of the size-1 configuration with 
#' largest Bayes factor, minus the threshold.
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

  # Estimate heritability of trait using truncated eigendecomposition.
  # Used as prior for the variance of effect sizes.
  her <- Map(estimate_her, ld_matrices, z_scores)
  
  # Store a sparse array of relevant configurations: i-th row
  # contains the variants which are causal for i-th config.
  # Separate numeric vector gives log Bayes factors of configs.
  configs <- Map(enumerate_configs, z_scores, ld_matrices,
                 her, rep(d, length(z_scores)), thresh=thresh)
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


# Estimate heritability of trait based on variants in a locus.
# The observed effects z are related to true effects z0 by
# z=ld %*% z0. Heritability is t(z0) %*% ld %*% z0 =
# t(z) %*% ld^{-1} %*% z. Since ld is typically singular,
# using a truncated eigendecomposition.
estimate_her <- function(ld, z, prop_eig = 0.95) {
  n <- nrow(ld)
  eigen <- eigen(ld, symmetric=TRUE)
  
  trunc_inv <- numeric(n)
  total_eig <- 0
  max_k <- 1
  cutoff <- prop_eig * n
  
  for (i in seq(n)) {
    total_eig <- total_eig + eigen$values[i]
    if (total_eig < cutoff | i==1) {
      max_k <- i
      trunc_inv[i] <- 1 / eigen$values[i]
    }
  }
  
  if (n==1) {
    trunc_ld_inv <- trunc_inv * eigen$vectors %*% t(eigen$vectors)
  } else {
    trunc_ld_inv <- eigen$vectors %*% diag(trunc_inv) %*% t(eigen$vectors)
  }
  her <- t(z) %*% trunc_ld_inv %*% z - max_k
  
  if (her <= 0) {
    # Clunky: find a better fix
    return(estimate_her(ld,z,prop_eig=prop_eig+0.01))
  }
  
  return(min(her[1,1],500))
}



