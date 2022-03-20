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
#' @param in_path A path.
#' @param ld_files A vector of files.
#' @param loci A vector of locus names.
#' @param d Maximal number of independent signals considered.
#' @param thresh Numeric. Configurations are not stored if their log 
#' Bayes factor is smaller than that of the size-1 configuration with 
#' largest Bayes factor, minus the threshold.
#' @export
compute_likelihoods <- function(in_path, ld_files, loci, d=2, thresh=12) {
  if (thresh < 10)
    warning("Using a threshold smaller than 10 may lead to inexact results.")
  
  # Read annotations and other inputs
  ld_matrices <- lapply(ld_files, FUN=read_ld_matrix, path=in_path, reg=1e-5)
  names(ld_matrices) <- loci

  z_scores <- lapply(loci, FUN=read_z_scores, path=in_path)
  names(z_scores) <- loci

  # clunky, change later
  z_scores <- Map(remove_na_vars, z_scores, ld_matrices, loci)

  prior_variances <- Map(estimate_prior_variance, z_scores, ld_matrices) %>%
    unlist()

  configs <- Map(enumerate_configs, z_scores, ld_matrices,
                 prior_variances, rep(d, length(z_scores)),
                 thresh=thresh)
  return(configs)
}

remove_na_vars <- function(z, ld, locus) {
  sel <- which(abs(diag(ld)-1) > 1e-4)
  if (length(sel)>0) {
    z[sel] <- 0
    print(length(sel))
    warning(paste("Diagonal elements not 1 for locus:", locus))
  }

  return(z)
}


estimate_prior_variance <- function(z, ld_matrix,
                                    prop_ld_eigenvalues = 0.95) {
  
  n <- nrow(ld_matrix)
  eigen <- eigen(ld_matrix, symmetric=TRUE)
  
  trunc_inv <- numeric(n)
  total_eig <- 0
  max_k <- 1
  cutoff <- prop_ld_eigenvalues * n
  
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
  prior_variance <- t(z) %*% trunc_ld_inv %*% z - max_k
  
  return(min(prior_variance[1,1],500))
}



