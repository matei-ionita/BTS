

#' @title Fit a model using EM algorithm.
#' @description Given likelihoods and annotations for all loci,
#' find optimal values for variant posterior probabilities to
#' be causal, as well as enrichment of annotations in causal
#' variants. Uses Expectation-Maximization to optimize total
#' data likelihood.
#' @param configs Configurations of causal variants and their
#' log Bayes factors.
#' @param annot For each locus, a matrix of annotations.
#' @param max_iter Positive integer, maximum number of 
#' EM iterations after which algorithm stops regardless of
#' convergence.
#' @param tol Tolerance for convergence.
#' @export
run_em <- function(configs, annot,
                   max_iter=20, tol=1e-2) {
  message(colnames(annot[[1]])[-1])
  enrich <- initialize_enrich(annot)
  lower <- rep(-20,length(enrich))
  upper <- rep( 20,length(enrich))
  iter <- 1
  prev_log_post <- -Inf

  while(iter <= max_iter) {
    # Expectation
    e_step <- Map(compute_post_locus, configs,
                  rep(list(enrich), length(annot)), annot)

    posteriors <- lapply(e_step, function(x) x$var_post)
    log_post_locus <- lapply(e_step, function(x) x$total_post)
    total_log_post <- Reduce(log_post_locus, f="+")

    # Maximization
    optim_res <- optim(par=enrich, fn=data_lik, gr=grad_lik,
                       posteriors, annot,
                       method = "L-BFGS-B",
                       lower=lower, upper=upper)
    enrich <- optim_res$par
    optim_val <- -1 * optim_res$value

    # Convergence criterion
    if (abs((total_log_post-prev_log_post)) < tol) {
      message(paste("Converged after", iter, "iterations."))
      break
    }

    iter <- iter+1
    prev_log_post <- total_log_post
  }

  if(iter > max_iter)
    message(paste("Did not converge after", max_iter, "iterations."))

  print(total_log_post)

  return(list(posteriors=posteriors,
              enrich=enrich,
              log_post=total_log_post,
              log_post_locus=log_post_locus))
}


# Initialize enrichment parameters to give
# uniform prior over all variants.
initialize_enrich <- function(annot) {
  loci_size <- vapply(annot, nrow, integer(1))
  mean_size <- sum(loci_size)/length(loci_size)

  enrich <- rep(0,ncol(annot[[1]]))
  names(enrich) <- colnames(annot[[1]])
  enrich["baseline"] <- log(mean_size-1)
  return(enrich)
}


# Function to optimize
data_lik <- function(enrich, posteriors, annot) {
  enrich_rep <- rep(list(enrich), length(annot))
  lik <- Map(data_lik_locus, posteriors, enrich_rep, annot) %>%
    Reduce(f="+")
  return(-lik)
}

data_lik_locus <- function(locus_post, enrich, locus_annot) {
  exponent <- locus_annot %*% matrix(enrich, ncol=1)
  lik_locus <- -locus_post * log(1 + exp(exponent)) -
    (1-locus_post) * log(1+exp(-exponent))
  return(sum(lik_locus))
}

# Gradient of function to optimize
grad_lik <- function(enrich, posteriors, annot) {
  enrich_rep <- rep(list(enrich), length(annot))
  grad <- Map(grad_lik_locus, posteriors, enrich_rep, annot) %>%
    Reduce(f="+")
  return(-grad)
}

grad_lik_locus <- function(locus_post, enrich, locus_annot) {
  exponent <- t(locus_annot %*% matrix(enrich, ncol=1))
  locus_post <- matrix(locus_post, nrow=1)
  grad_lik_locus <- -(locus_post / (1 + exp(-exponent))) %*% locus_annot +
    ((1-locus_post) / (1 + exp(exponent))) %*% locus_annot
  return(grad_lik_locus)
}



#' @title Summary and comparison of many models.
#' @description Returns a data frame of models, with enrichment
#' parameters, log likelihood, number of annotated variants, sorted
#' by p-value.
#' @param results List of results for models.
#' @param annot_inputs List of annotations for all models and loci.
#' @export
get_model_summary <- function(results, annot_inputs) {
  n <- length(results)
  result_null <- results[[1]]

  model_summary <- lapply(results[seq(2,n)], function(result) {
    annot_name <- names(result$enrich)[2]
    tibble(annotation=annot_name,
           enrich_baseline=result$enrich[1],
           enrich_annot=result$enrich[2],
           prior_odds=(1+exp( result$enrich[1] )) /
             (1+exp( result$enrich[1]+result$enrich[2] )),
           log_post=result$log_post)
  }) %>% do.call(what=rbind)

  L0 <- result_null$log_post
  row.names(model_summary) <- NULL

  model_summary$p <- (2*model_summary$log_post - 2*L0) %>%
    pchisq(df=2, lower.tail = FALSE)
  model_summary$FDR <- p.adjust(model_summary$p, method="fdr")

  model_summary$n_annot <- vapply(annot_inputs[seq(2,n)], function(annot) {
    lapply(annot, function(x) sum(x[,2])) %>%
      do.call(what=sum)
  }, numeric(1))
  n_total <- vapply(annot_inputs[[1]], nrow, integer(1)) %>% sum()
  model_summary$frac_annot <- model_summary$n_annot / n_total


  return(model_summary %>% arrange(p))
}

