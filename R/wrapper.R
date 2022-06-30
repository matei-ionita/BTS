
#' @title Wrapper for main BTS functionality.
#' @description Read input, compute and write output.
#' @param in_path path containing input data, in individual
#' directories for each locus.
#' @param out_path path for writing output.
#' @param d Positive integer, maximum size of causal sets
#' considered.
#' @export
BTS <- function(in_path, out_path, d=2) {
  if (as.integer(d)!=d || d < 1)
    stop("d must be a positive integer.")
  
  loci <- list.files(in_path)
  
  configs <- compute_configs(in_path, loci, d=d)
  annot_inputs <- get_annot_inputs(in_path, loci)
  results <- lapply(annot_inputs, run_em, configs=configs)
  model_summary <- get_model_summary(results, annot_inputs)
  
  write_output(results, model_summary, loci, out_path)
}

