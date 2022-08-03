
#' @title Reads and organizes functional annotations.
#' @description Loci must have corresponding directories
#' in the provided path. Assuming each locus directory
#' contains annotation files in the format
#' annotations_[group].tsv. [group] could be a tissue
#' category, e.g. Brain or Blood.
#' @param path Path to annotations.
#' @param loci A vector of loci.
#' @export
get_annot_inputs <- function(path, loci) {
  path <- paste0(path,"/")
  groups <- get_annot_groups(path, loci[1])
  annot_inputs <- lapply(groups, parse_annot_inputs,
                         loci=loci, path=path) %>%
    do.call(what=c)
  
  # Add a null model: uniform priors, no annotation
  annot_inputs <- add_null_input(annot_inputs)
  return(annot_inputs)
}


get_variant_names <- function(path, loci) {
  path <- paste0(path,"/")
  groups <- get_annot_groups(path, loci[1])
  variant_names <- lapply(loci, var_names_locus, group = groups[1], path=path)
  return(variant_names)
}


var_names_locus <- function(locus, group, path) {
  file <- paste0(path, locus, "/annotations_", group, ".tsv")
  info <- read_delim(file, delim=" ", col_types = cols(.default="i"), progress=FALSE)
  ##### FIX THIS #####
}

read_ld_matrix <- function(path, locus, reg=0) {
  file <- paste0(path, "/", locus, "/ld.txt")
  ld_matrix <- read_delim(file, delim=" ", col_types = cols(), 
          col_names = FALSE, progress = FALSE) %>% as.matrix()

  # Regularization: not necessary in our framework,
  # can set reg=0
  ld_matrix <- ld_matrix + diag(reg, nrow=nrow(ld_matrix))
  return(ld_matrix)
}

read_z_scores <- function(path, locus) {
  file <- paste0(path, "/", locus, "/z.txt")
  z_score <- read_delim(file, delim=" ", col_types = cols(),
                        progress=FALSE)[[1]]

  return(z_score)
}


get_annot_groups <- function(path, locus) {
  files <- list.files(path=paste0(path, locus),
                      pattern="annotations")
  
  keep <- vapply(files, function(file) {
    x <- read_delim(paste0(path, locus, "/", file), col_types = cols(), progress=FALSE)
    all(vapply(x, is.numeric, logical(1)))
  }, logical(1))

  if(!all(keep))
    warning("Skipping annotation files containing non-numeric values.")

  files[keep] %>%  
  str_remove("annotations_") %>%
    str_remove(".tsv")
}
  


read_annot <- function(path, locus, group) {
  file <- paste0(path, locus, "/annotations_", group, ".tsv")
  read_delim(file, delim=" ", col_types = cols(.default="i"), progress=FALSE)
}


parse_annot_inputs <- function(path, loci, group) {
  message(group)
  
  annots <- lapply(loci, FUN=read_annot, path=path, group=group)
  names(annots) <- loci

  # Count variants which have each annotation.
  # Remove annotations with count of 0: equivalent to null model
  annot_names <- names( annots[[1]] )
  annot_count <- vapply(annot_names, function(annot_name) {
    lapply(annots, function(annot) sum(annot[,annot_name])) %>%
      do.call(what=sum)
  }, numeric(1))
  annot_names <- annot_names[which(annot_count>0)]

  annot_in <- lapply(annot_names, function(annot_name) {
    lapply(annots, select_and_pad, annot_names=annot_name)
  })
  return(annot_in)
}

# add a baseline annotation which is 1 for all variants
select_and_pad <- function(input_locus, annot_names) {
  baseline <- matrix(1L, nrow=nrow(input_locus), ncol=1)
  annot_locus <- cbind(baseline, as.integer(as.matrix(input_locus[,annot_names,drop=FALSE])))
  colnames(annot_locus) <- c("baseline", annot_names)
  return(annot_locus)
}


add_null_input <- function(annot_inputs) {
  null_input <- annot_inputs[[1]]
  for(locus in names(null_input)) {
    null_input[[locus]] <- null_input[[locus]][,1,drop=FALSE]
  }
  annot_inputs <- c(list(null_input), annot_inputs)
  return(annot_inputs)
}



write_output <- function(results, model_summary, loci, out_path, pval_cutoff = 0.5) {
  # Write model summary
  out_fn <- paste0(out_path, "/model_summary.tsv")
  write_tsv(model_summary, out_fn, progress=FALSE)
  
  # Write locus likelihood by model
  model_names <- c("Null", sapply(results[seq(2,length(results))], 
                                  function(result) names(result$enrich)[2] ))
  
  log_lik <- lapply(results, function(result) unlist(result$log_post_locus)) %>%
    do.call(what=cbind)
  colnames(log_lik) <- model_names

  log_lik_diff <- sweep( log_lik , 1, log_lik[,"Null"])
  pvals <- pchisq(2*log_lik_diff, df=2, lower.tail = FALSE) %>% 
    as_tibble() %>%
    mutate(locus = loci) %>%
    relocate(locus)
  
  log_lik <- log_lik %>%
    as_tibble() %>%
    mutate(locus = loci) %>%
    relocate(locus)
  
  out_fn <- paste0(out_path, "/locus_log_lik_by_model.tsv")
  write_tsv(log_lik, out_fn, progress=FALSE)
  out_fn <- paste0(out_path, "/locus_pval_by_model.tsv")
  write_tsv(pvals, out_fn, progress=FALSE)
  
  # Write variant posterior probabilities
  post_tall <- lapply(seq_along(loci), function(i) {
    
    post_locus <- lapply(results, function(result) {
      return(result$posteriors[[i]])
    }) %>% do.call(what=cbind)
    
    colnames(post_locus) <- model_names
    post_locus <- as_tibble(post_locus)
    
    out_fn <- paste0(out_path, "/variant_posteriors_", loci[i], ".tsv")
    write_tsv(post_locus, out_fn, progress=FALSE)
    
    post_locus_tall <- lapply(model_names, function(mod) get_credible_set(pull(post_locus, mod))) %>%
      do.call(what=rbind)
    post_locus_tall$locus <- loci[i]
    post_locus_tall$locus_size <- nrow(post_locus)
    post_locus_tall$model <- model_names

    return(post_locus_tall)
  }) %>% do.call(what=rbind)
  
  # Write locus prioritization
  p_tall <- pvals %>%
    pivot_longer(-locus, names_to="model", values_to="pval") %>%
    inner_join(post_tall, by=c("locus", "model")) %>%
    arrange(pval, model)
  
  out_fn <- paste0(out_path, "/locus_prioritization.tsv")
  write_tsv(p_tall, out_fn, progress=FALSE)
  
}


get_credible_set <- function(probs, d, frac=0.8) {
  s <- sum(probs)
  or <- order(probs, decreasing = TRUE)
  
  cumsum <- 0
  i <- 0
  while(cumsum < frac * s) {
    i <- i+1
    cumsum <- cumsum + probs[or[i]]
  }
  
  return(tibble(expected_n_causal = s, credible_set_size = i, 
                credible_set = paste(or[seq_len(i)], collapse = ",")))
}



