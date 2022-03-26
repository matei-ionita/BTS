
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
                        col_names=FALSE, progress=FALSE)[[1]]
  return(z_score)
}


get_annot_groups <- function(path, locus) {
  list.files(path=paste0(path, locus),
             pattern="annotations") %>%
    str_remove("annotations_") %>%
    str_remove(".tsv")
}


read_annot <- function(path, locus, group) {
  file <- paste0(path, locus, "/annotations_", group, ".tsv")
  read_tsv(file, col_types = cols(.default="i"), progress=FALSE)
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



