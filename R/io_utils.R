
#' @title Reads and organizes functional annotations.
#' @description Currently assumes annotations are organized by tissue.
#' @param in_path Path to annotations.
#' @param loci A vector of loci.
#' @export
get_annot_inputs <- function(in_path, loci) {
  tissues <- get_tissues(in_path, loci)
  annot_inputs <- lapply(tissues, get_tissue_inputs,
                         loci=loci, in_path=in_path) %>%
    do.call(what=c)
  annot_inputs <- add_null_input(annot_inputs)
  return(annot_inputs)
}


read_ld_matrix <- function(ld_file, path, reg=1e-5) {
  # file <- paste0(path, "IBD.", locus, ".ld")
  file <- paste0(path, ld_file)
  ld_matrix <- read_delim(file, delim=" ", col_types = cols(), col_names = FALSE) %>%
    as.matrix()
  
  # regularization as done by paintor: default 1e-5
  ld_matrix <- ld_matrix + diag(reg, nrow=nrow(ld_matrix))
  return(ld_matrix)
}

read_z_scores <- function(locus, path, tissue="Tongue") {
  inputs <- read_inputs(locus, path, tissue=tissue)
  z_score <- inputs$normalized_beta / inputs$SE
  return(z_score)
}

# read_z_scores <- function(locus, path) {
#   z_score <- read_tsv(paste0(path, locus))$Zscore
#   return(z_score)
# }

read_inputs <- function(locus, path, tissue) {
  tissue_rep <- str_replace(tissue, "_", " ")
  locus <- str_remove(locus, "kunkle.")
  file <- paste0(path, "filer_overlap/", locus, "/annot_matrix/annot_matrix.", tissue_rep, ".annot")
  inputs <- read_tsv(file, col_types = cols())
  return(inputs)
}



get_tissues <- function(in_path, loci) {
  loci <- str_remove(loci, "kunkle.")
  tissues <- list.files(path=paste0(in_path, "filer_overlap/",loci[1]),
                        pattern=".annot", recursive=TRUE) %>%
    str_remove("annot_matrix/annot_matrix.") %>%
    str_remove(".annot")
  return(setdiff(tissues, c("Not applicable", "Tissue category")))
}


get_tissue_inputs <- function(tissue, loci, in_path) {
  message(tissue)
  inputs <- lapply(loci, FUN=read_inputs, path=in_path, tissue=tissue)
  names(inputs) <- loci

  annot_tissue <- grep(tissue, names(inputs[[1]]), value=TRUE)
  annot_count <- vapply(annot_tissue, function(annot_name) {
    lapply(inputs, function(input) sum(input[,annot_name])) %>%
      do.call(what=sum)
  }, numeric(1))
  annot_tissue <- annot_tissue[which(annot_count>0)]

  annots <- lapply(annot_tissue, function(annot_name) {
    lapply(inputs, select_and_pad, annot_names=annot_name)
  })
  return(annots)
}


select_and_pad <- function(input_locus, annot_names) {
  baseline <- matrix(1, nrow=nrow(input_locus), ncol=1)
  annot_locus <- cbind(baseline, input_locus[,annot_names,drop=FALSE])
  return(as.matrix(annot_locus))
}



get_annot_inputs_example <- function(loci, in_path) {
  inputs <- lapply(loci, function(locus) {
    read_delim(paste0(in_path, locus, ".annotations"), delim=" ")
  })
  names(inputs) <- loci

  annots <- lapply(inputs, function(input) {
    baseline <- matrix(1, nrow=nrow(input), ncol=1)
    as.matrix(cbind(baseline, input[,"DHS"]))
  })
  return(list(annots))
}


add_null_input <- function(annot_inputs) {
  null_input <- annot_inputs[[1]]
  for(locus in names(null_input)) {
    null_input[[locus]] <- null_input[[locus]][,1,drop=FALSE]
  }
  annot_inputs <- c(list(null_input), annot_inputs)
  return(annot_inputs)
}



