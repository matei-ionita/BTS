# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

enumerate_configs <- function(z0, ld, pr_var, d, thresh) {
    .Call(`_BTS_enumerate_configs`, z0, ld, pr_var, d, thresh)
}

compute_post_locus <- function(locus_configs, enrich0, locus_annot) {
    .Call(`_BTS_compute_post_locus`, locus_configs, enrich0, locus_annot)
}

