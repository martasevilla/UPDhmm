#' Convert a collapsed VCF into a structured data.frame with inferred states
#'
#' This function transforms a large collapsed VCF object into a structured `data.frame` 
#' that contains predicted states for each variant, along with genomic coordinates, 
#' genotype codings, and optional depth/quality metrics derived from DP or AD fields.
#'
#' @param largeCollapsedVcf Name of the large collapsed VCF file.
#' @param add_ratios Logical; if `TRUE`, the function extracts depth/quality information 
#'   for the proband, mother, and father, using DP or AD fields from the VCF.
#' @param field_DP Optional string specifying the name of an alternative DP-like field 
#'   to use instead of the standard DP field.
#'
#' @return data.frame

asDfVcf <- function(largeCollapsedVcf, add_ratios = FALSE, field_DP = NULL) {

  ## Extract column metadata and variant metadata
  colData_vcf <- SummarizedExperiment::colData(largeCollapsedVcf)
  mcols_vcf <- S4Vectors::mcols(largeCollapsedVcf)

  ## Extract genomic coordinates
  start_pos <- GenomicRanges::start(largeCollapsedVcf) 
  end_pos <- GenomicRanges::end(largeCollapsedVcf)
  seqnames_chr <- as.character(GenomicRanges::seqnames(largeCollapsedVcf))
  geno_coded <- mcols_vcf$geno_coded
  
  ## --------------------------------------------------------------
  ## Build the base data.frame with variant info and inferred states
  ## --------------------------------------------------------------
  dt <- data.frame(
    ID         = colData_vcf$ID[1],
    start      = start_pos,
    end        = end_pos,
    group      = mcols_vcf$states,
    seqnames   = seqnames_chr,
    geno_coded = geno_coded
  )
  
  ## --------------------------------------------------------------
  ## Optional: incorporate depth or allele-depth–derived metrics
  ## --------------------------------------------------------------
  if (add_ratios) {
    geno_list <- VariantAnnotation::geno(largeCollapsedVcf)
    expected_samples <- c("proband", "mother", "father")
    
    ## Preferentially use user-specified DP-like field if present
    if (!is.null(field_DP) && field_DP %in% names(geno_list)) {
      quality_matrix <- as.matrix(geno_list[[field_DP]])[, expected_samples, drop = FALSE]

    ## Otherwise fallback to DP if available
    } else if ("DP" %in% names(geno_list)) {
      quality_matrix <- as.matrix(geno_list$DP)[, expected_samples, drop = FALSE]
    
    ## If no DP is present, reconstruct per-sample depth from AD
    } else if ("AD" %in% names(geno_list)) {
      ad_array <- geno_list$AD
      n_variants <- dim(ad_array)[1]
      n_samples  <- dim(ad_array)[2]

      ## Collapse allele depths across alleles for each sample
      quality_matrix <- matrix(
        colSums(matrix(ad_array, nrow = n_variants * dim(ad_array)[3])), 
        nrow = n_variants, ncol = n_samples
      )[, expected_samples, drop = FALSE]

    ## No suitable quality metric found → fill with NA
    } else {
      quality_matrix <- matrix(
        NA_real_, 
        nrow = length(start_pos), 
        ncol = length(expected_samples),
        dimnames = list(NULL, expected_samples)
      )
      warning("No DP or AD field found in VCF. Quality columns set to NA.")
    }
    
    ## Convert quality matrix to data.frame and merge with base 
    quality_dt <- data.frame(
      quality_proband = as.numeric(quality_matrix[, "proband"]),
      quality_mother  = as.numeric(quality_matrix[, "mother"]),
      quality_father  = as.numeric(quality_matrix[, "father"]),
      stringsAsFactors = FALSE
    )

    dt <- cbind(dt, quality_dt)
  }

  return(dt)
}
