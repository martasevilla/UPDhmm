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
    
    ## Decide which field to use
    dp_field <- if (!is.null(field_DP) && field_DP %in% names(geno_list)) { field_DP 
                } else if ("DP" %in% names(geno_list)) { "DP" 
                } else if ("AD" %in% names(geno_list)) { "AD" 
                } else { NULL }
    
    if (!is.null(dp_field)) {
      if (dp_field == "AD") {
        quality_matrix <- apply(geno_list$AD, 2, function(col) {
          vapply(col, function(x) {if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)}, numeric(1))
        })
        
      } else {
        quality_matrix <- as.matrix(geno_list[[dp_field]])
      }
      
      ## Ensure order proband, mother, father
      quality_matrix <- quality_matrix[, expected_samples, drop = FALSE]
      
      ## Convert to numeric safely
      quality_dt <- data.frame(
        quality_proband = as.numeric(quality_matrix[, "proband"]),
        quality_mother  = as.numeric(quality_matrix[, "mother"]),
        quality_father  = as.numeric(quality_matrix[, "father"]),
        stringsAsFactors = FALSE
      )
      
      ## Bind to base dataframe
      dt <- cbind(dt, quality_dt)
      
    } else {
      ## No DP or AD field
      warning("No DP or AD field found in VCF. Quality columns set to NA.")
      dt$quality_proband <- dt$quality_mother <- dt$quality_father <- NA_real_
    }
  }
  
  return(dt)
}
