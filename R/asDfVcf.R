#' Convert a collapsed VCF into a structured data.frame with inferred states
#'
#' This function transforms a CollapsedVCF object into a
#' structured data.frame containing genomic coordinates, inferred hidden
#' states, genotype codings, and optional read depth metrics.
#'
#'
#' @param largeCollapsedVcf A CollapsedVCF object containing variants
#'   already processed with \code{vcfCheck()} and annotated with inferred states
#'   via \code{applyViterbi()}.
#'
#' @param add_ratios Logical; default FALSE.  
#'   If \code{TRUE}, the function extracts depth or allele-depth–derived metrics
#'   for the trio using the field specified in field_DP
#'
#' @param field_DP Optional character string specifying which VCF FORMAT field to use for depth metrics (e.g., DP, AD, or a custom field). 
#'
#' @return A data.frame with one row per variant, including:
#' \itemize{
#'   \item ID: sample identifier
#'   \item start, end: genomic coordinates of the variant
#'   \item seqnames: chromosome name
#'   \item group: inferred HMM hidden state
#'   \item geno_coded: numeric genotype coding for the trio
#' }
#'
#' If add_ratios = TRUE, the following additional columns are included:
#' \itemize{
#'   \item quality_proband, quality_mother, quality_father:
#'         per-variant depth values or summed allele-depth metrics.
#' }
#' If the requested depth field is unavailable, these columns are returned as \code{NA}.
#' 
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
