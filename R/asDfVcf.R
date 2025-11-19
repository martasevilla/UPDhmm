#' Convert a collapsed VCF into a structured data.frame with inferred states
#'
#' This internal helper function transforms a \code{CollapsedVCF} object into a
#' structured \code{data.frame} containing genomic coordinates, inferred hidden
#' states, genotype codings, and optional read depth or allele depth metrics.
#'
#' The output is suitable for downstream processing, such as collapsing contiguous
#' variants into blocks or calculating per-block metrics.
#'
#' @param largeCollapsedVcf A \code{CollapsedVCF} object containing variants
#'   already processed with \code{vcfCheck()} and annotated with inferred states
#'   (e.g., via \code{applyViterbi()}).
#'
#' @param add_ratios Logical; default \code{FALSE}.  
#'   If \code{TRUE}, the function extracts depth or allele-depth–derived metrics
#'   for the trio (proband, mother, father), using fields \code{DP} or \code{AD} 
#'   from the VCF.
#'
#' @param field_DP Optional character string specifying the name of a non-standard
#'   depth field in the VCF to use instead of \code{DP} or \code{AD}.
#'
#' @details
#'
#' The resulting data.frame includes the following columns:
#' \itemize{
#'   \item \code{ID}: original sample ID
#'   \item \code{start}, \code{end}: genomic coordinates of the variant
#'   \item \code{seqnames}: chromosome
#'   \item \code{group}: inferred hidden state from the HMM
#'   \item \code{geno_coded}: numeric genotype coding for the trio
#' }
#'
#' If \code{add_ratios = TRUE}, additional columns are added:
#' \itemize{
#'   \item \code{quality_proband}, \code{quality_mother}, \code{quality_father}:
#'         per-variant depth or summed allele depth metrics.
#' }
#' If the requested depth field is not found, these columns are filled with \code{NA}.
#'
#' @return A \code{data.frame} with one row per variant and columns as described above.
#'
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
