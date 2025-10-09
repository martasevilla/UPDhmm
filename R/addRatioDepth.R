#' Add depth ratios (inside vs outside events) per sample
#'
#' This function computes the average depth inside each block and compares it
#' to the average depth outside the block, generating a ratio for each sample
#' (proband, mother, father).
#'
#' @param largeCollapsedVcf Input VCF file
#' @param filtered_def_blocks_states data.frame object containing the block
#' @param field Optional character. Read Depth field to use (e.g. "DP", "DS").
#'   If NULL, the function will try DP first, then AD as fallback.
#'
#' @return data.frame with added depth and inside/outside ratios per sample
addRatioDepth <- function(
    filtered_def_blocks_states,
    largeCollapsedVcf,
    field = NULL
) {
  # Extract depth info
  if (is.null(field)) {
    if ("DP" %in% names(VariantAnnotation::geno(largeCollapsedVcf))) {
      dp_matrix <- VariantAnnotation::geno(largeCollapsedVcf)$DP
    } else if ("AD" %in% names(VariantAnnotation::geno(largeCollapsedVcf))) {
      ad_matrix <- VariantAnnotation::geno(largeCollapsedVcf)$AD
      dp_matrix <- rowSums(ad_matrix)
    } else {
      stop("No DP or AD field found in VCF")
    }
  } else {
    if (!(field %in% names(VariantAnnotation::geno(largeCollapsedVcf)))) {
      stop(paste("Field", field, "not found in VCF"))
    }
    dp_matrix <- VariantAnnotation::geno(largeCollapsedVcf)[[field]]
  }
  
  dp_matrix <- as.matrix(dp_matrix)
  samples <- colnames(dp_matrix)
  
  if (!all(c("proband", "mother", "father") %in% samples)) {
    stop("Samples proband/mother/father not found in VCF colnames")
  }
  
  # Convert VCF to GRanges
  vcf_gr <- GenomicRanges::granges(largeCollapsedVcf)
  
  block <- filtered_def_blocks_states
  gr_block <- GenomicRanges::GRanges(
    seqnames = block$seqnames,
    ranges   = IRanges::IRanges(start = block$start, end = block$end)
  )
  
  inside_idx  <- S4Vectors::queryHits(IRanges::findOverlaps(vcf_gr, gr_block))
  outside_idx <- setdiff(seq_along(vcf_gr), inside_idx)
  
  if (length(inside_idx) == 0) {
    warning(sprintf(
      "Block %s:%s-%s has no variants inside. Ratios set to NA.",
      block$seqnames, block$start, block$end
    ))
    ratios_df <- data.frame(
      ratio_proband = NA,
      ratio_mother  = NA,
      ratio_father  = NA,
      stringsAsFactors = FALSE
    )
    
  } else if (length(outside_idx) == 0) {
    Biobase::note(sprintf(
      "Block %s:%s-%s spans entire chromosome. Ratios set to 1.",
      block$seqnames, block$start, block$end
    ))
    ratios_df <- data.frame(
      ratio_proband = 1,
      ratio_mother  = 1,
      ratio_father  = 1,
      stringsAsFactors = FALSE
    )
    
  } else {
    inside_dp  <- dp_matrix[inside_idx, , drop = FALSE]
    outside_dp <- dp_matrix[outside_idx, , drop = FALSE]
    
    mean_inside  <- colMeans(inside_dp, na.rm = TRUE)
    mean_outside <- colMeans(outside_dp, na.rm = TRUE)
    
    ratios <- mean_inside / mean_outside
    ratios_df <- data.frame(
      ratio_proband = ratios["proband"],
      ratio_mother  = ratios["mother"],
      ratio_father  = ratios["father"],
      stringsAsFactors = FALSE
    )
  }
  
  # Combine block info with ratios
  df <- cbind(block, ratios_df)
  rownames(df) <- NULL
  return(df)
}
