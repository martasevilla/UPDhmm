#' Process a single chromosome for UPD detection
#'
#' Internal helper function to run the full pipeline on one chromosome:
#' It performs the following steps:
#' 1. Applies the Viterbi algorithm to infer hidden states using the provided HMM.
#' 2. Converts the VCF to a data.frame including optional read depth/quality metrics.
#' 3. Collapses contiguous variants with the same inferred state into blocks.
#' 4. Calculates Mendelian error counts per block.
#' 5. Optionally computes per-block read depth ratios relative to the rest of the chromosome.
#'
#' @param vcf_chr CollapsedVCF object for one chromosome
#' @param hmm Hidden Markov Model object
#' @param add_ratios Logical; if `TRUE`, read depth/quality ratios will be included in the analysis.
#' @param field_DP Optional character specifying the FORMAT field in the VCF for depth metrics.
#' @param total_sum Numeric vector of total read depths per sample for the chromosome.
#' @param total_valid Numeric vector of total valid positions per sample for the chromosome.
#' @param quality_cols Character vector of sample names corresponding to columns used for read depth/quality metrics (default `c("proband", "mother", "father")`).
#' @param sum_cols Character vector of column names in `blk` containing per-block summed depth/quality.
#' @param count_cols Character vector of column names in `blk` containing per-block counts of valid positions.
#' @param ratio_cols Character vector of column names to store inside/outside block ratios.
#' @param mendelian_error_values Character vector of genotype codes considered Mendelian errors, based on the HMM emissions.
#' 
#' @return A data.frame of detected blocks for the chromosome, or NULL if error

processChromosome <- function(vcf_chr, hmm, add_ratios = FALSE, field_DP = NULL, total_sum = NULL, total_valid = NULL,
                              quality_cols = c("proband", "mother", "father"), sum_cols = c("total_sum_quality_proband", "total_sum_quality_mother", "total_sum_quality_father"), 
                              count_cols = c("total_count_quality_proband", "total_count_quality_mother", "total_count_quality_father") , ratio_cols = c("ratio_proband", "ratio_mother", "ratio_father"), 
                              mendelian_error_values) {
  tryCatch({
    
    chr_name <- as.character(GenomeInfoDb::seqnames(vcf_chr)[1])
    
    #################################################
    # 1. Run Viterbi to infer hidden states
    #################################################
    vcf_vit <- applyViterbi(vcf_chr, hmm)
      
    if (!inherits(vcf_vit, "CollapsedVCF")) {
      stop(sprintf("[Chromosome %s] applyViterbi did not return CollapsedVCF.", chr_name))
    }
    
    #################################################
    # 2. Convert VCF to data.frame with optional depth/quality metrics
    #################################################
    df_vit <- asDfVcf(vcf_vit, add_ratios, field_DP)

    if (!inherits(df_vit, "data.frame")) {
      stop(sprintf("[Chromosome %s] asDfVcf did not return a data.frame.", chr_name))
    }
    
    #################################################
    # 3. Collapse contiguous variants with the same inferred state into blocks
    #################################################
    blk <- blocksVcf(df_vit)
    
    if (!inherits(blk, "data.frame")) {
      stop(sprintf("[Chromosome %s] blocksVcf did not return a data.frame.", chr_name))
    }

    #################################################
    # 4. Count Mendelian-inconsistent genotypes per block
    #################################################
    if (!is.null(hmm)) {
      
      # Extract the pre-computed numeric genotype strings from the VCF
      geno_coded <- S4Vectors::mcols(vcf_chr)$geno_coded
      positions <- GenomicRanges::start(vcf_chr)

      # Create IRanges for positions with Mendelian errors
      snp_error_gr <- IRanges::IRanges(
        start = positions[geno_coded %in% mendelian_error_values],
        end   = positions[geno_coded %in% mendelian_error_values]
      )
      
      # Create IRanges for block boundaries
      blk_gr <- IRanges::IRanges(start = blk$start, end = blk$end)
      
      # Count number of Mendelian errors overlapping each block
      blk$n_mendelian_error <- IRanges::countOverlaps(blk_gr, snp_error_gr)
      
      blk$geno_coded <- NULL
      
    } else {
      # If no HMM provided, fill with NA
      blk$n_mendelian_error <- NA_integer_
    }

    #################################################
    # 5. Optional: compute per-block read depth ratios
    #################################################
    if (!is.null(total_sum) & !is.null(total_valid)) {

      # Extract inside-block sums and counts
      inside_sum   <- as.matrix(blk[, sum_cols, drop = FALSE])
      inside_count <- as.matrix(blk[, count_cols, drop = FALSE])

      # Compute outside-block sums and counts: chromosome total minus inside-block values
      outside_sum   <- matrix(total_sum, nrow = nrow(blk), ncol = length(quality_cols), byrow = TRUE) - inside_sum
      outside_count <- matrix(total_valid, nrow = nrow(blk), ncol = length(quality_cols), byrow = TRUE) - inside_count

      # Compute mean depth inside and outside blocks (avoid division by zero)
      inside_mean  <- inside_sum / pmax(inside_count, 1)
      outside_mean <- outside_sum / pmax(outside_count, 1)

      # Compute ratio: inside / outside
      ratio_mat <- inside_mean / outside_mean
      colnames(ratio_mat) <- ratio_cols

      # Add ratios to blocks data.frame
      blk[ratio_cols] <- as.data.frame(ratio_mat)
      
      blk <- blk[, setdiff(colnames(blk), c(sum_cols, count_cols)), drop = FALSE]
    }
    rownames(blk) <- NULL
    return(blk)
   
  }, error = function(e) {
    message("processChromosome failed: ", conditionMessage(e))
    return(NULL)
  })
}

