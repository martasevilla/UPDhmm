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
#'
#' @return A data.frame of detected blocks for the chromosome, or NULL if error

processChromosome <- function(vcf_chr, hmm, add_ratios = FALSE, field_DP = NULL, total_sum = NULL, total_valid = NULL) {
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
    # 3. Collapse contiguous variants into blocks
    #################################################
    blk <- blocksVcf(df_vit)
    
    if (!inherits(blk, "data.frame")) {
      stop(sprintf("[Chromosome %s] blocksVcf did not return a data.frame.", chr_name))
    }

    #################################################
    # 4. Calculate Mendelian error counts per block
    #################################################
    if (!is.null(hmm)) {
      
      # Identify genotype values considered Mendelian errors based on emission probabilities
      emission_probs <- hmm$emissionProbs["normal", ]
      mendelian_error_values <- names(emission_probs[emission_probs == min(emission_probs)])

      # Use the pre-existing genotype codings from the VCF
      geno_coded <- S4Vectors::mcols(vcf_chr)$geno_coded
      positions <- GenomicRanges::start(vcf_chr)

      # Count number of Mendelian-inconsistent genotypes within each block
      blk$n_mendelian_error <- vapply(seq_len(nrow(blk)), function(i) {
        idx <- positions >= blk$start[i] & positions <= blk$end[i]
        sum(geno_coded[idx] %in% mendelian_error_values)
      }, integer(1))

      blk$geno_coded <- NULL

    } else {
      # If no HMM provided, fill with NA
      blk$n_mendelian_error <- NA_integer_
    }

    #################################################
    # 5. Optional: compute per-block read depth ratios
    #################################################
    quality_cols <- c("proband", "mother", "father")
    if (!is.null(total_sum) & !is.null(total_valid)) {
      
      # Columns in blk with sum and count of quality/depth values per block
      sum_cols   <- paste0("total_sum_quality_", quality_cols)
      count_cols <- paste0("total_count_quality_", quality_cols)

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
      colnames(ratio_mat) <- paste0("ratio_", quality_cols)

      # Add ratios to blocks data.frame
      blk <- cbind(blk, as.data.frame(ratio_mat))
      
      keep_cols <- setdiff(colnames(blk), c(sum_cols, count_cols))
      blk <- blk[, keep_cols, drop = FALSE]
    }
    rownames(blk) <- NULL
    return(blk)
   
  }, error = function(e) {
    message("processChromosome failed: ", conditionMessage(e))
    return(NULL)
  })
}

