#' Process a single chromosome for UPD detection
#'
#' This internal helper function runs the full UPD-calling pipeline for a single
#' chromosome. It applies the Viterbi algorithm to infer hidden states from trio
#' genotype observations, converts the resulting Viterbi-annotated VCF into a
#' structured per-variant data.frame, summarizes consecutive variants 
#' with identical inferred states into blocks, counts Mendelian-inconsistent
#' sites per block, and optionally computes block-level read-depth ratios.
#'
#'
#' @param vcf_chr A CollapsedVCF object representing a single chromosome.
#'
#' @param hmm A Hidden Markov Model object.
#' 
#' @param add_ratios Logical; default = FALSE.
#' 
#' If TRUE, computes normalized per-block read depth ratios for each individual based on total mean depth.
#'
#' @param field_DP Optional character string specifying which VCF FORMAT field to use for depth metrics (e.g., DP, AD, or a custom field). 
#'
#' @param total_sum Optional numeric vector containing the total read depth per
#'   sample across the entire VCF, computed by \code{computeTrioTotals()} in 
#'   \code{calculateEvents()}. Used in combination with 
#'   \code{total_valid} to normalize per-block depth ratios.
#'
#' @param total_valid Optional numeric vector containing the total number of 
#'   valid (non-missing) depth observations per sample across the entire VCF, 
#'   also computed by \code{computeTrioTotals()}. Used in combination with 
#'   \code{total_sum} to normalize per-block depth ratios.
#'   
#' @param sum_cols Character vector with names of columns in the block data.frame
#'   that store the sum of depth metrics inside each block.  
#'
#' @param count_cols Character vector with names of columns in the block 
#'   data.frame storing the count of valid depth positions inside each block.  
#'   
#' @param ratio_cols Character vector containing the names of the output columns
#'   in which to store block-level inside/outside depth ratios.  
#'
#' @param mendelian_error_values Character vector of genotype codes considered
#'   Mendelian errors (i.e., observations with minimal emission probability in 
#'   the "normal" state).  
#'   Provided by \code{calculateEvents()}.
#'
#' @return A data.frame summarizing blocks detected on the chromosome. Columns include:
#' \itemize{
#'   \item `seqnames` – chromosome name
#'   \item `start`, `end` – genomic coordinates of the block
#'   \item `group` – inferred HMM state
#'   \item `n_snps` – number of SNPs in the block
#'   \item `n_mendelian_error` – number of Mendelian-inconsistent genotypes in the block
#'   \item depth-ratio metrics (if add_ratios = TRUE)
#' }
#'
#' If an error occurs during processing, a message is printed and `NULL` is returned.
#'
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

