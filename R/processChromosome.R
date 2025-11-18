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
#' @param total_mean   
#' @param mendelian_error_values Character vector of genotype codes considered Mendelian errors, based on the HMM emissions.
#' 
#' @return A data.frame of detected blocks for the chromosome, or NULL if error

processChromosome <- function(vcf_chr, hmm, add_ratios = FALSE, field_DP = NULL, total_mean = NULL, mendelian_error_values) {
  
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
    # 2. Convert VCF to blocks and optionally compute ratios directly
    #################################################
    blk <- blocksVcfNew(vcf_vit, add_ratios, field_DP, total_mean)

    if (!inherits(blk, "data.frame")) {
      stop(sprintf("[Chromosome %s] blocksVcfNew did not return a data.frame.", chr_name))
    }
    
    #################################################
    # 3. Count Mendelian-inconsistent genotypes per block
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

    rownames(blk) <- NULL
    return(blk)
   
  }, error = function(e) {
    message("processChromosome failed: ", conditionMessage(e))
    return(NULL)
  })
}
