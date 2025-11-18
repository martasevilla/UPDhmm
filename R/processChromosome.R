#' Process a single chromosome for UPD detection
#'
#' Internal helper function to run the full pipeline on one chromosome:
#' It performs the following steps:
#' 1. Runs the Viterbi algorithm to assign inferred hidden states to all variants.
#' 2. Converts the Viterbi-annotated VCF into a block-level data.frame, optionally adding depth-based ratios.
#' 3. Counts Mendelian-inconsistent genotypes within each block using the pre-computed `geno_coded` values.
#' 
#' @param vcf_chr CollapsedVCF object for one chromosome
#' @param hmm Hidden Markov Model object
#' @param add_ratios Logical; if `TRUE`, read depth/quality ratios will be included in the analysis.
#' @param field_DP Optional character specifying the FORMAT field in the VCF for depth metrics.
#' @param total_mean Optional numeric value representing the chromosome-wide mean depth for ratio calculations.
#' @param mendelian_error_values Character vector of genotype codes considered Mendelian errors, based on the HMM emissions.
#' 
#' @return A data.frame containing the inferred blocks and associated metrics for the chromosome,
#'         or `NULL` if any processing step fails.
processChromosome <- function(vcf_chr, hmm, add_ratios = FALSE, field_DP = NULL, total_mean = NULL, mendelian_error_values) {
  
  tryCatch({
    
    chr_name <- as.character(GenomeInfoDb::seqnames(vcf_chr)[1])
    
    #################################################
    # 1. Run Viterbi to infer hidden states for all variants
    #################################################
    vcf_vit <- applyViterbi(vcf_chr, hmm)
      
    if (!inherits(vcf_vit, "CollapsedVCF")) {
      stop(sprintf("[Chromosome %s] applyViterbi did not return a CollapsedVCF object.", chr_name))
    }
    
    #################################################
    # 2. Build block-level representation and optionally compute depth ratios
    #################################################
    blk <- blocksVcfNew(vcf_vit, add_ratios, field_DP, total_mean)

    if (!inherits(blk, "data.frame")) {
      stop(sprintf("[Chromosome %s] blocksVcfNew did not return a data.frame.", chr_name))
    }
    
    #################################################
    # 3. Count Mendelian-inconsistent genotypes per block
    #################################################
    if (!is.null(hmm)) {
      
      # Pre-computed genotype codes (character representation) for each variant
      geno_coded <- S4Vectors::mcols(vcf_chr)$geno_coded
      
      # Genomic positions of all variants
      positions <- GenomicRanges::start(vcf_chr)

      # IRanges marking each variant that is a Mendelian error
      snp_error_gr <- IRanges::IRanges(
        start = positions[geno_coded %in% mendelian_error_values],
        end   = positions[geno_coded %in% mendelian_error_values]
      )
      
      # IRanges corresponding to block boundaries
      blk_gr <- IRanges::IRanges(start = blk$start, end = blk$end)
      
      # Count overlapping Mendelian-error positions within each block
      blk$n_mendelian_error <- IRanges::countOverlaps(blk_gr, snp_error_gr)
      
      # Clean up internal field if present
      blk$geno_coded <- NULL
      
    } else {
      # No HMM provided → no Mendelian error model available
      blk$n_mendelian_error <- NA_integer_
    }

    rownames(blk) <- NULL
    return(blk)
   
  }, error = function(e) {
    message("processChromosome failed: ", conditionMessage(e))
    return(NULL)
  })
}