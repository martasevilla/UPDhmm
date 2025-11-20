#' Process a single chromosome for UPD detection
#'
#' Internal helper function to run the full pipeline on one chromosome:
#' - applyViterbi
#' - asDfVcf
#' - blocksVcf
#'
#' @param vcf_chr CollapsedVCF object for one chromosome
#' @param hmm Hidden Markov Model object
#' @param genotypes Named vector mapping genotype strings to numeric states
#'
#' @return A data.frame of detected blocks for the chromosome, or NULL if error
processChromosome <- function(vcf_chr, hmm, add_ratios = FALSE, field_DP = NULL, total_mean = NULL, mendelian_error_values) {
  
  tryCatch({
    
    chr_name <- as.character(GenomeInfoDb::seqnames(vcf_chr)[1])
    
    #################################################
    # 1) Run Viterbi
    #################################################
    vcf_vit <- tryCatch(
      applyViterbi(largeCollapsedVcf = vcf_chr,
                   hmm = hmm),
      error = function(e) {
        stop(sprintf("[Chromosome %s] Error in applyViterbi: %s",
                     chr_name, conditionMessage(e)))
      }
    )
    if (!inherits(vcf_vit, "CollapsedVCF")) {
      stop(sprintf("[Chromosome %s] applyViterbi did not return CollapsedVCF.", chr_name))
    }
    
    #################################################
    # 2) Create blocks and optionally compute depth ratios
    #################################################
    blk <- tryCatch(
      blocksVcf(vcf_vit, add_ratios, field_DP, total_mean),
      error = function(e) {
        stop(sprintf("[Chromosome %s] Error in blocksVcf: %s",
                     chr_name, conditionMessage(e)))
      }
    )
    if (!inherits(blk, "data.frame")) {
      stop(sprintf("[Chromosome %s] blocksVcf did not return a data.frame.", chr_name))
    }
    
    #################################################
    # 3) Count Mendelian-inconsistent genotypes per block
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