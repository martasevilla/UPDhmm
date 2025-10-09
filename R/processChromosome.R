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
processChromosome <- function(vcf_chr, hmm, genotypes) {
  tryCatch({
    
    chr_name <- as.character(GenomeInfoDb::seqnames(vcf_chr)[1])
    
    #################################################
    # 1) Run Viterbi
    #################################################
    vcf_vit <- tryCatch(
      applyViterbi(largeCollapsedVcf = vcf_chr,
                   hmm = hmm, genotypes = genotypes),
      error = function(e) {
        stop(sprintf("[Chromosome %s] Error in applyViterbi: %s",
                     chr_name, conditionMessage(e)))
      }
    )
    if (!inherits(vcf_vit, "CollapsedVCF")) {
      stop(sprintf("[Chromosome %s] applyViterbi did not return CollapsedVCF.", chr_name))
    }
    
    #################################################
    # 2) Convert to dataframe
    #################################################
    df_vit <- tryCatch(
      asDfVcf(largeCollapsedVcf = vcf_vit, genotypes = genotypes),
      error = function(e) {
        stop(sprintf("[Chromosome %s] Error in asDfVcf: %s",
                     chr_name, conditionMessage(e)))
      }
    )
    if (!inherits(df_vit, "data.frame")) {
      stop(sprintf("[Chromosome %s] asDfVcf did not return a data.frame.", chr_name))
    }
    if (ncol(df_vit) != 6) {
      stop(sprintf("[Chromosome %s] asDfVcf dataframe has %d columns, expected 6.",
                   chr_name, ncol(df_vit)))
    }
    
    #################################################
    # 3) Create blocks
    #################################################
    blk <- tryCatch(
      blocksVcf(df_vit),
      error = function(e) {
        stop(sprintf("[Chromosome %s] Error in blocksVcf: %s",
                     chr_name, conditionMessage(e)))
      }
    )
    if (!inherits(blk, "data.frame")) {
      stop(sprintf("[Chromosome %s] blocksVcf did not return a data.frame.", chr_name))
    }
    if (ncol(blk) != 6) {
      stop(sprintf("[Chromosome %s] blocksVcf dataframe has %d columns, expected 6.",
                   chr_name, ncol(blk)))
    }
    
    #################################################
    # If everything worked, return blocks
    #################################################
    return(blk)
    
  }, error = function(e) {
    message("processChromosome failed: ", conditionMessage(e))
    return(NULL)
  })
}
