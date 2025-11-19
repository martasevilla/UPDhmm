#' Process a single chromosome for UPD detection
#'
#' This internal helper function runs the full UPD-calling pipeline for a single
#' chromosome. It applies the Viterbi algorithm to infer hidden states from trio
#' genotype observations, summarizes consecutive variants with identical inferred
#' states into blocks, counts Mendelian-inconsistent sites per block, and
#' optionally computes block-level read-depth ratios.
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
#' @param total_mean Optional numeric vector of per-sample mean read depths across the entire VCF, used to normalize per-block depth ratios computed via \code{computeTrioTotals()} in \code{calculateEvents()}.
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