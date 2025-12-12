#' Process a single chromosome for UPD detection
#'
#' Internal helper function to run the full pipeline on one chromosome:
#' - applyViterbi
#' - blocksVcf
#'
#' @param vcf_chr CollapsedVCF object for one chromosome
#' 
#' @param hmm Hidden Markov Model object
#' 
#' @param add_ratios Logical; default = FALSE.
#' 
#' @param field_DP Default = `NULL`. Character string specifying which FORMAT field in the VCF
#' contains the read depth information to use in `addRatioDepth()`.
#' If `NULL` (default), the function will automatically try `"DP"` (standard depth)
#' or `"AD"` (allelic depths, summed across alleles).
#' Use this parameter if your VCF uses a non-standard field name for depth,
#' e.g. `field = "NR"` or `"field_DP"`.
#' 
#' If TRUE, computes normalized per-block read depth ratios for each individual based on total mean depth.
#' 
#' @param total_mean Optional numeric vector of per-sample mean read depths across the entire VCF, used to normalize per-block depth ratios computed via \code{computeTrioTotals()} in \code{calculateEvents()}.
#' 
#' @param mendelian_error_values Character vector of genotype codes considered
#'   Mendelian errors (i.e., observations with minimal emission probability in 
#'   the "normal" state).  
#'   Provided by \code{calculateEvents()}.
#'
#' @return A data.frame of detected blocks for the chromosome, or NULL if error
#' Columns include:
#' \itemize{
#'   \item `chromosome` – chromosome name
#'   \item `start`, `end` – genomic coordinates of the block
#'   \item `group` – inferred HMM state
#'   \item `n_snps` – number of SNPs in the block
#'   \item `n_mendelian_error` – number of Mendelian-inconsistent genotypes in the block
#'   \item depth-ratio metrics (if add_ratios = TRUE)
#' }
#'
#' @keywords internal
#' 
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