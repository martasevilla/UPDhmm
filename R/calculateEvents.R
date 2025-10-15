#' Calculate UPD events in trio VCFs.
#'
#' This function predicts the hidden states by applying the Viterbi algorithm
#' using the Hidden Markov Model (HMM) from the UPDhmm package. It takes the
#' genotypes of the trio as input and includes a final step to simplify the
#' results into blocks.
#'
#' @param largeCollapsedVcf The VCF file in the general format 
#' (largeCollapsedVcf) with VariantAnnotation package. Previously edited with 
#' `vcfCheck()` function from UPDhmm package.
#'
#' @param hmm Default = `NULL`. If no arguments are added, the package 
#' will use the default HMM already implemented, based on Mendelian 
#' inheritance. If an optional HMM is desired, it should adhere to the 
#' general HMM format from `HMM` package with the following elements inside 
#' a list:
#'   1. The hidden state names in the "States" vector.
#'   2. All possible observations in the "Symbols" vector.
#'   3. Start probabilities of every hidden state in the "startProbs" vector.
#'   4. Transition probabilities matrix between states in "transProbs".
#'   5. Probabilities associated between every hidden state and all possible 
#'      observations in the "emissionProbs" matrix.
#'
#' @param field_DP Default = `NULL`. Character string specifying which FORMAT field in the VCF
#' contains the read depth information to use in `addRatioDepth()`.
#' If `NULL` (default), the function will automatically try `"DP"` (standard depth)
#' or `"AD"` (allelic depths, summed across alleles).
#' Use this parameter if your VCF uses a non-standard field name for depth,
#' e.g. `field = "NR"` or `"field_DP"`.
#'
#' @param BPPARAM Parallelization settings, passed to
#'   \link[BiocParallel]{bplapply}.
#'   By default `BiocParallel::SerialParam()` (serial execution).
#'   To enable parallelization, provide a BiocParallel backend, e.g.
#'   `BiocParallel::MulticoreParam(workers = min(2, parallel::detectCores()))`
#'   or `BiocParallel::SnowParam(workers = 2)`.
#'   Note: when running under R CMD check or Bioconductor build systems,
#'   the number of workers may be automatically limited (usually less or equal to 2).
#'
#' @param verbose Logical, default = `FALSE`. 
#'   If `TRUE`, progress messages will be printed during processing.
#'
#' @return A `data.frame` object containing all detected events in the provided trio. 
#' If no events are found, the function will return an empty `data.frame`.
#'
#' @export
#'
#' @examples
#' file <- system.file(package = "UPDhmm", "extdata", "test_het_mat.vcf.gz")
#' vcf <- VariantAnnotation::readVcf(file)
#' processedVcf <- vcfCheck(vcf,
#'     proband = "NA19675", 
#'     mother = "NA19678",
#'     father = "NA19679"
#' )
#'
#' # Run in serial mode (default)
#' res <- calculateEvents(processedVcf)
#'
#' # Run in parallel with 2 cores
#' library(BiocParallel)
#' param <- MulticoreParam(workers = 2)
#' res_parallel <- calculateEvents(processedVcf, BPPARAM = param)
calculateEvents <- function(largeCollapsedVcf,
                            hmm = NULL,
                            field_DP = NULL,
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose = FALSE) {
  # 0. Check input
  if (!inherits(largeCollapsedVcf, "CollapsedVCF")) {
    stop("Argument 'largeCollapsedVcf' must be a CollapsedVCF object.")
  }
  
  if (is.null(hmm)) {
    utils::data("hmm", package = "UPDhmm", envir = environment())
  }
  
  genotypes <- c(
    "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
    "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
  )
  
  # 1. Split VCF into chromosomes
  split_vcf_raw <- split(largeCollapsedVcf,
                         f = GenomicRanges::seqnames(largeCollapsedVcf))
  split_vcf_raw <- split_vcf_raw[lengths(split_vcf_raw) > 0]
  
  if (length(split_vcf_raw) == 0) {
    if (verbose) message("No chromosomes found in VCF.")
    return(data.frame())
  }
  
  if (verbose) message("Processing ", length(split_vcf_raw), " chromosomes...")
  
  # 2. Run pipeline per chromosome (serial or parallel)
  if (inherits(BPPARAM, "SerialParam")) {
    blocks_state <- lapply(split_vcf_raw, processChromosome,
                           hmm = hmm, genotypes = genotypes)
  } else {
    blocks_state <- BiocParallel::bplapply(split_vcf_raw, processChromosome,
                                           hmm = hmm, genotypes = genotypes,
                                           BPPARAM = BPPARAM)
  }
  
  
  
  # Drop NULLs  
  blocks_state <- Filter(Negate(is.null), blocks_state)
  if (length(blocks_state) == 0) {
    stop("calculateEvents failed: no valid chromosomes processed. 
       Likely cause: applyViterbi does not recognized all GT or 
       check your VCF formatting and trio sample IDs.")
  }
  
  
  # 3. Clean results
  def_blocks_states <- do.call(rbind, blocks_state)
  
  # 4. Filter events (skip normal state, low SNPs, sex chromosomes)
  filtered_def_blocks_states <- def_blocks_states[
    def_blocks_states$n_snps > 1 &
      def_blocks_states$group != "normal" &
      !(def_blocks_states$seqnames %in% c("chrX", "X")), ]
  
  if (nrow(filtered_def_blocks_states) == 0) {
    if (verbose) message("No non-normal events found.")
    return(data.frame())
  }
  
  if (verbose) {
    message("Found ", nrow(filtered_def_blocks_states), " candidate events.")
  }
  
  # 5. Add OR + depth ratios
  blocks_list <- lapply(seq_len(nrow(filtered_def_blocks_states)), function(i) {
    df_or <- addOr(filtered_def_blocks_states[i, , drop = FALSE],
                   largeCollapsedVcf, hmm, genotypes)
    df_ratio <- addRatioDepth(filtered_def_blocks_states[i, , drop = FALSE],
                              largeCollapsedVcf, field = field_DP)
    cbind(df_or, df_ratio[, setdiff(names(df_ratio), names(df_or))])
  })
  
  block_def <- do.call(rbind, blocks_list)
  rownames(block_def) <- NULL
  return(block_def)
}
