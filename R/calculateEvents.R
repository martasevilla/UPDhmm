#' Calculate UPD events in trio VCFs.
#'
#' This function predicts hidden states for a trio (proband, mother, father) VCF by 
#' applying the Viterbi algorithm using a Hidden Markov Model (HMM) from the UPDhmm 
#' package. It optionally computes per-sample read depth ratios and summarizes 
#' contiguous variants with the same inferred state into blocks suitable for 
#' downstream analysis.
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
#' @param field_DP Optional character specifying the FORMAT field in the VCF to use 
#'   for read depth metrics. If `NULL` (default), the function will try `"DP"` (standard depth) 
#'   or `"AD"` (allele depths summed across alleles). Use this if your VCF uses a 
#'   non-standard depth field, e.g., `"NR"`.
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
#' @return A `data.table` object containing all detected events in the provided trio. 
#' If no events are found, the function will return an empty `data.table`.
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
                            add_ratios = TRUE,
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose = FALSE) {
  ## --------------------------------------------------------------
  ## 0. Input validation
  ## --------------------------------------------------------------
  if (!inherits(largeCollapsedVcf, "CollapsedVCF")) {
    stop("Argument 'largeCollapsedVcf' must be a CollapsedVCF object.")
  }

  if (is.null(hmm)) {
    utils::data("hmm", package = "UPDhmm", envir = environment())
  }

  ## --------------------------------------------------------------
  ## 1. Optional: compute per-sample depth/quality ratios
  ## --------------------------------------------------------------
  total_sum_per_individual <- total_valid_per_individual <- NULL
  if (add_ratios) {
    geno_list <- VariantAnnotation::geno(largeCollapsedVcf)
    expected_samples <- c("proband", "mother", "father")
  
    dp_field <- if (!is.null(field_DP) && field_DP %in% names(geno_list)) { field_DP } 
                else if ("DP" %in% names(geno_list)) { "DP" } 
                else if ("AD" %in% names(geno_list)) { "AD" } 
                else { NULL }

    if (!is.null(dp_field)) {
      # Handle AD separately (sum across alleles)
      if (dp_field == "AD") {
        quality_matrix <- apply(geno_list$AD, c(1,2), sum)
      } else {
        quality_matrix <- as.matrix(geno_list[[dp_field]])
      }
      # Keep only the expected trio samples if present
      present <- intersect(expected_samples, colnames(quality_matrix))
      if (length(present) > 0L) {

        quality_matrix <- quality_matrix[, present, drop = FALSE]
        total_sum_per_individual <- colSums(quality_matrix, na.rm = TRUE)
        total_valid_per_individual <- colSums(!is.na(quality_matrix))

        # Ensure order proband,mother,father
        total_sum_per_individual <- total_sum_per_individual[expected_samples]
        total_valid_per_individual <- total_valid_per_individual[expected_samples]
      }
    } else {
      warning("No DP or AD field found in VCF.")
    }
  }

  ## --------------------------------------------------------------
  ## 2. Split VCF by chromosome
  ## --------------------------------------------------------------
  split_vcf_raw <- S4Vectors::splitAsList(
    largeCollapsedVcf, GenomicRanges::seqnames(largeCollapsedVcf)
  )
  split_vcf_raw <- split_vcf_raw[lengths(split_vcf_raw) > 0L]

  if (length(split_vcf_raw) == 0L) {
    if (verbose) message("No chromosomes found in VCF.")
    return(data.table::data.table())
  }

  if (verbose) message("Processing ", length(split_vcf_raw), " chromosomes...")

  ## --------------------------------------------------------------
  ## 3. Run HMM pipeline per chromosome
  ## --------------------------------------------------------------
  blocks_state <- if (inherits(BPPARAM, "SerialParam")) {
    lapply(split_vcf_raw, processChromosome,
           total_sum = total_sum_per_individual,
           total_valid = total_valid_per_individual,
           field_DP = field_DP,
           add_ratios = add_ratios,
           hmm = hmm)
  } else {
    BiocParallel::bplapply(split_vcf_raw, processChromosome,
                           total_sum = total_sum_per_individual,
                           total_valid = total_valid_per_individual,
                           field_DP = field_DP,
                           add_ratios = add_ratios,
                           hmm = hmm, BPPARAM = BPPARAM)
  }

  blocks_state <- Filter(Negate(is.null), blocks_state)
  if (length(blocks_state) == 0L) {
    stop("calculateEvents failed: no valid chromosomes processed. 
       Likely cause: applyViterbi does not recognize all GT or 
       check your VCF formatting and trio sample IDs.")
  }

  ## --------------------------------------------------------------
  ## 4. Merge chromosome-level results and filter events
  ## --------------------------------------------------------------
  def_blocks_states <- data.table::rbindlist(blocks_state, use.names = TRUE, fill = TRUE)

  filtered_def_blocks_states <- def_blocks_states[def_blocks_states$n_snps > 1 & def_blocks_states$group != "normal" & !(def_blocks_states$seqnames %in% c("chrX","X")),]

  if (nrow(filtered_def_blocks_states) == 0L) {
    if (verbose) message("No non-normal events found.")
    return(data.table::data.table())
  }

  if (verbose) {
    message("Found ", nrow(filtered_def_blocks_states), " candidate events.")
  }

  return(filtered_def_blocks_states[])
}
