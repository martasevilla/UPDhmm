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
#' @param add_ratios Logical, default = `FALSE`. If `TRUE`, the function computes 
#'   per-sample totals of read depth and counts of valid calls for the trio, which 
#'   can be used for normalization or downstream metrics.
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
                            add_ratios = FALSE,
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose = FALSE) {
  ## --------------------------------------------------------------
  ## 0. Input validation
  ## --------------------------------------------------------------
  if (!inherits(largeCollapsedVcf, "CollapsedVCF")) {
    stop("Argument 'largeCollapsedVcf' must be a CollapsedVCF object.")
  }

  # Load the default HMM from the UPDhmm package if no custom HMM is provided
  if (is.null(hmm)) {
    utils::data("hmm", package = "UPDhmm", envir = environment())
  }

  ## --------------------------------------------------------------
  ## 1. Optional: compute per-sample depth/quality ratios
  ## --------------------------------------------------------------
  
  mean_depth_per_individual <- NULL
  if (add_ratios) {
    mean_depth_per_individual <- computeTrioTotals(vcf = largeCollapsedVcf, field_DP = field_DP)
  }
  
  ## --------------------------------------------------------------
  ## 2. Split VCF by chromosome
  ## --------------------------------------------------------------
  split_vcf_raw <- split(
    largeCollapsedVcf, GenomicRanges::seqnames(largeCollapsedVcf)
  )
  split_vcf_raw <- split_vcf_raw[lengths(split_vcf_raw) > 0L]

  if (length(split_vcf_raw) == 0L) {
    if (verbose) message("No chromosomes found in VCF.")
    return(data.frame())
  }

  if (verbose) message("Processing ", length(split_vcf_raw), " chromosomes...")
  
  # Determine which genotype codes correspond to Mendelian errors (lowest emission probability for 'normal' state)
  emission_probs <- hmm$emissionProbs["normal", ]
  mendelian_error_values <- names(emission_probs[emission_probs == min(emission_probs)])

  ## --------------------------------------------------------------
  ## 3. Run HMM pipeline per chromosome
  ## --------------------------------------------------------------
  blocks_state <- if (inherits(BPPARAM, "SerialParam")) {
    lapply(split_vcf_raw, processChromosome,
           total_mean = mean_depth_per_individual,
           field_DP = field_DP,
           add_ratios = add_ratios,
           hmm = hmm, 
           mendelian_error_values = mendelian_error_values)
  } else {
    BiocParallel::bplapply(split_vcf_raw, processChromosome,
                           total_mean = mean_depth_per_individual,
                           field_DP = field_DP,
                           add_ratios = add_ratios,
                           hmm = hmm, BPPARAM = BPPARAM, 
                           mendelian_error_values = mendelian_error_values)
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
  def_blocks_states <- do.call(rbind, blocks_state)
  
  # Keep only blocks with >1 SNP, non-normal state, and autosomal chromosomes
  filtered_def_blocks_states <- def_blocks_states[def_blocks_states$n_snps > 1 & def_blocks_states$group != "normal" & !(def_blocks_states$seqnames %in% c("chrX","X","chrY","Y")),]

  if (nrow(filtered_def_blocks_states) == 0L) {
    if (verbose) message("No non-normal events found.")
    return(data.frame())
  }

  if (verbose) {
    message("Found ", nrow(filtered_def_blocks_states), " candidate events.")
  }
  rownames(filtered_def_blocks_states) <- NULL
  return(filtered_def_blocks_states)
}


computeTrioTotals <- function(vcf, expected_samples = c("proband","mother","father"), field_DP = NULL) {
  mean_depth <- NULL
  geno_list <- VariantAnnotation::geno(vcf)
  
  # ---------------------------------------------------------------
  # Determine which depth/coverage field to use for calculations
  # Priority:
  # 1) Use 'field_DP' if specified and exists in VCF
  # 2) Use standard 'DP' field if present
  # 3) Use 'AD' (allele depths) if present
  # 4) Otherwise, no depth field available
  # ---------------------------------------------------------------
  
  dp_field <- if (!is.null(field_DP) && field_DP %in% names(geno_list)) {
    field_DP 
  } else if ("DP" %in% names(geno_list)) {
    "DP" 
  } else if ("AD" %in% names(geno_list)) {
    "AD" 
  } else {
    NULL
  }
  
  if (!is.null(dp_field)) {
    if (dp_field == "AD") {
      # If using allele depths (AD), sum across all alleles for each sample
      # Handle cases where all values are NA by returning NA
      depth_matrix <- apply(geno_list$AD, 2, function(col) {
        vapply(col, function(x) {if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)}, numeric(1))
      })
    } else {
      depth_matrix <- as.matrix(geno_list[[dp_field]])
    }
    
    # Compute mean depth per individual, ignoring NA values
    mean_depth <- colMeans(depth_matrix, na.rm = TRUE)
    
    # Ensure the order proband, mother, father
    mean_depth <- mean_depth[expected_samples]
    
  } else {
    warning("No DP or AD field found in VCF.")
  }

  return(mean_depth)
}
