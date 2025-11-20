#' Detect UPD events in trio VCFs using a HMM
#'
#' This function detects regions of uniparental disomy (UPD) in trio genotypes 
#' by applying a Hidden Markov Model (HMM). It runs the Viterbi algorithm 
#' chromosome by chromosome, optionally computes per-sample read depth ratios, 
#' identifies consecutive variants sharing the same inferred state, and summarizes 
#' them into blocks. Blocks are then filtered to retain only those with more than 
#' one SNP, non-normal HMM states, and located on autosomes. The resulting table 
#' includes genomic coordinates, HMM state, number of SNPs, number of Mendelian 
#' errors, and optionally per-block depth-ratio metrics for downstream analysis.
#'
#' @param largeCollapsedVcf A VCF previously processed 
#'   with \code{vcfCheck()} function from UPDhmm package.
#'
#' @param hmm Optional custom Hidden Markov Model object Default = NULL. 
#'  
#'   If NULL, the function uses the default Mendelian HMM included in the UPDhmm package.  
#'   A custom HMM must be a list following the structure of the HMM package, containing:
#'   
#'   \itemize{
#'     \item States – character vector of hidden state names
#'     \item Symbols – vector of allowed observation symbols (genotype codes)
#'     \item startProbs – named vector of initial state probabilities
#'     \item transProbs – state transition probability matrix
#'     \item emissionProbs – matrix of emission probabilities for each state × symbol
#'   }
#'   
#' @param field_DP Default = NULL. Character string specifying which FORMAT field in the VCF
#' contains the read depth information.
#' 
#' If NULL (default), the function will automatically try *DP* (standard depth)
#' or *AD* (allelic depths, summed across alleles).
#' Use this parameter if your VCF uses a non-standard field name for depth,
#' e.g. *field = "NR"* or *field_DP*.
#'
#' @param add_ratios Logical; default = FALSE.
#'   
#'   If TRUE, per-sample read depth sums and counts of valid calls are computed; 
#'   these values are then used to calculate the genome-wide mean per individual, 
#'   which is used to normalize the per-block depth ratios.
#'
#' @param BPPARAM Parallelization settings, passed to
#'   \link[BiocParallel]{bplapply}.
#'   By default *BiocParallel::SerialParam()* (serial execution).
#'   To enable parallelization, provide a BiocParallel backend, e.g.
#'   *BiocParallel::MulticoreParam(workers = min(2, parallel::detectCores()))*
#'   or *BiocParallel::SnowParam(workers = 2)*.
#'   
#'   Note: when running under R CMD check or Bioconductor build systems,
#'   the number of workers may be automatically limited (usually less or equal to 2).
#'
#' @param verbose Logical, default = FALSE. 
#' 
#'   If TRUE, progress messages will be printed during processing.
#'
#' @details
#'
#' The function performs the following major steps:
#'
#' *1. Optional per-sample read depth totals*
#' 
#' If add_ratios = TRUE, the function computes the total read depth and the number of valid calls per individual across the entire VCF, using the field specified in field_DP or, if unavailable, DP or AD. These totals are later used to normalize per-block depth ratios.
#' 
#' *2. Chromosomal splitting and per-chromosome HMM processing*
#' 
#' The VCF is split by chromosome and \code{processChromosome()} is applied to each, which runs the Viterbi algorithm to infer hidden states and groups consecutive variants with the same state into blocks, generating summary metrics for each.  
#' This step can be executed in series or in parallel depending on the BPPARAM parameter.
#' 
#' *3. Consolidation and filtering of detected UPD events*
#' 
#' All blocks from all chromosomes are combined into a single data.frame and filtered to retain only those with more than one SNP, a state different from normal, and located on autosomes. The final output summarizes detected UPD events, including genomic coordinates, HMM state, number of SNPs, number of Mendelian errors per block, and, if calculated, per-block depth ratios.
#' 
#' @return A data.frame summarizing all detected UPD-like events.  
#' Columns include:
#' \itemize{
#'   \item seqnames – chromosome name  
#'   \item start, end – genomic coordinates  
#'   \item group – inferred HMM state  
#'   \item n_snps – number of SNPs in the block  
#'   \item n_mendelian_error – number of Mendelian errors in the block  
#'   \item depth-ratio metrics (if add_ratios = TRUE)  
#' }
#'
#' If no events are detected, returns an empty data.frame.
#'
#' @examples
#' file <- system.file(package = "UPDhmm", "extdata", "test_het_mat.vcf.gz")
#' vcf <- VariantAnnotation::readVcf(file)
#'
#' processedVcf <- vcfCheck(
#'     vcf,
#'     proband = "NA19675",
#'     mother  = "NA19678",
#'     father  = "NA19679"
#' )
#'
#' # Run in serial mode (default)
#' res <- calculateEvents(processedVcf)
#'
#' # Run in parallel with 2 cores
#' library(BiocParallel)
#' param <- MulticoreParam(workers = 2)
#' res_parallel <- calculateEvents(processedVcf, BPPARAM = param)
#'
#' @export
#' 
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
  
  total_sum_per_individual <- total_valid_per_individual <- NULL
  if (add_ratios) {
    trio_totals <- computeTrioTotals(largeCollapsedVcf, field_DP = field_DP)
    total_sum_per_individual <- trio_totals$total_sum
    total_valid_per_individual <- trio_totals$total_valid
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
           total_sum = total_sum_per_individual,
           total_valid = total_valid_per_individual,
           field_DP = field_DP,
           add_ratios = add_ratios,
           hmm = hmm, 
           mendelian_error_values = mendelian_error_values)
  } else {
    BiocParallel::bplapply(split_vcf_raw, processChromosome,
                           total_sum = total_sum_per_individual,
                           total_valid = total_valid_per_individual,
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
  total_sum <- total_valid <- NULL
  geno_list <- VariantAnnotation::geno(vcf)
  
  # ---------------------------------------------------------------
  # Determine which depth/coverage field to use for calculations
  # Priority:
  # 1) Use 'field_DP' if specified and exists in VCF
  # 2) Use standard 'DP' field if present
  # 3) Use 'AD' (allele depths) if present
  # 4) Otherwise, no depth field available
  # ---------------------------------------------------------------
  
  dp_field <- if (!is.null(field_DP) && field_DP %in% names(geno_list)) { field_DP 
              } else if ("DP" %in% names(geno_list)) { "DP" 
              } else if ("AD" %in% names(geno_list)) { "AD" 
              } else { NULL }
  
  if (!is.null(dp_field)) {
    if (dp_field == "AD") {
      # If using allele depths (AD), sum across all alleles for each sample
      # Handle cases where all values are NA by returning NA
      quality_matrix <- apply(geno_list$AD, 2, function(col) {
        vapply(col, function(x) {if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)}, numeric(1))
      })
      
    } else {
      quality_matrix <- as.matrix(geno_list[[dp_field]])
    }
      
    # Compute total read depth per individual, ignoring NAs
    total_sum <- colSums(quality_matrix, na.rm = TRUE)
    
    # Compute number of valid (non-NA) calls per individual
    total_valid <- colSums(!is.na(quality_matrix))
    
    # Ensure order proband, mother, father
    total_sum <- total_sum[expected_samples]
    total_valid <- total_valid[expected_samples]
    
  } else {
    warning("No DP or AD field found in VCF.")
  }
  
  list(total_sum = total_sum, total_valid = total_valid)
}
