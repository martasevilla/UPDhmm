#' Calculate UPD events in trio VCFs.
#'
#' This function predicts the hidden states by applying the Viterbi algorithm
#' using the Hidden Markov Model (HMM) from the UPDhmm package. It takes the
#' genotypes of the trio as input and includes a final step to simplify the
#' results into blocks.
#'
#' @param largeCollapsedVcf The VCF file in the general format 
#' (largeCollapsedVcf) with VariantAnnotation package. Previously edited with 
#' `vcfCheck()` function from UPDhmm package
#' @param hmm Default = NULL. If no arguments are added, the package 
#' will use the default HMM already implemented, based on Mendelian 
#' inheritance. If an optional HMM is desired, it should adhere to the 
#' general HMM format from `HMM` package with the following elements inside 
#' a list:
#' 1. The hidden state names in the "States" vector.
#' 2. All possible observations in the "Symbols" vector.
#' 3. Start probabilities of every hidden state in the "startProbs" vector.
#' 4. Transition probabilities matrix between states in "transProbs".
#' 5. Probabilities associated between every hidden state and all possible 
#' observations in the "emissionProbs" matrix. 
#' @export
#' @examples
#' file <- system.file(package = "UPDhmm", "extdata", "test_het_mat.vcf.gz")
#' vcf <- VariantAnnotation::readVcf(file)
#' processedVcf <- vcfCheck(vcf,
#'     proband = "NA19675", mother = "NA19678",
#'     father = "NA19679"
#' )
#'
#'@return A data.frame object containing all detected events in the provided trio. 
#'If no events are found, the function  will return an empty data.frame. 
#'
#'
calculateEvents <-
    function(largeCollapsedVcf, 
        hmm = NULL 
        ) {
        
  # Check if `largeCollapsedVcf` is a VCF object
    if (!inherits(largeCollapsedVcf, "CollapsedVCF")) {
            stop("Argument 'largeCollapsedVcf' must be a VCF object.")
        }

  # Check if `hmm` is a list or assign default hmm
    if (base::is.null(hmm)) {
            utils::data("hmm", package = "UPDhmm", envir = environment())
        }

  # Create genotypes vector    
    genotypes <- c(
                "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
                "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
            )
#################################################
# 1 split the vcf into chromosomes
#################################################
  
split_vcf_raw <- split(largeCollapsedVcf, f = GenomicRanges::seqnames(largeCollapsedVcf))
    
# Check if elements inside `split_vcf_raw` are VCF object
names_split_vcf_raw <- names(split_vcf_raw)
# Iterate through the elements using their names
for (name in names_split_vcf_raw) {
  element <- split_vcf_raw[[name]]
   if (!inherits(element, "CollapsedVCF")) {
     stop("Something went wrong during splitting vcf into chromosomes.")
          }
        }


#################################################
# 2 apply viterbi (Return NULL if an error occurs)
#################################################

split_vcf <- lapply(split_vcf_raw, function(x) {
  tryCatch(applyViterbi(largeCollapsedVcf = x, hmm = hmm, genotypes = genotypes),
          error = function(e) NULL)
        })


# Check if `split_vcf` is a list
if (!inherits(split_vcf, "list")) {
  stop("Something went wrong during applyViterbi.")
}

# Check elements inside `split_vcf`
for (gcf in split_vcf) {
  if (!inherits(gcf, "CollapsedVCF")) {
    stop("Something went wrong during applyViterbi.")
            }
        }

#################################################
# 3 as_df
#################################################
        
split_vcf_df <- lapply(split_vcf, function(x) {
  tryCatch(asDfVcf(largeCollapsedVcf = x, genotypes = genotypes),
                error = function(e) NULL)
        })

# Check if `split_vcf_df` is a list
    if (!inherits(split_vcf_df, "list")) {
      stop("Something went wrong while creating dataframes of detected events.")
        }

# Check elements inside `split_vcf_df`
    for (df in split_vcf_df) {
      if (!inherits(df, "data.frame")) {
        stop("Something went wrong while creating dataframes of detected events.")
            }
        }

    for (df in split_vcf_df) {
      if (ncol(df) != 5) {
        stop("Something went wrong while creating dataframes of detected events.")
            }
        }
#################################################
# 4 Create blocks of contiguous positions with same state
#################################################

blocks_state <- lapply(split_vcf_df, function(df) {
  tryCatch(blocksVcf(df), error = function(e) NULL)
        })

# Check if `blocks_state` is a list
if (!inherits(blocks_state, "list")) {
    stop("Something went wrong during the transformation into blocks of coordinates.")
        }

# Check elements inside `blocks_state`
for (df in blocks_state) {
  if (!inherits(df, "data.frame")) {
  stop("Something went wrong during the transformation 
                     into blocks of coordinates.")
            }
        }

for (df in blocks_state) {
  if (ncol(df) != 5) {
  stop("Something went wrong during the transformation 
                     into blocks of coordinates.")
    }
        }

#################################################
# 5 simplify all chr objects into one data.frame
#################################################
       
def_blocks_states <- base::Reduce(rbind, blocks_state)

# Check if `def_blocks_state` is a data.frame and with corrected format
if (!inherits(def_blocks_states, "data.frame")) {
  stop("Something went wrong during the reduction of blocks.")
        }

if (ncol(def_blocks_states) != 5) {
  stop("Something went wrong during the reduction of blocks.")
        }

#########################################################################
# 6 Filter normal state blocks , sexual chromosomes and isolated variants
#########################################################################
filtered_def_blocks_states <-
  def_blocks_states[def_blocks_states$n_snps > 1 &
  def_blocks_states$group != "normal" &
  !(def_blocks_states$seqnames %in% c("chrX", "X")), ]

# Check if `filtered_def_blocks_states` is a data.frame and with corrected format
if (!inherits(filtered_def_blocks_states, "data.frame")) {
  stop("Something went wrong during filtering of events.")
        }
if (ncol(filtered_def_blocks_states) != 5) {
  stop("Something went wrong during filtering of events.")
        }


####################################################
# 7 Calculate statistics parameters
####################################################
if (nrow(filtered_def_blocks_states) > 0) {
  blocks_list <- lapply(seq_len(nrow(filtered_def_blocks_states)),
                function(i) {
                  addOr(filtered_def_blocks_states = 
                  filtered_def_blocks_states[i, , drop = FALSE],
                  largeCollapsedVcf = largeCollapsedVcf,
                  hmm = hmm,
                  genotypes = genotypes
                    )
                }
            )
        } else {
            block_def <- data.frame()
        }


# Check if `blocks_list` is in corrected format and new columns have been added.
for (df in blocks_list) {
  if (!inherits(df, "data.frame")) {
  stop("Something went wrong during the calculation of statistics parameters.")
            }
        }

for (df in blocks_list) {
  if (ncol(df) != 8) {
  stop("Something went wrong during the calculation 
                     of statistics parameters.")
            }
        }


########################################################
# 8 Transform final output to data.frame
########################################################

block_def <- base::Reduce(rbind, blocks_list)

# Check if `block_def` is in corrected format
 if (!inherits(block_def, "data.frame")) {
 stop("Something went wrong during final reduction and creation of final output.")
   }

if (ncol(block_def) != 8) {
  stop("Something went wrong during final reduction and creation of final output.")
        }


return(block_def)
    }

