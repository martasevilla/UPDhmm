#' Calculate UPD events in trio VCFs.
#'
#' This function predicts the hidden states by applying the Viterbi algorithm
#' using the Hidden Markov Model (HMM) from the UPDhmm package. It takes the
#' genotypes of the trio as input and includes a final step to simplify the
#' results into blocks.
#'
#' @param largecollapsedVcf The VCF file in the general format
#' @param hmm Hidden Markov Model used to infer the events
#' @param genotypes Possible GT formats and its correspondency with the hmm
#' (largecollapsedVcf) with VariantAnnotation package.
#' @return Dataframe object containing blocks of predicted events.
#' @export
#' @examples
#' fl <- system.file("extdata", "test.vcf.gz", package = "UPDhmm")
#' vcf <- VariantAnnotation::readVcf(fl)
#' largecollapsedVcf <- vcf_check(vcf,
#'   proband = "Sample1", mother = "Sample3",
#'   father = "Sample2"
#' )
#' calculate_events(largecollapsedVcf)
calculate_events <- function(largecollapsedVcf=NULL,genotypes=NULL) {

  utils::data("hmm")
  hmm <- hmm
  genotypes <-  c("0/0" = "1", "0/1" = "2","1/0" = "2", "1/1" = "3",
                  "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3" )

  # 1 split the vcf into chromosomes
  split_vcf_raw <- split(largecollapsedVcf,
                     f = GenomicRanges::seqnames(largecollapsedVcf))

  #2 apply viterbi
  split_vcf <- lapply(split_vcf_raw, function(x) {
   tryCatch(
     UPDhmm:::apply_viterbi(largecollapsedVcf = x, hmm = hmm, genotypes = genotypes),
     error = function(e) NULL  # Return NULL if an error occurs
   )
 })



  #3 as_df
  split_vcf_df <- lapply(split_vcf, function(x)
  { tryCatch(UPDhmm:::as_df_vcf(largecollapsedVcf=x,
                       genotypes = genotypes),
             error = function(e) NULL) })


  # #4 Create blocks of contiguous positions with same state
   blocks_state <- lapply(split_vcf_df, function(df) {
    tryCatch(UPDhmm:::blocks_vcf(df),error = function(e) NULL)})

   #5 simplify all chr objects into one data.frame
   def_blocks_states <- data.table::rbindlist(blocks_state[!sapply(blocks_state, is.null)])


  #6 Filter normal state blocks , sexual chromosomes and isolated
  #variants of a certain state
   filtered_def_blocks_states <-
   def_blocks_states[def_blocks_states$n_snps > 1 &
                       def_blocks_states$group != "normal" &
                       !(def_blocks_states$seqnames %in% c("chrX", "X")), ]

  # 7 Calculate statistics parameters
   if (nrow(filtered_def_blocks_states) > 0) {
   blocks_list <- lapply(seq_len(nrow(filtered_def_blocks_states)), function(i) {
     UPDhmm:::add_or(
       filtered_def_blocks_states = filtered_def_blocks_states[i, , drop = FALSE],
       largecollapsedVcf = largecollapsedVcf,
       hmm = hmm,
       genotypes = genotypes
     )

     })

   }
   else{

     print("No events found")
     block_def<-data.frame()
   }

  # 8 Transform final output to data.frame
  block_def <- data.table::rbindlist(blocks_list)
  return(block_def)
}


