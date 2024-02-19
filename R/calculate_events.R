#' Calculate UPD events in trio vcfs
#'
#' This function is for appliying the hidden-markov-model using the viterbi
#' algorithm after that it simplify into blocks and finally export into a
#' bed file
#'
#' @param largecollapsedVcf The general format of vcf with VariantAnnotation
#' package
#' @return optional save object (raw and simplify ) and export bed file
#'
#' @export
#' @examples
#' fl <- system.file("extdata", "test.vcf.gz", package = "UPDhmm")
#' vcf <- VariantAnnotation::readVcf(fl)
#' largecollapsedVcf <- vcf_check(vcf,
#'   proband = "Sample1", mother = "Sample3",
#'   father = "Sample2"
#' )
#' calculate_events(largecollapsedVcf)
calculate_events <- function(largecollapsedVcf) {
  # split the vcf into chromosomes
  split_vcf <- split(largecollapsedVcf,
    f = GenomicRanges::seqnames(largecollapsedVcf)
  )

  #apply viterbi
   split_vcf <- lapply(
     split_vcf,
     purrr::possibly(apply_viterbi, otherwise = NULL)
   )

  # as_df
  split_vcf_df <- sapply(split_vcf, function(x) {
    tryCatch(
      as_df_vcf(x),
      error = function(e) NULL
    )
  }, simplify = FALSE)

  blocks_state <- sapply(split_vcf_df, function(df) {
    tryCatch(
      blocks_vcf(df),
      error = function(e) NULL
    )
  }, simplify = FALSE)
  # simplify all chr into one datafram
  def_blocks_states <- data.table::rbindlist(blocks_state[!sapply(
    blocks_state,
    is.null
  )])

  colnames(def_blocks_states) <- c("start", "end", "group", "n_snps", "seqnames")
  filtered_def_blocks_states <- subset(
    def_blocks_states,
    def_blocks_states$n_snps > 1 &
      def_blocks_states$group != "normal" &
      !(def_blocks_states$seqnames %in%
        c("chrX", "X"))
  )
  filtered_def_blocks_states$log_OR <- numeric(nrow(filtered_def_blocks_states))
  filtered_def_blocks_states$p_value <- numeric(nrow(filtered_def_blocks_states))
  if (nrow(filtered_def_blocks_states) > 0) {
    for (i in seq_len(nrow(filtered_def_blocks_states))) {
      filtered_def_blocks_states$log_OR[[i]] <-
        as.numeric(add_or(
          filtered_def_blocks_states$start[[i]],
          filtered_def_blocks_states$end[[i]],
          filtered_def_blocks_states$seqnames[[i]],
          filtered_def_blocks_states$group[[i]],
          split_vcf_df
        )["log_OR"])
      filtered_def_blocks_states$p_value[[i]] <-
        as.numeric(add_or(
          filtered_def_blocks_states$start[[i]],
          filtered_def_blocks_states$end[[i]],
          filtered_def_blocks_states$seqnames[[i]],
          filtered_def_blocks_states$group[[i]],
          split_vcf_df
        )["p_value"])
    }
  } else {
    methods::show("No events found")
  }
  return(filtered_def_blocks_states)
}
