#' Function to simplify contiguous variants with the same state into
#' blocks.
#'
#' @param df data.frame resulting from the `as_df_vcf` function.
#' @return data.frame containing information on the chromosome, start
#' #' position of the block, end position of the block, and predicted state.


blocksVcf <- function(df) {
    # Add a column for n_snps_raw
    df$n_snps_raw <-
        base::cumsum(c(TRUE, df$group[-1L] != df$group[-length(df$group)]))
    # Create a vector with unique ids
    unique_n_snps_raw <- base::unique(df$n_snps_raw)
    
    # Create a vector with the Sample ID
    Sample_ID <- unique(df$ID)
    
    # Iterate through data.frame to create the ultimate data.frame
    result_list <- lapply(unique_n_snps_raw, function(n) {
        subset_df <- df[df$n_snps_raw == n, ]
        data.frame(
            ID = Sample_ID,
            start = min(subset_df$start),
            end = max(subset_df$end),
            group = unique(subset_df$group),
            seqnames = unique(subset_df$seqnames),
            n_snps = nrow(subset_df),
            stringsAsFactors = FALSE
        )
    })

    simplified_df <- base::Reduce(rbind, result_list)
    
    return(simplified_df)
}
