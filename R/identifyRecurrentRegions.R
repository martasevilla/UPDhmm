#' Identify recurrent genomic regions across samples
#'
#' This function identifies recurrent genomic regions supported by multiple
#' samples based on overlapping genomic intervals. UPD events are clustered based
#' on their reciprocal overlap. Recurrent regions are defined as clusters supported
#' by a minimum number of samples.
#' 
#' @param df A data.frame with columns:
#'   - ID: Sample identifier.
#'   - chromosome: Chromosome name.
#'   - start: Start coordinate of the region.
#'   - end: End coordinate of the region.
#'   - n_mendelian_error or total_mendelian_error: Number of Mendelian errors in the region.
#' @param error_threshold Numeric, default = 100.
#'   Maximum number of Mendelian errors allowed for a region to be considered.
#' @param min_support Integer, default = 3.
#'   Minimum number of unique samples required to call a region recurrent.
#' @param ID_col Character string indicating the column name containing
#'   sample identifiers. Default is `"ID"`.
#' @param max_dist Numeric, default = 0.3.
#'   Maximum distance allowed to group events in the same cluster when cutting
#'   the hierarchical clustering dendrogram.
#' @param linkage Character string specifying the linkage method for the
#'   agglomerative hierarchical clustering. Default `"complete"`, where the
#'   distance between two clusters is defined as the largest distance between
#'   all possible pairs of events belonging to the two clusters.
#'
#' @return A GRanges object containing the recurrent regions that meet
#'   the minimum support threshold. Metadata columns include:
#'   \itemize{
#'     \item \code{n_samples} – number of distinct samples supporting the region.
#'     \item \code{supporting_events} – a GRangesList with the individual events
#'       that compose the recurrent region.
#'   }
#'   
#' The \code{supporting_events} correspond exactly to the original events that
#'   were grouped into the cluster and define the boundaries of the resulting
#'   recurrent region.
#'   
#' If no clusters meet the filtering criteria, the function returns \code{NULL}.
#' 
#' @details
#' The function operates chromosome-wise and proceeds in several steps. First,
#' events with a high number of Mendelian errors are filtered out. Then, a
#' pairwise distance matrix is computed for overlapping events, where the
#' distance between two events e1 and e2 is defined as:
#'
#' \deqn{1 - (overlap\_width / max(width(e1), width(e2)))}
#'
#' Event pairs that do not overlap are not explicitly evaluated and are 
#' assigned a distance of 1. In addition, only the upper triangular part of the
#' distance matrix is computed, exploiting its symmetry to avoid redundant
#' calculations and reduce computational cost.
#'
#' Agglomerative hierarchical clustering is applied to this distance matrix. At
#' each step, the most similar clusters are merged according to the specified
#' linkage method (\code{linkage}), with the default \code{"complete"} linkage
#' using the maximum distance between all possible pairs of events from the two
#' clusters. The resulting dendrogram is cut at a height defined by
#' \code{max_dist}, which represents the maximum distance allowed for events to
#' be grouped in the same cluster. Smaller values enforce stricter overlap
#' similarity.
#'
#' Only clusters supported by at least \code{min_support} distinct samples are
#' retained. Each valid cluster is then collapsed into a single recurrent region
#' whose coordinates span all contributing events, and the original events are
#' stored as \code{supporting_events} with coverage metrics.
#' 
#' These steps are performed internally by a helper function \code{identifyRecurrentRegionsByChr}.
#' Results from all chromosomes are then combined into a single GRanges object returned by this function.
#' 
#' @export
#' @examples
#' df <- data.frame(
#'   ID = c("S1", "S2", "S3", "S4"),
#'   chromosome = c("chr1", "chr1", "chr1", "chr1"),
#'   start = c(100, 120, 500, 510),
#'   end = c(150, 170, 550, 560),
#'   n_mendelian_error = c(10, 20, 5, 5)
#' )
#' identifyRecurrentRegions(df, ID_col = "ID", error_threshold = 50, min_support = 2)
#' 
identifyRecurrentRegions <- function(df,
                                     ID_col = "ID",
                                     error_threshold = 100,
                                     min_support = 3,
                                     max_dist = 0.3,
                                     linkage = "complete") {
  
  #---------------------------------------------------------------
  # 1. Validate inputs
  #---------------------------------------------------------------
  if (!"chromosome" %in% names(df)) {
    stop("Input must contain column 'chromosome'.")
  }
  
  #---------------------------------------------------------------
  # 2. Split data by chromosome
  #---------------------------------------------------------------
  data_by_chr <- split(df, df$chromosome)
  
  #---------------------------------------------------------------
  # 3. Apply recurrent region detection per chromosome
  #---------------------------------------------------------------
  recurrent_by_chr <- lapply(data_by_chr, function(x) {
    identifyRecurrentRegionsByChr(
      df = x,
      ID_col = ID_col,
      error_threshold = error_threshold,
      min_support = min_support,
      max_dist = max_dist,
      linkage = linkage
    )
  })
  
  #---------------------------------------------------------------
  # 4. Remove chromosomes without recurrent regions
  #---------------------------------------------------------------
  recurrent_by_chr <- recurrent_by_chr[!sapply(recurrent_by_chr, is.null)]
  
  if (length(recurrent_by_chr) == 0) {
    return(NULL)
  }
  
  #---------------------------------------------------------------
  # 5. Combine results into a single GRanges object
  #---------------------------------------------------------------
  grl <- GenomicRanges::GRangesList(recurrent_by_chr)
  
  return (unlist(grl, use.names = FALSE))
}


identifyRecurrentRegionsByChr  <- function(df,
                                           ID_col = "ID",
                                           error_threshold = 100,
                                           min_support = 3,
                                           max_dist = 0.3,
                                           linkage = "complete") {
  
  #---------------------------------------------------------------
  # 1. Validate inputs
  #---------------------------------------------------------------
  if (!ID_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in input dataframe.", ID_col))
  }
  
  if (!all(c("start", "end") %in% names(df))) {
    stop("Input must contain columns 'start' and 'end'.")
  }
  
  err_col <- if ("n_mendelian_error" %in% names(df)) {
    "n_mendelian_error"
  } else if ("total_mendelian_error" %in% names(df)) {
    "total_mendelian_error"
  } else {
    stop("Input must contain either 'n_mendelian_error' or 'total_mendelian_error'.")
  }
  
  #---------------------------------------------------------------
  # 2. Convert to GRanges
  #---------------------------------------------------------------
  gr <- GenomicRanges::makeGRangesFromDataFrame(
    df,
    seqnames.field = "chromosome",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE
  )
  
  
  #---------------------------------------------------------------
  # 3. Filter by error threshold
  #---------------------------------------------------------------
  keep <- S4Vectors::mcols(gr)$n_mendelian_error < error_threshold
  if (!any(keep)) return(NULL)
  
  gr <- gr[keep]
  n  <- length(gr)
  
  # Check if there are enough events to cluster
  if (n < 2) {
    return(NULL)
  }
  
  #---------------------------------------------------------------
  # 4. Distance matrix
  #---------------------------------------------------------------
  
  # Initialize distance matrix with maximum distance (1 = no overlap)
  dist_mat <- matrix(1, n, n)
  diag(dist_mat) <- 0
  
  # Identify overlapping event pairs (distance is computed only for these to reduce computational cost)
  hits <- GenomicRanges::findOverlaps(gr, gr)
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  
  # Retain only the upper triangular matrix to avoid redundant computations
  upper <- qh < sh 
  qh <- qh[upper]
  sh <- sh[upper]
  
  # Compute distance as 1 - (overlap / maximum event width)
  ov <- GenomicRanges::pintersect(gr[qh], gr[sh])
  max_len <- pmax(IRanges::width(gr[qh]), IRanges::width(gr[sh]))
  d <- 1 - (IRanges::width(ov) / max_len)
  
  # Fill symmetric distance matrix
  dist_mat[cbind(qh, sh)] <- d
  dist_mat[cbind(sh, qh)] <- d
  
  # Perform agglomerative hierarchical clustering on the distance matrix
  hc <- stats::hclust(stats::as.dist(dist_mat), method = linkage) 
  
  # Cut the dendrogram at height max_dist to define clusters
  cluster <- stats::cutree(hc, h = max_dist) 
  S4Vectors::mcols(gr)$cluster <- cluster
  
  #---------------------------------------------------------------
  # 5. Filter clusters by minimum support
  #---------------------------------------------------------------
  meta <- S4Vectors::mcols(gr)
  
  n_samples <- tapply(meta$ID, meta$cluster, function(x) length(unique(x)))
  valid_clusters <- names(n_samples)[n_samples >= min_support]
  
  if (is.integer(cluster)) {
    valid_clusters <- as.integer(valid_clusters)
  } else if (is.numeric(cluster)) {
    valid_clusters <- as.numeric(valid_clusters)
  }
  
  if (length(valid_clusters) == 0) return(NULL)
  
  #---------------------------------------------------------------
  # 6. Build recurrent regions
  #---------------------------------------------------------------
  regions_gr <- lapply(valid_clusters, function(cl) {
    
    events <- gr[cluster == cl]
    
    # Define region boundaries spanning all supporting events
    region_start <- min(GenomicRanges::start(events))
    region_end   <- max(GenomicRanges::end(events))
    
    GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(events)[1],
      ranges   = IRanges::IRanges(region_start, region_end),
      n_samples = length(unique(S4Vectors::mcols(events)$ID)),
      supporting_events = GenomicRanges::GRangesList(events)
    )
  })
  
  do.call(c, regions_gr)
}
