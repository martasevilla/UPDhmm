df <- data.frame(
  ID = c("S1", "S2", "S3", "S3"),
  seqnames = c("chr1", "chr1", "chr1", "chr1"),
  start = c(100, 120, 500, 510),
  end = c(150, 170, 550, 560),
  n_mendelian_error = c(10, 20, 5, 5)
)


expected <- GenomicRanges::GRanges(
  seqnames = "chr1",
  ranges = IRanges::IRanges(
    start = 100,
    end = 170
  ),
  n_samples = 2
)



test_that("identifyRecurrentRegions", {
  # Run function
  out <- identifyRecurrentRegions(df, ID_col = "ID", error_threshold = 50, min_support = 2)
  expect_equal(out, expected)
})
