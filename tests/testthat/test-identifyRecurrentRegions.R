df <- data.frame(
  ID = c("S1", "S2", "S3", "S4", "S5", "S6"),
  chromosome = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
  start = c(100, 120, 130, 150, 400, 420),
  end   = c(300, 320, 800, 350, 500, 520),
  n_mendelian_error = c(10, 20, 5, 20, 5, 5)
)

events <- GenomicRanges::GRanges(
  seqnames = "chr1",
  ranges = IRanges::IRanges(start = c(100, 120, 150),
                   end   = c(300, 320, 350)),
  ID = c("S1", "S2", "S4"),
  n_mendelian_error = c(10, 20, 20),
  cluster = c(1,1,1)
)

expected <- GenomicRanges::GRanges(
  seqnames = "chr1",
  ranges = IRanges::IRanges(
    start = 100,
    end = 350
  ),
  n_samples = 3,
  supporting_events = GenomicRanges::GRangesList(events)
)

test_that("identifyRecurrentRegions", {
  # Run function
  out <- identifyRecurrentRegions(df)
  expect_equal(out, expected)
})
