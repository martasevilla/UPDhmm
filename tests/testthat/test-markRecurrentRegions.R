input <- data.frame(
  ID = c("S1", "S2", "S3", "S4"),
  chromosome = c("chr1", "chr1", "chr1", "chr1"),
  start = c(100, 150, 182, 400),
  end   = c(300, 350, 382, 500)
)

recurrent_gr <- GenomicRanges::GRanges(
  seqnames = "chr1",
  ranges = IRanges::IRanges(
    start = 120, 
    end = 320
  ),
  n_samples = 3
)

expected <- data.frame(
  ID = c("S1", "S2", "S3", "S4"),
  chromosome = c("chr1", "chr1", "chr1", "chr1"),
  start = c(100, 150, 182, 400),
  end   = c(300, 350, 382, 500),
  Recurrent = c("Yes", "Yes", "No", "No"),
  n_samples = c(3, 3, 1, 1)
)

test_that("Test if markRecurrentRegions works", {
  out <- markRecurrentRegions(input, recurrent_gr)
  expect_equal(expected,out)
})
