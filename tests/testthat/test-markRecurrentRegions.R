input <- data.frame(
  ID = c("S1", "S2", "S3", "S4"),
  seqnames = c("chr1", "chr1", "chr1", "chr2"),
  start = c(100, 120, 500, 100),
  end = c(150, 170, 550, 150),
  n_mendelian_error = c(10, 20, 5, 200)
)

recurrent_gr <- GenomicRanges::GRanges(
  seqnames = "chr1",
  ranges = IRanges::IRanges(
    start = 100,
    end = 170
  ),
  n_samples = 2
)

expected <- data.frame(
  ID = c("S1", "S2", "S3", "S4"),
  seqnames = c("chr1", "chr1", "chr1", "chr2"),
  start = c(100, 120, 500, 100),
  end = c(150, 170, 550, 150),
  n_mendelian_error = c(10, 20, 5, 200),
  Recurrent = c("Yes", "Yes", "No", "No"),
  n_samples = c(2,2,1,1)
)

test_that("Test if markRecurrentRegions works", {
  out <- markRecurrentRegions(input, recurrent_gr)
  expect_equal(expected,out)
})
