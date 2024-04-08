expected_df <- data.frame(
    start = 32489853,
    end = 32489876,
    group = as.character("iso_mat"),
    seqnames = "6",
    n_snps = 3
)


expected_df$group <- as.character(expected_df$group)

input <- data.frame(
  start = c(32489853, 32489856, 32489876),
  end = c(32489853, 32489856, 32489876),
  group = c("iso_mat", "iso_mat", "iso_mat"),
  seqnames = c("6", "6", "6")
)
test_that("Test if simplification into blocks works", {
    out <- blocksVcf(input)
    out <- as.data.frame(out)
    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})
