# test-blocksVcf.R

# Expected UPD event blocks 
expected_df <- data.frame(
  ID = "NA19675",
  seqnames = "6",
  start = 32489853,
  end = 32489876,
  group = as.character("iso_mat"),
  n_snps = 3L,
  geno_coded = I(setNames(list(c("133","133","121")), "1")),
  total_sum_quality_proband = 185,
  total_sum_quality_mother = 178,
  total_sum_quality_father = 184,
  total_count_quality_proband = 3,
  total_count_quality_mother = 3,
  total_count_quality_father = 3
)

input <- data.frame(
  ID = c("NA19675","NA19675","NA19675"),
  start = c(32489853, 32489856, 32489876),
  end = c(32489853, 32489856, 32489876),
  group = c("iso_mat", "iso_mat", "iso_mat"),
  seqnames = c("6", "6", "6"),
  geno_coded = c("133", "133", "121"),
  quality_proband = c(60, 60, 65),
  quality_mother = c(56, 58, 64),
  quality_father = c(60, 64, 60)
)

# ------------------------------------------------------------------------- #
# Test basic block merging with full quality information
# ------------------------------------------------------------------------- #

test_that("Test if simplification into blocks works", {
    out <- blocksVcf(input)
    
    out <- as.data.frame(out)
    
    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Test basic block merging without quality metrics
# ------------------------------------------------------------------------- #

test_that("Test if simplification into blocks works", {
  out <- blocksVcf(input[, !(names(input) %in% c("quality_proband", "quality_mother", "quality_father"))])
  
  out <- as.data.frame(out)
  
  expected_no_quality <- expected_df[, !(names(expected_df) %in% c(
    "total_sum_quality_proband", "total_sum_quality_mother", "total_sum_quality_father",
    "total_count_quality_proband", "total_count_quality_mother", "total_count_quality_father"
  ))]
  out_no_quality <- out[, names(expected_no_quality), drop = FALSE]
  
  expect_equal(out_no_quality, out_no_quality)
  expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Introduce NA values to test robustness of block-level aggregation
# ------------------------------------------------------------------------- #

input$quality_proband[1] <- NA
expected_df$total_sum_quality_proband <- expected_df$total_sum_quality_proband - 60
expected_df$total_count_quality_proband <- expected_df$total_count_quality_proband - 1

# ------------------------------------------------------------------------- #
# Test block merging with NA values in quality
# ------------------------------------------------------------------------- #

test_that("Test if simplification into blocks works", {
  out <- blocksVcf(input)
  
  out <- as.data.frame(out)
  
  expect_equal(out, expected_df)
  expect_s3_class(out, "data.frame")
})
