# test-asDfVcf.R

# Expected output reference data.frame 
expect_df <- data.frame(
  ID = rep("NA19685", 15),
  start = c(
    32489853, 32589856, 32599876, 33489910, 33499925,
    22368862, 22368905, 22369426, 22382655, 22382828,
    22925851, 22929905, 22929906, 42109223, 42109975
  ),
  end = c(
    32489853, 32589856, 32599876, 33489910, 33499925,
    22368862, 22368905, 22369426, 22382655, 22382828,
    22925851, 22929905, 22929906, 42109223, 42109975
  ),
  group = rep("iso_mat", 15),
  seqnames = c(rep("6", 5), rep("15", 10)),
  geno_coded = c(
    "133", "133", "121", "122", "133",
    "123", "122", "133", "321", "312",
    "231", "212", "323", "231", "123"
  ),
  quality_proband = c(
    60, 60, 65, 60, 50,
    60, 60, 60, 63, 60,
    62, 60, 60, 64, 60
  ),
  quality_mother = c(
    56, 58, 64, 58, 50,
    58, 60, 60, 58, 62,
    60, 64, 60, 60, 58
  ),
  quality_father = c(
    60, 64, 60, 60, 52,
    63, 62, 60, 60, 60,
    60, 58, 60, 60, 63
  )
)

file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)

input <- vcfCheck(
  largeCollapsedVcf = input,
  father = "NA19689", mother = "NA19688",
  proband = "NA19685", check_quality = TRUE
)

S4Vectors::mcols(input)$states <- c("iso_mat", "iso_mat", "iso_mat")

# ------------------------------------------------------------------------- #
# Test basic conversion with add_ratios = FALSE
# ------------------------------------------------------------------------- #

test_that("Test if transformation of vcf with states to basic dataframe works", {
    out <- asDfVcf(largeCollapsedVcf = input)
    
    out$seqnames <- as.character(out$seqnames)
    
    expect_false(any(c("quality_proband", "quality_mother", "quality_father") %in% names(out)))
    expected_no_ratio <- expect_df[, !(names(expect_df) %in% c("quality_proband", "quality_mother", "quality_father"))]
    
    expect_equal(out[, names(expected_no_ratio)], expected_no_ratio)
    expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Test conversion with add_ratios = TRUE using explicit DP field
# ------------------------------------------------------------------------- #

test_that("Test if transformation of vcf with states to basic dataframe works", {
  out <- asDfVcf(largeCollapsedVcf = input, add_ratios = TRUE, field_DP = "DP")
  
  expect_s3_class(out, "data.frame")
  expect_equal(expect_df, out)
})

# ------------------------------------------------------------------------- #
# Test conversion with add_ratios = TRUE without specifying a field
# ------------------------------------------------------------------------- #

test_that("Test if transformation of vcf with states to basic dataframe works", {
  out <- asDfVcf(largeCollapsedVcf = input, add_ratios = TRUE)
  
  expect_s3_class(out, "data.frame")
  expect_equal(expect_df, out)
})

# ------------------------------------------------------------------------- #
# Test conversion with add_ratios = TRUE using explicit AD field
# ------------------------------------------------------------------------- #

test_that("Test if transformation of vcf with states to basic dataframe works", {
  out <- asDfVcf(largeCollapsedVcf = input, add_ratios = TRUE, field_DP = "AD")
  
  expect_s3_class(out, "data.frame")
  expect_equal(expect_df, out)
})

# ------------------------------------------------------------------------- #
# Modify DP and AD to introduce NA values
# ------------------------------------------------------------------------- #
g_dp <- VariantAnnotation::geno(input)$DP
g_ad <- VariantAnnotation::geno(input)$AD

# Introduce NA in proband DP at the first variant
g_dp[1, "proband"] <- NA
VariantAnnotation::geno(input)$DP <- g_dp

# Introduce NA in proband AD (both alleles NA) at the first variant
g_ad[1, "proband"][[1]] <- c(NA,NA)
VariantAnnotation::geno(input)$AD <- g_ad

expect_df$quality_proband[1] <- NA

# ------------------------------------------------------------------------- #
# Repeat the three tests under conditions where NA values are present
# ------------------------------------------------------------------------- #

test_that("Test if transformation of vcf with states to basic dataframe works", {
  out <- asDfVcf(largeCollapsedVcf = input, add_ratios = TRUE, field_DP = "DP")
  
  expect_s3_class(out, "data.frame")
  expect_equal(expect_df, out)
})

test_that("Test if transformation of vcf with states to basic dataframe works", {
  out <- asDfVcf(largeCollapsedVcf = input, add_ratios = TRUE)
  
  expect_s3_class(out, "data.frame")
  expect_equal(expect_df, out)
})

test_that("Test if transformation of vcf with states to basic dataframe works", {
  out <- asDfVcf(largeCollapsedVcf = input, add_ratios = TRUE, field_DP = "AD")
  
  expect_s3_class(out, "data.frame")
  expect_equal(expect_df, out)
})

