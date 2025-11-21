# test-calculateEvents.R

# Expected UPD event blocks 
expected_def_blocks <- data.frame(
  ID = c("NA19685", "NA19685"),
  seqnames = c("6", "15"),
  start = c(32489853, 22368862),
  end = c(33499925, 42109975),
  group = c("het_mat", "iso_mat"),
  n_snps = c(5, 10),
  ratio_proband = c(0.974526, 1.010190),
  ratio_mother = c(1.002257, 1.025959),
  ratio_father = c(0.951220, 0.997783),
  n_mendelian_error = c(3, 6)
)

file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)

input <- vcfCheck(
  largeCollapsedVcf = input,
  father = "NA19689", mother = "NA19688",
  proband = "NA19685", check_quality = TRUE
)

expected_mean <- c(proband = 904/15, mother = 886/15, father = 902/15)

# ------------------------------------------------------------------------- #
# Test computeTrioTotals using DP explicitly
# ------------------------------------------------------------------------- #
test_that("computeTrioTotals calculates mean read depths correctly with DP", {
  mean_dp <- computeTrioTotals(vcf = input, field_DP = "DP")

  expect_true(is.numeric(mean_dp))
  expect_named(mean_dp, c("proband", "mother", "father"))
  expect_equal(mean_dp, expected_mean)
})

# ------------------------------------------------------------------------- #
# Test computeTrioTotals using AD explicitly
# ------------------------------------------------------------------------- #
test_that("computeTrioTotals calculates mean read depths correctly with AD", {
  mean_dp <- computeTrioTotals(vcf = input, field_DP = "AD")

  expect_true(is.numeric(mean_dp))
  expect_named(mean_dp, c("proband", "mother", "father"))
  expect_equal(mean_dp, expected_mean)
})

# ------------------------------------------------------------------------- #
# Test computeTrioTotals without specifying a field
# The function should automatically fall back to DP → AD
# ------------------------------------------------------------------------- #
test_that("computeTrioTotals calculates mean read depths correctly with default field", {
  mean_dp <- computeTrioTotals(vcf = input)

  expect_true(is.numeric(mean_dp))
  expect_named(mean_dp, c("proband", "mother", "father"))
  expect_equal(mean_dp, expected_mean)
})

# ------------------------------------------------------------------------ #
# Modify DP and AD to introduce NA values
# This tests whether computeTrioTotals properly ignores NA values
# ------------------------------------------------------------------------ #

g_dp <- VariantAnnotation::geno(input)$DP
g_ad <- VariantAnnotation::geno(input)$AD

# Introduce NA in proband DP at the first variant
g_dp[1, "proband"] <- NA
VariantAnnotation::geno(input)$DP <- g_dp

# Introduce NA in proband AD (both alleles NA) at the first variant
g_ad[1, "proband"][[1]] <- c(NA,NA)
VariantAnnotation::geno(input)$AD <- g_ad

# Expected sums and valid counts after introducing NA values
expected_mean <- c(proband = 844/14, mother = 886/15, father = 902/15)

# ------------------------------------------------------------------------- #
# Repeat the three tests under conditions where NA values are present
# ------------------------------------------------------------------------- #
test_that("computeTrioTotals calculates mean read depths correctly with DP", {
  mean_dp <- computeTrioTotals(vcf = input, field_DP = "DP")

  expect_true(is.numeric(mean_dp))
  expect_named(mean_dp, c("proband", "mother", "father"))
  expect_equal(mean_dp, expected_mean)
})

test_that("computeTrioTotals calculates mean read depths correctly with AD", {
  mean_dp <- computeTrioTotals(vcf = input, field_DP = "AD")

  expect_true(is.numeric(mean_dp))
  expect_named(mean_dp, c("proband", "mother", "father"))
  expect_equal(mean_dp, expected_mean)
})

test_that("computeTrioTotals calculates mean read depths correctly with default field", {
  mean_dp <- computeTrioTotals(vcf = input)

  expect_true(is.numeric(mean_dp))
  expect_named(mean_dp, c("proband", "mother", "father"))
  expect_equal(mean_dp, expected_mean)
})

# ------------------------------------------------------------------------- #
# Test calculateEvents() using default HMM (add_ratios = FALSE)
# ------------------------------------------------------------------------- #
test_that("Test if the general function works (default HMM, add_ratios = FALSE)", {
  out <- calculateEvents(largeCollapsedVcf = input)
  
  out$seqnames <- as.character(out$seqnames)
  out <- as.data.frame(out)
  
  # Should not contain ratio columns when add_ratios = FALSE
  expect_false(any(c("ratio_proband", "ratio_mother", "ratio_father") %in% names(out)))
  expected_no_ratio <- expected_def_blocks[, !(names(expected_def_blocks) %in% c("ratio_proband", "ratio_mother", "ratio_father"))]
  
  # Compare structural output against expected UPD blocks
  expect_equal(out[, names(expected_no_ratio)], expected_no_ratio)
  expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Test calculateEvents() with default HMM and add_ratios = TRUE
# ------------------------------------------------------------------------- #
test_that("Test if the general function works (default HMM, add_ratios = TRUE)", {
    out <- calculateEvents(largeCollapsedVcf = input, add_ratios = TRUE, field_DP = "DP")
    
    out$seqnames <- as.character(out$seqnames)
    
    out <- as.data.frame(out)
    expect_equal(expected_def_blocks, out, tolerance = 1e-6)
    expect_s3_class(out, "data.frame")
})

# Custom HMM definition to test flexibility of calculateEvents()
new_hmm<-list(
  States = c("normal", "iso_fat", "iso_mat", "het_fat", "het_mat")
  ,
  Symbols = c("111", "112", "113", "121", "122", "123", "131", "132", "133", 
              "211", "212", "213", "221", "222", "223", "231", "232", "233", 
              "311", "312", "313", "321", "322", "323", "331", "332", "333", 
              "000")
  ,
  startProbs = c(normal = 0.996, iso_fat = 0.001, iso_mat = 0.001, 
                 het_fat = 0.001, het_mat = 0.001)
  ,
  transProbs = matrix(
    c(
      0.99996, 0.00001, 0.00001, 0.00001, 0.00001,
      0.00001, 0.99996, 0.00001, 0.00001, 0.00001,
      0.00001, 0.00001, 0.99996, 0.00001, 0.00001,
      0.00001, 0.00001, 0.00001, 0.99996, 0.00001,
      0.00001, 0.00001, 0.00001, 0.00001, 0.99996
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      from = c("normal", "iso_fat", "iso_mat", "het_fat", "het_mat"),
      to = c("normal", "iso_fat", "iso_mat", "het_fat", "het_mat")
    )
  )
  ,
  emissionProbs = matrix(
    c(0.09250,0.00001,0.00001,0.09250,0.09250,0.00001,0.00001,0.09250,0.00001,0.09250,0.09250,0.00001,0.09250,0.12500,0.09250,0.00001,0.09250,0.09250,0.00001,0.09250,0.00001,0.00001,0.09250,0.09250,0.00001,0.00001,0.09250,0.00006,0.09250,0.00001,0.00001,0.12500,0.00001,0.00001,0.09250,0.00001,0.00001,0.09250,0.00001,0.09250,0.12500,0.00001,0.12500,0.09250,0.00001,0.09250,0.00001,0.00001,0.09250,0.00001,0.00001,0.12500,0.00001,0.00001,0.09250,0.00003,0.09250,0.00001,0.00001,0.09250,0.00001,0.09250,0.00001,0.00001,0.09250,0.12500,0.00001,0.00001,0.12500,0.00001,0.12500,0.00001,0.00001,0.12500,0.09250,0.00001,0.00001,0.09250,0.00001,0.09250,0.00001,0.00001,0.09250,0.00003,0.09250,0.00001,0.00001,0.12500,0.00001,0.00001,0.09250,0.00001,0.00001,0.00001,0.12500,0.00001,0.00001,0.25000,0.00001,0.00001,0.12500,0.00001,0.00001,0.00001,0.09250,0.00001,0.00001,0.12500,0.00001,0.00001,0.09250,0.00000,0.09250,0.00001,0.00001,0.00001,0.12500,0.00001,0.00001,0.00001,0.09250,0.12500,0.00001,0.00001,0.00001,0.25000,0.00001,0.00001,0.00001,0.12500,0.09250,0.00001,0.00001,0.00001,0.12500,0.00001,0.00001,0.00001,0.09250,0.00000),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      states = c("normal", "iso_fat", "iso_mat", "het_fat", "het_mat"),
      symbols = c(
        "111", "112", "113", "121", "122", "123", "131", "132", "133",
        "211", "212", "213", "221", "222", "223", "231", "232", "233",
        "311", "312", "313", "321", "322", "323", "331", "332", "333",
        "000"
      )
    )
  )
)

expected_def_blocks <- data.frame(
  ID = c("NA19685", "NA19685"),
  seqnames = c("6", "15"),
  start = c(32489853, 22368862),
  end = c(33499925, 42109975),
  group = c("het_mat", "iso_mat"),
  n_snps = c(5, 10),
  ratio_proband = c(0.978982, 1.010509),
  ratio_mother = c(1.002257, 1.025959),
  ratio_father = c(0.951220, 0.997783),
  n_mendelian_error = c(3, 6)
)

  
file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)

input <- vcfCheck(
  largeCollapsedVcf = input,
  father = "NA19689", mother = "NA19688",
  proband = "NA19685", check_quality = TRUE
)

# ------------------------------------------------------------------------- #
# Test calculateEvents() using custom HMM, add_ratios = FALSE
# ------------------------------------------------------------------------- #

test_that("Test if the general function works (custom HMM, add_ratios = FALSE)", {
  out <- calculateEvents(largeCollapsedVcf = input, hmm = new_hmm, field_DP = "DP")
  out$seqnames <- as.character(out$seqnames)
  
  expect_false(any(c("ratio_proband", "ratio_mother", "ratio_father") %in% names(out)))
  expected_no_ratio <- expected_def_blocks[, !(names(expected_def_blocks) %in% c("ratio_proband", "ratio_mother", "ratio_father"))]
  expect_equal(out[, names(expected_no_ratio)], expected_no_ratio)
  
  out <- as.data.frame(out)
  expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Test calculateEvents() using custom HMM, add_ratios = TRUE
# ------------------------------------------------------------------------- #

test_that("Test if the general function works (custom HMM, add_ratios = TRUE)", {
  out <- calculateEvents(largeCollapsedVcf = input, hmm = new_hmm, field_DP = "DP", add_ratios = TRUE)
  out$seqnames <- as.character(out$seqnames)
  
  out <- as.data.frame(out)
  expect_equal(expected_def_blocks, out, tolerance = 1e-6)
  expect_s3_class(out, "data.frame")
})
