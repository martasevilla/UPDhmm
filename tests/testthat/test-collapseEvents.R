# Input test dataframe
test_df <- data.frame(
    ID = c("S1", "S1", "S1", "S1", "S1", "S2", "S2"),
    seqnames = c("1", "1", "1", "1", "1", "2", "2"),
    start = c(50, 55, 100, 150, 300, 500, 550),
    end = c(70, 65, 120, 180, 320, 520, 580),
    group = c("iso_mat", "iso_mat", "iso_mat", "iso_mat", "het_pat", "iso_mat", "iso_mat"),
    n_mendelian_error = c(1, 5, 5, 10, 2, 50, 30),
    stringsAsFactors = FALSE
  )
  
# Expected result after collapsing
expected_result <- data.frame(
    ID = c("S1", "S1", "S2"),
    seqnames = c("1", "1", "2"),
    group = c("het_pat","iso_mat", "iso_mat"),
    n_events = c(1, 2, 2),
    total_mendelian_error = c(2, 15, 80),
    total_size = c(20, 50, 50),
    collapsed_events = c( "1:300-320","1:100-120,1:150-180", "2:500-520,2:550-580"),
    min_start = c(300,100, 500),
    max_end = c(320, 180, 580),
    stringsAsFactors = FALSE
  )
  



test_that("Test if calculation collapseEvents works correctly", {
  # Run function
  out <- collapseEvents(subset_df = test_df, min_ME = 2, min_size = 20)
  # Test equality
  expect_equal(out, expected_result)
})


expected_result <- data.frame(
  ID = c("NA19685", "NA19685"),
  seqnames = c("6", "15"),
  group = c("het_mat", "iso_mat"),
  n_events = c(1, 1),
  total_mendelian_error = c(3, 6),
  total_size = c(1010072, 19741113),
  collapsed_events = c("6:32489853-33499925", "15:22368862-42109975"),
  min_start = c(32489853,22368862),
  max_end = c(33499925, 42109975),
  ratio_proband = c(0.978982, 1.010509),
  ratio_mother = c(1.002257, 1.025959),
  ratio_father = c(0.951220, 0.997783)
)

file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)

input <- vcfCheck(
  largeCollapsedVcf = input,
  father = "NA19689", mother = "NA19688",
  proband = "NA19685", check_quality = TRUE
)

test_df <- calculateEvents(largeCollapsedVcf = input)

test_that("Test if calculation collapseEvents calculates mean read depths correctly with DP", {
  out <- collapseEvents(
    subset_df = test_df, 
    largeCollapsedVcf = input, 
    field_DP = "DP",
    min_ME = 2, 
    min_size = 20
  )
  
  expect_equal(out, expected_result, tolerance = 1e-6)
})

test_that("Test if calculation collapseEvents calculates mean read depths correctly with AD", {
  out <- collapseEvents(
    subset_df = test_df, 
    largeCollapsedVcf = input, 
    field_DP = "AD",
    min_ME = 2, 
    min_size = 20
  )
  
  expect_equal(out, expected_result, tolerance = 1e-6)
})

test_that("Test if calculation collapseEvents calculates mean read depths correctly with default field", {
  out <- collapseEvents(
    subset_df = test_df, 
    largeCollapsedVcf = input, 
    min_ME = 2, 
    min_size = 20
  )
  expect_equal(out, expected_result, tolerance = 1e-6)
})


g_dp <- VariantAnnotation::geno(input)$DP
g_ad <- VariantAnnotation::geno(input)$AD

# Introduce NA in proband DP at the first variant
g_dp[1, "proband"] <- NA
VariantAnnotation::geno(input)$DP <- g_dp

# Introduce NA in proband AD (both alleles NA) at the first variant
g_ad[1, "proband"][[1]] <- c(NA,NA)
VariantAnnotation::geno(input)$AD <- g_ad

expected_result <- data.frame(
  ID = c("NA19685", "NA19685"),
  seqnames = c("6", "15"),
  group = c("het_mat", "iso_mat"),
  n_events = c(1, 1),
  total_mendelian_error = c(3, 6),
  total_size = c(1010072, 19741113),
  collapsed_events = c("6:32489853-33499925", "15:22368862-42109975"),
  min_start = c(32489853,22368862),
  max_end = c(33499925, 42109975),
  ratio_proband = c(0.974526, 1.010190),
  ratio_mother  = c(1.002257, 1.025959),
  ratio_father  = c(0.951220, 0.997783)
)

test_that("Test if calculation collapseEvents calculates mean read depths correctly with DP", {
  out <- collapseEvents(
    subset_df = test_df, 
    largeCollapsedVcf = input, 
    field_DP = "DP",
    min_ME = 2, 
    min_size = 20
  )
  
  expect_equal(out, expected_result, tolerance = 1e-6)
})

test_that("Test if calculation collapseEvents calculates mean read depths correctly with AD", {
  out <- collapseEvents(
    subset_df = test_df, 
    largeCollapsedVcf = input, 
    field_DP = "AD",
    min_ME = 2, 
    min_size = 20
  )
  
  expect_equal(out, expected_result, tolerance = 1e-6)
})

test_that("Test if calculation collapseEvents calculates mean read depths correctly with default field", {
  out <- collapseEvents(
    subset_df = test_df, 
    largeCollapsedVcf = input, 
    min_ME = 2, 
    min_size = 20
  )
  expect_equal(out, expected_result, tolerance = 1e-6)
})