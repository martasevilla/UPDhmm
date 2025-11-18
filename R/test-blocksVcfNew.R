# test-blocksVcfNew.R

expected_df <- data.frame(
  ID = "NA19685",
  seqnames = "6",
  start = 32489853,
  end = 42109975,
  group = as.character("iso_mat"),
  n_snps = 15L,
  geno_coded = I(setNames(list(c("133", "133", "121", "122", "133", "123", "122", "133", "321", "312", "231", "212", "323", "231", "123")), "1")),
  mean_quality_proband = 60.26667,
  mean_quality_mother = 60.13333,
  mean_quality_father = 59.06667
)


file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)

input <- vcfCheck(
  largeCollapsedVcf = input,
  father = "NA19689", mother = "NA19688",
  proband = "NA19685", check_quality = TRUE
)

S4Vectors::mcols(input)$states <- c("iso_mat", "iso_mat", "iso_mat")



test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(input)

    out <- as.data.frame(out)
    
    expected_no_ratio <- expected_df[, !(names(expected_df) %in% c("mean_quality_proband", "mean_quality_mother", "mean_quality_father"))]
    out_no_ratio <- out[, names(expected_no_ratio), drop = FALSE]
  
    expect_equal(out_no_ratio, expected_no_ratio)
    expect_s3_class(out, "data.frame")
})


test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(input, TRUE, "DP")

    out <- as.data.frame(out)
    
    out$mean_quality_proband <- round(out$mean_quality_proband,5)
    out$mean_quality_mother <- round(out$mean_quality_mother,5)
    out$mean_quality_father <- round(out$mean_quality_father,5)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})


test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(input, TRUE, "AD")

    out <- as.data.frame(out)
    
    out$mean_quality_proband <- round(out$mean_quality_proband,5)
    out$mean_quality_mother <- round(out$mean_quality_mother,5)
    out$mean_quality_father <- round(out$mean_quality_father,5)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(input, TRUE)

    out <- as.data.frame(out)
    
    out$mean_quality_proband <- round(out$mean_quality_proband,5)
    out$mean_quality_mother <- round(out$mean_quality_mother,5)
    out$mean_quality_father <- round(out$mean_quality_father,5)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
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

expected_df$mean_quality_proband <- 60.28571

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(input, TRUE, "DP")

    out <- as.data.frame(out)
    
    out$mean_quality_proband <- round(out$mean_quality_proband,5)
    out$mean_quality_mother <- round(out$mean_quality_mother,5)
    out$mean_quality_father <- round(out$mean_quality_father,5)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})


test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(input, TRUE, "AD")

    out <- as.data.frame(out)
    
    out$mean_quality_proband <- round(out$mean_quality_proband,5)
    out$mean_quality_mother <- round(out$mean_quality_mother,5)
    out$mean_quality_father <- round(out$mean_quality_father,5)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(input, TRUE)

    out <- as.data.frame(out)
    
    out$mean_quality_proband <- round(out$mean_quality_proband,5)
    out$mean_quality_mother <- round(out$mean_quality_mother,5)
    out$mean_quality_father <- round(out$mean_quality_father,5)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})
