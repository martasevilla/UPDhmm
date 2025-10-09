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
  genotype = c(
    "133", "133", "121", "122", "133",
    "123", "122", "133", "321", "312",
    "231", "212", "323", "231", "123"
  )
)


file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)
SummarizedExperiment::colData(input)$ID <- colnames(input)
colnames(input) <- c("proband", "father", "mother")
S4Vectors::mcols(input)$states <- c("iso_mat", "iso_mat", "iso_mat")

genotypes <- c(
    "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
    "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
)



test_that("Test if transformation of vcf with states to basic dataframe works", {
    out <- asDfVcf(largeCollapsedVcf = input, genotypes = genotypes)
    expect_s3_class(out, "data.frame")
    expect_equal(expect_df, out)
})
