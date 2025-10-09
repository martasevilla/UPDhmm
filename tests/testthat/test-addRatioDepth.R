expected_result <- data.frame(
     ID = "NA19675",
     start = 32489853,
     end = 32489876,
     group = "iso_mat",
     seqnames = "6",
     n_snps = 3,
     ratio_proband = 0.99,
     ratio_mother = 1.2 ,
     ratio_father = 0.93
   )


file <- system.file(package = "UPDhmm", "extdata", "test_het_mat.vcf.gz")
input <- VariantAnnotation::readVcf(file)
colnames(input) <- c("proband", "father", "mother")

  
position <- data.frame(
   ID = "NA19675",
   start = 32489853,
   end = 32489876,
   group = "iso_mat",
   seqnames = "6",
   n_snps = 3)
   

test_that("Test if calculation of statistic parameters works", {
  out <- addRatioDepth(
    filtered_def_blocks_states = position,
    largeCollapsedVcf = input
    
  )
  # round for errors in expect_equal
  out$ratio_proband <- round(out$ratio_proband, digits = 2)
  out$ratio_mother <- round(out$ratio_mother, digits = 2)
  out$ratio_father <- round(out$ratio_father, digits = 2)
  expect_equal(out, expected_result)
})
