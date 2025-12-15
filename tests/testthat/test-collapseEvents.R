# Input test dataframe
test_df <- data.frame(
    ID = c("S1", "S1", "S1", "S1", "S1", "S2", "S2"),
    chromosome = c("1", "1", "1", "1", "1", "2", "2"),
    start = c(50, 75, 100, 150, 300, 500, 550),
    end = c(70, 85, 120, 180, 320, 520, 580),
    n_snps = c(5, 3, 8, 10, 6, 12, 7),
    group = c("iso_mat", "iso_mat", "iso_mat", "iso_mat", "het_pat", "iso_mat", "iso_mat"),
    n_mendelian_error = c(1, 5, 5, 10, 2, 50, 30),
    ratio_proband = rep(NA_real_, 7),
    ratio_mother  = rep(NA_real_, 7),
    ratio_father  = rep(NA_real_, 7),
    stringsAsFactors = FALSE
  )


# Expected result after collapsing
expected_result <- data.frame(
    ID = c("S1", "S1", "S2"),
    chromosome = c("1", "1", "2"),
    start = c(300,100, 500),
    end = c(320, 180, 580),
    group = c("het_pat","iso_mat", "iso_mat"),
    n_events = c(1, 2, 2),
    total_mendelian_error = c(2, 15, 80),
    total_size = c(20, 50, 50),
    total_snps = c(6, 18, 19),
    prop_covered = c(1, 0.625, 0.625),
    ratio_proband = rep(NA_real_, 3),
    ratio_mother  = rep(NA_real_, 3),
    ratio_father  = rep(NA_real_, 3),
    collapsed_events = c( "1:300-320","1:100-120,1:150-180", "2:500-520,2:550-580"),
    stringsAsFactors = FALSE
  )
  

expected_cols <- c(
  "ID", "chromosome", "start", "end", "group", 
  "n_events", "total_mendelian_error", "total_size",
  "total_snps", "prop_covered", "ratio_proband", 
  "ratio_mother", "ratio_father", "collapsed_events"
)

test_that("Test if collapseEvents returns empty df with correct structure when all events are filtered out", {
  # Run function
  out <- collapseEvents(subset_df = test_df, min_ME = 2, min_size = 200)

  # Must be empty
  expect_equal(nrow(out), 0)
  
  # Must have correct column structure
  expect_equal(colnames(out), expected_cols)
  
}) 

test_that("Test if calculation collapseEvents works correctly (no ratios)", {
  # Run function
  out <- collapseEvents(subset_df = test_df, min_ME = 2, min_size = 20)
  # Test equality
  expect_equal(out, expected_result)
})


# Input test dataframe
test_df <- data.frame(
  ID = c("S1", "S1", "S1", "S1", "S1", "S2", "S2"),
  chromosome = c("1", "1", "1", "1", "1", "2", "2"),
  start = c(50, 75, 100, 150, 300, 500, 550),
  end = c(70, 85, 120, 180, 320, 520, 580),
  n_snps = c(5, 3, 8, 10, 6, 12, 7),
  group = c("iso_mat", "iso_mat", "iso_mat", "iso_mat", "het_pat", "iso_mat", "iso_mat"),
  n_mendelian_error = c(1, 5, 5, 10, 2, 50, 30),
  ratio_proband = c(0.98, 1.01, 0.99, 1.02, 0.97, 1.03, 1.01),
  ratio_mother  = c(1.00, 1.03, 1.01, 1.04, 0.98, 1.05, 1.02),
  ratio_father  = c(0.95, 0.97, 0.96, 0.98, 0.99, 1.00, 0.99),
  stringsAsFactors = FALSE
)


# Expected result after collapsing
expected_result <- data.frame(
  ID = c("S1", "S1", "S2"),
  chromosome = c("1", "1", "2"),
  start = c(300,100, 500),
  end = c(320, 180, 580),
  group = c("het_pat","iso_mat", "iso_mat"),
  n_events = c(1, 2, 2),
  total_mendelian_error = c(2, 15, 80),
  total_size = c(20, 50, 50),
  total_snps = c(6, 18, 19),
  prop_covered = c(1, 0.625, 0.625),
  ratio_proband = c(0.97, 1.01, 1.02),
  ratio_mother  = c(0.98, 1.03, 1.04),
  ratio_father  = c(0.99, 0.97, 1.00),
  collapsed_events = c( "1:300-320","1:100-120,1:150-180", "2:500-520,2:550-580"),
  stringsAsFactors = FALSE
)

test_that("Test if collapseEvents returns empty df with correct structure when all events are filtered out and ratios are present", {
  # Run function
  out <- collapseEvents(subset_df = test_df, min_ME = 2, min_size = 200)
  
  # Must be empty
  expect_equal(nrow(out), 0)
  
  # Must have correct column structure
  expect_equal(colnames(out), expected_cols)
  
})   
  
test_that("Test if calculation collapseEvents works correctly (with ratios)", {
  # Run function
  out <- collapseEvents(subset_df = test_df, min_ME = 2, min_size = 20)
  # Test equality
  expect_equal(out, expected_result, tolerance = 1e-2)
})

