file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

test_that("residue_regions returns correctly-shaped, internally consistent output", {
  a <- output_tp(file_nm, percent = FALSE)
  res <- residue_regions(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")

  expect_s3_class(res, "data.frame")
  expect_true(all(c("direction", "region_start", "region_end", "n_residues",
                     "n_peptides", "mean_per_residue_delta", "pooled_CI_Da") %in% colnames(res)))
  expect_true(all(res$direction %in% c("destabilized", "protected")))
  expect_true(all(res$region_end >= res$region_start))
  expect_equal(res$n_residues, res$region_end - res$region_start + 1)
  expect_true(all(res$n_peptides >= 2))  # default min_peptides = 2
  ## "destabilized" regions have positive per-residue delta, "protected" negative
  expect_true(all(res$mean_per_residue_delta[res$direction == "destabilized"] > 0))
  expect_true(all(res$mean_per_residue_delta[res$direction == "protected"] < 0))
})

test_that("residue_regions errors on an unmatched Deut.Time", {
  a <- output_tp(file_nm, percent = FALSE)
  expect_error(residue_regions(a, state1 = "Unbound", state2 = "bound", times = "not_a_time"))
})

test_that("raising min_peptides only shrinks or leaves regions unchanged", {
  a <- output_tp(file_nm, percent = FALSE)
  res_lenient <- residue_regions(a, state1 = "Unbound", state2 = "bound", times = "1800.00s", min_peptides = 2)
  res_strict  <- residue_regions(a, state1 = "Unbound", state2 = "bound", times = "1800.00s", min_peptides = 4)
  expect_true(sum(res_strict$n_residues) <= sum(res_lenient$n_residues))
})

test_that("residue_regions() correctly matches state names containing a space or a plus", {
  ## regression test: output_tp()'s wide-format columns embed each state
  ## name with spaces stripped and any other invalid character (e.g. "+")
  ## silently mangled to "." by data.frame()'s own make.names() - comparing
  ## resolve_raw_cols()/resolve_col() against the RAW state string directly
  ## used to error ("Expected 3 raw replicate columns... found 0") for any
  ## such state name. A realistic case: "Protein State 1" is HDX Examiner's
  ## own default, unedited state name.
  raw <- read.csv(file_nm)
  b5 <- raw[raw$Protein.State == "bound", ]
  b5$Protein.State <- "tubulin+b5"
  ps1 <- raw[raw$Protein.State == "Unbound", ]
  ps1$Protein.State <- "Protein State 1"
  raw2 <- rbind(ps1, b5)
  tmp <- tempfile(fileext = ".csv")
  write.csv(raw2, tmp, row.names = FALSE)

  a <- output_tp(tmp, states = c("Protein State 1", "tubulin+b5"), percent = FALSE)
  expect_no_error(res <- residue_regions(a, state1 = "Protein State 1", state2 = "tubulin+b5", times = "1800.00s"))
  expect_true(all(!is.na(res$n_residues)))
})

test_that("residue_regions_subset_summary() rolls up derived regions per named subset", {
  a <- output_tp(file_nm, percent = FALSE)
  res <- residue_regions(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
  summ <- residue_regions_subset_summary(res, regions = list(nterm = c(1, 100), cterm = c(101, 225)))

  expect_true(all(c("subset", "n_residues", "n_residues_significant", "frac_significant",
                     "n_regions_overlapping", "directions") %in% colnames(summ)))
  expect_equal(summ$n_residues, c(100, 125))
  expect_true(all(summ$n_residues_significant <= summ$n_residues))
  expect_true(all(summ$frac_significant >= 0 & summ$frac_significant <= 1))

  expect_error(residue_regions_subset_summary(data.frame(x = 1), regions = list(a = c(1, 2))))
  expect_error(residue_regions_subset_summary(res, regions = NULL))
})
