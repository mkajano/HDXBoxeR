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
