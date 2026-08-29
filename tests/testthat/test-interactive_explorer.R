file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

test_that(".hdx_comparison_table returns internally consistent, correctly-shaped rows", {
  a <- output_tp(file_nm)
  tbl <- HDXBoxeR:::.hdx_comparison_table(a, "Unbound", "bound", "1800.00s", 3, 0.01, 0.01)

  expect_s3_class(tbl, "data.frame")
  expect_true(all(c("peptide_id", "Start", "End", "Sequence", "state1_val", "state2_val",
                     "delta", "p_value", "rx", "ry", "ci", "is_sig") %in% colnames(tbl)))
  ## the join key must match the row it was built from - regression test
  ## for a real bug where sorting a display table separately from its key
  ## column desynced peptide_id from Sequence/Start/End/delta/p_value
  expect_equal(tbl$peptide_id, paste(tbl$Start, tbl$End, tbl$Sequence, sep = "_"))
  expect_equal(tbl$delta, tbl$state2_val - tbl$state1_val)
})

test_that("interactive functions require plotly/crosstalk (skipped if not installed)", {
  testthat::skip_if_not_installed("plotly")
  testthat::skip_if_not_installed("crosstalk")
  a <- output_tp(file_nm)

  v <- interactive_volcano(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
  expect_s3_class(v, "plotly")

  f <- interactive_fissure_vent(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
  expect_s3_class(f, "plotly")

  w <- interactive_woods(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
  expect_s3_class(w, "plotly")
})

test_that("significant_peptides_table's displayed rows match their own peptide_id key", {
  testthat::skip_if_not_installed("plotly")
  testthat::skip_if_not_installed("crosstalk")
  testthat::skip_if_not_installed("DT")
  a <- output_tp(file_nm)

  sd_obj <- HDXBoxeR:::.hdx_shared(a, "Unbound", "bound", "1800.00s", 3, 0.01, 0.01, NULL, "test")
  full <- sd_obj$data(withSelection = FALSE)

  t <- significant_peptides_table(a, state1 = "Unbound", state2 = "bound", times = "1800.00s",
                                   shared_data = sd_obj)
  expect_s3_class(t, "datatables")

  ## the underlying crosstalk-wrapped data (not just the DT widget object)
  ## must have Sequence/Start/End aligned to peptide_id, matching the
  ## non-table source-of-truth table built from the same shared_data -
  ## regression test for a real bug where the display table was sorted
  ## separately from the key column, desyncing every row
  view <- t$x$data
  expect_equal(nrow(view), nrow(full))
  expected_start <- full$Start[match(view$peptide_id, full$peptide_id)]
  expected_end <- full$End[match(view$peptide_id, full$peptide_id)]
  expected_seq <- full$Sequence[match(view$peptide_id, full$peptide_id)]
  expect_equal(view$Start, expected_start)
  expect_equal(view$End, expected_end)
  expect_equal(view$Sequence, expected_seq)
})
