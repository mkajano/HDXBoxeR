## compare_timecourse(): summed differences and per-time counts

file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

test_that("summed delta, its SE and the per-time calls match direct computation", {
  tc <- compare_timecourse(file_nm, "Unbound", "bound")
  expect_true(all(c("summed_delta", "se", "p_value", "ci_sum", "sig_summed", "n_timepoints_sig",
                    "sig_min_timepoints") %in% colnames(tc)))
  a <- suppressMessages(output_tp(file_nm))
  pep <- a[a$Start == 98 & a$End == 115, ]
  u <- as.matrix(pep[, paste0("Unbound_X..Deut_", 1:3)]); b <- as.matrix(pep[, paste0("bound_X..Deut_", 1:3)])
  row <- tc[tc$Start == 98 & tc$End == 115, ]
  expect_equal(row$summed_delta, sum(rowMeans(b) - rowMeans(u)), tolerance = 1e-8)
  expect_equal(row$se, sqrt(sum(apply(u, 1, var) / 3 + apply(b, 1, var) / 3)), tolerance = 1e-8)

  ## per-time calls agree with the two-state rule used by plot_fissure_vent()
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  fv <- suppressMessages(plot_fissure_vent(a, times = "60.00s", state1 = "Unbound", state2 = "bound"))
  grDevices::dev.off()
  m <- merge(fv[, c("Start", "End", "Sequence", "is_sig")], tc[, c("Start", "End", "Sequence", "sig_60.00s")])
  expect_equal(m$is_sig, m$`sig_60.00s`)
  expect_equal(tc$sig_min_timepoints, tc$n_timepoints_sig >= 2)
})

test_that("compare_timecourse() p_adjust only makes calls stricter", {
  raw <- compare_timecourse(file_nm, "Unbound", "bound")
  bh <- compare_timecourse(file_nm, "Unbound", "bound", p_adjust = "BH")
  expect_true(all(bh$p_adj >= raw$p_adj, na.rm = TRUE))
  expect_lte(sum(bh$sig_summed), sum(raw$sig_summed))
})
