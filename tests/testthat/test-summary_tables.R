file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

test_that("backHX_calculations() uses measured centroid x charge", {
  ## regression test: back-exchange used to be computed from the theoretical
  ## centroid m/z difference with no charge factor (~66-70% on this file)
  res <- backHX_calculations(file_nm)
  expect_match(attr(res, "method"), "Exp Cent")

  raw <- read.csv(file_nm)
  keys <- c("Protein.State", "Start", "End", "Sequence", "Charge")
  m <- merge(raw[raw$Deut.Time == "0s", c(keys, "Exp.Cent")],
             raw[raw$Deut.Time == "FD", c(keys, "Exp.Cent")], by = keys)
  m <- m[m$Protein.State == "Unbound", ]
  expected <- mean((1 - (m$Exp.Cent.y - m$Exp.Cent.x) * m$Charge / nb_exch_deut(m) / 0.85) * 100)
  expect_equal(as.numeric(res["Unbound", "<BackHX>"]), round(expected, 2))
  expect_true(all(as.numeric(res[, "<BackHX>"]) < 50))
})

test_that("general_info() reports the replicate count output_tp() detects", {
  ## regression test: was n(Experiment) / n(Deut.Time) = 0.75 on this file
  gi <- general_info(file_nm)
  expect_equal(as.numeric(gi[["# replicates"]]), c(3, 3))
})

test_that("summary_sd_CI() honours replicates", {
  ## regression test: output_tp() was called without replicates, giving a
  ## bogus extra "state" row when replicates != 3
  res <- summary_sd_CI(file_nm, replicates = 2)
  expect_equal(nrow(res), 2)
})
