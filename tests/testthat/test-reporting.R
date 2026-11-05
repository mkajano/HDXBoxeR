## Phase 5: reporting and provenance

file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

test_that("hdx_results_table() agrees with plot_fissure_vent() and records its parameters", {
  res <- hdx_results_table(file_nm)
  expect_true(all(c("Deut.Time", "state1", "state2", "mean1", "sd1", "n1", "delta", "p_value",
                    "p_adj", "ci", "is_sig") %in% colnames(res)))
  a <- suppressMessages(output_tp(file_nm))
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  fv <- suppressMessages(plot_fissure_vent(a, times = "1800.00s", state1 = "Unbound", state2 = "bound"))
  grDevices::dev.off()
  r <- res[res$Deut.Time == "1800.00s", ]
  m <- merge(r[, c("Start", "End", "Sequence", "is_sig", "ci")], fv[, c("Start", "End", "Sequence", "is_sig", "ci")],
             by = c("Start", "End", "Sequence"))
  expect_equal(m$is_sig.x, m$is_sig.y)
  expect_equal(m$ci.x, m$ci.y)

  st <- attr(res, "hdxboxer")
  expect_equal(st$parameters$p_adjust, "none")
  expect_equal(st$package_version, as.character(utils::packageVersion("HDXBoxeR")))
})

test_that("hdx_results_table(file =) writes the csv and a parameters file", {
  f <- tempfile(fileext = ".csv")
  hdx_results_table(file_nm, file = f, p_adjust = "BH")
  expect_true(file.exists(f))
  pf <- sub("\\.csv$", "_parameters.txt", f)
  expect_true(file.exists(pf))
  txt <- readLines(pf)
  expect_true(any(grepl("p_adjust: BH", txt)))
  expect_true(any(grepl("sessionInfo", txt)))
  expect_equal(nrow(read.csv(f)), nrow(hdx_results_table(file_nm)))
})

test_that("hdx_summary_table() puts metadata first and computes the Masson items", {
  s <- hdx_summary_table(file_nm, metadata = list("HDX reaction details" = "PBS pD 7.4",
                                                  "Temperature" = c("20 C", "21 C")))
  expect_equal(colnames(s), c("Item", "Unbound", "bound"))
  expect_equal(s$Item[1:2], c("HDX reaction details", "Temperature"))
  expect_equal(s$bound[2], "21 C")
  expect_true(all(c("Number of peptides", "Sequence coverage (%)", "Back-exchange (%): mean (range)",
                    "Significant differences in HDX") %in% s$Item))
  expect_error(hdx_summary_table(file_nm, metadata = list(x = c("a", "b", "c"))), "1 or 2 values")
  md <- tempfile(fileext = ".md")
  hdx_summary_table(file_nm, file = md)
  expect_match(readLines(md)[1], "^\\| Item \\|")
})

test_that("general_info() reports peptide length as the number of residues", {
  gi <- general_info(file_nm)
  a <- suppressMessages(output_tp(file_nm, states = "Unbound", times = "3.00s"))
  expect_equal(as.numeric(gi[["<Peptide Length>"]][1]), round(mean(a$End - a$Start + 1), 2))
})

test_that("compare_states_unmatched(seed =) is reproducible and leaves the caller's RNG alone", {
  set.seed(99); before <- runif(1)
  set.seed(99)
  r1 <- compare_states_unmatched(file_nm, times = "1800.00s", regions = list(all = c(1, 225)),
                                 method = "resample", n_perm = 50, n_boot = 50, seed = 1)
  after <- runif(1)
  r2 <- compare_states_unmatched(file_nm, times = "1800.00s", regions = list(all = c(1, 225)),
                                 method = "resample", n_perm = 50, n_boot = 50, seed = 1)
  expect_equal(r1$summary$p_value, r2$summary$p_value)
  expect_equal(r1$summary$ci_low, r2$summary$ci_low)
  expect_equal(after, before)
  expect_equal(attr(r1, "hdxboxer")$parameters$seed, 1)
})

test_that("render_report() produces an HTML report and the results files", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  out <- file.path(tempdir(), "test_report.html")
  suppressMessages(render_report(file_nm, output_file = out, times = c("60.00s", "1800.00s"),
                                 regions = list(nterm = c(1, 100))))
  expect_true(file.exists(out))
  expect_true(file.exists(file.path(tempdir(), "test_report_results.csv")))
  h <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(h, "Data quality")
  expect_match(h, "Regions")
})
