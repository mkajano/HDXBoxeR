file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

## compare_states_tp() needs 3+ states; the shipped example only has two
## ("Unbound", "bound"), so build a synthetic third by jittering "bound"
raw3_path <- local({
  raw <- read.csv(file_nm)
  synth <- raw[raw$Protein.State == "bound", ]
  synth$Protein.State <- "synthetic"
  set.seed(42)
  synth$X..Deut <- as.numeric(synth$X..Deut) + rnorm(nrow(synth), sd = 0.3)
  synth$Experiment <- gsub("^bound", "synthetic", synth$Experiment)
  raw3 <- rbind(raw, synth)
  tmp <- tempfile(fileext = ".csv")
  write.csv(raw3, tmp, row.names = FALSE)
  tmp
})

test_that("compare_states_tp() errors on fewer than three states", {
  expect_error(compare_states_tp(raw3_path, times = "1800.00s", states = c("Unbound", "bound")))
})

test_that("compare_states_tp() returns one row per peptide-by-pair, all three pairs by default", {
  res <- compare_states_tp(raw3_path, times = "1800.00s")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("state1", "state2", "delta", "p_value", "ci", "p_adj",
                     "p_omnibus", "p_omnibus_adj", "is_sig") %in% colnames(res)))
  pairs_seen <- unique(paste(res$state1, res$state2))
  expect_equal(length(pairs_seen), 3)
})

test_that("compare_states_tp() reference= only tests reference vs. each other state", {
  res <- compare_states_tp(raw3_path, times = "1800.00s", reference = "Unbound")
  expect_true(all(res$state1 == "Unbound"))
  expect_equal(sort(unique(res$state2)), c("bound", "synthetic"))
})

test_that("compare_states_tp() regions= splits output and tags a region column", {
  res <- compare_states_tp(raw3_path, times = "1800.00s",
                            regions = list(low = c(1, 100), high = c(101, 225)))
  expect_true("region" %in% colnames(res))
  expect_setequal(unique(res$region), c("low", "high"))
  expect_true(all(res$Start[res$region == "low"] >= 1 & res$End[res$region == "low"] <= 100))
  expect_true(all(res$Start[res$region == "high"] >= 101 & res$End[res$region == "high"] <= 225))
})

test_that("compare_states_tp() reproduces plot_fissure_vent()'s delta/p-value for a known peptide", {
  ## regression check: the Unbound-vs-bound pair inside the 3-state call
  ## should exactly match the two-state pipeline's own numbers (peptide
  ## 98-115 at 1800.00s, cross-checked by hand against interactive_explorer.R)
  res <- compare_states_tp(raw3_path, times = "1800.00s")
  row <- res[res$state1 == "Unbound" & res$state2 == "bound" & res$Start == 98 & res$End == 115, ]
  expect_equal(nrow(row), 1)
  expect_equal(row$delta, -0.495, tolerance = 1e-6)
  expect_equal(row$p_value, 0.003094462, tolerance = 1e-6)
})

test_that("compare_states_tp() by = \"position\" errors as not-yet-implemented", {
  expect_error(compare_states_tp(raw3_path, times = "1800.00s", by = "position"))
})

test_that("compare_states_unmatched() requires regions", {
  expect_error(compare_states_unmatched(file_nm, times = "1800.00s"))
})

test_that("compare_states_unmatched() returns summary and coverage with expected columns", {
  res <- compare_states_unmatched(file_nm, times = "1800.00s", regions = list(all = c(1, 225)))
  expect_type(res, "list")
  expect_true(all(c("summary", "coverage") %in% names(res)))
  expect_true(all(c("region", "state1", "state2", "shift_hl", "p_value", "p_omnibus", "p_adj") %in%
                    colnames(res$summary)))
  expect_true(all(c("region", "state", "n_peptides", "frac_shared_all_states") %in% colnames(res$coverage)))
  expect_true(res$coverage$frac_shared_all_states[1] > 0 && res$coverage$frac_shared_all_states[1] <= 1)
})

test_that("compare_states_unmatched() method = \"classic\" warns about pseudo-replication", {
  expect_warning(compare_states_unmatched(file_nm, times = "1800.00s", regions = list(all = c(1, 225)),
                                           method = "classic"),
                  "pseudo-replicat|correlated|inflates")
})

test_that("compare_states_unmatched() method = \"resample\" returns a bootstrap CI", {
  set.seed(1)
  res <- compare_states_unmatched(file_nm, times = "1800.00s", regions = list(all = c(1, 225)),
                                   method = "resample", n_perm = 100, n_boot = 100)
  expect_false(is.na(res$summary$ci_low[1]))
  expect_false(is.na(res$summary$ci_high[1]))
  expect_true(res$summary$ci_low[1] <= res$summary$ci_high[1])
})

test_that("plot_coverage_diagnostics() runs without error and errors on bad input", {
  res <- compare_states_unmatched(file_nm, times = "1800.00s",
                                   regions = list(nterm = c(1, 100), cterm = c(101, 225)))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_no_error(plot_coverage_diagnostics(res))
  grDevices::dev.off()

  expect_error(plot_coverage_diagnostics(list(summary = res$summary)))
})
