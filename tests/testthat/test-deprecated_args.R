## renamed arguments: pv_cutoff -> pvalue, cola / colors_initial -> col

file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

test_that("old argument names still work, with a deprecation warning, and give the same result", {
  a <- suppressMessages(output_tp(file_nm))
  ## as in 0.0.2, pv_cutoff sets both the p-value cutoff and alpha
  new <- diff_peptides_prep(a, pvalue = 0.05, alpha = 0.05)
  expect_warning(old <- diff_peptides_prep(a, pv_cutoff = 0.05), "'pv_cutoff' is deprecated; use 'pvalue'")
  expect_equal(old, new)
  expect_warning(old2 <- diff_peptides_prep(a, pv_cutoff = 0.05, alpha = 0.01), "deprecated")
  expect_equal(old2, diff_peptides_prep(a, pvalue = 0.05, alpha = 0.01))

  expect_equal(suppressWarnings(color_ranges_Spectral(c(-Inf, 0, Inf), colors_initial = "white")),
               color_ranges_Spectral(c(-Inf, 0, Inf), col = "white"))
  expect_warning(color_ranges_Blue_Red_heat_map(c(-Inf, 0, Inf), colors_initial = "white"), "deprecated")

  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  expect_warning(plots_av_tp(a, cola = c(1, 2)), "'cola' is deprecated; use 'col'")
  expect_no_warning(plots_av_tp(a, col = c(1, 2)))
  grDevices::dev.off()
})

test_that("positional calls are unaffected by the renames", {
  a <- suppressMessages(output_tp(file_nm))
  av <- ave_timepoint(a); pv <- pv_timepoint(a); s <- sd_timepoint(a)
  expect_equal(significant_peptide_uptake(av, pv, s, 0.01, 0.01, 3),
               significant_peptide_uptake(av, pv, s, alpha = 0.01, pvalue = 0.01, replicates = 3))
})

test_that("PyMOL scripts shift residue numbers by resno_offset", {
  a <- suppressMessages(output_tp(file_nm, times = "60.00s"))
  d0 <- file.path(tempdir(), "pym0"); d1 <- file.path(tempdir(), "pym1")
  dir.create(d0, showWarnings = FALSE); dir.create(d1, showWarnings = FALSE)
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  suppressMessages(pymol_script_significant_peptide(a, path = d0, pvalue = 0.05))
  suppressMessages(pymol_script_significant_peptide(a, path = d1, pvalue = 0.05, resno_offset = -10))
  grDevices::dev.off()
  f0 <- list.files(d0, full.names = TRUE)[1]; f1 <- list.files(d1, full.names = TRUE)[1]
  r0 <- regmatches(readLines(f0), regexpr("resi [0-9-]+", readLines(f0)))
  r1 <- regmatches(readLines(f1), regexpr("resi [0-9-]+", readLines(f1)))
  start0 <- as.numeric(sub("resi ([0-9]+)-.*", "\\1", r0[1])); start1 <- as.numeric(sub("resi ([0-9]+)-.*", "\\1", r1[1]))
  expect_equal(start1, start0 - 10)
})

test_that("three_states example data reads, with uneven replicates and mutated sequences", {
  f3 <- system.file("extdata", "three_states.csv", package = "HDXBoxeR")
  expect_true(nzchar(f3))
  a <- suppressMessages(output_tp(f3))
  expect_setequal(names(attr(a, "replicates")), c("WT_pH7.5", "D109H_pH7.5", "R120G_ph7.5"))
  pos <- suppressMessages(output_tp(f3, seq_match = TRUE))
  expect_true(any(grepl("\\.", pos$Sequence)))
  expect_gt(nrow(pos), nrow(a))
})
