## fixes from the pre-release review

file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
f3 <- system.file("extdata", "three_states.csv", package = "HDXBoxeR")
key <- function(d) paste(d$Start, d$End, d$Sequence, d$Charge)

test_that("pv_timecourse() and prep_timecourse_plot_ave() rows describe the same peptides", {
  for (f in c(file_nm, f3)) {
    th <- suppressMessages(output_tc(f))
    st <- unique(th$Protein.State)
    c1 <- th[th$Protein.State == st[1], ]; v1 <- th[th$Protein.State == st[2], ]
    pv <- pv_timecourse(c1, v1); av <- prep_timecourse_plot_ave(c1, v1)
    expect_identical(key(pv), key(av[[1]]))
    expect_identical(key(av[[1]]), key(av[[2]]))
    expect_false(is.unsorted(av[[1]]$Start))
  }
})

test_that("robot significance matches a by-key recomputation", {
  th <- suppressMessages(output_tc(file_nm)); thP <- suppressMessages(output_tc(file_nm, percent = TRUE))
  c1 <- th[th$Protein.State == "Unbound", ]; v1 <- th[th$Protein.State == "bound", ]
  pv <- pv_timecourse(c1, v1)
  av <- prep_timecourse_plot_ave(c1, v1)
  ci <- prep_timecourse_plot_sd(c1, v1)
  m <- match(key(av[[1]]), key(pv))
  sig <- sapply(7:ncol(av[[1]]), function(i) {
    s <- pv[m, i] < 0.01 & abs(av[[1]][, i] - av[[2]][, i]) > ci[i - 6]; s & !is.na(s)
  })
  expected <- av[[1]][rowSums(sig) >= 1, c("Start", "End", "Sequence", "Charge")]
  got <- robot_indexes_df(thP = thP, th = th)
  expect_setequal(key(got), key(expected))
})

test_that("compare_states_tp() takes one time point", {
  expect_error(compare_states_tp(f3, times = c("60.00s", "1800.00s")), "one time point")
})

test_that("deprecated pv_cutoff sets alpha too, as in 0.0.2", {
  a <- suppressMessages(output_tp(file_nm))
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp); on.exit(grDevices::dev.off())
  expect_warning(plot_heat_map_tp(a, pv_cutoff = 0.05), "deprecated")
  ## same result as giving both
  sd1 <- sd_timepoint(a); av1 <- ave_timepoint(a); pv1 <- pv_timepoint(a)
  expect_warning(old <- significant_peptide_uptake(av1, pv1, sd1, pv_cutoff = 0.05), "deprecated")
  expect_equal(old, significant_peptide_uptake(av1, pv1, sd1, alpha = 0.05, pvalue = 0.05))
})

test_that("per-residue functions stay within Start-End with joined or modified sequences", {
  expect_equal(.hdx_clean_seq(c("ACDP.ACEP", "PEPTIDE+Phospho")), c("ACDP", "PEPTIDE"))
  r <- HDXBoxeR:::peptide_exch_residues(10, "ACDEF.ACDEG", end = 14)
  expect_true(all(r >= 12 & r <= 14))
  a <- suppressMessages(output_tp(f3, seq_match = TRUE, times = "60.00s"))
  expect_no_error(res <- residue_level_average(suppressMessages(output_tp(f3, seq_match = TRUE, times = "60.00s", percent = TRUE)),
                                               protein_length = 175))
})

test_that(".hdx_strip_suffix() keeps full state names", {
  expect_equal(.hdx_strip_suffix(c("av_Unbound_X..Deut", "av_bound_Deut..", "Unbound_Deut.._3", "sd_WT_X..Deut.x")),
               c("Unbound", "bound", "Unbound", "WT"))
})

test_that("output_prep() works with uneven replicates", {
  out <- suppressMessages(output_prep(f3, output_name = tempfile(fileext = ".csv")))
  expect_true(nrow(out) > 0)
  expect_setequal(unique(out$Protein.State), c("WT_pH7.5", "D109H_pH7.5", "R120G_ph7.5"))
  expect_false(any(is.nan(out[["Mean.Deut.Uptake_[Da]"]])))
})

test_that("output_prep() reports replicates per time point; hdx_summary_table() keeps the states order", {
  expect_message(output_prep(file_nm, output_name = tempfile(fileext = ".csv")), "Unbound=3")
  s <- suppressMessages(hdx_summary_table(file_nm, states = c("bound", "Unbound")))
  expect_equal(colnames(s)[2:3], c("bound", "Unbound"))
  expect_equal(s[s$Item == "Significant differences in HDX", "bound"], "reference state")
})

test_that("output_UD() keeps a missing value of a measured 0.00s control as NA", {
  r <- read.csv(file_nm, check.names = FALSE)
  i <- which(r$`Deut Time` == "0.00s")[1]
  r$`# Deut`[i] <- NA
  f <- tempfile(fileext = ".csv"); utils::write.csv(r, f, row.names = FALSE)
  ud <- output_UD(f)
  ## that peptide is not filled with 0 (it has no value, so it is left out)
  hit <- ud[ud$Protein.State == r$`Protein State`[i] & ud$Start == r$Start[i] &
              ud$End == r$End[i] & ud$Charge == r$Charge[i], 7:ncol(ud), drop = FALSE]
  expect_false(any(as.matrix(hit) == 0, na.rm = TRUE))
  expect_equal(nrow(ud), nrow(output_UD(file_nm)) - 1)
})

test_that("PyMOL scripts work with the default path", {
  a <- suppressMessages(output_tp(file_nm))
  old <- setwd(tempdir()); on.exit(setwd(old))
  expect_no_error(suppressMessages(pymol_script_significant_peptide(a)))
})

test_that("backHX_calculations() gives one value range per state", {
  b <- backHX_calculations(file_nm)
  expect_equal(nrow(b), 2)
  expect_true(all(b[, 2] <= b[, 1] & b[, 1] <= b[, 3]))
})
