## p_adjust / min_timepoints / shared significance rule

file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
tp <- function(...) suppressMessages(output_tp(...))

test_that("pv_timepoint(p_adjust) adjusts within each time point, leaving the reference column", {
  a <- tp(file_nm)
  raw <- pv_timepoint(a); bh <- pv_timepoint(a, p_adjust = "BH")
  t <- a$Deut.Time == "60.00s"
  expect_equal(bh[t, 8], p.adjust(raw[t, 8], "BH"))
  expect_equal(bh[, 7], raw[, 7])
  expect_error(pv_timepoint(a, p_adjust = "nonsense"))
})

test_that("plot_fissure_vent() returns what it drew, with p_adjust applied to the outline", {
  a <- tp(file_nm)
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  raw <- suppressMessages(plot_fissure_vent(a, times = "60.00s", state1 = "Unbound", state2 = "bound"))
  bh <- suppressMessages(plot_fissure_vent(a, times = "60.00s", state1 = "Unbound", state2 = "bound",
                                           p_adjust = "BH"))
  grDevices::dev.off()
  expect_true(all(c("delta", "p_value", "p_adj", "ci", "is_sig") %in% colnames(raw)))
  expect_equal(raw$p_adj, raw$p_value)
  expect_equal(bh$p_adj, p.adjust(bh$p_value, "BH"))
  expect_equal(raw$is_sig, !is.na(raw$p_value) & abs(raw$delta) > raw$ci & raw$p_value < 0.01)
  expect_lte(sum(bh$is_sig), sum(raw$is_sig))
})

test_that("robot_indexes(min_timepoints) = 1 keeps the old selection; higher values select a subset", {
  th <- suppressMessages(output_tc(file_nm)); thP <- suppressMessages(output_tc(file_nm, percent = TRUE))
  any_t <- robot_indexes(thP = thP, th = th)
  two_t <- robot_indexes(thP = thP, th = th, min_timepoints = 2)
  expect_true(length(any_t) > 0)
  expect_true(all(two_t %in% any_t))
  expect_lte(length(two_t), length(any_t))
  df1 <- robot_indexes_df(thP = thP, th = th, min_timepoints = 2)
  expect_equal(nrow(df1), length(two_t))
})

test_that("residue_regions(): package-wide rule by default; fdr_cutoff reproduces the old rule with a warning", {
  a <- tp(file_nm)
  new <- residue_regions(a, "Unbound", "bound", times = "1800.00s")
  expect_warning(old <- residue_regions(a, "Unbound", "bound", times = "1800.00s", fdr_cutoff = 0.05),
                 "deprecated")
  expect_true(nrow(old) >= nrow(new))
  expect_true(all(new$n_peptides >= 2))
})

test_that("residue_regions() counts charge states of one peptide once", {
  f <- make_hdx_csv(states = c("A", "B"))
  x <- read.csv(f, check.names = FALSE)
  ## a second charge state of peptide 12-19 must not corroborate itself
  extra <- x[x$Start == 12, ]; extra$Charge <- extra$Charge + 1
  utils::write.csv(rbind(x, extra), f, row.names = FALSE)
  a <- tp(f)
  ## only peptide 12-19 (two charge states) overlaps residues 15-19 alone
  rr <- residue_regions(a, "A", "B", times = "100.00s", min_peptides = 2, pvalue = 1, alpha = 0.99)
  covered_by_one <- 15:19
  hit <- any(sapply(seq_len(nrow(rr)), function(i) any(covered_by_one %in% rr$region_start[i]:rr$region_end[i] &
                                                       !covered_by_one %in% 5:14)))
  expect_false(hit)
})

test_that("compare_states_tp(): raw p by default, Welch ANOVA omnibus, deprecated p_adjust_method", {
  raw <- local({
    r <- read.csv(file_nm)
    s <- r[r$Protein.State == "bound", ]; s$Protein.State <- "third"
    set.seed(3); s$X..Deut <- as.numeric(s$X..Deut) + rnorm(nrow(s), sd = 0.2)
    s$Experiment <- gsub("^bound", "third", s$Experiment)
    f <- tempfile(fileext = ".csv"); write.csv(rbind(r, s), f, row.names = FALSE); f
  })
  res <- suppressMessages(compare_states_tp(raw, times = "60.00s"))
  expect_equal(res$p_adj, res$p_value)
  kw <- suppressMessages(compare_states_tp(raw, times = "60.00s", omnibus = "kruskal"))
  expect_false(isTRUE(all.equal(res$p_omnibus, kw$p_omnibus)))
  expect_warning(old <- suppressMessages(compare_states_tp(raw, times = "60.00s", p_adjust_method = "BH")),
                 "deprecated")
  expect_equal(old$p_adj, p.adjust(old$p_value, "BH"))
})
