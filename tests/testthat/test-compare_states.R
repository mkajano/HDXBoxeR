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

test_that("compare_states_tp() reads the input file exactly once regardless of pairs/regions", {
  ## regression test for a real performance bug: output_tp() re-parsed the
  ## whole file from disk (via arg_df()) once per state pair PLUS once more
  ## per region for each pair (regions didn't restrict what got computed,
  ## only filtered it afterward) - for N states x R regions that's up to
  ## R * choose(N, 2) + 1 full file parses of an identical file. Fixed by
  ## parsing once (arg_df(filepath)) and threading it through every
  ## output_tp() call via its .arg_df argument, and by computing each
  ## pair's result once and reusing it across regions instead of
  ## recomputing per region.
  ns <- asNamespace("HDXBoxeR")
  orig_arg_df <- ns$arg_df
  n_calls <- 0
  counting_arg_df <- function(filepath) { n_calls <<- n_calls + 1; orig_arg_df(filepath) }
  assignInNamespace("arg_df", counting_arg_df, ns = "HDXBoxeR")
  on.exit(assignInNamespace("arg_df", orig_arg_df, ns = "HDXBoxeR"))

  res <- compare_states_tp(raw3_path, times = "1800.00s",
                            regions = list(low = c(1, 100), high = c(101, 225)))
  expect_equal(n_calls, 1)
  expect_true(nrow(res) > 0)
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

test_that("compare_states_tp() correctly matches a state name containing a space", {
  ## regression test: output_tp()'s wide-format columns strip spaces from a
  ## state name (gsub(' ', '', state)) but then run that string through
  ## data.frame()'s own make.names() too, which additionally mangles any
  ## OTHER invalid character - matching state_of_col against the raw state
  ## string (with its space still in it) silently matched nothing and every
  ## delta/p_value/p_adj/is_sig for that state came back NA/FALSE with no
  ## warning at all. "Protein State 1" is HDX Examiner's own default,
  ## unedited state name - a very realistic case, not a contrived one.
  raw <- read.csv(file_nm)
  synth <- raw[raw$Protein.State == "bound", ]
  synth$Protein.State <- "Protein State 1"
  set.seed(42)
  synth$X..Deut <- as.numeric(synth$X..Deut) + rnorm(nrow(synth), sd = 0.3)
  synth$Experiment <- gsub("^bound", "syn", synth$Experiment)
  raw3_space <- rbind(raw, synth)
  tmp <- tempfile(fileext = ".csv")
  write.csv(raw3_space, tmp, row.names = FALSE)

  res <- compare_states_tp(tmp, times = "1800.00s", states = c("Unbound", "bound", "Protein State 1"))
  sub <- res[res$state1 == "Protein State 1" | res$state2 == "Protein State 1", ]
  expect_true(nrow(sub) > 0)
  expect_true(all(!is.na(sub$delta)))
  expect_true(all(!is.na(sub$p_value)))
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

test_that("compare_states_unmatched() defaults to raw # Deut (Da), not a per-peptide normalization", {
  ## regression test: this function used to always pool the package's own
  ## derived frac_deut (Da / exchangeable-proton count), silently applying
  ## a per-peptide size correction nobody asked for - raw # Deut is now the
  ## default (value = "raw"), with frac_deut/percent as explicit opt-ins.
  raw_file <- read.csv(file_nm)
  raw_val <- raw_file$X..Deut[raw_file$Protein.State == "Unbound" & raw_file$Deut.Time == "1800.00s" &
                                raw_file$Start == 40 & raw_file$End == 58]

  res_raw <- compare_states_unmatched(file_nm, times = "1800.00s", regions = list(nterm = c(1, 100)))
  expect_equal(res_raw$value, "raw")
  d <- res_raw$data[res_raw$data$Start == 40 & res_raw$data$End == 58 & res_raw$data$Protein.State == "Unbound", ]
  expect_equal(d$deut_value, mean(raw_val), tolerance = 1e-6)

  res_frac <- compare_states_unmatched(file_nm, times = "1800.00s", regions = list(nterm = c(1, 100)),
                                        value = "frac_deut")
  expect_equal(res_frac$value, "frac_deut")
  d2 <- res_frac$data[res_frac$data$Start == 40 & res_frac$data$End == 58 & res_frac$data$Protein.State == "Unbound", ]
  expect_false(isTRUE(all.equal(d2$deut_value, d$deut_value)))  # genuinely different scale

  ## the per-peptide diagnostic is scale-invariant: dividing both of one
  ## peptide's states by the same per-peptide constant cannot change a
  ## paired t-test's p-value, so this must be identical under raw/frac_deut
  expect_equal(res_raw$summary$n_peptides_significant, res_frac$summary$n_peptides_significant)
  expect_equal(res_raw$summary$frac_peptides_significant, res_frac$summary$frac_peptides_significant)
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

test_that("compare_states_unmatched() flags dilution: pooled non-sig but peptides individually differ", {
  ## regression check for the exact failure mode that motivated this
  ## diagnostic: nterm/cterm here pool a mix of peptides that individually
  ## differ (by the package's own Welch-t criterion) with ones that don't,
  ## so the pooled Wilcoxon test washes out to non-significant while a real
  ## fraction of the pooled peptides are individually significant.
  res <- compare_states_unmatched(file_nm, times = "1800.00s",
                                   regions = list(nterm = c(1, 100), cterm = c(101, 225)))
  expect_true(all(c("n_peptides_tested", "n_peptides_significant", "frac_peptides_significant") %in%
                    colnames(res$summary)))
  expect_true(all(res$summary$p_adj > 0.01))
  expect_true(all(res$summary$frac_peptides_significant > 0))
  expect_equal(res$summary$n_peptides_significant, c(5L, 8L))
  expect_equal(res$summary$n_peptides_tested, c(9L, 19L))
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

test_that("plot_peptide_coverage_map() runs for all states and a single state, errors on unknown state", {
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_no_error(plot_peptide_coverage_map(file_nm))
  expect_no_error(plot_peptide_coverage_map(file_nm, states = "Unbound"))
  grDevices::dev.off()

  expect_error(plot_peptide_coverage_map(file_nm, states = "not_a_state"))
})

test_that(".hdx_palette() distinguishes an explicit color from a brewer palette name", {
  ## regression test: a single explicit color like "steelblue" is also
  ## length-1, same as a brewer palette name like "Set2" - without checking
  ## membership in brewer.pal.info, "steelblue" was passed straight to
  ## brewer.pal() and errored ("steelblue is not a valid palette name").
  expect_equal(HDXBoxeR:::.hdx_palette(1, "steelblue"), "steelblue")
  expect_equal(HDXBoxeR:::.hdx_palette(3, "steelblue"), rep("steelblue", 3))
  expect_length(HDXBoxeR:::.hdx_palette(3, "Set2"), 3)
  expect_equal(HDXBoxeR:::.hdx_palette(2, c("red", "blue")), c("red", "blue"))
  expect_error(HDXBoxeR:::.hdx_palette(3, c("red", "blue")))
})

test_that("styling parameters on the ANOVA-workflow plots run without error", {
  res <- compare_states_unmatched(file_nm, times = "1800.00s",
                                   regions = list(nterm = c(1, 100), cterm = c(101, 225)))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_no_error(plot_coverage_diagnostics(res, palette = "Dark2", bg = "white", show_grid = FALSE, legend = FALSE))
  expect_no_error(plot_compare_states_unmatched(res, plot_type = "jitter", palette = "steelblue",
                                                  bg = "white", show_grid = FALSE))
  expect_no_error(plot_compare_states_unmatched(res, plot_type = "box", notch = TRUE))
  expect_no_error(plot_compare_states_unmatched(res, plot_type = "violin"))
  expect_no_error(plot_peptide_coverage_map(file_nm, states = "Unbound", palette = "steelblue",
                                              bg = "white", show_grid = FALSE, block_height = 0.9))
  grDevices::dev.off()
})

test_that(".hdx_tint() lightens toward white and always returns an opaque color", {
  ## regression test for the reported "I cannot see the points" issue:
  ## fills/points must never rely on alpha compositing (some graphics
  ## devices render semi-transparent fills as invisible or fully solid) -
  ## .hdx_tint() blends in RGB space instead, so its output has no alpha
  ## channel at all.
  white <- HDXBoxeR:::.hdx_tint("steelblue", 1)
  same  <- HDXBoxeR:::.hdx_tint("steelblue", 0)
  mid   <- HDXBoxeR:::.hdx_tint("steelblue", 0.5)
  expect_equal(toupper(white), "#FFFFFF")
  expect_equal(nchar(mid), 7)  # "#RRGGBB", no alpha channel appended
  expect_true(all(grDevices::col2rgb(mid) >= grDevices::col2rgb("steelblue")))
})

test_that("plot_compare_states_tp() draws every peptide, not only significant ones, for 3+ states", {
  ## build a synthetic third state, same trick as the compare_states_tp() tests above
  raw <- read.csv(file_nm)
  synth <- raw[raw$Protein.State == "bound", ]
  synth$Protein.State <- "synthetic"
  set.seed(42)
  synth$X..Deut <- as.numeric(synth$X..Deut) + rnorm(nrow(synth), sd = 0.3)
  synth$Experiment <- gsub("^bound", "synthetic", synth$Experiment)
  raw3 <- rbind(raw, synth)
  tmp_csv <- tempfile(fileext = ".csv")
  write.csv(raw3, tmp_csv, row.names = FALSE)

  res_tp <- compare_states_tp(tmp_csv, times = "1800.00s",
                               regions = list(nterm = c(1, 100), cterm = c(101, 225)))
  tmp_png <- tempfile(fileext = ".png")
  grDevices::png(tmp_png)
  expect_no_error(plot_compare_states_tp(res_tp, plot_type = "jitter"))
  expect_no_error(plot_compare_states_tp(res_tp, plot_type = "box"))
  expect_no_error(plot_compare_states_tp(res_tp, plot_type = "violin"))
  grDevices::dev.off()

  expect_error(plot_compare_states_tp(data.frame(x = 1)))
})
