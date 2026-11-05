## compare_regions(): region-level ANOVA with peptides as blocks

file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

test_that("with two states the region test is a paired t-test of per-peptide differences", {
  res <- compare_regions(file_nm, times = "1800.00s", regions = list(nterm = c(1, 100)))
  a <- suppressMessages(output_tp(file_nm, times = "1800.00s"))
  a <- a[a$Start >= 1 & a$End <= 100, ]
  u <- rowMeans(a[, paste0("Unbound_X..Deut_", 1:3)]); b <- rowMeans(a[, paste0("bound_X..Deut_", 1:3)])
  tt <- t.test(b, u, paired = TRUE)
  expect_equal(res$summary$p_state, tt$p.value, tolerance = 1e-8)
  expect_equal(res$pairs$estimate, mean(b - u), tolerance = 1e-8)
  expect_equal(res$summary$n_peptides, nrow(a))
  expect_equal(res$summary$df_resid, nrow(a) - 1)
})

test_that("compare_regions() returns one row per region and time, and handles p_adjust", {
  reg <- list(nterm = c(1, 100), cterm = c(101, 225))
  res <- compare_regions(file_nm, regions = reg)
  expect_equal(nrow(res$summary), 2 * 4)
  expect_true(all(c("F_state", "p_state", "p_state_adj", "p_state_x_peptide") %in% colnames(res$summary)))
  expect_equal(res$summary$p_state_adj, res$summary$p_state)
  bh <- compare_regions(file_nm, regions = reg, p_adjust = "BH")
  expect_equal(bh$summary$p_state_adj, p.adjust(bh$summary$p_state, "BH"))
})

test_that("time_model = TRUE fits one model per region with a state x time test", {
  res <- compare_regions(file_nm, regions = list(cterm = c(101, 225)), time_model = TRUE)
  expect_equal(res$summary$Deut.Time, "all")
  expect_true("p_state_x_time" %in% colnames(res$summary))
  expect_false(is.na(res$summary$p_state_x_time))
  expect_error(compare_regions(file_nm, times = "60.00s", regions = list(a = c(1, 225)), time_model = TRUE),
               "two or more time points")
})

test_that("compare_regions() works for 3+ states with pairwise contrasts and a reference", {
  f <- make_hdx_csv(states = c("A", "B", "C"), replicates = c(A = 3, B = 4, C = 3))
  ## add small noise so the additive model does not fit exactly
  x <- read.csv(f, check.names = FALSE)
  set.seed(7)
  lab <- !x[["Deut Time"]] %in% c("0s", "FD")
  x[["# Deut"]][lab] <- as.numeric(x[["# Deut"]][lab]) + rnorm(sum(lab), sd = 0.03)
  utils::write.csv(x, f, row.names = FALSE)
  res <- compare_regions(f, times = "100.00s", regions = list(all = c(1, 30)))
  expect_equal(res$summary$df_state, 2)
  expect_equal(nrow(res$pairs), 3)
  ## synthetic data: state shift 0.5 Da per state step
  expect_equal(res$pairs$estimate[res$pairs$state1 == "A" & res$pairs$state2 == "C"], 1.0, tolerance = 0.05)
  ref <- compare_regions(f, times = "100.00s", regions = list(all = c(1, 30)), reference = "B")
  expect_equal(nrow(ref$pairs), 2)
  expect_true(all(ref$pairs$state1 == "B"))
})

test_that("compare_regions() errors clearly when no region can be tested", {
  expect_error(suppressMessages(compare_regions(file_nm, times = "1800.00s", regions = list(one = c(98, 115)))),
               "two or more peptides")
  expect_error(compare_regions(file_nm, times = "1800.00s"), "regions is required")
})
