## End-to-end: the comparison functions accept states with different
## replicate counts and peptides with a missing replicate.

uneven_file <- function() {
  make_hdx_csv(states = c("WT", "WT_mut", "B"), replicates = c(WT = 3, WT_mut = 4, B = 3),
               drop = data.frame(state = "B", time = "100.00s", replicate = 1, Start = 5))
}

test_that("plot_fissure_vent() and the interactive table work with uneven replicates", {
  f <- uneven_file()
  a <- suppressMessages(output_tp(f, states = c("WT", "WT_mut")))
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  expect_no_error(suppressMessages(plot_fissure_vent(a, times = "100.00s", state1 = "WT", state2 = "WT_mut")))
  grDevices::dev.off()

  tbl <- HDXBoxeR:::.hdx_comparison_table(a, "WT", "WT_mut", "100.00s", 3, 0.01, 0.01)
  x1 <- as.numeric(a[a$Deut.Time == "100.00s" & a$Start == 1, paste0("WT_X..Deut_", 1:3)])
  x2 <- as.numeric(a[a$Deut.Time == "100.00s" & a$Start == 1, paste0("WT_mut_X..Deut_", 1:4)])
  expect_equal(tbl$p_value[tbl$Start == 1], t.test(x1, x2, conf.level = 0.99)$p.value)
})

test_that("state lookup is exact first: 'WT' does not pick up 'WT_mut'", {
  cols <- c("Deut.Time", "Start", "End", "Sequence", "Search.RT", "Charge",
            paste0("WT_X..Deut_", 1:3), paste0("WT_mut_X..Deut_", 1:4))
  expect_equal(HDXBoxeR:::.hdx_state_cols("WT", cols), paste0("WT_X..Deut_", 1:3))
  expect_equal(HDXBoxeR:::.hdx_state_cols("WT_m", cols), paste0("WT_mut_X..Deut_", 1:4))
  expect_equal(HDXBoxeR:::.hdx_summary_col("WT", c(cols[1:6], "av_WT_X..Deut", "av_WT_mut_X..Deut"), "av_"),
               "av_WT_X..Deut")
})

test_that("compare_states_tp() handles uneven replicates and a missing replicate", {
  f <- uneven_file()
  res <- suppressMessages(compare_states_tp(f, times = "100.00s"))
  expect_equal(length(unique(paste(res$state1, res$state2))), 3)
  expect_false(anyNA(res$p_omnibus))
  row <- res[res$state1 == "WT_mut" & res$state2 == "B" & res$Start == 5, ]
  expect_equal(nrow(row), 1)
  expect_false(is.na(row$p_value))
})

test_that("woods_CI_plot() / robot plots run on uneven replicate time courses", {
  f <- make_hdx_csv(replicates = c(A = 3, B = 4))
  th <- suppressMessages(output_tc(f)); thP <- suppressMessages(output_tc(f, percent = TRUE))
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  expect_no_error(suppressMessages(robot_indexes(thP = thP, th = th, replicates = 3)))
  grDevices::dev.off()
})
