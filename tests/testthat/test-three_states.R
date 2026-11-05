## checks on the 3-state example (uneven replicates, mutant peptides, "0s" control)

f3 <- system.file("extdata", "three_states.csv", package = "HDXBoxeR")

test_that("robot_indexes_df() returns peptide keys and robot_indexes() row numbers that point to them", {
  th <- suppressMessages(output_tc(f3)); thP <- suppressMessages(output_tc(f3, percent = TRUE))
  df1 <- robot_indexes_df(thP = thP, th = th)
  expect_true(all(c("Start", "End", "Sequence", "Charge") %in% colnames(df1)))
  expect_false(anyDuplicated(df1[, c("Start", "End", "Sequence", "Charge")]) > 0)
  idx <- robot_indexes(thP = thP, th = th)
  ## the rows point to significant peptides of the WT vs D109H table; peptides
  ## that cover the D109H site have another sequence there and are left out
  st <- unique(thP$Protein.State)
  tab <- prep_timecourse_plot_ave(thP[thP$Protein.State == st[1], ], thP[thP$Protein.State == st[2], ], 3)[[1]]
  key <- function(d) paste(d$Start, d$End, d$Sequence, d$Charge)
  expect_true(all(key(tab[idx, ]) %in% key(df1)))
  expect_lte(length(idx), nrow(df1))
  missing <- df1[!key(df1) %in% key(tab), ]
  expect_true(all(missing$Start <= 109 & missing$End >= 109))
})

test_that("hdx_summary_table() handles uneven replicates and names the reference state", {
  s <- suppressMessages(hdx_summary_table(f3))
  expect_equal(ncol(s), 4)
  rep_row <- s[s$Item == "Number of replicates", ]
  expect_equal(unname(unlist(rep_row[1, -1])), c("4", "3-4", "3-4"))
  sig_row <- s[s$Item == "Significant differences in HDX", ]
  expect_equal(sig_row[1, 2], "reference state")
  expect_match(sig_row[1, 3], "critical interval")
  expect_match(s[s$Item == "Sequence coverage (%)", 2], "residues 1-175")
})

test_that("output_FD() returns every state; output_UD() finds a '0s' control", {
  fd <- output_FD(f3)
  expect_setequal(unique(fd$Protein.State), c("WT_pH7.5", "D109H_pH7.5", "R120G_ph7.5"))
  expect_equal(colnames(fd)[1:6], c("Protein.State", "Start", "End", "Sequence", "Search.RT", "Charge"))
  ## the "0s" rows are found; HDXExaminer leaves their uptake empty, so it is 0
  ud3 <- output_UD(f3)
  expect_setequal(unique(ud3$Protein.State), unique(fd$Protein.State))
  expect_true(all(as.matrix(ud3[, 7:ncol(ud3)]) == 0, na.rm = TRUE))
  expect_no_error(suppressMessages(average_timecourse(f3)))
  expect_no_error(suppressMessages(sd_timecourse_proc(f3)))
  ## the 2-state example has an empty "0s" and a "0.00s" with values
  ud <- output_UD(system.file("extdata", "All_results_table.csv", package = "HDXBoxeR"))
  expect_true(nrow(ud) > 0)
  expect_setequal(unique(ud$Protein.State), c("Unbound", "bound"))
})

test_that("_proc functions stop when the two tables are not aligned", {
  a <- suppressMessages(output_tp(f3)); ap <- suppressMessages(output_tp(f3, percent = TRUE))
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp); on.exit(grDevices::dev.off())
  expect_error(plot_heat_map_tp_proc(input_proc = ap[-1, ], input_up = a))
})

test_that("peptides are matched without retention time", {
  r <- read.csv(f3, check.names = FALSE)
  base <- suppressMessages(output_tp(f3))
  ## shift the RT of one state; the wide table must not change in size
  r$Search.RT[r$Protein.State == "R120G_ph7.5"] <- as.numeric(r$Search.RT[r$Protein.State == "R120G_ph7.5"]) + 0.3
  f <- tempfile(fileext = ".csv"); utils::write.csv(r, f, row.names = FALSE)
  shifted <- suppressMessages(output_tp(f))
  expect_equal(nrow(shifted), nrow(base))
})
