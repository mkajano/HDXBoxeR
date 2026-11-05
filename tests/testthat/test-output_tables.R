## Core parsers: output_tp(), output_tc(), output_prep(), select_indices()

tp <- function(...) suppressMessages(output_tp(...))
tc <- function(...) suppressMessages(output_tc(...))

test_that("output_tp() lays out one row per peptide x time and one column per state x replicate", {
  f <- make_hdx_csv()
  a <- tp(f)
  expect_equal(nrow(a), 2 * nrow(default_peptides))
  expect_equal(colnames(a)[1:6], c("Deut.Time", "Start", "End", "Sequence", "Search.RT", "Charge"))
  expect_equal(colnames(a)[7:12], c(paste0("A_X..Deut_", 1:3), paste0("B_X..Deut_", 1:3)))
  expect_equal(unique(a$Deut.Time), c("10.00s", "100.00s"))

  ## every value lands in the right cell
  for (t_idx in 1:2) for (p in seq_len(nrow(default_peptides))) for (r in 1:3) {
    row <- a$Deut.Time == c("10.00s", "100.00s")[t_idx] & a$Start == default_peptides$Start[p]
    expect_equal(a[row, paste0("A_X..Deut_", r)], expected_uptake(1, t_idx, p, r))
    expect_equal(a[row, paste0("B_X..Deut_", r)], expected_uptake(2, t_idx, p, r))
  }
})

test_that("output_tp() percent = TRUE returns Deut % columns", {
  f <- make_hdx_csv()
  a <- tp(f, percent = TRUE)
  expect_true(all(grepl("_Deut\\.\\._[0-9]+$", colnames(a)[7:ncol(a)])))
  nb <- nb_exch_deut(default_peptides)
  row1 <- a$Deut.Time == "10.00s" & a$Start == 1
  expect_equal(a[row1, "A_Deut.._1"], 100 * expected_uptake(1, 1, 1, 1) / nb[1])
})

test_that("output_tp() honours replicates, states and times", {
  f <- make_hdx_csv(states = c("A", "B", "C"), times = c("10.00s", "100.00s", "1000.00s"))
  a <- tp(f, replicates = 2, states = c("C", "A"), times = "100.00s")
  expect_equal(colnames(a)[7:10], c(paste0("C_X..Deut_", 1:2), paste0("A_X..Deut_", 1:2)))
  expect_equal(unique(a$Deut.Time), "100.00s")
})

test_that("output_tp() uses every replicate; states may differ; replicates= caps", {
  f <- make_hdx_csv(replicates = c(A = 3, B = 4))
  a <- tp(f)
  expect_equal(colnames(a)[7:13], c(paste0("A_X..Deut_", 1:3), paste0("B_X..Deut_", 1:4)))
  expect_equal(attr(a, "replicates"), c(A = 3L, B = 4L))
  row <- a$Deut.Time == "100.00s" & a$Start == 20
  expect_equal(a[row, "B_X..Deut_4"], expected_uptake(2, 2, 4, 4))
  capped <- tp(f, replicates = 3)
  expect_equal(attr(capped, "replicates"), c(A = 3L, B = 3L))
})

test_that("output_tc() keeps every replicate per time point", {
  f <- make_hdx_csv(replicates = c(A = 3, B = 4))
  b <- suppressMessages(output_tc(f))
  expect_equal(attr(b, "replicates"), c(`10.00s` = 4L, `100.00s` = 4L))
  expect_true(all(is.na(b[b$Protein.State == "A", "t10.00s_X..Deut_4"])))
  expect_false(anyNA(b[b$Protein.State == "B", "t10.00s_X..Deut_4"]))
})

test_that("output_tp() seq_match: FALSE requires identical sequence, TRUE matches by position", {
  f <- make_hdx_csv()
  x <- read.csv(f, check.names = FALSE)
  hit <- x[["Protein State"]] == "B" & x$Start == 1
  x$Sequence[hit] <- "MKTAYIAKQA"
  utils::write.csv(x, f, row.names = FALSE)

  exact <- tp(f, seq_match = FALSE)
  by_pos <- tp(f, seq_match = TRUE)
  expect_false(any(exact$Start == 1))
  expect_true(all(c(1, 5, 12, 20) %in% by_pos$Start))
  expect_equal(unique(by_pos$Sequence[by_pos$Start == 1]), "MKTAYIAKQR.MKTAYIAKQA")
})

test_that("output_tp() keeps a peptide with a missing replicate as NA; min_replicates sets the floor", {
  f <- make_hdx_csv(drop = data.frame(state = "B", time = "10.00s", replicate = 2, Start = 12))
  a <- tp(f)
  row <- a$Deut.Time == "10.00s" & a$Start == 12
  expect_equal(sum(row), 1)
  expect_true(is.na(a[row, "B_X..Deut_2"]))
  expect_equal(a[row, "B_X..Deut_3"], expected_uptake(2, 1, 3, 3))

  strict <- tp(f, min_replicates = 3)
  expect_false(any(strict$Deut.Time == "10.00s" & strict$Start == 12))
  expect_true(any(strict$Deut.Time == "100.00s" & strict$Start == 12))

  ## two of three missing -> below the default floor of 2 -> dropped
  f2 <- make_hdx_csv(drop = data.frame(state = "B", time = "10.00s", replicate = 1:2, Start = 12))
  expect_message(a2 <- output_tp(f2), "1 peptide row\\(s\\) dropped")
  expect_false(any(a2$Deut.Time == "10.00s" & a2$Start == 12))
})

test_that("output_tc() lays out one row per peptide x state and one column per time x replicate", {
  f <- make_hdx_csv()
  b <- tc(f)
  expect_equal(nrow(b), 2 * nrow(default_peptides))
  expect_equal(colnames(b)[1:6], c("Sequence", "Protein.State", "Start", "End", "Search.RT", "Charge"))
  expect_equal(colnames(b)[7:12], c(paste0("t10.00s_X..Deut_", 1:3), paste0("t100.00s_X..Deut_", 1:3)))
  row <- b$Protein.State == "B" & b$Start == 20
  expect_equal(b[row, "t100.00s_X..Deut_2"], expected_uptake(2, 2, 4, 2))
})

test_that("output_prep() writes the Masson-style table with mean/SD and controls", {
  f <- make_hdx_csv()
  out_file <- tempfile(fileext = ".csv")
  suppressMessages(output_prep(f, output_name = out_file))
  o <- read.csv(out_file, check.names = FALSE)
  up <- "Mean.Deut.Uptake_.Da."
  expect_true(all(c("Protein.State", "Deut.Time", "Start", "End", "Sequence", "Charge",
                    up, "st.dev_Deut_Uptake") %in% colnames(o)))
  r <- o$Protein.State == "A" & o$Deut.Time == "10.00s" & o$Start == 1
  expect_equal(o[r, up], round(mean(sapply(1:3, function(k) expected_uptake(1, 1, 1, k))), 4))
  expect_equal(o[r, "st.dev_Deut_Uptake"], 0.01, tolerance = 1e-6)
  ## FD rows carry the measured FD uptake and SD 0; 0s rows are 0
  fd <- o$Protein.State == "A" & o$Deut.Time == "FD" & o$Start == 1
  expect_equal(o[fd, up], 0.7 * nb_exch_deut(default_peptides)[1])
  expect_equal(o[fd, "st.dev_Deut_Uptake"], 0)
  ud <- o$Protein.State == "A" & o$Deut.Time == "0s" & o$Start == 1
  expect_equal(o[ud, up], 0)
})

test_that("select_indices() filters by start, end, length and time", {
  a <- tp(make_hdx_csv())
  expect_equal(a$Start[select_indices(a, start = 10)], c(12, 20, 12, 20))
  expect_equal(a$End[select_indices(a, end = 14)], c(10, 14, 10, 14))
  expect_true(all((a$End - a$Start)[select_indices(a, length = 8)] <= 8))
  expect_true(all(a$Deut.Time[select_indices(a, times = "100.00s")] == "100.00s"))
  expect_equal(length(select_indices(a, start = 5, end = 19, times = "10.00s")), 2)
})
