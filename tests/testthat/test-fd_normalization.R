## percent = "FD" and peptide_qc()

tp <- function(...) suppressMessages(output_tp(...))

test_that("percent = 'FD' is 100 x uptake / mean FD uptake of the same peptide and state", {
  f <- make_hdx_csv()   # FD uptake = 0.7 * nb_exch_deut
  a <- tp(f, percent = "FD")
  expect_equal(attr(a, "percent"), "FD")
  expect_true(all(grepl("_Deut\\.\\._[0-9]+$", colnames(a)[7:ncol(a)])))
  nb <- nb_exch_deut(default_peptides)
  row <- a$Deut.Time == "100.00s" & a$Start == 20
  expect_equal(a[row, "B_Deut.._2"], 100 * expected_uptake(2, 2, 4, 2) / (0.7 * nb[4]))

  tc <- suppressMessages(output_tc(f, percent = "FD"))
  r2 <- tc$Protein.State == "A" & tc$Start == 1
  expect_equal(tc[r2, "t10.00s_Deut.._1"], 100 * expected_uptake(1, 1, 1, 1) / (0.7 * nb[1]))
})

test_that("percent = 'FD' leaves peptides without a usable FD value empty instead of dropping them", {
  f <- make_hdx_csv()
  x <- read.csv(f, check.names = FALSE)
  x <- x[!(x[["Deut Time"]] == "FD" & x[["Protein State"]] == "B" & x$Start == 12), ]   # no FD
  x[["# Deut"]][x[["Deut Time"]] == "FD" & x$Start == 20] <- 0                             # FD <= 0
  utils::write.csv(x, f, row.names = FALSE)

  expect_message(a <- output_tp(f, percent = "FD"), "no usable fully deuterated")
  expect_equal(nrow(a), 2 * nrow(default_peptides))
  r12 <- a$Start == 12
  expect_true(all(is.na(a[r12, paste0("B_Deut.._", 1:3)])))
  expect_false(anyNA(a[r12, paste0("A_Deut.._", 1:3)]))
  expect_true(all(is.na(a[a$Start == 20, 7:12])))
})

test_that("percent = TRUE / FALSE are unchanged by the FD option", {
  f <- make_hdx_csv()
  expect_equal(attr(tp(f, percent = TRUE), "percent"), TRUE)
  expect_equal(attr(tp(f), "percent"), FALSE)
})

test_that("peptide_qc() flags missing references, all-zero uptake and bad FD values", {
  f <- make_hdx_csv()
  x <- read.csv(f, check.names = FALSE)
  ## peptide 12-19 in state A: no 0s control and uptake reported as 0 everywhere
  hit <- x[["Protein State"]] == "A" & x$Start == 12
  x <- x[!(hit & x[["Deut Time"]] == "0s"), ]
  hit <- x[["Protein State"]] == "A" & x$Start == 12
  x[["# Deut"]][hit] <- 0
  ## peptide 20-30 in B: FD above the theoretical maximum
  x[["# Deut"]][x[["Deut Time"]] == "FD" & x[["Protein State"]] == "B" & x$Start == 20] <- 50
  utils::write.csv(x, f, row.names = FALSE)

  qc <- peptide_qc(f)
  expect_equal(nrow(qc), 2 * nrow(default_peptides))
  q1 <- qc[qc$Protein.State == "A" & qc$Start == 12, ]
  expect_false(q1$has_undeut); expect_true(q1$all_zero_uptake); expect_true(q1$fd_nonpositive)
  expect_match(q1$issues, "no undeuterated control")
  q2 <- qc[qc$Protein.State == "B" & qc$Start == 20, ]
  expect_true(q2$fd_above_max)
  ok <- qc[qc$Protein.State == "B" & qc$Start == 1, ]
  expect_equal(ok$issues, "")
  expect_equal(ok$max_uptake, 9)          # MKTAYIAKQR: 10 - 1 - 0 prolines
  expect_equal(ok$fd_fraction, round(0.7 * 8 / 9, 3))
})

test_that("theoretical maximum follows HDXExaminer: length - 1 - non-N-terminal prolines", {
  expect_equal(HDXBoxeR:::.hdx_theoretical_max(c("GPLGSKAVVPGPAEHPLQY", "PAAA", "AAPA", "AA+Ox")),
               c(14, 3, 2, 1))
})
