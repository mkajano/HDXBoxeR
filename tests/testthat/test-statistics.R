## Core statistics: sd/ave/pv_timepoint(), dif_ave(), CI functions, nb_exch_deut()

tp <- function(...) suppressMessages(output_tp(...))

test_that("ave_timepoint(), sd_timepoint(), dif_ave() match direct computation", {
  a <- tp(make_hdx_csv())
  av <- ave_timepoint(a, 3)
  s <- sd_timepoint(a, 3)
  expect_equal(colnames(av)[7:8], c("av_A_X..Deut", "av_B_X..Deut"))
  expect_equal(colnames(s)[7:8], c("sd_A_X..Deut", "sd_B_X..Deut"))
  expect_equal(av$av_A_X..Deut, rowMeans(a[, 7:9]))
  expect_equal(av$av_B_X..Deut, rowMeans(a[, 10:12]))
  expect_equal(s$sd_A_X..Deut, apply(a[, 7:9], 1, sd))
  expect_equal(s$sd_B_X..Deut, rep(0.01, nrow(a)), tolerance = 1e-9)
  d <- dif_ave(av)
  expect_equal(d[, 8], av[, 8] - av[, 7])
  expect_equal(d[, 7], rep(0, nrow(a)))
})

test_that("pv_timepoint() is a Welch t-test of every state against the first", {
  a <- tp(make_hdx_csv(states = c("A", "B", "C")))
  pv <- pv_timepoint(a, 3)
  for (i in seq_len(nrow(a))) {
    x1 <- as.numeric(a[i, 7:9])
    expect_equal(pv[i, 8], t.test(x1, as.numeric(a[i, 10:12]))$p.value)
    expect_equal(pv[i, 9], t.test(x1, as.numeric(a[i, 13:15]))$p.value)
  }
})

test_that("CI_2pts() follows the Hageman & Weis pooled-SD formula with qt(alpha, 2n-2)", {
  s1 <- c(0.10, 0.12, 0.13, 0.09, 0.11, 0.10)
  s2 <- c(0.18, 0.11, 0.13, 0.08, 0.11, 0.06)
  for (n in c(2, 3, 4)) for (alpha in c(0.05, 0.01)) {
    ## pooled SD per state = sqrt(mean(s^2)) for equal n; SE of difference =
    ## sqrt((sp1^2 + sp2^2) / n); threshold = SE * |t(alpha, 2n - 2)|
    expected <- sqrt((mean(s1^2) + mean(s2^2)) / n) * abs(qt(alpha, 2 * n - 2))
    expect_equal(CI_2pts(s1, s2, replicates = n, alpha = alpha), expected)
  }
})

test_that("CI_single(), CI_tp() and CI_tc() agree with CI_2pts()", {
  s1 <- c(0.10, 0.12, 0.13, 0.09)
  expect_equal(CI_single(s1, 3, 0.01), sqrt(mean(s1^2) / 3) * abs(qt(0.01, 4)))

  a <- tp(make_hdx_csv(states = c("A", "B", "C")))
  s <- sd_timepoint(a, 3)
  ci <- CI_tp(s, 3, 0.01)
  expect_length(ci, 2)
  expect_equal(ci[1], CI_2pts(s[, 7], s[, 8], 3, 0.01))
  expect_equal(ci[2], CI_2pts(s[, 7], s[, 9], 3, 0.01))

  sc <- data.frame(matrix(0, 4, 6), x = c(0.1, 0.2, 0.1, 0.3))
  sv <- data.frame(matrix(0, 4, 6), x = c(0.2, 0.1, 0.2, 0.1))
  expect_equal(CI_tc(sc, sv, 3, 0.01), CI_2pts(sc$x, sv$x, 3, 0.01))
})

test_that("sd_ci_timepoint() recovers SD of the difference from the two-sided t interval", {
  a <- tp(make_hdx_csv())
  s <- sd_ci_timepoint(a, replicates = 3, alpha = 0.01)
  x1 <- as.numeric(a[1, 7:9]); x2 <- as.numeric(a[1, 10:12])
  tt <- t.test(x1, x2, conf.level = 0.99)
  expect_equal(s[1, 8], diff(tt$conf.int) / (2 * abs(qt(0.005, 4))) * sqrt(3))
})

test_that("nb_exch_deut(): length - 2 - non-N-terminal prolines", {
  d <- data.frame(Start = c(1, 1, 1, 1, 1), End = c(10, 10, 10, 2, 3),
                  Sequence = c("AAAAAAAAAA",   # no Pro: 10 - 2 = 8
                               "AAAPAAAAAA",   # internal Pro: 7
                               "PAAAAAAAAA",   # N-terminal Pro is not subtracted again: 8
                               "AA",           # 2 residues: 0
                               "APA"),         # Pro at position 2: 3 - 2 - 1 = 0
                  stringsAsFactors = FALSE)
  expect_equal(nb_exch_deut(d), c(8, 7, 8, 0, 0))
  expect_equal(nb_exch_deut(default_peptides), c(8, 7, 6, 9))
})
