## Input layer: read_hdx(), read_hdxexaminer(), read_dynamx(), as_hdx_long()

file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
p <- 1.00727646688

## DynamX 3.x cluster export with known uptake: neutral mass of a peptide is
## m0 + uptake; Center (m/z) = mass / z + proton mass. Uptake(state, minutes,
## file) = base * f(t) + shift(state) + 0.01 * replicate.
make_dynamx_cluster <- function(version = 3, fd_state = TRUE, two_charges = FALSE) {
  peps <- data.frame(Start = c(1, 5), End = c(10, 14), Sequence = c("MKTAYIAKQR", "YIAKQRPISF"),
                     MaxUptake = c(9, 8), MHP = c(1200.5, 1300.6), m0 = c(1199.5, 1299.6),
                     base = c(3, 4), stringsAsFactors = FALSE)
  states <- c("apo", "holo"); shift <- c(apo = 0, holo = 0.5)
  rows <- list()
  add <- function(state, expo, file, k, upt, z = 2, inten = 1e5) {
    mass <- peps$m0[k] + upt
    rows[[length(rows) + 1]] <<- data.frame(
      Protein = "Prot", Start = peps$Start[k], End = peps$End[k], Sequence = peps$Sequence[k],
      Modification = "", Fragment = "", MaxUptake = peps$MaxUptake[k], MHP = peps$MHP[k],
      State = state, Exposure = expo, File = file, z = z, RT = 5 + k + 0.01 * nchar(file),
      Inten = inten, Center = mass / z + p, stringsAsFactors = FALSE)
  }
  for (k in 1:2) {
    for (s in states) {
      add(s, 0, "undeut_1", k, 0); add(s, 0, "undeut_2", k, 0)
      for (expo in c(0.5, 5)) for (r in 1:3) {
        u <- peps$base[k] * c(`0.5` = 0.3, `5` = 0.9)[[as.character(expo)]] + shift[[s]] + 0.01 * r
        add(s, expo, paste0(s, "_", expo, "_", r), k, u)
        if (two_charges) add(s, expo, paste0(s, "_", expo, "_", r), k, u, z = 3, inten = 3e5)
      }
      if (two_charges) { add(s, 0, "undeut_1", k, 0, z = 3); add(s, 0, "undeut_2", k, 0, z = 3) }
    }
    if (fd_state) for (r in 1:2) {
      add("FD control", 0, "undeut_1", k, 0)
      add("FD control", 0.167, paste0("fd_", r), k, 0.8 * peps$MaxUptake[k] + 0.01 * r)
    }
  }
  d <- do.call(rbind, rows)
  if (version == 2) {
    d$Fragment <- NULL
    colnames(d)[colnames(d) == "MaxUptake"] <- "Max Exchangers"
  }
  f <- tempfile(fileext = ".csv")
  utils::write.csv(d, f, row.names = FALSE, quote = FALSE)
  f
}

test_that("read_hdxexaminer()/read_hdx() give the same results as a file path", {
  d <- read_hdx(file_nm)
  expect_s3_class(d, "hdx_data")
  a1 <- suppressMessages(output_tp(file_nm)); a2 <- suppressMessages(output_tp(d))
  expect_equal(a2, a1)
  expect_equal(nm_states(d), c("Unbound", "bound"))
  expect_equal(backHX_calculations(d), backHX_calculations(file_nm))
})

test_that("read_dynamx() cluster data: files are replicates, uptake from centroid mass", {
  f <- make_dynamx_cluster()
  d <- read_dynamx(f, fd = "FD control")
  expect_s3_class(d, "hdx_data")
  expect_setequal(unique(d$Protein.State), c("apo", "holo"))
  expect_setequal(unique(d$Deut.Time), c("0s", "30.00s", "300.00s", "FD"))

  a <- suppressMessages(output_tp(d))
  expect_equal(attr(a, "replicates"), c(apo = 3L, holo = 3L))
  row <- a$Deut.Time == "300.00s" & a$Start == 5
  expect_equal(a[row, "holo_X..Deut_2"], 4 * 0.9 + 0.5 + 0.02, tolerance = 1e-6)
  expect_equal(a[row, "apo_X..Deut_1"], 4 * 0.9 + 0.01, tolerance = 1e-6)

  pct <- suppressMessages(output_tp(d, percent = TRUE))
  expect_equal(pct[row, "apo_Deut.._1"], 100 * (4 * 0.9 + 0.01) / 8, tolerance = 1e-6)

  ## FD control copied to every state; back-exchange from centroids
  ## FD uptake = 0.8 * MaxUptake + 0.01 * r (mean r = 1.5); N = nb_exch_deut = 8 and 7
  bx <- backHX_calculations(d, Dfact = 1)
  expected <- mean(c(1 - (0.8 * 9 + 0.015) / 8, 1 - (0.8 * 8 + 0.015) / 7)) * 100
  expect_equal(as.numeric(bx[, "<BackHX>"]), rep(round(expected, 2), 2), tolerance = 1e-6)
})

test_that("read_dynamx() reads DynamX 2.x headers and passes a path straight to output_tp()", {
  f2 <- make_dynamx_cluster(version = 2, fd_state = FALSE)
  a <- suppressMessages(output_tp(f2))
  expect_equal(nrow(a), 4)
  expect_true(all(grepl("^(apo|holo)_X\\.\\.Deut_[1-3]$", colnames(a)[7:12])))
})

test_that("read_dynamx(): charge states separate by default, combined on request", {
  f <- make_dynamx_cluster(fd_state = FALSE, two_charges = TRUE)
  sep <- suppressMessages(output_tp(read_dynamx(f)))
  expect_setequal(unique(sep$Charge), c(2, 3))
  comb <- suppressMessages(output_tp(read_dynamx(f, combine_charges = TRUE)))
  expect_equal(nrow(comb), nrow(sep) / 2)
  row <- comb$Deut.Time == "30.00s" & comb$Start == 1
  expect_equal(comb[row, "apo_X..Deut_3"], 3 * 0.3 + 0.03, tolerance = 1e-6)
})

test_that("read_dynamx() fd = 'last' turns each state's longest exposure into the FD control", {
  d <- read_dynamx(make_dynamx_cluster(fd_state = FALSE), fd = "last")
  expect_setequal(unique(d$Deut.Time), c("0s", "30.00s", "FD"))
})

test_that("read_dynamx() state data warns that it has no replicates", {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Protein = "P", Start = 1, End = 10, Sequence = "MKTAYIAKQR",
                              Modification = "", Fragment = "", MaxUptake = 9, MHP = 1200.5,
                              State = c("apo", "apo"), Exposure = c(0, 1), Center = c(1200.5, 1203.5),
                              `Center SD` = 0.1, Uptake = c(0, 3), `Uptake SD` = 0.1, RT = 5, `RT SD` = 0.01,
                              check.names = FALSE), f, row.names = FALSE)
  expect_warning(d <- read_hdx(f), "no replicates")
  expect_equal(d$X..Deut[d$Deut.Time == "60.00s"], 3)
})

test_that("as_hdx_long() builds input from any long table", {
  long <- data.frame(state = rep(c("apo", "holo"), each = 6),
                     t_min = rep(rep(c(1, 10), each = 3), 2),
                     run = rep(paste0("r", 1:3), 4),
                     from = 5, to = 14, seq = "YIAKQRPISF",
                     D = c(1.1, 1.2, 1.15, 2.4, 2.5, 2.45, 0.9, 1.0, 0.95, 2.0, 2.1, 2.05))
  d <- as_hdx_long(long, state = "state", time = "t_min", replicate = "run", start = "from",
                   end = "to", sequence = "seq", uptake = "D", time_unit = "min")
  a <- suppressMessages(output_tp(d))
  expect_equal(unique(a$Deut.Time), c("60.00s", "600.00s"))
  expect_equal(a$holo_X..Deut_2[a$Deut.Time == "600.00s"], 2.1)
  expect_equal(a$apo_X..Deut_1[1], 1.1)
  expect_error(as_hdx_long(long, state = "nope", time = "t_min", replicate = "run", start = "from",
                           end = "to", sequence = "seq", uptake = "D"), "not found")
})

test_that("an unrecognised file is reported clearly", {
  f <- tempfile(fileext = ".csv"); writeLines(c("a,b", "1,2"), f)
  expect_error(output_tp(f), "unrecognised input format")
})
