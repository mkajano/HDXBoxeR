file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
protein_length <- 225

test_that("residue_deconvolution returns bounded, correctly-shaped output", {
  a <- output_tp(file_nm, percent = FALSE)
  res <- residue_deconvolution(a, protein_length = protein_length)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), protein_length * length(unique(a$Deut.Time)))
  expect_true(all(c("residue", "Deut.Time") %in% colnames(res)))

  value_cols <- grep("_X\\.\\.Deut_[0-9]+$", colnames(res), value = TRUE)
  expect_true(length(value_cols) > 0)
  vals <- unlist(res[, value_cols])
  ## tiny numerical tolerance for optim's floating-point precision at the
  ## box constraints, not a real bound violation
  expect_true(all(vals >= -1e-6 & vals <= 1 + 1e-6, na.rm = TRUE))
})

test_that("residue_deconvolution errors on percent-deuteration input", {
  a <- output_tp(file_nm, percent = TRUE)
  expect_error(residue_deconvolution(a, protein_length = protein_length))
})

test_that("residue_level_average returns bounded, correctly-shaped output", {
  a <- output_tp(file_nm, percent = TRUE)
  res <- residue_level_average(a, protein_length = protein_length)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), protein_length * length(unique(a$Deut.Time)))

  value_cols <- grep("_Deut\\.\\._[0-9]+$", colnames(res), value = TRUE)
  expect_true(length(value_cols) > 0)
  vals <- unlist(res[, value_cols])
  ## real peptide-level percent deuteration can run slightly outside 0-100
  ## due to experimental/back-exchange-correction noise; this is a sanity
  ## check against gross errors, not a strict physical bound
  expect_true(all(vals >= -20 & vals <= 120, na.rm = TRUE))
})

test_that("residue_level_average errors on uptake (non-percent) input", {
  a <- output_tp(file_nm, percent = FALSE)
  expect_error(residue_level_average(a, protein_length = protein_length))
})

test_that("peptide_exch_residues() count matches nb_exch_deut() for every peptide", {
  ## regression test for the off-by-one where peptide_exch_residues()
  ## started at local index 2 (excluding only the first residue) instead of
  ## 3 (excluding the first two), one too many per peptide.
  ##
  ## nb_exch_deut() has its own pre-existing quirk (not touched here): its
  ## formula is length-2-nP_nonfirst, which subtracts a Proline at local
  ## position 2 TWICE - once via the generic "-2" and again because it is a
  ## "non-first" Proline - so it reports one fewer exchangeable residue than
  ## peptide_exch_residues() specifically when the peptide's second residue
  ## is a Proline. Accounted for explicitly below rather than silently
  ## tolerated, so a real future regression still fails this test.
  a <- output_tp(file_nm, percent = FALSE)
  one_per_peptide <- a[!duplicated(paste(a$Start, a$End, a$Sequence)), ]

  n_from_helper <- vapply(seq_len(nrow(one_per_peptide)), function(i) {
    length(peptide_exch_residues(one_per_peptide$Start[i], one_per_peptide$Sequence[i]))
  }, integer(1))
  second_residue_is_proline <- substr(one_per_peptide$Sequence, 2, 2) == "P"

  expect_equal(n_from_helper, nb_exch_deut(one_per_peptide) + as.integer(second_residue_is_proline))
})

test_that("lambda_smooth actually changes the deconvolution result", {
  a <- output_tp(file_nm, percent = FALSE)
  res0 <- residue_deconvolution(a, protein_length = protein_length, lambda_smooth = 0)
  res1 <- residue_deconvolution(a, protein_length = protein_length, lambda_smooth = 1)

  value_cols <- grep("_X\\.\\.Deut_[0-9]+$", colnames(res0), value = TRUE)
  vals0 <- unlist(res0[, value_cols])
  vals1 <- unlist(res1[, value_cols])
  expect_false(isTRUE(all.equal(vals0, vals1)))
})
