#' Per-peptide quality checks of the input data
#'
#' Lists, for every peptide (start, end, sequence, charge) in every protein
#' state, simple checks of the controls and uptake values. Nothing is
#' removed from the data; use the table to decide what to exclude.
#'
#' Checks:
#' \describe{
#'   \item{`has_undeut`}{an undeuterated (`"0s"`) control is present. Without
#'     it HDXExaminer cannot compute uptake and reports 0.}
#'   \item{`all_zero_uptake`}{uptake is exactly 0 at every labeling time and
#'     in the FD control - in practice a missing reference, not a real
#'     measurement of no exchange.}
#'   \item{`negative_uptake`}{at least one replicate has uptake below 0.}
#'   \item{`fd_uptake`}{mean uptake (Da) of the fully deuterated (`"FD"`)
#'     control, `NA` if there is none.}
#'   \item{`max_uptake`}{theoretical maximum: length - 1 - prolines after
#'     the first residue (as HDXExaminer and DynamX count it).}
#'   \item{`fd_fraction`}{`fd_uptake / max_uptake`; 1 minus this is the
#'     peptide's back-exchange (without correction for the D2O fraction).}
#'   \item{`fd_missing`, `fd_nonpositive`, `fd_above_max`}{no FD value; FD
#'     uptake <= 0; FD uptake above the theoretical maximum.}
#'   \item{`issues`}{the checks that failed, as text (empty if none).}
#' }
#'
#' @param filepath input file or table, as for [output_tp()].
#' @return data frame with one row per state and peptide.
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' qc <- peptide_qc(file_nm)
#' qc[qc$issues != "", ]
#' @export
peptide_qc <- function(filepath) {
  raw <- .hdx_raw(filepath)
  raw$X..Deut <- suppressWarnings(as.numeric(raw$X..Deut))
  pep_key <- function(d) paste(d$Protein.State, d$Start, d$End, d$Sequence, d$Charge, sep = "\r")

  peps <- unique(raw[, c("Protein.State", "Start", "End", "Sequence", "Charge")])
  k <- pep_key(peps)
  kr <- pep_key(raw)

  undeut <- raw$Deut.Time == "0s"
  fd_rows <- raw$Deut.Time == "FD"
  labeled <- !undeut & !fd_rows & raw$Deut.Time != "0.00s"

  peps$has_undeut <- k %in% kr[undeut]
  vals <- raw$X..Deut[labeled | fd_rows]; kv <- kr[labeled | fd_rows]
  peps$all_zero_uptake <- vapply(k, function(x) {
    v <- vals[kv == x]; v <- v[!is.na(v)]
    length(v) > 0 && all(v == 0)
  }, logical(1), USE.NAMES = FALSE)
  peps$negative_uptake <- k %in% kr[labeled & !is.na(raw$X..Deut) & raw$X..Deut < 0]

  fd <- tapply(raw$X..Deut[fd_rows], kr[fd_rows], mean, na.rm = TRUE)
  peps$fd_uptake <- as.numeric(fd[k])
  peps$fd_uptake[is.nan(peps$fd_uptake)] <- NA_real_
  peps$max_uptake <- .hdx_theoretical_max(peps$Sequence)
  peps$fd_fraction <- round(peps$fd_uptake / peps$max_uptake, 3)
  peps$fd_missing <- is.na(peps$fd_uptake)
  peps$fd_nonpositive <- !peps$fd_missing & peps$fd_uptake <= 0
  peps$fd_above_max <- !peps$fd_missing & peps$fd_uptake > peps$max_uptake

  flags <- c(has_undeut = "no undeuterated control", all_zero_uptake = "uptake 0 at every time",
             negative_uptake = "negative uptake", fd_missing = "no FD control",
             fd_nonpositive = "FD uptake <= 0", fd_above_max = "FD above theoretical maximum")
  bad <- cbind(!peps$has_undeut, peps$all_zero_uptake, peps$negative_uptake,
               peps$fd_missing, peps$fd_nonpositive, peps$fd_above_max)
  peps$issues <- apply(bad, 1, function(b) paste(flags[b], collapse = "; "))

  peps <- peps[order(peps$Protein.State, peps$Start, peps$End, peps$Charge), ]
  rownames(peps) <- NULL
  peps
}
