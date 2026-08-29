#' Per-residue uptake as a length-weighted mean across covering peptides
#'
#' A fast, simple companion to [residue_deconvolution()]: for each residue,
#' at each Protein.State/replicate/timepoint, computes a weighted mean of
#' percent deuteration across every peptide that covers that residue, with
#' each peptide weighted by `1 / nb_exch_deut(peptide)` - i.e. shorter, more
#' localized peptides contribute more than long ones. This is the same
#' peptide-to-residue weighting scheme used by PyHDX's residue mapping (its
#' "Z matrix", see Smit et al., PyHDX, <https://pyhdx.readthedocs.io/>)
#' before PyHDX goes on to fit a kinetic model per residue. This function
#' stops at the mapping step and does not fit any kinetic model or
#' cross-residue regularization - see Details.
#'
#' Unlike [residue_deconvolution()], this does not attempt to separate the
#' contribution of overlapping peptides algebraically, so a residue near the
#' edge of a long peptide will still carry some signal from its neighbors.
#' Provided as a quick, PyHDX-consistent baseline to compare against
#' [residue_deconvolution()].
#'
#' @param df data.frame, output of [output_tp()] called with `percent = TRUE`.
#'   Percent deuteration is required here (not uptake in Daltons): it is
#'   already normalized per peptide length, so a weighted average across
#'   peptides of different lengths is meaningful. Raw uptake is an extensive
#'   quantity (scales with peptide length) and averaging it directly across
#'   peptides of different lengths would not be - use
#'   [residue_deconvolution()] for uptake-based, length-aware reconstruction.
#' @param protein_length integer, total number of residues in the protein.
#' @return data.frame with one row per residue and `Deut.Time`, and one
#'   column per Protein.State/replicate, giving the weighted mean percent
#'   deuteration across all peptides covering that residue. `NA` where no
#'   peptide covers a residue.
#' @seealso [residue_deconvolution()], [output_tp()], [nb_exch_deut()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm, percent = TRUE)
#' res <- residue_level_average(a, protein_length = 225)
#' @export
residue_level_average <- function(df, protein_length) {

  value_cols <- grep("_Deut\\.\\._[0-9]+$", colnames(df), value = TRUE)
  if (length(value_cols) == 0) {
    stop("No percent-deuteration columns found - df must be the output of output_tp(..., percent = TRUE).")
  }

  residues <- seq_len(protein_length)
  timepoints <- unique(df$Deut.Time)

  out_list <- vector("list", length(timepoints))
  for (t in seq_along(timepoints)) {
    sub <- df[df$Deut.Time == timepoints[t], ]

    nb <- nb_exch_deut(sub)
    valid_w <- !is.na(nb) & nb > 0  ## a peptide with <=0 exchangeable residues
                                     ## (e.g. very short or Proline-rich) has no
                                     ## meaningful 1/nb_exch_deut weight
    w <- rep(NA_real_, length(nb))
    w[valid_w] <- 1 / nb[valid_w]

    ## which residues each peptide covers - same convention as nb_exch_deut()
    ## (peptide_exch_residues(), shared with residue_deconvolution())
    cov_mat <- matrix(FALSE, nrow(sub), protein_length)
    for (i in seq_len(nrow(sub))) {
      keep <- peptide_exch_residues(sub$Start[i], sub$Sequence[i])
      keep <- keep[keep >= 1 & keep <= protein_length]
      cov_mat[i, keep] <- TRUE
    }

    res_mat <- matrix(NA_real_, protein_length, length(value_cols),
                       dimnames = list(NULL, value_cols))
    for (r in residues) {
      covering <- cov_mat[, r] & valid_w
      if (!any(covering)) next
      wi <- w[covering]
      vals <- as.matrix(sub[covering, value_cols, drop = FALSE])
      ok <- !is.na(vals) & !is.na(wi)
      num <- colSums(vals * wi * ok, na.rm = TRUE)
      den <- colSums(wi * ok, na.rm = TRUE)
      res_mat[r, ] <- ifelse(den > 0, num / den, NA_real_)
    }
    out_list[[t]] <- data.frame(residue = residues, Deut.Time = timepoints[t],
                                 res_mat, check.names = FALSE)
  }

  do.call(rbind, out_list)
}
