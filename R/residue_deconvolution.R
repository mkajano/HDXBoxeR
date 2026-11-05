#' Per-residue deconvolution of HDX-MS uptake from overlapping peptides
#'
#' Estimates a single-residue-resolution deuterium uptake profile from
#' overlapping peptide-level data, by solving a constrained least-squares
#' problem rather than simply averaging peptides that cover a residue.
#' Each peptide's observed uptake (in Daltons) is modelled as the sum of the
#' per-residue uptake of every exchangeable residue it covers (the first two
#' residues of each peptide and all non-first Prolines are excluded, matching
#' the convention used in [nb_exch_deut()]). The resulting system is solved per
#' Protein State, replicate, and timepoint with residue-level uptake
#' constrained to `\[0, 1\]` (a residue cannot be less than fully protected or
#' more than fully deuterated). With typical overlapping-peptide coverage
#' (redundancy > ~5x), this recovers sharper, more localized residue-level
#' resolution than a plain per-residue mean across covering peptides,
#' particularly at region boundaries where long peptides would otherwise
#' blend signal from neighboring residues into a region's apparent value.
#'
#' @section What it does not do:
#' Each time point is solved on its own. There is no kinetic model and no
#' link between time points. `lambda_smooth > 0` adds a penalty on the
#' difference between neighboring residues within one time point.
#'
#' @param df data.frame, output of [output_tp()] called with `percent = FALSE`.
#'   Uptake in Daltons is required because it is additive across residues;
#'   percent deuteration is not (it depends on each peptide's own number of
#'   exchangeable residues) and cannot be used directly for deconvolution.
#' @param protein_length integer, total number of residues in the protein.
#'   Used to set the residue index range for the deconvolution.
#' @param lambda_smooth non-negative numeric, default `0.1`. Strength of an L2
#'   penalty on the difference between adjacent covered residues (only
#'   applied between residues that are actually consecutive in sequence -
#'   a gap in coverage breaks the chain). `0` reproduces the original,
#'   unregularized per-timepoint fit. Larger values pull neighboring
#'   residues toward the same value unless the peptide data specifically
#'   distinguish them. The default of `0.1` is a modest starting point, not
#'   a universally correct value - the right strength depends on peptide
#'   redundancy and protein size. On the example data, values between 0 and
#'   1 predicted left-out peptides equally well.
#' @return data.frame with one row per residue and `Deut.Time`, and one column
#'   per Protein.State/replicate (same naming as the `X..Deut` columns in
#'   `df`), giving the deconvolved per-residue fractional uptake (0-1).
#'   Residues not covered by any peptide at a given timepoint are `NA`.
#' @seealso [nb_exch_deut()], [output_tp()], [coverage_residue()],
#'   [residue_level_average()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm, percent = FALSE)
#' res <- residue_deconvolution(a, protein_length = 225)
#' res_unsmoothed <- residue_deconvolution(a, protein_length = 225, lambda_smooth = 0)
#' @export
residue_deconvolution <- function(df, protein_length, lambda_smooth = 0.1) {

  value_cols <- grep("_X\\.\\.Deut_[0-9]+$", colnames(df), value = TRUE)
  if (length(value_cols) == 0) {
    stop("No uptake columns found - df must be the output of output_tp(..., percent = FALSE).")
  }

  residues <- seq_len(protein_length)
  timepoints <- unique(df$Deut.Time)

  ## peptide_exch_residues() (shared with residue_level_average(), defined in
  ## nb_exchangeable_protons.R) gives the exchangeable residue indices for one
  ## peptide, matching nb_exch_deut()'s convention.
  build_design_matrix <- function(sub) {
    A <- matrix(0, nrow(sub), protein_length)
    for (i in seq_len(nrow(sub))) {
      keep <- peptide_exch_residues(sub$Start[i], sub$Sequence[i], sub$End[i])
      keep <- keep[keep >= 1 & keep <= protein_length]
      A[i, keep] <- 1
    }
    A
  }

  ## adjacency (first-difference) matrix D restricted to residues that are
  ## both "covered" (in the optimization) and consecutive in raw sequence
  ## index - a coverage gap breaks the chain rather than being bridged
  adjacency_pairs <- function(covered_idx) {
    keep <- which(diff(covered_idx) == 1)
    if (length(keep) == 0) return(NULL)
    cbind(keep, keep + 1L)
  }

  deconvolve_column <- function(A, b, covered_idx) {
    ok <- !is.na(b)
    covered <- colSums(A[ok, , drop = FALSE]) > 0
    out <- rep(NA_real_, ncol(A))
    n <- sum(covered)
    if (n == 0) return(out)

    Ai <- A[ok, covered, drop = FALSE]
    bi <- b[ok]
    pairs <- if (lambda_smooth > 0) adjacency_pairs(covered_idx[covered]) else NULL

    obj <- function(x) {
      loss <- sum((Ai %*% x - bi)^2)
      if (!is.null(pairs)) loss <- loss + lambda_smooth * sum((x[pairs[,1]] - x[pairs[,2]])^2)
      loss
    }
    grd <- function(x) {
      g <- as.numeric(2 * crossprod(Ai, (Ai %*% x - bi)))
      if (!is.null(pairs)) {
        d <- x[pairs[,1]] - x[pairs[,2]]
        g[pairs[,1]] <- g[pairs[,1]] + 2 * lambda_smooth * d
        g[pairs[,2]] <- g[pairs[,2]] - 2 * lambda_smooth * d
      }
      g
    }

    fit <- stats::optim(rep(0.5, n), obj, grd, method = "L-BFGS-B",
                         lower = 0, upper = 1)
    out[covered] <- fit$par
    out
  }

  out_list <- vector("list", length(timepoints))
  for (t in seq_along(timepoints)) {
    sub <- df[df$Deut.Time == timepoints[t], ]
    A <- build_design_matrix(sub)

    res_mat <- matrix(NA_real_, protein_length, length(value_cols),
                       dimnames = list(NULL, value_cols))
    for (vc in value_cols) {
      res_mat[, vc] <- deconvolve_column(A, sub[[vc]], residues)
    }
    out_list[[t]] <- data.frame(residue = residues, Deut.Time = timepoints[t],
                                 res_mat, check.names = FALSE)
  }

  do.call(rbind, out_list)
}
