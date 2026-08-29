#' Number of exchangeable protons
#'
#' Provides a vector with number of exchangeable protons, calculated from the input table.
#' Number of protons calculated as peptide_length - 2 - number of Prolines in the peptide
#' that are not in the first position. Note: if the peptide's *second* residue is a
#' Proline, it is excluded both by the generic "-2" and by the Proline count, so the
#' result is one lower than `peptide_exch_residues()`'s residue count for that peptide.
#'
#' @param df standard deviation from one sample
#' @return vector with number of exchangeable protons
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' nb_exch_deut(a)
#' @export
nb_exch_deut<-function(df){
  nb_exch_deut<-c(c(df$End-df$Start)-str_count(df$Sequence,'P') )+
  as.numeric(substr(df$Sequence, 1, 1)=="P")-1
  return(nb_exch_deut)
}

#' Exchangeable residue indices for one peptide (internal)
#'
#' Residue indices considered exchangeable for a single peptide - the first
#' two residues of the peptide and any Proline other than the first residue
#' are excluded. Matches the convention used by [nb_exch_deut()]; shared by
#' [residue_deconvolution()] and [residue_level_average()] so the two stay
#' consistent about which residues a peptide covers.
#'
#' @param start integer, peptide Start position (1-based residue index).
#' @param seq_string character, peptide Sequence.
#' @return integer vector of residue indices (relative to the full protein).
#' @noRd
peptide_exch_residues <- function(start, seq_string) {
  chars <- strsplit(seq_string, "")[[1]]
  if (length(chars) < 3) return(integer(0))
  local_idx <- 3:length(chars)
  resid_idx <- start + local_idx - 1
  resid_idx[chars[local_idx] != "P"]
}
