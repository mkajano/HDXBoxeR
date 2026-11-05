######
#' Provides summary table with Critical interval and standard deviation within the set.
#'
#' For each protein state: the pooled replicate standard deviation of all
#' peptides and time points (each peptide's SD weighted by its replicates
#' minus one), and the critical interval of that state alone ([CI_single()])
#' using the state's average number of replicates.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @param replicates optional maximum number of replicates per state.
#'   Default `NULL` uses every replicate.
#' @param alpha critical interval, to have more restrictive use lower values, default=0.01
#' @return Returns summary table.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- summary_sd_CI(file_nm, alpha=0.01)
#' @export
summary_sd_CI<-function(filepath, replicates=NULL, alpha=0.01){
  a <- suppressMessages(output_tp(filepath, replicates = replicates))
  sd1 <- sd_timepoint(a)
  n <- attr(sd1, "n", exact = TRUE)
  roots <- colnames(n)
  sds <- vapply(seq_along(roots), function(j) .hdx_pooled_sd(sd1[[j + 6]], n[[j]]), numeric(1))
  cis <- vapply(seq_along(roots), function(j) {
    ok <- !is.na(sd1[[j + 6]]) & n[[j]] >= 2
    N <- mean(n[[j]][ok])
    sds[j] / sqrt(N) * abs(qt(alpha / 2, 2 * N - 2))
  }, numeric(1))
  out <- cbind(st.dev = round(sds, 3), CI = round(cis, 3))
  rownames(out) <- sub("_(X\\.\\.Deut|Deut\\.\\.)$", "", roots)
  out
}
