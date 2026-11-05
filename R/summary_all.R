######
#' Returns full summary table.
#'
#' Returns summary data.
#' Function returns: Protein states, timepoints, number of replicates,  # peptides, % coveregae, average peptide length and redundancy.
#' backexchange calculations (average and range), Critical interval and standard deviation.
#' Function requires undeuterated and Fully deuterated sets marked in Deut.time as 0s and FD respectively.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @param replicates optional maximum number of replicates per state.
#'   Default `NULL` uses every replicate.
#' @param Dfact Dfact is the fraction of D/H in the labeling buffer used. Default set up to 0.85
#' @param alpha Critical interval, level of stringency
#' @param protein_range residue range for the coverage, as in [general_info()].
#' @return Returns summary table. Back-exchange columns are `NA` when the
#'   input has no usable undeuterated/fully deuterated controls.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- all_summary(file_nm, Dfact=0.85)
#' @export
all_summary<-function(filepath, replicates=NULL, Dfact=0.85, alpha=0.01, protein_range=NULL){
  raw <- .hdx_raw(filepath)
  gi <- general_info(raw, protein_range = protein_range)
  sdci <- summary_sd_CI(raw, replicates, alpha = alpha)
  bx <- tryCatch(backHX_calculations(raw, Dfact), error = function(e) NULL)
  states <- gi$Protein.State
  bx_df <- data.frame(`<BackHX>` = NA_real_, `BackHX Range From` = NA_real_, To = NA_real_,
                      check.names = FALSE)[rep(1, length(states)), ]
  if (!is.null(bx)) {
    m <- match(make.names(gsub(" ", "", states)), make.names(gsub(" ", "", rownames(bx))))
    ok <- !is.na(m) & is.finite(bx[pmax(m, 1), 1])
    bx_df[ok, ] <- bx[m[ok], , drop = FALSE]
  }
  m2 <- match(.state_column_prefix(states), rownames(sdci))
  sum <- data.frame(gi, st.dev = sdci[m2, "st.dev"], CI = sdci[m2, "CI"], bx_df,
                    check.names = FALSE, row.names = states)
  return(sum)
}


