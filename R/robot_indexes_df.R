## rows of a time-course comparison significant (p < pvalue AND |delta| > CI)
## at >= min_timepoints time points; columns 7+ are the time points
.robot_sig_rows <- function(pv, av_c, av_v, CI_all, pvalue, min_timepoints = 1) {
  cols <- 7:ncol(av_c)
  sig <- sapply(cols, function(i) {
    s <- pv[, i] < pvalue & abs(av_c[, i] - av_v[, i]) > CI_all[i - 6]
    s & !is.na(s)
  })
  if (is.null(dim(sig))) sig <- matrix(sig, nrow = nrow(av_c))
  which(rowSums(sig) >= min_timepoints)
}

## significant peptides (Start, End, Sequence, Charge) of every state compared
## with states[1]; collected by peptide identity, not by row number, because
## each comparison is merged separately and row numbers differ between them
.robot_sig_peptides <- function(thP, th, states, replicates, pvalue, alpha, p_adjust, min_timepoints) {
  keys <- list()
  for (state in states[2:length(states)]) {
    control_df_up <- th[th$Protein.State == states[1], ]
    variant_df_up <- th[th$Protein.State == state, ]
    pv1 <- pv_timecourse(df_c = control_df_up, df_v = variant_df_up, replicates, p_adjust = p_adjust)
    lav.proc_up <- prep_timecourse_plot_ave(control_df_up, variant_df_up, replicates)
    sh_avc_up <- lav.proc_up[[1]]; sh_avv_up <- lav.proc_up[[2]]
    CI_all <- prep_timecourse_plot_sd(control_df_up, variant_df_up, replicates = replicates, alpha = alpha)
    rows <- .robot_sig_rows(pv1, sh_avc_up, sh_avv_up, CI_all, pvalue, min_timepoints)
    keys[[state]] <- sh_avc_up[rows, c("Start", "End", "Sequence", "Charge")]
  }
  out <- unique(do.call(rbind, keys))
  out <- out[order(out$Start, out$End, out$Charge), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Returns dataframe with peptides which exhibit significant difference between two sets
#'
#' Function to help decide which peptides will be drawn on Robot plots.
#'
#'
#' @param thP output of output_tc(percent = TRUE). Raw data for procent deuteration for time courses
#' @param th output of output_tc(). Raw data for uptake deuteration for time courses
#' @param pvalue p-value cutoff. Default set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param states Protein states from the set. As default all states are chosen.
#' @param alpha cutoff for Critical interval. Default=0.01
#' @param p_adjust multiple-testing correction of the Welch p-values, passed to [stats::p.adjust()]: `"none"` (default, raw p-values) or e.g. `"BH"`. A peptide is significant when its (adjusted) p-value is below the cutoff AND its difference exceeds the critical interval.
#' @param min_timepoints a peptide is selected when it is significant at at least this many time points. Default `1` (any time point, as before); `2` or more ignores one-off calls at a single time point.
#'@return Returns dataframe listing peptides that are significantly different between sets.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' tm_df<-output_tc(filepath=file_nm)
#' tmP_df<-output_tc(filepath=file_nm, percent=TRUE)
#'
#' # more restictive peptide selection
#' robot_indexes_df(thP = tmP_df, th=tm_df, pvalue=0.001, alpha=0.01)
#' @export

robot_indexes_df<-function(thP, th, replicates=3,
                           pvalue=0.01, alpha=0.01, states, p_adjust="none", min_timepoints=1){
  if(missing(states)) states=unique(thP$Protein.State)
  sig <- .robot_sig_peptides(thP, th, states, replicates, pvalue, alpha, p_adjust, min_timepoints)
  sig_pep <- data.frame(nb = seq_len(nrow(sig)), sig)
  return(sig_pep)}

#' Returns indexes for peptides with significant difference between two sets
#'
#' Function to help decide which peptides will be drawn on Robot plots.
#'
#'
#' @param thP output of output_tc(percent = TRUE). Raw data for procent deuteration for time courses
#' @param th output of output_tc(). Raw data for uptake deuteration for time courses
#' @param pvalue p-value cutoff. Default set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param states Protein states from the set. As default all states are chosen.
#' @param alpha critical interval, to have more restrictive use lower values, default=0.01
#' @param p_adjust multiple-testing correction of the Welch p-values, passed to [stats::p.adjust()]: `"none"` (default, raw p-values) or e.g. `"BH"`. A peptide is significant when its (adjusted) p-value is below the cutoff AND its difference exceeds the critical interval.
#' @param min_timepoints a peptide is selected when it is significant at at least this many time points. Default `1` (any time point, as before).
#' @return Returns indexes of significant peptides: row numbers in the
#'   states\[1\] vs states\[2\] table used by the robot plots. A significant
#'   peptide that is not in that table (e.g. its sequence differs in a
#'   mutant) is left out; [robot_indexes_df()] lists every significant peptide.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' tm_df<-output_tc(filepath=file_nm)
#' tmP_df<-output_tc(filepath=file_nm, percent=TRUE)
#'
#'  # more restictive peptide selection
#' robot_indexes(thP = tmP_df, th=tm_df, pvalue=0.001, alpha=0.01)
#' @export
robot_indexes<-function(thP, th, replicates=3,
                        pvalue=0.01, states, alpha=0.01, p_adjust="none", min_timepoints=1){
  if(missing(states)) states=unique(thP$Protein.State)
  sig <- .robot_sig_peptides(thP, th, states, replicates, pvalue, alpha, p_adjust, min_timepoints)
  ## row numbers in the states[1] vs states[2] table drawn by robot_2states_indexes()
  sh_avc <- prep_timecourse_plot_ave(thP[thP$Protein.State == states[1], ],
                                     thP[thP$Protein.State == states[2], ], replicates)[[1]]
  key <- function(d) paste(d$Start, d$End, d$Sequence, d$Charge, sep = "\r")
  peptide_all <- sort(which(key(sh_avc) %in% key(sig)))
  return(peptide_all)}

