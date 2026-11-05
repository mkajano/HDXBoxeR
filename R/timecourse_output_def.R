####################
#' Prepares output for HDX-MS for the deuteration uptake or percent deuteration for the time courses.
#'
#' Returns a data frame organized for additional analysis.
#' In columns are deuteration uptake or percent deuteration data for the given protein states.
#' Function allows for writing csv with data, matching sequences of peptide.
#' Protein.States, Deut.times, or number of replicates can be specified.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @param replicates optional maximum number of replicates per state and time. By default every replicate is used; missing measurements are `NA`. Replicates are numbered in the order their experiments appear in the file.
#' @param min_replicates a peptide is kept for a state only if it has at least this many replicate values at every time. Default `2`.
#' @param states function allows to choose what states should be used for analysis. Default all states are used.
#' @param times lists the deuteration times to be used in analysis. Default all states used.
#' @param seq_match How peptides are matched between time points, as in [output_tp()]. `FALSE` (default): same Sequence, Start, End, Search.RT and Charge. `TRUE`: matched by position only, ignoring Sequence; differing sequences are joined with `"."`.
#' @param csv Flag allowing saving the output as csv. With default csv="NA", data is not saved. If csv output is desided, provide output name.
#' @param percent What the values are. `FALSE` (default): deuterium uptake in Da. `TRUE`: percent deuteration as reported by the input (HDXExaminer's `Deut %`, i.e. relative to the theoretical maximum; for DynamX 100 x uptake / `MaxUptake`). `"FD"`: percent of the measured fully deuterated control, 100 x uptake / FD uptake of the same peptide and state (FD replicates averaged, charge states pooled); this corrects for back-exchange and D2O content. Peptides without a usable FD value (missing, or FD uptake <= 0) keep their row with empty (`NA`) values. `"per_residue"`: uptake per exchangeable amide, Da / [nb_exch_deut()] (length - 2 - prolines after the first residue), so peptide length does not affect the critical-interval half of the significance rule (`TRUE` and `"FD"` are length-normalized too). The p-values are the same on every scale; only the critical interval and the differences change.
#' @return data frame with reorganized data where in columns is the deuteration uptake for the time points. Attributes `"replicates"` (replicate columns per time point), `"states"`, `"times"` and `"percent"` record how it was built.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(filepath=file_nm) ###all default parameters used
#'
#' # all possible flags listed & percent deuteration output,
#' #with sequences matching for protein states.
#'
#' a<-output_tc(filepath=file_nm, replicates=3, states=c("bound", "Unbound"),
#' times=c("3.00s", "72000.00s"), seq_match=TRUE, csv="NA", percent=TRUE)
#' @export
output_tc<- function(filepath, replicates, states, times, seq_match=FALSE, csv="NA", percent=FALSE,
                     min_replicates=2){
  a<-arg_df(filepath)
  if(missing(states)) { states=unique(a$Protein.State); .hdx_msg("Protein.States used:", states)}
  if(missing(times)) { times=arguments_call2(filepath, states); .hdx_msg("Deut.times used:", times)}
  if(missing(replicates)) replicates <- NULL

  if (!(isTRUE(seq_match) || identical(seq_match, FALSE))) {
    stop("seq_match must be TRUE or FALSE", call. = FALSE)
  }
  all1 <- .hdx_wide(a, states = states, times = times, replicates = replicates,
                    min_replicates = min_replicates, seq_match = seq_match,
                    percent = percent, layout = "tc",
                    fd = if (identical(percent, "FD")) .hdx_fd_table(.hdx_raw(filepath)) else NULL)
  .hdx_msg("Number of replicates used:",
           paste0(names(attr(all1, "replicates")), "=", attr(all1, "replicates")))

  if (csv=="NA"){
    message("no csv file written")
  } else {
    write.csv(all1, file=csv)}

  return(all1)}
