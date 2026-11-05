####################
#' Prepares output for HDX-MS for the deuteration uptake or percent deuteration for the time points.
#'
#' Returns a data frame organized for additional analysis.
#' In columns are deuteration uptake or percent deuteration data for the given protein states.
#' Function allows for writing csv with data, matching sequences of peptide.
#' Protein.States, Deut.times, or number of replicates can be specified.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @param replicates optional maximum number of replicates per state and time. By default every replicate is used, so states may have different numbers of replicates (e.g. 3 and 4); missing measurements are `NA`. Replicates are numbered in the order their experiments appear in the file.
#' @param min_replicates a peptide is kept only if it has at least this many replicate values in every state at that time. Default `2`, the minimum for a standard deviation and a t-test.
#' @param states function allows to choose what states should be used for analysis. Default all states are used.
#' @param times lists the deuteration times to be used in analysis. Default all states used.
#' @param seq_match How peptides are matched between states. `FALSE` (default): a peptide is kept only if the same Sequence, Start, End, Search.RT and Charge are found in every state (exact sequence identity). `TRUE`: peptides are matched by position only (Start, End, Search.RT and Charge), ignoring Sequence, so peptides covering e.g. a point mutation are kept; where the sequences differ between states they are joined with "." in the `Sequence` column.
#' @param csv Flag allowing saving the output as csv. With default csv="NA", data is not saved. If csv output is desided, provide output name.
#' @param percent What the values are. `FALSE` (default): deuterium uptake in Da. `TRUE`: percent deuteration as reported by the input (HDXExaminer's `Deut %`, i.e. relative to the theoretical maximum; for DynamX 100 x uptake / `MaxUptake`). `"FD"`: percent of the measured fully deuterated control, 100 x uptake / FD uptake of the same peptide and state (FD replicates averaged, charge states pooled); this corrects for back-exchange and D2O content. Peptides without a usable FD value (missing, or FD uptake <= 0) keep their row with empty (`NA`) values.
#' @return data frame with reorganized data where in columns is the deuteration uptake for Protein States. Attributes `"replicates"` (number of replicate columns per state), `"states"`, `"times"` and `"percent"` record how it was built.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(filepath=file_nm) ###all default parameters used
#'
#'
#' # all possible flags listed & percent deuteration output,
#' # with sequences matching for protein states.
#'
#' a<-output_tp(filepath=file_nm, replicates=3, states=c("bound", "Unbound"),
#' times=c("3.00s", "72000.00s"), seq_match=TRUE, csv="NA", percent=TRUE)
#' @param .arg_df advanced/internal: a pre-computed [arg_df()] data.frame for
#'   this same `filepath`, to skip re-reading and re-parsing the file from
#'   disk. Default `NULL` reads `filepath` as usual - only pass this when
#'   calling `output_tp()` many times for the same file (e.g. once per pair
#'   of states) and `arg_df(filepath)` has already been computed once by the
#'   caller; passing a data.frame that doesn't match `filepath` will
#'   silently produce wrong results, since nothing here re-checks it.
#' @export
output_tp<- function(filepath, replicates, states, times, seq_match=FALSE, csv="NA", percent=FALSE,
                     min_replicates=2, .arg_df=NULL){
  a<-if (is.null(.arg_df)) arg_df(filepath) else .arg_df
  if(missing(states)) { states=unique(a$Protein.State); .hdx_msg("Protein.States used:", states)}
  if(missing(times)) { times=arguments_call2(filepath, states); .hdx_msg("Deut.times used:", times)}
  if(missing(replicates)) replicates <- NULL

  if (!(isTRUE(seq_match) || identical(seq_match, FALSE))) {
    stop("seq_match must be TRUE or FALSE", call. = FALSE)
  }
  all1 <- .hdx_wide(a, states = states, times = times, replicates = replicates,
                    min_replicates = min_replicates, seq_match = seq_match,
                    percent = percent, layout = "tp",
                    fd = if (identical(percent, "FD")) .hdx_fd_table(.hdx_raw(filepath)) else NULL)
  .hdx_msg("Number of replicates used:",
           paste0(names(attr(all1, "replicates")), "=", attr(all1, "replicates")))

  if (csv=="NA"){
    message("no csv file written")
  } else {
    write.csv(all1, file=csv)}

  return(all1)}
