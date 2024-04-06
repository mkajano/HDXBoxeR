#' Returns csv with averages from analysis for procent deuteration file, standard deviation for time courses.
#'
#' Returns information from analysis and save it as csv file.
#' Sets are compared to the first state in the input file.
#'
#' @param filepath path to All.Data.csv input from HDX-Examiner.
#' @param output_name name of the output in csv format.
#' @param replicates number of replicates used
#' @param ... other variables for output_tc
#' @return csv with analysis for procent deuteration: standard deviation, for all protein states for time courses.
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' verbose_timecourse_output(file_nm,tempfile(), replicates=3)
#' names_states<- nm_states(file_nm)
#' verbose_timecourse_output(file_nm, tempfile(), seq_match=TRUE, percent=TRUE,
#' states=names_states, replicates=3, times="3.00s")
#' }
#' @export
verbose_timecourse_output<-function(filepath, output_name, replicates=3, ...){

  df<-output_tc(filepath, ...)

  ##average for timepoints
  s1<-sd_timepoint(df, replicates)
  av1<-ave_timepoint(df, replicates)

  df_List<-list(av1, s1)
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  write.csv(bp, output_name)
  return(bp)
}
