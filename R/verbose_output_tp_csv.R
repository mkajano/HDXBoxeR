
#' Returns csv with averages from analysis for uptake file, standard deviation, p-values.
#'
#' Returns information from analysis and save it as csv file.
#' Sets are compared to the first state in the input file.
#'
#' @param filepath path to All.Data.csv input from HDX-Examiner.
#' @param output_name name of the output in csv format.
#' @param replicates number of replicates used
#' @param ... other variables for output_tp
#' @return csv with analysis for uptake file, standard deviation, p-values for all protein states.
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' verbose_timepoint_output(file_nm, tempfile())
#' names_states<- nm_states(file_nm)
#' verbose_timepoint_output(file_nm, tempfile(), seq_match=TRUE, percent=TRUE,
#' states=names_states, replicates=3, times="3.00s")
#' }
#' @export
verbose_timepoint_output<-function(filepath, output_name, replicates=3, ...){
  df<-output_tp(filepath, ...)
  pv1<-pv_timepoint(df, replicates)
  s1<-sd_timepoint(df, replicates)
  av1<-ave_timepoint(df, replicates)
  df_List<-list(av1, s1, pv1)
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  write.csv(bp, output_name)
  return(bp)}
