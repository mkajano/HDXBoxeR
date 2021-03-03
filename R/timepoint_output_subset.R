#' Returns data.frame with averages from analysis for uptake file, standard deviation, p-values for selected protein states.
#'
#' Returns information from analysis data.frame for selected protein states.
#' Sets are compared to the first state in the input file.
#'
#' @param filepath path to All.Data.csv input from HDX-Examiner.
#' @param states list of states that will be chosen for analysis
#' @return data frame with analysis for uptake file, standard deviation, p-values for selected protein states.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' names_states<- nm_states(file_nm)
#' summary<-timepoint_output_subset(file_nm, states=names_states)
#' @export
timepoint_output_subset<-function(filepath, states){
  df<-output_tp_states(filepath, states)
  pv1<-pv_timepoint(df)
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  df_List<-list(av1, s1, pv1)
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  return(bp)
}

#' Returns data.frame with averages from analysis for procent deuteration, standard deviation, p-values for selected protein states.
#'
#' Returns information from analysis data.frame for selected protein states.
#' Sets are compared to the first state in the input file.
#'
#' @param filepath path to All.Data.csv input from HDX-Examiner.
#' @param states list of states that will be chosen for analysis
#' @return data.frame analysis for procent deuteration, standard deviation, p-values for selected protein states.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' names_states<- nm_states(file_nm)
#' summary<-timepoint_output_proc_subset(file_nm, states=names_states)
#' @export
timepoint_output_proc_subset<-function(filepath, states){
  df<-output_tp_proc_states(filepath, states)
  pv1<-pv_timepoint(df)
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  df_List<-list(av1, s1, pv1)
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  return(bp)
}


