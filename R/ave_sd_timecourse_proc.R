
average_timecourse_proc<-function(filepath){
  Start<-c()
  End<-c()
  Charge<-c()
  df<-output_tc(filepath, percent=T)
  fp<-output_FD_proc(filepath)
  ud<-output_UD_proc(filepath)

  ##average for timepoints
  av1<-ave_timepoint(df)

  ### get averages or non-averages values for FD and undeuterated.
  ###if function to check for number of replicates. if replicates =1 or smaller,
  if ((dim(fp)[2]-6) >1){
    av_fp<-ave_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) ==1) {
    av_fp=fp
  } else if ((dim(fp)[2]-6) == 0){
    message("Full deuteration sample not provided")
  }

  if ((dim(ud)[2]-6) >1){
    av_ud<-ave_timepoint(ud, replicates=(dim(ud)[2]-6))
  } else if ((dim(ud)[2]-6) ==1) {
    av_ud<-ud
  }  else if ((dim(ud)[1]+dim(ud)[2]) == 0){
    message("Non-deuterated sample not provided")
  }

  if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av_ud,av1, av_fp)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av1, av_fp)
  } else if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av_ud,av1)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av1)
  }

  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  return(bp)
}

#' Returns standard deviation for percent deuteration data for timecourses.
#'
#' Calculates standard deviation for time course data.
#'
#' @param filepath filepath to the All_results input file.
#' @return Data.frame with standard deviation.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' sd_timecourse(filepath=file_nm)
#' @export
sd_timecourse_proc<-function(filepath){
  Start<-c()
  End<-c()
  Charge<-c()
  df<-output_tc(filepath, percent=T)
  fp<-output_FD_proc(filepath)
  ud<-output_UD_proc(filepath)

  ##average for timepoints
  s1<-sd_timepoint(df)


  ### get averages or non-averages values for FD and undeuterated.
  ###if function to check for number of replicates. if replicates =1 or smaller,
  if ((dim(fp)[2]-6) >1){
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) ==1) {
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) == 0){
    message("Full deuteration sample not provided")
  }

  if ((dim(ud)[2]-6) >1){
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))

  } else if ((dim(ud)[2]-6) ==1) {
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
  }  else if ((dim(ud)[2]-6) == 0){
    message("Non-deuterated sample not provided")
  }

  if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list( sd_ud, s1, sd_fp)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list( s1, sd_fp)
  } else if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list( sd_ud, s1)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list( s1)
  }

  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  return(bp)
}
