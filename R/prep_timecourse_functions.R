#' Prepares function for plotting averages in timecourse
#'
#' Preparatory function
#'
#' @param control_df dataframe of control
#' @param variant_df dataframe for variant
#' @param replicates number of replicates. Default set to 3.
#' @return dataframes with matched peptides in time course
#' @export
prep_timecourse_plot_ave<-function(control_df, variant_df, replicates=3){
  av_c<-  ave_timepoint(control_df, replicates)
  av_v<-  ave_timepoint(variant_df, replicates)
  comb_av<-merge(av_c, av_v, by = c('Start','End', 'Sequence', 'Search.RT','Charge'))
  comb_av<-arrange(comb_av, c(Start))
  sh_avc<-comb_av[,c(1:dim(av_c)[2])]
  sh_avv<-comb_av[,c(1:5, (dim(av_c)[2]+1):(dim(av_c)[2]*2-5))]
  return(list(sh_avc, sh_avv))
}

#' Prepares function for Critical interval for timecourses
#'
#' Preparatory function
#'
#' @param control_df_up dataframe of control
#' @param variant_df_up dataframe for variant
#' @param replicates number of replicates. Default set to 3.
#' @param pv_cutoff cut off of pvalue used in calculation of critical interval. Default set to 0.01
#' @return Critial interval for all sets
#' @export
prep_timecourse_plot_sd<-function(control_df_up, variant_df_up, replicates=3, pv_cutoff=0.01){
  sd_c<-  sd_timepoint(control_df_up, replicates)
  sd_v<-  sd_timepoint(variant_df_up, replicates)
  CI_all<-CI_tc(sd_c, sd_v, replicates, pv_cutoff)
  return(CI_all)}
