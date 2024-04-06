#' Returns dataframe with peptides which exhibit significant difference between two sets
#'
#' Function to help decide which peptides will be drawn on Robot plots.
#'
#'
#' @param thP output of output_tcourse_proc() function. Raw data for procent deuteration for time courses
#' @param th output of output_tcourse() function. Raw data for uptake deuteration for time courses
#' @param pvalue p-value cutoff. Default set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param states Protein states from the set. As default all states are chosen.
#' @param CI_factor Multiplication factor for Critical Interval. Allows for more restrictive selection of Critial interval.
#' @return Returns dataframe listing peptides that are significantly different between sets.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' tm_df<-output_tc(filepath=file_nm)
#' tmP_df<-output_tc(filepath=file_nm, percent=TRUE)
#'
#' # more restictive peptide selection
#' robot_indexes_df(thP = tmP_df, th=tm_df, pvalue=0.001, CI_factor=3)
#' @export

robot_indexes_df<-function(thP, th, replicates=3,
                           pvalue=0.01, states, CI_factor=1){
  if(missing(states)) states=unique(thP$Protein.State)

  for ( state in states[2:length(states)]) {

    control_df<- thP[thP$Protein.State==states[1],]
    variant_df<- thP[thP$Protein.State==state,]

    control_df_up<- th[th$Protein.State==states[1],]
    variant_df_up<- th[th$Protein.State==state,]

    pv1<-pv_timecourse(df_c = control_df_up, df_v=variant_df_up, replicates)
    lav.proc<-prep_timecourse_plot_ave(control_df, variant_df, replicates)
    lav.proc_up<-prep_timecourse_plot_ave(control_df, variant_df, replicates)

    sh_avc<-lav.proc[[1]]
    sh_avv<-lav.proc[[2]]
    sh_avc_up<-lav.proc_up[[1]]
    sh_avv_up<-lav.proc_up[[2]]

    CI_all<-prep_timecourse_plot_sd(control_df_up, variant_df_up, replicates=3, pv_cutoff=pvalue)
    CI_all<-CI_all*CI_factor
    peptide_all<-c()
    for (i in 7:dim(sh_avc)[2]) {
      peptide_all<-c(peptide_all, which(pv1[, i]<pvalue &
                                          abs(sh_avc_up[, i]-sh_avv_up[, i]) > CI_all[i-6]))
    }
    peptide_all<-sort(unique(peptide_all))
  }

  sig_pep<-data.frame(1:length(peptide_all),sh_avc[peptide_all,c(1:3,5)])
  colnames(sig_pep)[1]<-"nb"
  return(sig_pep)}

#' Returns indexes for peptides with significant difference between two sets
#'
#' Function to help decide which peptides will be drawn on Robot plots.
#'
#'
#' @param thP output of output_tcourse_proc() function. Raw data for procent deuteration for time courses
#' @param th output of output_tcourse() function. Raw data for uptake deuteration for time courses
#' @param pvalue p-value cutoff. Default set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param states Protein states from the set. As default all states are chosen.
#' @param CI_factor Multiplication factor for Critical Interval. Allows for more restrictive selection of Critial interval.
#' @return Returns indexes of significant peptides
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' tm_df<-output_tc(filepath=file_nm)
#' tmP_df<-output_tc(filepath=file_nm, percent=TRUE)
#'
#'  # more restictive peptide selection
#' robot_indexes(thP = tmP_df, th=tm_df, pvalue=0.001, CI_factor=3)
#' @export
robot_indexes<-function(thP, th, replicates=3,
                        pvalue=0.01, states, CI_factor=1){
  if(missing(states)) states=unique(thP$Protein.State)

  for ( state in states[2:length(states)]) {

    control_df<- thP[thP$Protein.State==states[1],]
    variant_df<- thP[thP$Protein.State==state,]

    control_df_up<- th[th$Protein.State==states[1],]
    variant_df_up<- th[th$Protein.State==state,]

    pv1<-pv_timecourse(df_c = control_df_up, df_v=variant_df_up, replicates)
    lav.proc<-prep_timecourse_plot_ave(control_df, variant_df, replicates)
    lav.proc_up<-prep_timecourse_plot_ave(control_df, variant_df, replicates)

    sh_avc<-lav.proc[[1]]
    sh_avv<-lav.proc[[2]]
    sh_avc_up<-lav.proc_up[[1]]
    sh_avv_up<-lav.proc_up[[2]]

    CI_all<-prep_timecourse_plot_sd(control_df_up, variant_df_up, replicates=3, pv_cutoff = pvalue)
    CI_all<-CI_all*CI_factor
    peptide_all<-c()
    for (i in 7:dim(sh_avc)[2]) {
      peptide_all<-c(peptide_all, which(pv1[, i]<pvalue &
                                          abs(sh_avc_up[, i]-sh_avv_up[, i]) > CI_all[i-6]))
    }
    peptide_all<-sort(unique(peptide_all))

  }
  return(peptide_all)}

