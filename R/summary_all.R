######
#' Returns full summary table.
#'
#' Returns summary data.
#' Function returns: Protein states, timepoints, number of replicates,  # peptides, % coveregae, average peptide length and redundancy.
#' backexchange calculations (average and range), Critical interval and standard deviation.
#' Function requires undeuterated and Fully deuterated sets marked in Deut.time as 0s and FD respectively.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param replicates number of replicates. Default set to 3.
#' @param Dfact Dfact is the fraction of D/H in the labeling buffer used. Default set up to 0.85
#' @return Returns summary table.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- all_summary(file_nm, replicates=3, Dfact=0.85)
#' @export
all_summary<-function(filepath, replicates=3, Dfact=0.85){
  sum<-cbind(general_info(filepath), summary_sd_CI(filepath, replicates), backHX_calculations(filepath, Dfact))
  return(sum)
}


