
#' Prepares output for deuteration uptake of HDX-MS for the time points
#'
#' In reprocessed data frame columns have time points data for uptake data.
#' Function takes raw data from HDXexaminer and organize data into columns with uptake of deuterons
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param output name of the output csv file
#' @return Returns reprocessed data. Data.frame with reorganized data where in columns are deuteration uptake for different timepoints in the analysis.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp_csv(file_nm, "output.csv")
#' @export
output_tp_csv<-function(filepath, output){
  write.csv(output_tp(filepath), output)}

#' Prepares output for procent deuteration of HDX-MS for the time points
#'
#' In reprocessed data frame columns have time points data for procent deuteration.
#' Function takes raw data from HDXexaminer and organize data into columns with procent deuteration.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param output name of the output csv file
#' @return Returns reprocessed data. Data.frame with reorganized data where in columns are procent deuteration for different timepoints in the analysis.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp_proc_csv(file_nm, "output.csv")
#' @export
output_tp_proc_csv<-function(filepath, output){
  write.csv(output_tp_proc(filepath), output)}

