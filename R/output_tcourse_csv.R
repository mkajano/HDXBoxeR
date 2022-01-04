######
#' Prepares output for deuteration uptake of HDX-MS for the timecourses.
#'
#' In reprocessed data frame columns have time course data for procent deuteration.
#' Function takes raw data from HDXexaminer and organize data into columns with % deuteration
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param output name of the output csv file
#' @return Returns reprocessed data. Data.frame with reorganized data where in columns are deuterium uptake for different timepoints in the analysis.
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tcourse_csv(file_nm, "output.csv")
#' }
#' @export
output_tcourse_csv<-function(filepath, output){
  write.csv(output_tcourse(filepath=filepath), output)}


#' Prepares output for deuteration uptake of HDX-MS for the timecourses.
#'
#' In reprocessed data frame columns have time course data for procent deuteration.
#' Function takes raw data from HDXexaminer and organize data into columns with % deuteration
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param output name of the output csv file
#' @return Returns reprocessed data. Data.frame with reorganized data where in columns are procent deuteration for different timepoints in the analysis.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tcourse_proc_csv(file_nm, "output.csv")
#' @export
output_tcourse_proc_csv<-function(filepath, output){
  write.csv(output_tcourse_proc(filepath), output)}

