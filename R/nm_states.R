####################
#' Lists names of states in data sets
#'
#' Returns vector with name of states used for choosing states for input functions generation.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return list of Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' names_states<- nm_states(file_nm)
#' @export
nm_states<- function(filepath){
  nm_states<-c()
  a<-read.csv(file=filepath,  header = FALSE, skip = 1)### load the file witout headers
  nm_states<-unique(a[,1])
  return(nm_states)
}
