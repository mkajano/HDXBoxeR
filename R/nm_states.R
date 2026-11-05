####################
#' Lists names of states in data sets
#'
#' Returns vector with name of states used for choosing states for input functions generation.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @return list of Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' names_states<- nm_states(file_nm)
#' @export
nm_states<- function(filepath){
  a <- .hdx_raw(filepath)
  nm_states<-unique(a$Protein.State)
  return(nm_states)
}
