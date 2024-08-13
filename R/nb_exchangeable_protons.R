#' Number of exchangeable protons
#'
#' Provides a vector with number of exchangeable protons, calculated from the input table.
#' Number of protons calculated as peptide_length - 2 - number of Prolines in the peptide
#' that are not in the first position
#'
#' @param df standard deviation from one sample
#' @return vector with number of exchangeable protons
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' nb_exch_deut(a)
#' @export
nb_exch_deut<-function(df){
  nb_exch_deut<-c(c(df$End-df$Start)-str_count(df$Sequence,'P') )+
  as.numeric(substr(df$Sequence, 1, 1)=="P")-1
  return(nb_exch_deut)
}
