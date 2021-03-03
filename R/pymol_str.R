#' Preparatory function writing pymol scripts
#'
#' Function rearrange vector to string by adding + sign between the numbers.
#'
#' @param ind1 vector of numbers (residues)
#' @return string with + as a separator.
#' @examples
#' res<-c(1,5, 19, 100, 109)
#' pymol_str(res)
#' @export
  pymol_str<-function(ind1) {ind2<-paste(as.character(ind1), sep="' '", collapse="+")
  return(ind2)}
