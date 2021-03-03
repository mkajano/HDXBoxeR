#' Global confidence interval treshold from experimental standard deviation for 1 sample
#'
#' Calculation of global confidence interval using approach by:
#' Reliable Identification of Significant Differences in Differential Hydrogen Exchange-Mass Spectrometry Measurements Using a Hybrid Significance Testing Approach
#' Tyler S. Hageman and David D. Weis
#' Analytical Chemistry 2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325
#' calculations for alpha 0.99
#'
#' @param s1 standard deviation from one sample
#' @param replicates number of replicates. Default set to 3.
#' @return treshold for determining significance.
#' @examples
#' sd1<-data.frame(c(0.1, 0.12, 0.13, 0.09, 0.11, 0.10))
#' CI_single(s1=sd1, replicates=3)
#' @export
CI_single<-function(s1, replicates=3){
  tvalue=abs(qt(0.01/2, replicates*2-2))
  sp1<-sqrt(sum(s1^2*(replicates-1))/((replicates-1)*length(s1)))
  spa<-sqrt((sp1^2)/replicates)
  CI<-spa*tvalue ###
  return(CI)
}

