#' Global confidence interval treshold from experimental standard deviation for 2 samples.
#'
#' Calculation of global confidence interval using approach by:
#' Reliable Identification of Significant Differences in Differential Hydrogen Exchange-Mass Spectrometry Measurements Using a Hybrid Significance Testing Approach
#' Tyler S. Hageman and David D. Weis
#' Analytical Chemistry 2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325
#' calculations for alpha 0.99
#'
#' @param s1 standard deviation from one sample
#' @param s2 standard deviation from seconda sample
#' @param replicates number of replicates. Default set to 3.
#' @return treshold for determining significance.
#' @examples
#' sd1<-data.frame(c(0.1, 0.12, 0.13, 0.09, 0.11, 0.10))
#' sd2<-data.frame(c(0.18, 0.11, 0.13, 0.08, 0.11, 0.06))
#' CI_2pts(s1=sd1, s2=sd2, replicates=3)
#' @export
CI_2pts<-function(s1, s2, replicates=3){
  tvalue=abs(qt(0.01/2, replicates*2-2))
  sp1<-sqrt(sum(s1^2*(replicates-1))/((replicates-1)*length(s1)))
  sp2<-sqrt(sum(s2^2*(replicates-1))/((replicates-1)*length(s2)[1]))
  spa<-sqrt((sp1^2+sp2^2)/replicates)
  CI<-spa*tvalue ###if number of replicates changes this number has to change!!! Loook in the table!!!
  return(CI)
}

#' Global confidence interval treshold from experimental standard deviation
#'
#' Calculation of global confidence interval using approach by for all protein states compared to first state in the data.frame.
#' Reliable Identification of Significant Differences in Differential Hydrogen Exchange-Mass Spectrometry Measurements Using a Hybrid Significance Testing Approach
#' Tyler S. Hageman and David D. Weis
#' Analytical Chemistry 2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325
#'
#'
#' @param df standard deviation dataframe.
#' @param replicates number of replicates. Default set to 3.
#' @param alpha significance level. Set as default to 0.01
#' @return treshold for determining significance.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(file_nm, seq_match=FALSE)
#' sd<-sd_timepoint(df=a, replicates=3)
#' CI_tp(df=sd, replicates=3, alpha=0.01 )
#' CI_tp(sd)
#' @export
CI_tp<-function(df, replicates=3, alpha=0.01 ){
  CI_all<-c()
  tvalue=abs(qt(alpha/2, replicates*2-2))
  for ( i in 8:dim(df)[2]){
    sp1<-sqrt(sum(df[,7]^2*(replicates-1))/((replicates-1)*length(df[,7])))
    sp2<-sqrt(sum(df[,i]^2*(replicates-1))/((replicates-1)*length(df[,i])[1]))
    spa<-sqrt((sp1^2+sp2^2)/replicates)
    CI<-spa*tvalue
    CI_all<-c(CI_all, CI)}
  return(CI_all)
}
