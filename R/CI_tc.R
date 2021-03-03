
#' Critial interval calculation two sets of timecourses
#'
#' Preparatory function for calculation of pvalue between sets.
#'
#' @param sd_c dataframe of control
#' @param sd_v dataframe for variant
#' @param replicates number of replicates. Default set to 3.
#' @param pv_cutoff pvalue cutoff. Default set to 0.01
#' @return Critical interval for 2 sets
#' @export
CI_tc<-function(sd_c, sd_v, replicates=3, pv_cutoff=0.01 ){
  CI_all<-c()
  tvalue=abs(qt(pv_cutoff/2, replicates*2-2))

  for ( i in 7:dim(sd_c)[2]){
    sp1<-sqrt(sum(sd_c[,i]^2*(replicates-1))/((replicates-1)*length(sd_c[,i])))
    sp2<-sqrt(sum(sd_v[,i]^2*(replicates-1))/((replicates-1)*length(sd_v[,i])[1]))
    spa<-sqrt((sp1^2+sp2^2)/replicates)
    CI<-spa*tvalue
    CI_all<-c(CI_all, CI)}
  return(CI_all)
}
