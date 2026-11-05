
#' Critial interval calculation two sets of timecourses
#'
#' Preparatory function for calculation of pvalue between sets.
#'
#' @param sd_c dataframe of control
#' @param sd_v dataframe for variant
#' @param replicates number of replicates. Default set to 3.
#' @param alpha pvalue cutoff. Default set to 0.01
#' @return Critical interval for 2 sets
#' @export
CI_tc<-function(sd_c, sd_v, replicates=3, alpha=0.01 ){
  CI_all<-c()
  for ( i in 7:dim(sd_c)[2]){
    CI_all<-c(CI_all, CI_2pts(sd_c[,i], sd_v[,i], replicates, alpha,
                              n1 = .hdx_n_for(sd_c, colnames(sd_c)[i]),
                              n2 = .hdx_n_for(sd_v, colnames(sd_v)[i])))}
  return(CI_all)
}
