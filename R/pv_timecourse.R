#' pvalue calculation between two sets of the data at certain timepoint
#'
#' Preparatory function for calculation of pvalue between sets.
#'
#' @param df_c dataframe of control
#' @param df_v dataframe for variant
#' @param replicates number of replicates. Default set to 3.
#' @return pvalue comparisons between two sets.
#' @export
pv_timecourse<-function(df_c,df_v, replicates=3) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df_c)[2]-6)/replicates
  nm_root<-colnames(df_c[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  pv_nm<-paste("pv_", nm_root, sep="")

  pv1<-c(); for ( j in 1:nb_sets) {
    combined_df<-merge(df_c, df_v, by = c('Start','End', 'Sequence', 'Search.RT',
                                          'Charge'))
    every=(nb_sets*replicates)+7
    for (i in 1:dim(combined_df)[1]) {x1<-combined_df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)];
    x2<-combined_df[i,(every+(j-1)*replicates+1):(every+(j-1)*replicates+replicates)]
    tt<-c(); tt<-t.test(x1, x2)
    pv1<-c(pv1, tt$p.val)}
  }
  pv2<-data.frame(matrix(pv1, ncol=nb_sets , byrow = FALSE))
  colnames(pv2)<-pv_nm
  pv2<-data.frame(combined_df[,1:6], pv2)
  return(pv2)
}
