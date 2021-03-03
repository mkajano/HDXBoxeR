#' Returns coverage per residue
#'
#' returns vector with coverage information
#'
#' @param df1 output from functions output_tp or output_tp_proc.
#' @param start_col number of "Start" column in data.frame
#' @param end_col   number of "Start" column in data.frame
#' @return vector with coverage per residue
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' coverage_residue(df1=a,start_col=2, end_col=3 )
#' @export
coverage_residue<-function(df1, start_col, end_col){

  ##prep of coverage
  vc<-c()
  vc1<-c()
  for ( i in 1:dim(df1)[1]){
    vc<-rep(0, length=max(df1[,end_col]))
    vc[df1[i,start_col]:df1[i,end_col]]=1##make multiple vectors which have 1 at position which peptide covers
    vc1<-c(vc1, vc)}
  cov.pos=data.frame(matrix(vc1, ncol =max(df1[,end_col]) , byrow = TRUE))  ##make coverage dataframe
  coverage<-colSums(cov.pos) ##sumarized occurances of peptides
  return(coverage)
}


#' Function returns which peptides are significantly based of pv_cutoff and Critial interval
#'
#' Returns data frame with significant peptides.
#'
#' @param df_av data.frame with averages created using ave_timepoint() function
#' @param pv data.frame with pvalues created using pv_timepoint() function
#' @param sd data.frame with standard deviations created using sd_timepoint() function
#' @param replicates number of replicates as default set to 3.
#' @param pv_cutoff cuttoff for Critical interval. Default=0.01
#' @return ranges per set
#' @export
significant_peptide_uptake<-function(df_av, pv, sd, pv_cutoff=0.01, replicates=3){
  CI=CI_tp(sd, replicates,pv_cutoff)
  abs.a<-c()
  for ( i in 8:dim(df_av)[2]){
    abs1<-abs(df_av[,i])>CI[i-7]
    abs.a<-c(abs.a, abs1)
  }
  abs.a=data.frame(matrix(abs.a, ncol =length(CI) , byrow = FALSE))

  cl1a<-pv[,8:dim(df_av)[2]]<pv_cutoff
  cl1<-abs.a*cl1a
  return(cl1)
}
