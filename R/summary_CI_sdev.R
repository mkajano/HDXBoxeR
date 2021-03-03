######
#' Provides summary table with Critical interval and standard deviation within the set.
#'
#' Returns summary data.
#' Function returns: Protein states, timepoints, number of replicates,  # peptides, % coveregae, average peptide length and redundancy.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param replicates number of replicates. Default set to 3.
#' @return Returns summary table.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- summary_sd_CI(file_nm, replicates=3)
#' @export
summary_sd_CI<-function(filepath, replicates=3){
a<- output_tp(filepath)
sd1<-sd_timepoint(a,replicates)
sds<-c()
cis<-c()
for ( i in 7:dim(sd1)[2]){
  sds<-c(sds,(mean(sd1[,i])))
  cis<-c(cis, CI_single(sd1[,i], replicates) )
}
nm1<-names(sd1)[7:dim(sd1)[2]]
nm1<-str_sub(nm1, start=4, end=-9)
sds<-round(sds,3)
cis<-round(cis, 3)
stat1<-c(sds, cis)
sum1<-data.frame(matrix(stat1,nrow=2, byrow = TRUE))
names(sum1)<-nm1
row.names(sum1)<-c("st.dev", "CI")
return(t(sum1))}

