#' Returns standard deviation for dataframe.
#'
#' Calculates standard deviation for the number of replicates in the function.
#'
#' @param df output from functions output_tp or output_tp_proc.
#' @param replicates number of replicates used. Default is set to replicates=3
#' @return Data.frame with standard deviation.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' sd<-sd_timepoint(df=a, replicates=3)
#' @export
sd_timepoint<-function(df, replicates=3) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/replicates
  nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  sd_nm<-paste("sd_", nm_root, sep="")

  sd1<-c(); for ( j in 1:nb_sets) {
    for (i in 1:dim(df)[1]) {sd1<-c(sd1, sd(df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]))}
  }
  sd2<-data.frame(matrix(sd1, ncol=nb_sets , byrow = FALSE))
  colnames(sd2)<-sd_nm
  sd2<-data.frame(df[,1:6], sd2)
  return(sd2)
}

#' Returns average value for either uptake of procent data.
#'
#' Calculates average of uptake or procent data. Returns data frame with average values.
#' Default for the number of replicates is 3.
#'
#' @param df output from functions output_tp or output_tp_proc.
#' @param replicates number of replicates used. Default is set to replicates=3
#' @return Data.frame with average values
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' ave<-ave_timepoint(df=a) ##if number of replicates is equal 3
#' ave<-ave_timepoint(df=a, replicates=4) ##if number of replicates is equal 4
#' @export
ave_timepoint<-function(df, replicates=3) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/replicates
  nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  ave_nm<-paste("av_", nm_root, sep="")

  ave1<-c(); for ( j in 1:nb_sets) {
    for (i in 1:dim(df)[1]) {ave1<-c(ave1, rowMeans(df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]))}
  }
  ave2<-data.frame(matrix(ave1, ncol=nb_sets , byrow = FALSE))
  colnames(ave2)<-ave_nm
  ave2<-data.frame(df[,1:6], ave2)
  return(ave2)
}

#' Calculation of pvalue between first protein state and any other state from all_states file
#'
#' Compares means of sets of uptake data and return dataframe with pvalues.
#' Welch t.test is used for analysis.
#' Sets are compared to the first state in the input file.
#' If other order of the sets is required use
#' Default for the number of replicates is 3.
#'
#' @param df output from functions output_tp or output_tp_proc.
#' @param replicates number of replicates used. Default is set to replicates=3
#' @return Data.frame with p-values
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' pv<-pv_timepoint(df=a) ##if number of replicates is equal 3
#' # pv1<-pv_timepoint(df=a, replicates=4) ##if number of replicates is equal 4
#' #b<-output_tp_states(file_nm, states=c("State4", "State2", "State3" ))
#' #pv_states<-pv_timepoint(df=b) ### here means of State4, will be compared to State2 and State4
#' @export
pv_timepoint<-function(df, replicates=3) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/replicates
  nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  pv_nm<-paste("pv_", nm_root, sep="")

  pv1<-c(); for ( j in 1:nb_sets) {
    for (i in 1:dim(df)[1]) {x1<-df[i,7:(7+replicates-1)]; x2<-df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]
    tt<-c(); tt<-t.test(x1, x2)
    pv1<-c(pv1, tt$p.val)}
  }

  pv2<-data.frame(matrix(pv1, ncol=nb_sets , byrow = FALSE))
  colnames(pv2)<-pv_nm
  pv2<-data.frame(df[,1:6], pv2)
  return(pv2)}


#' Returns data frame with difference of averages between State1 and other states provided.
#'
#' Returns average difference data.frame.
#' Sets are compared to the first state in the input file.
#' If other order of the sets is required use
#' Default for the number of replicates is 3.
#'
#' @param df output from functions output_tp, output_tp_proc,  output_tp_states or output_tp_proc_states.
#' @return Data.frame with difference values btw control and other protein states.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' pv<-pv_timepoint(df=a) ##if number of replicates is equal 3
#' pv1<-pv_timepoint(df=a, replicates=3) ##if number of replicates is equal 4
#' #b<-output_tp_states(file_nm, states=c("4EHP", "State2", "State3" ))
#' #pv_states<-pv_timepoint(df=b) ### here means of State4, will be compared to State2 and State4
#' @export
dif_ave<-function(df){
  da1<-data.frame(df[,1:6],df[,7:dim(df)[2]]-df[,7])
  return(da1)
}
