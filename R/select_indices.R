#' Allows for selecting some peptide from input data
#'
#' Function allows for picking indices from the inputs based on:
#' peptide start or end residue, length, state or timepoint.
#' If parameters set to NA, condition is skipped.
#'
#' @param df input file (output of output_tc or output_tp)
#' @param start provide number for the staring residue, default NA
#' @param end  provide number for the end residue, default NA
#' @param length provide max length of the peptide
#' @param times timepoints, only for the output_tp functions
#' @param states states, only for the output_tc functions
#' @return Row indices of the peptides that are fulfilling the conditions required.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' indb<-select_indices(a,length=12, start=100, end=200)
#' smaller_df<-a[indb,]
#' @export
select_indices<-function(df, start=NA, end=NA, length=NA,
                         times=NA, states=NA){
  if(is.na(start)==FALSE){
    condition1<-c(df$Start >= start)
  } else (condition1<-1:dim(df)[1])
  if(is.na(end)==FALSE){
    condition2<-df$End <= end
  } else (condition2<-1:dim(df)[1])
  if(is.na(length)==FALSE){
    condition3<-c((df$End -df$Start) <= length)
  } else (condition3<-1:dim(df)[1])
  if(is.na(times)[1]==FALSE){
    condition4<- rowSums(sapply(times, function(x) df$Deut.Time==x))
  } else (condition4<-1:dim(df)[1])
  if(is.na(states)[1]==F){
    condition5<-rowSums(sapply(states, function(x) df$Protein.State==x))
  } else (condition5<-1:dim(df)[1])

  indb<-which(condition1 &condition2 &condition3 & condition4 & condition5)
  return(indb)}

