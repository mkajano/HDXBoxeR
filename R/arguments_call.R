
#' Returns default arguments for the output_tp functions. States
#'
#' Function used as internal function
#'
#' @param filepath input file location
#' @return The default arguments to output_tp functions.
#' @export
arguments_call1<-function (filepath) {
  a<-arg_df(filepath)
  states = unique(a$Protein.State)
  return(states)
}

#' Returns default arguments for the output_tp functions. Deut.Time
#'
#' Function used as internal function
#'
#' @param filepath input file location
#' @param states states used
#' @return The default arguments to output_tp functions.
#' @export
arguments_call2<-function (filepath, states) {

  a<-arg_df(filepath)
  un_times<-c()
  for ( i in states) {

    un_times<-c(un_times,unique(a[which(a$Protein.State == i),
                                  which(colnames(a) == "Deut.Time")])) }

  dt.df<-as.data.frame(table(un_times))

  vdtFT<-dt.df[,2]==length(states)

  if (length(unique(vdtFT))==1 & all(vdtFT) == FALSE) {
    stop("No common Deut.Times between the Protein.States, or Protein.States named incorrectly, program will halt")
  } else  {
    times<-as.vector(dt.df[which(dt.df[,2]==length(states)),1])
  }


  return(times)
}




#' Returns default arguments for the output_tp functions. # replicates
#'
#' Function used as internal function
#'
#' @param filepath input file location
#' @param states states used
#' @param times deuteration times
#' @return The default arguments to output_tp functions.
#' @export
arguments_call3<-function (filepath, states, times) {
  a<-arg_df(filepath)
  len_rep <- c()
  for (i in states) {
    for (j in times) {
      len_rep <- c(len_rep, length(unique(a[which(a$Protein.State ==
                                                    i & a$Deut.Time == j), which(colnames(a) == "Experiment")])))
    }
  }
  replicates = min(len_rep)
  return(replicates)
}


