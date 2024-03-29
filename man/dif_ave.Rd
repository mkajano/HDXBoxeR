% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sd_ave_diff_pv_tp.R
\name{dif_ave}
\alias{dif_ave}
\title{Returns data frame with difference of averages between State1 and other states provided.}
\usage{
dif_ave(df)
}
\arguments{
\item{df}{output from functions output_tp, output_tp_proc,  output_tp_states or output_tp_proc_states.}
}
\value{
Data.frame with difference values btw control and other protein states.
}
\description{
Returns average difference data.frame.
Sets are compared to the first state in the input file.
If other order of the sets is required use
Default for the number of replicates is 3.
}
\examples{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
pv<-pv_timepoint(df=a) ##if number of replicates is equal 3
pv1<-pv_timepoint(df=a, replicates=3) ##if number of replicates is equal 4
#b<-output_tp_states(file_nm, states=c("4EHP", "State2", "State3" ))
#pv_states<-pv_timepoint(df=b) ### here means of State4, will be compared to State2 and State4
}
