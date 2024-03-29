% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/av_timecourse_plots.R
\name{plots_av_tcourse}
\alias{plots_av_tcourse}
\title{Generates average deuteration plot for the time-course.}
\usage{
plots_av_tcourse(df, replicates = 3, cola)
}
\arguments{
\item{df}{output from functions output_tcourse or output_tcourse_proc.}

\item{replicates}{number of replicates in set as default set to 3.}

\item{cola}{color pallette for different Protein States. As default Paired pallette from RColorBrewer is used.}
}
\value{
average deuteration plots
}
\description{
Returns plots with average deuteration at each peptide.
}
\examples{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tc(file_nm)
plots_av_tcourse(df=a, replicates=3, cola=c(1:4))
plots_av_tcourse(df=a)
}
