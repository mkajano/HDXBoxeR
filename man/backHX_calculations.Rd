% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backexchange_calculation.R
\name{backHX_calculations}
\alias{backHX_calculations}
\title{Summary of backexchange summary}
\usage{
backHX_calculations(filepath, Dfact = 0.85)
}
\arguments{
\item{filepath}{filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.}

\item{Dfact}{is the fraction of D/H in the labeling buffer used. Default set up to 0.85}
}
\value{
Returns summary table for backexchange.
}
\description{
Returns average and ranges of backexchange. Function calculates as: 1- (m100\%-m0\%)/N/Dfact.
m0\% is the non-deuterated peptide centroid mass, m100\% is the maximally labeled peptide centroid mass, N is
the theoretical number of backbone amides in the peptide and Dfrac is the fraction of D/H in the labeling buffer used.
Function requires undeuterated and Fully deuterated sets marked in Deut.time as 0s and FD respectively.
}
\examples{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- backHX_calculations(filepath=file_nm, Dfact=0.85)
}
