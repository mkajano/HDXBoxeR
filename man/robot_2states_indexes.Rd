% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robot_2states_indexes.R
\name{robot_2states_indexes}
\alias{robot_2states_indexes}
\title{Returns a robot plot for selected peptides for 2 protein states.}
\usage{
robot_2states_indexes(
  thP,
  th,
  indexes,
  states,
  replicates = 3,
  pvalue = 0.01,
  ylim,
  xlim,
  CI_factor = 1
)
}
\arguments{
\item{thP}{output of output_tcourse_proc() function. Raw data for procent deuteration for time courses}

\item{th}{output of output_tcourse() function. Raw data for uptake deuteration for time courses}

\item{indexes}{indexes of peptides to be drawn.}

\item{states}{Need to choose only two protein states}

\item{replicates}{number of replicates in sample. Default set to 3.}

\item{pvalue}{p-value cutoff here set up to 0.01}

\item{ylim}{y-axis range}

\item{xlim}{x-axis range. Set as default from max and minimum residues for the protein}

\item{CI_factor}{Multiplication factor for Critical Interval. Allows for more restrictive selection of Critial interval.}
}
\value{
Robot maps for timecourses for 2 protein states and selected indexes.
}
\description{
Modification of butterfly plot. x axis residues.
y axis \% deuteration for one variant above the axis and for second peptide below the axis.
Peptides are compared between the sets for the significance change between sets.
If there is significant change beteween sets peptides are plotted for all timepoints.
Significanty different timepoints for the peptides are colored.
Peptides ranges are plotted as a line at corresponding \% deuteration values.
}
\examples{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
tm_df<-output_tc(filepath=file_nm)
tmP_df<-output_tc(filepath=file_nm, percent=TRUE)
names_states<- nm_states(file_nm) ### returns states names
ind1<-robot_indexes(thP = tmP_df, th=tm_df, pvalue=0.001, CI_factor=3, states=names_states[1:2])
robot_2states_indexes(thP = tmP_df, th=tm_df,
 states=names_states[1:2],indexes =ind1, pvalue=0.001, CI_factor=3)
}
