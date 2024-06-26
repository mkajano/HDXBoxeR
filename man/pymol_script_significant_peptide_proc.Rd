% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pymol_script_significant_peptide_proc.R
\name{pymol_script_significant_peptide_proc}
\alias{pymol_script_significant_peptide_proc}
\title{Writes a text files with pymol scripts to list significant peptides}
\usage{
pymol_script_significant_peptide_proc(
  input_proc,
  input_up,
  path = "",
  ranges = c(-Inf, seq(-30, 30, by = 10), Inf),
  pv_cutoff = 0.01,
  replicates = 3,
  order.pep = TRUE
)
}
\arguments{
\item{input_proc}{Dataframe with organized procent deuteration data. Input generated using output_tp(, percent=T) function.}

\item{input_up}{Dataframe with organized deuteration uptake. Input generated using output_tp() function.}

\item{path}{location where the Pymol scripts will be saved}

\item{ranges}{ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)}

\item{pv_cutoff}{p-value cutoff here set up to 0.01}

\item{replicates}{number of replicates in sample. Default set to 3.}

\item{order.pep}{flag allowing to either order peptide acccording to the peptide length (default), or to position in the protein sequence.}
}
\value{
pymol script with colors assigned per peptide
}
\description{
Function write a script that can be used in pymol to color structure.
Number of colors and corresponding to them ranges can be defined by user.
}
\examples{
\donttest{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a_up<- output_tp(file_nm)
a_proc<- output_tp(file_nm, percent=TRUE)
pymol_script_significant_peptide_proc(input_proc=a_proc,
input_up=a_up,  path=tempdir(),replicates=3, pv_cutoff=0.01,
ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf), order.pep=TRUE)
}
}
