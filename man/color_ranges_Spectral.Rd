% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color_Spectral_pallette.R
\name{color_ranges_Spectral}
\alias{color_ranges_Spectral}
\title{Returns Spectral pallette with colors matching defined ranges}
\usage{
color_ranges_Spectral(ranges, colors_initial)
}
\arguments{
\item{ranges}{vector of numbers. Should have the same mumber of positive and negative values and contain 0.}

\item{colors_initial}{additional color that should be first in the pallette.}
}
\value{
color scheme for number
}
\description{
Spectral pallette for timecourse data
}
\examples{
color_ranges_Spectral(ranges=c(-Inf, -100, -50, 0, 50, 100, Inf), colors_initial="white")
}
