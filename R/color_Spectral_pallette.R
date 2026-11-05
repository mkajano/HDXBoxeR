#' Returns Spectral pallette with colors matching defined ranges
#'
#' Spectral pallette for timecourse data
#'
#' @param ranges vector of numbers. Should have the same mumber of positive and negative values and contain 0.
#' @param col additional color that should be first in the pallette.
#' @param colors_initial deprecated, use `col`.
#' @return color scheme for number
#' @examples
#' color_ranges_Spectral(ranges=c(-Inf, -100, -50, 0, 50, 100, Inf), col="white")
#' @export
color_ranges_Spectral<-function(ranges, col, colors_initial) {
  if (!missing(colors_initial)) { .hdx_deprecated("colors_initial", "col"); col <- colors_initial }
  Spectral<-brewer.pal(n = 11, name = "Spectral")
  cbr1<-c(col, rev(colorRampPalette(Spectral)(floor(length(ranges))-1)))
  return(cbr1)
}
