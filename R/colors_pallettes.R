#' Returns color pallete from red to blue with number of colors for defined ranges
#'
#'
#'
#' @param ranges vector of numbers. Should have the same mumber of positive and negative values and contain 0.
#' @param col additional color that should be first in the pallette.
#' @param colors_initial deprecated, use `col`.
#' @return color scheme for number
#' @examples
#' color_ranges_Blue_Red_heat_map(ranges=c(-Inf, -100, -50, 0, 50, 100, Inf), col="white")
#' @export
color_ranges_Blue_Red_heat_map<-function(ranges, col, colors_initial) {
  if (!missing(colors_initial)) { .hdx_deprecated("colors_initial", "col"); col <- colors_initial }

  blue1<-brewer.pal(n = 9, name = "Blues")[4:8]
  red1<-brewer.pal(n = 9, name = "Reds")[4:8]
  cbr1<-c(col, rev(colorRampPalette(red1)(floor(length(ranges)/2))),
          colorRampPalette(blue1)(floor(length(ranges)/2) ))
  return(cbr1)
}
