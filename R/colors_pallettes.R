#' Returns color pallete from red to blue with number of colors for defined ranges
#'
#'
#'
#' @param ranges vector of numbers. Should have the same mumber of positive and negative values and contain 0.
#' @param colors_initial additional color that should be first in the pallette.
#' @return color scheme for number
#' @examples
#' color_ranges_Blue_Red_heat_map(ranges=c(-Inf, -100, -50, 0, 50, 100, Inf), colors_initial="white")
#' @export
color_ranges_Blue_Red_heat_map<-function(ranges, colors_initial) {

  blue1<-brewer.pal(n = 9, name = "Blues")[4:8]
  red1<-brewer.pal(n = 9, name = "Reds")[4:8]
  cbr1<-c(colors_initial, rev(colorRampPalette(red1)(floor(length(ranges)/2))),
          colorRampPalette(blue1)(floor(length(ranges)/2) ))
  return(cbr1)
}
