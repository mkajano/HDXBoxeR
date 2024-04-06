#' Preparatory function returns legends for the timecourses.
#'
#' Preparatory function
#'
#' @param df data frame from which names will be extracted
#' @param cols colors to be used in legend
#' @return legend at the bottom of the plot
#' @export

legend_tc_bottom<-function(df, cols){

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))

  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start = 5, end = -10)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", nm1, xpd = TRUE, horiz = TRUE, inset = c(0,   0), bty = "n",  fill = cols, cex = 0.8)
}
