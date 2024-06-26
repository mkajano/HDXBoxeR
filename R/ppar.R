
####################
#' Preparation of figure window.
#'
#' Prepares a plotting window with specified margins with specific number of figure row and columns.
#'
#' @param mfrow2 mfrow: number of Multiple Figures (use ROW-wise).
#' @return  modified par function with adjusted parameters
#' @examples
#' ppar(c(2,1))
#' @export
ppar<-function(mfrow2){
  par(mfrow=mfrow2,mar = c(1.5, 1.5, 1, 1), oma=c(4,4,1,1), cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  return()
}

####################
#' Preparation of figure window. small margins
#'
#' Prepares a plotting window with specified margins with specific number of figure row and columns.
#'
#' @param mfrow2 mfrow: number of Multiple Figures (use ROW-wise).
#' @return  modified par function with adjusted parameters
#' @examples
#' pparLM(c(2,1))
#' @export
pparLM<-function(mfrow2){
  par(mfrow=mfrow2,mar = c(1.5, 1.5, 1, 1), oma=c(4,4,2,2), cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  return()
}

#' Preparation of figure window with area for figure at the bottom.
#'
#' Prepares a plotting window with specified margins with specific number of figure row and columns.
#'
#' @param mfrow2 mfrow: number of Multiple Figures (use ROW-wise).
#' @return  modified par function with adjusted parameters
#' @examples
#' ppar_bottom_legend(c(2,3))
#' @export
ppar_bottom_legend<-function(mfrow2){
  par(mfrow=mfrow2,mar = c(1.5, 1.5, 1.5, 1.5), oma=c(5, 5,5,5),
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  return()}  ### function I use to make fast changes in numer of plots in plot area


#' Preparation of figure window with more area on west side of plot.
#'
#' Prepares a plotting window with specified margins with specific number of figure row and columns.
#'
#' @param mfrow2 mfrow: number of Multiple Figures (use ROW-wise).
#' @return default plotting window
#' @examples
#' ppar_wider(c(2,1))
#' @export
ppar_wider<-function(mfrow2){
  par(mfrow=mfrow2,mar = c(1.5, 1.5, 1.5, 1.5), oma=c(3,2.4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  return()}  ### function I use to make fast changes in numer of plots in plot area

#' Reset plotting window parameters to default
#'
#' function by Farid Cheraghi,
#' https://stackoverflow.com/questions/9292563/reset-the-graphical-parameters-back-to-default-values-without-use-of-dev-off
#' function resets plotting window parameters
#'
#' @examples
#' reset_par()
#' @return default plotting window parameters
#' @export
  reset_par <- function(){
  op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                       ask = FALSE, bg = "transparent", bty = "o", cex = 1,
                       cex.axis = 1, cex.lab = 1, cex.main = 1.2, cex.sub = 1,
                       col = "black", col.axis = "black", col.lab = "black",
                       col.main = "black", col.sub = "black", crt = 0, err = 0L,
                       family = "", fg = "black", fig = c(0, 1, 0, 1),
                       fin = c(6.99999895833333, 6.99999895833333), font = 1L,
                       font.axis = 1L, font.lab = 1L, font.main = 2L,
                       font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L,
                       lend = "round", lheight = 1, ljoin = "round", lmitre = 10,
                       lty = "solid", lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42),
                       mar = c(5.1, 4.1, 4.1, 2.1), mex = 1, mfcol = c(1L, 1L),
                       mfg = c(1L, 1L, 1L,1L), mfrow = c(1L, 1L),
                       mgp = c(3, 1, 0), mkh = 0.001, new = FALSE,
                       oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1),
                       omi = c(0, 0, 0,0), pch = 1L,
                       pin = c(5.75999895833333, 5.15999895833333),
                       plt = c(0.117142874574832, 0.939999991071427,
                               0.145714307397962, 0.882857125425167),
                       ps = 12L, pty = "m", smo = 1, srt = 0, tck = NA_real_,
                       tcl = -0.5, usr = c(0.568, 1.432, 0.568, 1.432),
                       xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s",
                       xpd = FALSE, yaxp = c(0.6, 1.4, 4), yaxs = "r",
                       yaxt = "s", ylbias = 0.2),
                  .Names = c("xlog", "ylog", "adj", "ann", "ask", "bg",
                             "bty", "cex", "cex.axis", "cex.lab", "cex.main", "cex.sub",
                             "col", "col.axis", "col.lab", "col.main", "col.sub", "crt",
                             "err", "family", "fg", "fig", "fin", "font", "font.axis",
                             "font.lab", "font.main", "font.sub", "lab", "las", "lend",
                             "lheight", "ljoin", "lmitre", "lty", "lwd", "mai", "mar",
                             "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                             "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo",
                             "srt", "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd",
                             "yaxp", "yaxs", "yaxt", "ylbias"))
  par(op)
}
