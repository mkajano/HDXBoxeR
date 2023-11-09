
#' Prepares the plot window for the woods functions
#'
#' Internal function
#'
#' @param df dataframe
#' @param ddlab label
#' @param ... other
#' @return Plot window
#' @export
pl_gen_ch2<-function(df,  ddlab=1, ...){
  plot(x = 1, type = "n",  xlim = c(min(df$Start),
                                    max(df$End)), ylab = "",
       xlab = "",   yaxt = "n", ...)


  axis(1, at = seq(0, 1000, by = 10), cex.axis = 1, labels = F,
       tcl = -0.2)
  axis(2, at = seq(-1000, 1000, by = 20), cex.axis = 1,
       labels = c(rev(seq(20, 1000, by = 20)), seq(0, 1000,
                                                   by = 20)))
  axis(2, at = seq(-1000, 1000, by = 10), cex.axis = 1,
       labels = F, tcl = -0.2)
  box(lwd=2)
  exp_ddu <- c(expression("% Deuteration"), expression(Delta * " %D"[2] * "O"))
  mtext(c("Residue"), c(SOUTH <- 1), line = 0.6,
        outer = TRUE, cex=1)
  mtext(exp_ddu[ddlab], c(WEST <- 2), line = 0.7, outer = TRUE)
  nb1 = 1

  abline(h=seq(-180, 180, by=10), col="grey70", lty=3, lwd=0.5)
  abline(h=0, lty=2)}


#' Legend, bottom of the plots
#'
#' Internal function
#'
#' @param df dataframe
#' @param cols colors
#' @return legend at the bottom of the plot
#' @export
legend_states_PerD_bottom<-function (df, cols)
{
  nm1 <- unique(str_sub(colnames(df[7:dim(df)[2]]), start = 1, end = -11))
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0,
                                                        0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", nm1, xpd = TRUE, horiz = TRUE, inset = c(0,
                                                            0), bty = "n", fill = cols, cex = 0.8)
}

#' Legend, bottom of the plots
#'
#' Internal function
#'
#' @param names labels
#' @param cols colors
#' @return legend at the bottom of the plot
#' @export
legend_nm_bottom<-function (names, cols)
{
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0,
                                                        0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", names, xpd = TRUE, horiz = TRUE, inset = c(0,
                                                              0), bty = "n", fill = cols, cex = 0.8)
}

#' Duplicate set function
#'
#' Internal function
#'
#' @param df dataframe
#' @return duplicate sets
#' @export
duplicate_sets <- function(df) {
  # assuming a data.frame
  xvec <- do.call("paste", c(df, sep = "\r"))
  uni_vec<-unique(xvec)
  lapply(uni_vec, function(x) which(xvec == x))
}


#' Prepares the plot window for the woods functions
#'
#' Internal function
#'
#' @param df dataframe
#' @param ddlab label
#' @param timepoints deuteration times used
#' @param ... other
#' @return Plot window
#' @export
pl_gen_uptake<-function(df, timepoints, ddlab=1, ...){
  plot(x = timepoints, type = "n",  ylab = "",
       xlab = "",   yaxt = "n", log="x", xlim=c(2, 100000), ...)


  axis(2, at = seq(-1000, 1000, by = 20), cex.axis = 1,
       labels = c(rev(seq(20, 1000, by = 20)), seq(0, 1000,
                                                   by = 20)))
  axis(2, at = seq(-1000, 1000, by = 10), cex.axis = 1,
       labels = F, tcl = -0.2)
  box(lwd=2)
  exp_ddu <- c(expression("% Deuteration"), expression(Delta * " %D"[2] * "O"))
  mtext(c("log(Time)"), c(SOUTH <- 1), line = 0.6,
        outer = TRUE, cex=1)
  mtext(exp_ddu[ddlab], c(WEST <- 2), line = 0.7, outer = TRUE)
  nb1 = 1

  abline(h=seq(-180, 180, by=10), col="grey70", lty=3, lwd=0.5)
  abline(h=0, lty=2)}

