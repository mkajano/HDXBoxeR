#' Checks for NaN is data.frame
#'
#' Function by Hong Ooi; https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame
#' @param x Data frame to be checked for NaN
#' @examples
#' ## this function will overwrite the is.nan function that works only on vectors and matrices
#' df<-data.frame(c(0,NaN), c(1, 2))
#' is.nan(df)
#' df[is.nan(df)]<- 0
#' @return logical. Returns info if data.frame contains NaNs.
#' @export
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
