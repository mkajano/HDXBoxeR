#' function from plotfunctions package
#'
#' Margin coordinates
#' @param pos position
#' @param side side of plot
#' @param input plot or figure position
#' @return  coordinates of margins
#' @export
getCoords1<-function (pos = 1.1, side = 1, input = "p")
{
  p <- par()
    x.width = p$usr[2] - p$usr[1]
    y.width = p$usr[4] - p$usr[3]
    out <- rep(NA, length(pos))
    if (length(side) == 1) {
      side <- rep(side, length(pos))
    }
    out[which(side %in% c(1, 3))] <- pos[which(side %in%
                                                 c(1, 3))] * x.width + p$usr[1]
    out[which(side %in% c(2, 4))] <- pos[which(side %in%
                                                 c(2, 4))] * y.width + p$usr[3]
    return(out)
}
