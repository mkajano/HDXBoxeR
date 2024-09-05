#' Correlation functions plots and prep data
#'
#' Correlation of two states, the st.dev are shown as ellipsoids.
#'
#' @param s1 standard deviation from one sample
#' @param replicates number of replicates. Default set to 3.
#' @return treshold for determining significance.
#' @examples
#' sd1<-data.frame(c(0.1, 0.12, 0.13, 0.09, 0.11, 0.10))
#' CI_single(s1=sd1, replicates=3)
#' @export
#'

plot_cor_pts_elps<-function(df, replicates=3, times, state1, state2,
                            position_col=NA, length1=NA,  xlim = NULL, ylim = NULL, ...){
  av1<-ave_timepoint(df, replicates)
  sd1<-sd_timepoint(df, replicates)


  plot(av1[,7], av1[,7], type="n", xlab="", ylab="", xlim, ylim, ...)
  mtext(c("Deuteration (state 1)", "Deuteration (state 2)"),  c(SOUTH<-1, WEST<-2),line=0.7, outer=TRUE)
  points(c(-10, 110), c(-10, 110), type="l")


  if (is.na(position_col[1])==TRUE){
    for (i in 1:length(times)){
      ind1<-select_indices(av1, times = times[i], length = length1)
      ci<-CI_2pts(sd1[ind1,state1],sd1[ind1,state2], replicates )
      message(ci)
      points(c(-10, 110), c(-10, 110)+ci/2, type="l", col="grey")
      points(c(-10, 110), c(-10, 110)-ci/2, type="l", col="grey")

      for ( i in ind1){
        DrawEllipse(av1[i,state1],av1[i,state2]/2, sd1[i,state1]/2,
                    sd1[i,state2], ...)
      }
    }
  } else {
    for (i in 1:length(times)){
      ind1<-select_indices(av1, times = times[i] , length =length1)
      ci<-CI_2pts(sd1[ind1,state1],sd1[ind1,state2], replicates )
      points(c(-10, 110), c(-10, 110)+ci/2, type="l", col="grey")
      points(c(-10, 110), c(-10, 110)-ci/2, type="l", col="grey")

      message(ci)
      for ( i in ind1){
        DrawEllipse(av1[i,state1],av1[i,state2], sd1[i,state1]/2, sd1[i,state2]/2,
                    col = position_col[i], border= position_col[i], ...)
      }
    }
  }
}



plot_correlation_functions<-function(df, replicates=3){

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  av1<-ave_timepoint(df, replicates)
  par(mfrow=c(length(unique(av1$Protein.State)),1),
      mar = c(1.5, mar_x, 1, 1.1), oma=c(3,2.4,1,1),
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( i in(unique(df$Protein.State))){
    a1=av1[av1$Protein.State==i,]
    colmp<-heat_map_tc(a1, ranges)
    legend_heat_map_tc(a1)
    mtext(i, side=3, outer=FALSE, line=0, cex=0.65)}
  mtext(c("Residues"),  c(NORTH<-1),line=0.7, outer=TRUE, cex=0.8)
  return()}


