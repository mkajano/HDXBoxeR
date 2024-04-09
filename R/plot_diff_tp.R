#' Preparatory function for difference plot
#'
#' Returns plots with difference deuteration at each peptide.
#'
#' @param df output from functions output_tp or output_tp_proc.
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return plots of difference in average
#' @export

dif_tp<-function(df, cola) {
  if(missing(cola)) cola=c(1, brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"));

  n1=max((dim(df)[2]-7), 3)
  #cola<-c(brewer.pal(n = n1, name = "Paired"))
  plot(df[,7], type="n", xlab="", ylab="", lwd=2, col=cola[1],
       ylim=c(min(df[,7:dim(df)[2]])-0.5, max(df[,7:dim(df)[2]])+0.5))
  abline(h=0)
  for ( i in 8:dim(df)[2]){
    points(df[,i], type="l", xlab="", ylab="", col=cola[i-7])}
  axis(1, at=seq(0, 1000, by=10), cex.axis=1, labels=FALSE,tcl=-0.2)
  axis(2, at=seq(-1000, 1000, by=1), cex.axis=1, labels=FALSE,tcl=-0.2)

}


#' Legend for difference in averages plot.
#'
#' Returns legend for difference in average plots. Preparatory function.
#'
#' @param df output from functions average difference
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return legend for difference in average plot for time points
#' @export
lab_dif<-function(df, cola){
  if(missing(cola)) cola=c(1, brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"));
  yloc =1
  yadj =0
  xloc =getCoords1(1, input = "p")
  xadj =0
  n1=max((dim(df)[2]-7), 3)
  #cola<-c(brewer.pal(n = n1, name = "Paired"))
  nmx<-str_sub(colnames(df[8:dim(df)[2]]), start = 4, end=-9)
  nm1<-c(paste(str_sub(colnames(df)[7], start = 4, end=-9), "-"), nmx)
  cola2<-c("white", cola)
  legend(x = xloc, y = yloc,  legend = nm1,xjust = xadj,
         yjust = yadj,
         fill=cola2,  bty="n", cex=0.6, xpd = TRUE, border = cola2 )

}



#' Returns difference in average plot for timepoints in the data frame
#'
#' Returns plots with difference in avarage for each peptide.
#'
#' @param df output from functions output_tp or output_tp_proc.
#' @param replicates number of replicates in set as default set to 3.
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return plots of difference of averages
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' plots_diff_tp(df=a, replicates=3, cola=c(1:4))
#' plots_diff_tp(df=a)
#' @export
plots_diff_tp<-function(df, replicates=3, cola){
  if(missing(cola)) cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df, replicates))[2]-7), 3), name = "Paired"));

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow=c(length(unique(df$Deut.Time)), 1),mar = c(1.4, 1, 1, 9), oma=c(4,4,1,0.1), cex.axis=1, cex.main=1, cex.lab=1.1,
      mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)

  av1<-ave_timepoint(df, replicates)
  da1<-dif_ave(av1)
  for ( i in(unique(df$Deut.Time))){
    dif_tp(da1[da1$Deut.Time==i,], cola)

    mtext(i,  c(NORTH<-3),line=-1, outer=FALSE, cex=0.5)}
  exp_ddu<-expression(Delta*' D'[2]*'O uptake [Da]')
  mtext(c("Index"),  c(SOUTH<-1),line=0.7, outer=TRUE, adj=0.36)
  mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE)
  par(mfrow=c(1, 1))
  lab_dif(da1, cola)
}

