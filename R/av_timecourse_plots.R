#' Generates average deuteration plot for the time-course.
#'
#' Returns plots with average deuteration at each peptide.
#'
#' @param df output from functions output_tcourse or output_tcourse_proc.
#' @param replicates number of replicates in set as default set to 3.
#' @param cola color pallette for different Protein States. As default Paired pallette from RColorBrewer is used.
#' @return average deuteration plots
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(file_nm)
#' plots_av_tcourse(df=a, replicates=3, cola=c(1:4))
#' plots_av_tcourse(df=a)
#' @export
plots_av_tcourse<-function(df, replicates=3,cola){

  if(missing(cola)) cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  av1<-ave_timepoint(df, replicates)
  ppar(c(1,1))
  par(mar = c(1,1,1,1), mfrow=c(length(unique(av1$Protein.State)), 1), oma = c(3.5, 2.75, 1, 1))

  for ( i in(unique(av1$Protein.State))){
    av_tc(av1[av1$Protein.State==i,], cola)
    mtext(i,  c(South<-1), line=-1, outer=FALSE, cex=0.6)
    }
  legend_raw_ave_tc(av1, cola)
}


#' Preparatory function for average plot for timecourses
#'
#' Returns plots with average deuteration at each peptide.
#'
#' @param df output from functions output_tp or output_tp or output_tp_proc.
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return plots of averages
#' @export
av_tc<-function(df, cola) {
  if(missing(cola)) cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))
  plot(df[,7], type="n", xlab="", ylab="", lwd=2, cex.axis=0.8,
       col=cola[1], ylim=c(-5, max(df[,7:dim(df)[2]])+10))
  for ( i in 7:dim(df)[2]){
    points(df[,i], type="l", xlab="", ylab="", col=cola[i-6])}
  axis(1, at=seq(0, 1000, by=10), cex.axis=0.75, labels=FALSE,tcl=-0.2)
  axis(2, at=seq(0, 1000, by=5), cex.axis=0.75, labels=FALSE,tcl=-0.2)}

#' Legend for average deuteration plot for timecourse.
#'
#' Returns legend with average plots. Preparatory function.
#'
#' @param df output from functions output_tp or output_tp_proc.
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return legend for average plot for time course
#' @export
legend_raw_ave_tc<-function(df,cola){
  ##draw boxplots ave and sd1
  if(missing(cola)) cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))

  mtext(c("Index"),  c(SOUTH<-1),line=0.7, outer=TRUE, cex=0.8)
  mtext("% Deuteration",  c(WEST<-2),line=0.7, outer=TRUE, cex=0.8)
  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start = 4, end = -9)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", nm1, xpd = TRUE, horiz = TRUE, inset = c(0,   0), bty = "n",  fill = cola, cex = 0.8)}
