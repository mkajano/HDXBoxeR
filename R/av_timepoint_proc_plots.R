#' Preparatory function to draw legends for average procent
#'
#' Returns legend with average procent deuteration at each peptide.
#'
#' @param df output from functions output_tp or output_tp_proc.
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return legend for average deuteration procent for timepoints
#' @export
legend_raw_ave_proc<-function(df, cola) {
  if(missing(cola)) cola=c(1, brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"));
  #cola<-c(1, brewer.pal(n = n1, name = "Paired")
  ##draw boxplots ave and sd1
  exp_du<-expression('% D'[2]*'O ')
  mtext(c("Index"),  c(SOUTH<-1),line=0.7, outer=TRUE)
  mtext(exp_du,  c(WEST<-2),line=0.7, outer=TRUE)
  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start = 4, end = -10)
  legend(c("right"), nm1,
         fill=cola,  bty="n", cex=0.57, inset=c(-0.335,0), xpd = TRUE )
}

#' Returns average procent deuteration plot for time points
#'
#' Returns plots with average procent deuteration at each peptide.
#'
#' @param df output from functions output_tp_proc.
#' @param replicates number of replicates in set as default set to 3.
#' @param cola color pallette for different Protein States. As default Paired pallette from RColorBrewer is used.
#' @return average deuteration plots
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp_proc(file_nm)
#' plots_av_tp_proc(df=a, replicates=3, cola=c(1:4))
#' plots_av_tp_proc(df=a)
#' @export
plots_av_tp_proc<-function(df,replicates=3, cola) {
  if(missing(cola)) cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df, replicates))[2]-7), 3), name = "Paired"));
  par(mfrow=c(length(unique(df$Deut.Time)), 1),mar = c(1, 1, 1, 7), oma=c(4,4,1,2), cex.axis=1, cex.main=1, cex.lab=1.1,
      mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  av1<-ave_timepoint(df, replicates)
  for ( i in(unique(df$Deut.Time))){
    av_tp(av1[av1$Deut.Time==i,], cola)
    legend_raw_ave_proc(av1, cola)
    mtext(i,  c(NORTH<-3),line=-1, outer=FALSE, cex=0.5, adj=0.36)}
  reset_par()
}



