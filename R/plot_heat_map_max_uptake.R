#' Preparatory function for heat map of maximum uptake per residue.
#'
#' Returns heat map
#'
#' @param df average data frame. Generated using ave_timepoint() function.
#' @param pv pvalues dataframes calculated using pv_timepoint() function
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @param sd standard deviation data.frame generated using sd_timepoint function
#' @return maxiumum uptake heat map for timepoints
#' @export
heat_map_tp_maxuptake<-function(df, pv, sd, ranges=c(-Inf, seq(-30, 30, by=10), Inf),
                                pv_cutoff=0.01, replicates=3){
  #####
  #preparation significance per residue & coverage
  cl1<-significant_peptide_uptake(df, pv, sd, pv_cutoff, replicates)

  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])/(df[,7])*10*cl1)###vector which has significant average
  fc.d[fc.d=="-Inf"]<-0
  max.ac1<-c()
  for ( j in 1:dim(fc.d)[2]){
    ac<-c()
    ac1<-c()
    ac2<-c()
    for ( i in 1:dim(df)[1]){
      ac<-rep(0, length=max(df[,end_col]))
      ac[df[i,start_col]:df[i,end_col]]=fc.d[i,j]
      ##make multiple vectors which have 1 at position which peptide covers
      ac1<-c(ac1, ac)}
    ac2=data.frame(matrix(ac1, nrow =dim(df)[1], byrow=TRUE))
    max.a<-c()
    for ( k in 1:dim(ac2)[2]){
      ind1<-which.max(abs(ac2[,k]))
      nb1<-(ac2[ind1,k])
      max.a<-c(max.a, nb1)}
    max.ac1<-c(max.ac1, max.a)}
  max.ac2=data.frame(matrix(max.ac1, ncol = dim(fc.d)[2]))



  ##prep of coverage
  coverage<-coverage_residue(df, start_col, end_col)

  ranges_function(df,max.ac2 ) ###print ranges of data

  #####preparation of average per residue data.frame, which will have
  xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    max.ac2[xli[i]< max.ac2 & max.ac2 < xli[i+1]] <- num_ass[i]
  }
  max.ac2[max.ac2==0]<- (-10000)
  cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)
  si_apc<-abs(max.ac2)-9999
  ###define lack of coverage
  si_apc<-si_apc*cv_mis+1

  ###define missing coverage
  ##pallette definition
  ##color set up

  cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c( "grey45", "white"))


  plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])-5),
       ylim=c(0, (dim(df)[2]-7)), yaxt="n") ## mock plot, just to have it drawn correct limits set up
  xl <- 1; yb <- (0); xr <- max(df[,end_col]); yt <- (1)
  for ( i in 0:(dim(df)[2]-8)){
    yb=i; yt=i+1 ##loop to have initial values for y postions in loop to use multiple postion
    rect(head(seq(xl,xr,(xr-xl)/xr),-1),yb,
         tail(seq(xl,xr,(xr-xl)/xr),-1), yt,col=cbr1[si_apc[,i+1]], border = NA)

    ###coverage
    # rect(head(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc,yb,
    #      tail(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc, yt,col=cc[col_cv_mis+1], border = NA)
  }
  axis(1, at=seq(0, 700, by=25),  tcl=-0.2, labels = F)
  abline(h=0:7, lwd=0.5, col="grey30")
  box(lwd=2)
  return()
}

#' Plots heat maps for maximum uptake per residue.
#'
#' Returns heat map with maximum uptake per residue.
#'
#' @param df average data frame. Generated using ave_timepoint() function.
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param mar_x margin x width. Default=3.5
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @return heat map for maximum uptake per residue
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' plot_heat_map_max_uptake_tp(df=a, replicates=3, pv_cutoff=0.01,
#' ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf) )
#' plot_heat_map_max_uptake_tp(df=a)
#' @export
plot_heat_map_max_uptake_tp<-function(df, replicates=3,
                                      mar_x=3.5, ranges=c(-Inf, seq(-30, 30, by=10), Inf), pv_cutoff=0.01){

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  pv1<-pv_timepoint(df, replicates)
  s1<-sd_timepoint(df, replicates)
  av1<-ave_timepoint(df, replicates)
  par(mfrow=c(length(unique(av1$Deut.Time)),1),
      mar = c(1.5, mar_x, 1, 1.1), oma=c(3,2.4,1,1),
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( i in(unique(av1$Deut.Time))){
    a1=av1[av1$Deut.Time==i,]
    p1=pv1[pv1$Deut.Time==i,]
    sd1=s1[s1$Deut.Time==i,]
    message(paste("For timepoint", i))
    colmp<-heat_map_tp_maxuptake(a1, p1, sd1, ranges, pv_cutoff, replicates)
    legend_heat_map_tp(av1)
    mtext(i, side=3, outer=FALSE, line=0, cex=0.65)}
  mtext(c("Residues"),  c(NORTH<-1),line=0.7, outer=TRUE, cex=0.8)
  return()}



