

#' Plots heat maps for time courses.
#'
#' Returns heat map on timecourses with raw data.
#'
#' @param df timecourse input
#' @param ranges ranges for coloring scheme. Default set to c(seq(0, 100, by=10), Inf)
#' @return heat map for timecourses
#' @export
heat_map_tc<-function(df, ranges=c(seq(0, 100, by=10), Inf)){
  #####
  #preparation significance per residue & coverage
  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-df[,7:dim(df)[2]]###vector which has significant average
  ac<-data.frame(matrix(ncol= dim(df)[1], nrow =max(df[,end_col]), rep(0, length.out=max(df[,end_col])*dim(df)[1])))##make mock matrix

  sum.nm=0
  ac1<-c()
  #ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df[,end_col]), rep(0, length.out=dim(max(df[,end_col]))*dim(fc.d)[2])))##make mock matrix
  for ( j in 1:dim(fc.d)[2]){
    for ( i in 1:dim(df)[1]){

      sum.nm<-fc.d[i,j]
      ac[df[i,start_col]:df[i,end_col],i]<-sum.nm[1]
      ###create a vector of which is a sum of significance at position #assign new value to position
    }
    ac1<-c(ac1, rowSums(ac))}
  ac2=data.frame(matrix(ac1, ncol =dim(fc.d)[2] , byrow = FALSE))
  ##prep of coverage
  coverage<-coverage_residue(df, start_col, end_col)
  ave.p.cov<-ac2/coverage ## sums of the significant avererages divided by coverage.
  ave.p.cov[ave.p.cov=="NaN"]<-0 ##remove NAN divisions introduced by division by no coverage
  ranges_function_tc(df,ave.p.cov)
  #####preparation of average per residue data.frame, which will have
  xli=ranges; num_ass<-c(-10000001:(-10000001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    ave.p.cov[xli[i]< ave.p.cov & ave.p.cov < xli[i+1]] <- num_ass[i]
  }

  ave.p.cov[ave.p.cov==0]<- (-10000000)

  si_apc<-abs(ave.p.cov)-9999999
  cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
  si_apc<-si_apc*cv_mis+1

  ###define missing coverage
  cbr1<-color_ranges_Spectral(ranges=xli, c("white", 1))

  plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])-5),
       ylim=c(0, (dim(si_apc)[2])), yaxt="n") ## mock plot, just to have it drawn correct limits set up
  xl <- 1; yb <- (0); xr <- max(df[,end_col]); yt <- (1)
  for ( i in 0:(dim(si_apc)[2]-1)){
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


#' Plots heat maps for time courses.
#'
#' Returns heat map on timecourses with raw data.
#'
#'
#' @param df output from function output_tcourse
#' @param replicates number of replicates in sample. Default set to 3.
#' @param mar_x margin x width. Default=3.5
#' @param ranges ranges for coloring scheme. Default set to c(seq(0, 100, by=10), Inf)
#' @return heat map for time courses
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(file_nm)
#' plot_heat_map_tc(df=a, replicates=3, ranges=c(seq(0, 100, by=5), Inf))
#' plot_heat_map_tc(df=a)
#' @export
plot_heat_map_tc<-function(df, replicates=3, mar_x=3.5, ranges=c(-Inf, seq(0, 100, by=10), Inf)){

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


