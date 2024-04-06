#' Preparatory function for heat map for procent deuteration
#'
#' Returns heat map
#'
#' @param df average data frame for procent deuteration. Generated using ave_timepoint() function.
#' @param dfup average data frame for deuteration uptake. Generated using ave_timepoint() function.
#' @param pv pvalues dataframes calculated using pv_timepoint() function
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @param sd standard deviation data.frame generated using sd_timepoint function
#' @return heat map for timepoints
#' @export
heat_map_tp_proc<-function(df,dfup, pv, sd, ranges=c(-Inf, seq(-30, 30, by=10), Inf),
                           pv_cutoff=0.01, replicates=3){
  #####
  #preparation significance per residue & coverage

  cl1<-significant_peptide_uptake(dfup, pv, sd, pv_cutoff, replicates)

  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])*cl1)###vector which has significant average

  ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df[,end_col]), rep(0, length.out=max(df[,end_col])*dim(fc.d)[2])))##make mock matrix
  #ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df[,end_col]), rep(0, length.out=dim(max(df[,end_col]))*dim(fc.d)[2])))##make mock matrix
  for ( j in 1:dim(ac)[2]){
    for ( i in 1:dim(df)[1]){

      sum.nm<-(ac[df[i,start_col]:df[i,end_col],j]+fc.d[i,j])
      ###create a vector of which is a sum of significance at position
      ac[df[i,start_col]:df[i,end_col],j]=sum.nm[1] ## assign new value to position
    }}
  ##prep of coverage
  coverage<-coverage_residue(df, start_col, end_col)

  ave.p.cov<-ac/coverage ## sums of the significant avererages divided by coverage.
  ave.p.cov[ave.p.cov=="NaN"]<-0 ##remove NAN divisions introduced by division by no coverage

  ranges_function(df,ave.p.cov )

  #####preparation of average per residue data.frame, which will have
  xli=ranges; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    ave.p.cov[xli[i]< ave.p.cov & ave.p.cov < xli[i+1]] <- num_ass[i]
  }
  ave.p.cov[ave.p.cov==0]<- (-10000)

  si_apc<-abs(ave.p.cov)-9999
  cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
  si_apc<-si_apc*cv_mis+1

  ###define missing coverage
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

#' Plots heat maps for significant peptides.
#'
#' Returns heat map with average values for significant uptake per residue.
#'
#' @param input_proc Dataframe with organized procent deuteration data. Input generated using output_tp_proc() function.
#' @param input_up Dataframe with organized deuteration uptake. Input generated using output_tp() function.
#' @param mar_x margin x width. Default=3.5
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @return heat map for average uptake per residue for significant peptides.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a_up<- output_tp(file_nm)
#' a_proc<- output_tp(file_nm, percent=TRUE)
#' plot_heat_map_tp_proc(input_proc=a_proc, input_up=a_up, replicates=3, pv_cutoff=0.01,
#' ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf) )
#' plot_heat_map_tp_proc(input_proc=a_proc, input_up=a_up)
#' @export
plot_heat_map_tp_proc<-function(input_proc, input_up, mar_x=3.5, ranges=c(-Inf, -3,-2,-1, 0,1, 2,3, Inf),
                                pv_cutoff=0.01, replicates=3){

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  pv1<-pv_timepoint(input_up, replicates)
  s1<-sd_timepoint(input_up, replicates)
  av1<-ave_timepoint(input_proc, replicates)
  avu<-ave_timepoint(input_up, replicates)
  par(mfrow=c(length(unique(av1$Deut.Time)),1),
      mar = c(1.5, mar_x, 1, 1.1), oma=c(3,2.4,1,1),
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( i in(unique(av1$Deut.Time))){
    message(paste("For time point", i))
    a1=av1[av1$Deut.Time==i,]
    au=avu[avu$Deut.Time==i,]
    p1=pv1[pv1$Deut.Time==i,]
    sd1=s1[s1$Deut.Time==i,]
    colmp<-heat_map_tp_proc(a1,au, p1, sd1, ranges, pv_cutoff, replicates)
    legend_heat_map_tp_proc(av1)
    mtext(i, side=3, outer=FALSE, line=0, cex=0.65)}
  mtext(c("Residues"),  c(NORTH<-1),line=0.7, outer=TRUE, cex=0.8)
  }
