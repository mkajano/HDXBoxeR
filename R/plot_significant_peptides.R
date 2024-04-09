#' Preparatory function for significant peptide plots
#'
#' Returns plot where significant peptides are colored in blue-red scheme.
#'
#' @param df average data frame. Generated using ave_timepoint() function.
#' @param pv pvalues dataframes calculated using pv_timepoint() function
#' @param sd standard deviation data.frame generated using sd_timepoint function
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @param nb_row number of peptides in each row. Plotting parameter.
#' @return plot with peptides which are significantly different between sets.
#' @export
peptide_pv_tp<-function(df, pv, sd, nb_row, ranges=c(-Inf, seq(-30, 30, by=10), Inf),
                        pv_cutoff=0.01, replicates=3){
  #####
  #preparation significance per residue & coverage
  cl1<-significant_peptide_uptake(df, pv, sd, pv_cutoff, replicates)
  cl1[is.na(cl1)]<- 0

  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')

  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])/df[,7]*10*cl1)
  fc.d[is.nan(fc.d)]<- 0
  fc.d[fc.d == "-Inf"]<- 0

  # fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])*10*cl1)###vector which has significant average
  #
  # for ( i in 1:dim(fc.d)[2]){
  #   fc.d[,i]<-(fc.d[,i]/((df[,7]+df[,i+7])*0.5)) }
  si.fv<-(fc.d)

  ranges_function(df,si.fv ) ###print ranges of data

  xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    si.fv[xli[i]< si.fv & si.fv < xli[i+1]] <- num_ass[i]}
  si.fv[si.fv==0]<- (-10000)
  si.fv<-abs(si.fv)-9999
  ######
  ##color set up

  cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("black"))


  for ( j in 1:dim(si.fv)[2]){
    y1=c(rep(1:nb_row, times=floor(dim(df)[1]/nb_row)), 1:(dim(df)[1]%%nb_row)) ##y values on the plot corresponding to peptides index
    plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])),
         ylim=c(nb_row, 0), yaxt="n", xaxt="n") ## mock plot, just to have it drawn correct limits set up

    for ( i in 1:dim(df)[1]){
      points(c(df[i,start_col], df[i,end_col]), c(y1[i], y1[i]), type="l",
             col=cbr1[si.fv[i,j]])}

    main_nm=str_sub(colnames(df)[j+7], end = -9, start=4)
    mtext(main_nm, side=1, outer=FALSE, line=0, cex=0.6)
    axis(3, at=seq(0, max(df[, end_col])+25, by=25),  tcl=-0.2, labels = F)
    axis(3, at=seq(0, max(df[, end_col])+25, by=50), labels = T, cex.axis=0.65)
    box(lwd=2)
  }}

#' Significant peptide plots.
#'
#' Returns plot where significant peptides are colored in blue-red scheme.
#'
#' @param df1 average data frame. Generated using ave_timepoint() function.
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @param nb_pep_row number of peptides in each row. Plotting parameter. Default set to 100.
#' @return plot with peptides which are significantly different between sets.
#' @export
plot_peptide_sig_tp<-function(df1, replicates=3, nb_pep_row=100,
                              ranges=c(-Inf, seq(-30, 30, by=10), Inf), pv_cutoff=0.01){
  pv1<-pv_timepoint(df1,  replicates)
  s1<-sd_timepoint(df1, replicates)
  av1<-ave_timepoint(df1, replicates)
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  par(mar = c(1.5, 1.5, 1.5, 1.5), oma=c(3,2.4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( k in(unique(av1$Deut.Time))){
    message(paste("For timepoint", k))
    a1=av1[av1$Deut.Time==k,]
    p1=pv1[pv1$Deut.Time==k,]
    sd1=s1[s1$Deut.Time==k,]
    peptide_pv_tp(a1, p1, sd1, nb_pep_row,ranges, pv_cutoff, replicates)
    #legend_heat_map_tp(av1)
    mtext(k, side=4, outer=FALSE, line=0.2, cex=0.7)}
  mtext(c("Residues"),  c(NORTH<-3),line=0, outer=TRUE, cex=0.7)

  }
