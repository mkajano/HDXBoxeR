#' Preparatory function for volcano plot legends
#'
#' Returns volcano plots
#'
#' @param df output from functions output_tp
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return legends for volcano plots
#' @export
lab_vol<-function(df, cola){
  if(missing(cola)) cola=c(brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"));
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
         fill=cola2,  bty="n", cex=0.6, xpd = TRUE, border = cola2 )}

#' Preparatory function for volcano plot
#'
#' Returns volcano plots
#'
#' @param df1 differences in averages data.frame calculated using diff_ave function
#' @param CI critical interval, here is multiple sets are using maximun CI is used.
#' @param pv pvalues dataframes calculated using pv_timepoint function
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return volcano plots
#' @export
vol_tp<-function(df1, pv, CI, pv_cutoff=0.01, cola){
  if(missing(cola)) cola=c(brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"));
  abs.a<-abs(df1[,8:dim(df1)[2]])>CI
  cl1<-pv[,8:dim(df1)[2]]< pv_cutoff
  col_res<-data.frame(abs.a*cl1)
  n1=max((dim(df1)[2]-7), 3)###either 3 or length of the variants compared to chosen state
  #cola<-c(brewer.pal(n = n1, name = "Paired"))
  xl_val<-max(abs(c(min(df1[,8:dim(df1)[2]]),max(df1[,8:dim(df1)[2]]))))
  plot(df1[,8], pv[,8], log="y", ylim=c(1, 10^-7), xlab="", ylab="", col=1,
       pch=20, yaxt="n", type="n", xlim=c(-xl_val, xl_val ))
  rect(xleft=c(min(df1[,8:dim(df1)[2]])-5),xright = -CI, ybottom = pv_cutoff, ytop=10^-50, col="grey95", lty=3)
  rect(xleft=CI,xright = max(df1[,8:dim(df1)[2]] )+5, ybottom = pv_cutoff, ytop=10^-50, col="grey95", lty=3)
  for ( i in 8:dim(df1)[2]){
    nm_of_state<-c(paste(str_sub(colnames(df1)[7], start = 4, end=-9)), " - ",
            str_sub(colnames(df1)[i], start = 4, end=-9))
    #message(paste(nm_of_state, collapse = ""))

    cola1<-c("black", cola[i-7])
    c1<-cola1[(col_res[,i-7])+1]
    points(df1[,i], pv[,i],col=c1, pch=20, type="p") }
  box(lwd=2)
  exp_ddu<-expression(Delta*' D'[2]*'O uptake [Da]')
  mtext(c(exp_ddu),  c(SOUTH<-1),line=0.9, outer=TRUE, adj=c(0.36), cex=0.9)
  mtext(c("p-value"),  c(WEST<-2),line=0.9, outer=TRUE, cex=0.9)

  axis(2, at=c(10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10), cex.axis=1, labels=T)
  at.y <- outer(1:9, 10^(2:(-7)))
  axis(2, at=at.y, labels=F,tcl=-0.2)}

#' Returns volcano plots for timepoints in the data frame
#'
#' Returns volcano plots for each peptide. Critical interval is calculated according to #' Reliable Identification of Significant Differences in Differential Hydrogen Exchange-Mass Spectrometry Measurements Using a Hybrid Significance Testing Approach
#' Tyler S. Hageman and David D. Weis
#' Analytical Chemistry 2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325
#' calculations for alpha 0.99
#' pvalues calculated using Welch t-test.
#'
#' @param df output from functions output_tp
#' @param replicates number of replicates in set as default set to 3.
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param cola color pallette for different Protein States. As default Paired pallette from color.Brewer is used.
#' @return volcano plots
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' plots_vol_tp(df=a, replicates=3, cola=c(1:4), pv_cutoff=0.01 )
#' plots_vol_tp(df=a, pv_cutoff=0.05)
#' @export
plots_vol_tp<-function(df,replicates=3, pv_cutoff=0.01, cola){
  if(missing(cola)) cola=c(brewer.pal(n = max((dim(ave_timepoint(df, replicates))[2]-7), 3), name = "Paired"));

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow=c(length(unique(df$Deut.Time)), 1),mar = c(1, 1, 1, 10), oma=c(2.5,2.5,0.1,0.1), cex.axis=1, cex.main=1, cex.lab=1.1,
      mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)

  av1<-ave_timepoint(df, replicates)
  da1<-dif_ave(av1)
  pv1<-pv_timepoint(df, replicates)
  s1<-sd_timepoint(df,  replicates)

  for ( i in(unique(df$Deut.Time))){
    CI<-max(CI_tp(s1[s1$Deut.Time==i,], replicates, pv_cutoff))
    vol_tp(da1[da1$Deut.Time==i,], pv1[pv1$Deut.Time==i,], CI, pv_cutoff, cola)
    message(paste("CI @", i, ":", round(CI, digits = 2)))

    mtext(i,  c(NORTH<-3),line=-1, outer=FALSE, cex=0.5)
  }
  par(mfrow=c(1, 1))

  lab_vol(da1, cola)
  }

