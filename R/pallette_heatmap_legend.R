#' Color scheme using heatmap. Legend Extracts names from data.frame
#'
#' Returns names for legend for the heatmaps
#'
#' @param col_pallette pallette to be used in the heat map
#' @return legend for the heatmap
#' @export
pallette_legend<-function(col_pallette){
  xl <- 1; yb <- 1; xr <- 1.5; yt <- 2 ## seq values
  plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")# mock plot
  rect(xl,head(seq(yb,yt,(yt-yb)/length(col_pallette)),-1),xr,tail(seq(yb,yt,(yt-yb)/length(col_pallette)),-1),col=col_pallette)# make rect
  return()}

#' Color scheme using heatmap. Legend extracts names from data frame
#'
#' Returns names for legend for the heatmaps
#'
#' @param pallette pallette to be used in the heat map
#' @param lab labels to be used in pallette
#' @return legend for the heatmap
#' @export
pallette_ll<-function(pallette, lab){
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))
  ppar(c(1,1))
  pallette_legend(pallette)
  yb=1; yt=2
  ypos=seq(yb,yt,(yt-yb)/length(pallette))
  ypos.a<-c()

  for ( i in 1:(length(ypos)-1)){
    ypos.a<-c(ypos.a, (ypos[i]+ypos[i+1])/2)}
  axis(4, at=ypos.a,labels =lab,las=2, cex.axis=0.6, tcl=0, pos=1.5)
}

#' Legend for the heatmaps.Extracts names from data.frame
#'
#' Returns names for legend for the heatmaps
#'
#' @param df average data frame. Generated using ave_timepoint() function.
#' @return legend for the heatmap
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' legend_heat_map_tp(df=a)
#' @export
legend_heat_map_tp<-function(df){
  nm1<-str_sub(colnames(df[8:dim(df)[2]]), start=4, end=-9)
  axis(2, at=0.5:(dim(df)[2]-7.5), labels=nm1, las=2, line = , cex.axis=0.7)
}

#' Legend for the heatmaps percent.Extracts names from data.frame
#'
#' Returns names for legend for the heatmaps
#'
#' @param df average data frame.
#' @return legend for the heatmap prercent
#' @export
legend_heat_map_tp_proc<-function(df){
  nm1<-str_sub(colnames(df[8:dim(df)[2]]), start=4, end=-8)
  axis(2, at=0.5:(dim(df)[2]-7.5), labels=nm1, las=2, line = , cex.axis=0.7)
}

#' Legend for the heatmaps prep function.
#'
#' Returns names for legend for the heatmaps
#'
#' @param ranges ranges that are to be colored in the legend. Default ranges=c(-Inf,seq(-30, 30, by=10), Inf )
#' @return legend for the heatmap
#' @export
legend_heat_map<-function(ranges=c(-Inf,seq(-30, 30, by=10), Inf )){
  ##color set up
  cbr1<-color_ranges_Blue_Red_heat_map(ranges, c( "grey45", "white"))

  leg_nm<-c("no coverage", "not-significant")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],"%:", ranges[i+1], "%", sep=""))}
  pallette_ll(cbr1, leg_nm)}

#' Legend for the significant peptides
#'
#' Returns names for legend for the significant peptides plots.
#'
#' @param ranges ranges that are to be colored in the legend. Default ranges=c(-Inf,seq(-30, 30, by=10), Inf )
#' @return legend for the heatmap
#' @export
legend_sig_peptides<-function(ranges=c(-Inf,seq(-30, 30, by=10), Inf )){
  ##color set up
  cbr1<-color_ranges_Blue_Red_heat_map(ranges, c("black"))

  leg_nm<-c("not-significant")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],"%:", ranges[i+1], "%", sep=""))}
  pallette_ll(cbr1, leg_nm)}

