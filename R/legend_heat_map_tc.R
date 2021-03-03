#' Legend for the heatmaps for timecourses.
#'
#' Returns names for legend for the heatmaps. Extracts names from data.frame
#'
#' @param df generated using output_tcourse()
#' @return legend for the heatmap
#' @export
legend_heat_map_tc<-function(df){
  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start=4, end=-9)
  axis(2, at=0.5:(dim(df)[2]-6.5), labels=nm1, las=2, line = , cex.axis=0.7)
}

#' Legend for the heatmaps prep function for timecourses.
#'
#' Returns names for legend for the heatmaps
#'
#' @param ranges ranges that are to be colored in the legend. Default ranges=c(-Inf,seq(-30, 30, by=10), Inf )
#' @return legend for the heatmap
#' @export
legend_heat_map_timecourse<-function(ranges=c(-Inf,seq(0, 100, by=10), Inf )){
  cbr1<-color_ranges_Spectral(ranges, c("white"))

  leg_nm<-c("no coverage")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],"%:", ranges[i+1], "%", sep=""))}

  pallette_ll(cbr1, leg_nm)}
