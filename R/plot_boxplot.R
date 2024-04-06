#' Plots boxplots for all the averages in the set
#'
#' Returns boxplots to compare sets between each other
#'
#' @param df average data frame. Generated using ave_timepoint() function.
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ... inherited boxplot parameters
#' @return boxplots for average deuterium uptake per set.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' boxplot_tp(df=a, replicates=3)
#' @export

boxplot_tp<-function(df, replicates=3, ...){
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))

  pparLM(c(length(unique(df$Deut.Time)),1))
  av1<-ave_timepoint(df, replicates)
  s1<-sd_timepoint(df, replicates)
  nm1<-str_sub(colnames(av1[c(6+(1:(dim(s1)[2]-6)))]), start = 4, end=-9)
  for ( i in(unique(df$Deut.Time))){
    boxplot(av1[av1$Deut.Time==i,7:dim(av1)[2]], names = nm1, ...)
    mtext(i, side=3, outer=FALSE, line=-0.9, cex=0.5 )}
  mtext("Average", side=2, outer=TRUE, line=0.8, cex=1 )
}

