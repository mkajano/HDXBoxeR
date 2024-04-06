#' Uptake plots
#'
#' Uptake plots per peptide
#'
#' @param input_data output from function output_tp(..., percent=T)
#' @param cola colors, default NA
#' @param timepoints the labeling times
#' @param replicates replicates
#' @param seq_match Flag TRUE or FALSE, default TRUE, match sequence of the protein states
#' @return  Uptake plots
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(file_nm, percent=TRUE)
#' x=c(3,60, 1800, 72000)
#' uptake_plots(a, x)
#' }
#' @export
uptake_plots<-function(input_data, timepoints,replicates=3,cola=NA, seq_match=TRUE) {
  indd<-c()
  a1<-ave_timepoint(input_data, replicates)
  s1<-sd_timepoint(input_data, replicates)

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  if (seq_match==TRUE){
  indd<-duplicate_sets(input_data[,c(3,4,6,1)])
  } else {
    indd<-duplicate_sets(input_data[,c(3,4,6)])
  }

  states=unique(a1$Protein.State)

  if (is.na(cola[1])==FALSE){

  } else if (length(7:dim(a1)[2])<9){
    cola<-c(1,brewer.pal(9,"Set1"))
  } else{
    cola<-c(1,colorRampPalette(brewer.pal(9,"Set1"))(length(7:dim(a1)[2])))}

  par(mar = c(1.5, 1.5, 1.5, 1.5), oma=c(4,4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2,
      bg="white", font.lab=2, font.axis=2)

  for ( i in 1:length(indd)){
    pl_gen_uptake(input_data, timepoints, ylim=c(5,115))
    coli=0
    mtext(paste("Res: ", a1$Start[indd[[i]][1]], "-",a1$End[indd[[i]][1]]),  c(North<-3), line=-1, outer=FALSE, cex=0.6)
    mtext( a1$Sequence[indd[[i]][1]],  c(South<-1),
           line=-1, outer=FALSE, cex=0.4)

    for ( j in 1:length(indd[[i]])){
      coli=coli+1
      arrows(timepoints, as.numeric(a1[indd[[i]][j], 7:dim(a1)[2]]+s1[indd[[i]][j], 7:dim(a1)[2]]),
             timepoints, as.numeric(a1[indd[[i]][j],7:dim(a1)[2]]-s1[indd[[i]][j], 7:dim(a1)[2]]),
             length=0.05,
             angle=90, code=3, col=cola[coli])
      points(timepoints, a1[indd[[i]][j], 7:dim(a1)[2]], type = "o", col=cola[coli], pch=20, lwd=4)
    }
  }


  legend_nm_bottom(states, cola[1:length(states)])
  reset_par()
}

