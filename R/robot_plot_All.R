#' Returns a robot plot for comparisons of the timepoints samples
#'
#' Modification of butterfly plot. x axis residues.
#' y axis % deuteration for one variant above the axis and for second peptide below the axis.
#' Peptides are compared between the sets for the significance change between sets.
#' If there is significant change beteween sets peptides are plotted for all timepoints.
#' Significanty different timepoints for the peptides are colored.
#' Peptides ranges are plotted as a line at corresponding % deuteration values.
#'
#'
#' @param thP output of output_tcourse_proc() function. Raw data for procent deuteration for time courses
#' @param th output of output_tcourse() function. Raw data for uptake deuteration for time courses
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param states Protein states from the set. As default all states are chosen.
#' @param CI_factor Multiplication factor for Critical Interval. Allows for more restrictive selection of Critial interval.
#' @return Robot maps for timecourses
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' tm_df<-output_tc(filepath=file_nm)
#' tmP_df<-output_tc(filepath=file_nm, percent=TRUE)
#' robot_plot_All(thP = tmP_df, th=tm_df, pv_cutoff=0.001)
#'
#' # more restrictive peptide selection
#' robot_plot_All(thP = tmP_df, th=tm_df, pv_cutoff=0.001, CI_factor=3)
#'
#' @export
robot_plot_All<-function(thP, th, replicates=3,
                         pv_cutoff=0.01, states, CI_factor=1){
  if(missing(states)) states=unique(thP$Protein.State)

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  par(mfrow=c(length(states)-1, 1), mar = c(1.5, 1.5, 1.5, 1.5), oma=c(4,4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)

  for ( state in states[2:length(states)]) {

    control_df<- thP[thP$Protein.State==states[1],]
    variant_df<- thP[thP$Protein.State==state,]

    control_df_up<- th[th$Protein.State==states[1],]
    variant_df_up<- th[th$Protein.State==state,]

    pv1<-pv_timecourse(df_c = control_df_up, df_v=variant_df_up, replicates)
    lav.proc<-prep_timecourse_plot_ave(control_df, variant_df, replicates)
    lav.proc_up<-prep_timecourse_plot_ave(control_df, variant_df, replicates)

    sh_avc<-lav.proc[[1]]
    sh_avv<-lav.proc[[2]]
    sh_avc_up<-lav.proc_up[[1]]
    sh_avv_up<-lav.proc_up[[2]]

    CI_all<-prep_timecourse_plot_sd(control_df_up, variant_df_up, replicates=3, pv_cutoff = pv_cutoff)
    CI_all<-CI_all*CI_factor

    cola<-(brewer.pal(n = length(7:dim(sh_avc)[2])+1, name = "Oranges"))


    plot(x=1, type = "n", ylim=c(-120, 120), xlim=c(min(thP$Start), max(thP$End)), ylab="",
         xlab="", yaxt="n")
    axis(1, at=seq(0, 1000, by=10), cex.axis=1, labels=FALSE,tcl=-0.2)
    axis(2, at=seq(-1000, 1000, by=50), cex.axis=1, labels=c(rev(seq(50,1000, by=50)), seq(0,1000, by=50)))
    axis(2, at=seq(-1000, 1000, by=10), cex.axis=1, labels=FALSE,tcl=-0.2)
    exp_ddu<-expression('% Deuteration')
    mtext(c("Residue"),  c(SOUTH<-1),line=0.3, outer=TRUE, cex=0.8)
    mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE, cex=0.85)

    nb1=1
    peptide_all<-c()
    for (i in 7:dim(sh_avc)[2]) {
      peptide_all<-c(peptide_all, which(pv1[, i]<pv_cutoff &
                                          abs(sh_avc_up[, i]-sh_avv_up[, i]) > CI_all[i-6]))
    }
    peptide_all<-sort(unique(peptide_all))

    colg<-(brewer.pal(n = length(7:dim(sh_avc)[2])+2, name = "Blues"))

    for ( i in dim(sh_avc)[2]:7){
      xpoly<-c((sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2,
               rev((sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2))
      ypoly<-c(sh_avc[peptide_all,i], rev(sh_avv[peptide_all,i]*(-1)))
      polygon(x =xpoly,                           # X-Coordinates of polygon
              y = ypoly,                             # Y-Coordinates of polygon
              col = colg[i-6])}
    abline(h=0)

    for ( j in 7:dim(sh_avc)[2]){
      peptide_index<-which(pv1[,j]<pv_cutoff & abs(sh_avc_up[,j]-sh_avv_up[,j]) > CI_all[j-6])
      nb1=nb1+1
      for ( i in peptide_all){
        points(c(sh_avc$Start[i], sh_avc$End[i]), c(sh_avc[i,j],sh_avc[i,j] ), type="l", col="grey45")
        points(c(sh_avv$Start[i], sh_avv$End[i]), c(sh_avv[i,j],sh_avv[i,j])*(-1), type="l", col="grey45")
      }
      # for ( pinx in peptide_index){
      #   points(c(sh_avc$Start[pinx], sh_avc$End[pinx]), c(sh_avc[pinx,j],sh_avc[pinx,j] ), type="l", col=cola[nb1],  lwd=2)
      #   points(c(sh_avv$Start[pinx], sh_avv$End[pinx]), c(sh_avv[pinx,j],sh_avv[pinx,j])*(-1), type="l", col=cola[nb1], lwd=2)
      # }

      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avc[peptide_all,j] ), type="p", col="grey45", pch=20, lwd=2)
      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avv[peptide_all,j])*(-1), type="p", col="grey45",pch=20, lwd=2)
      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avc[peptide_all,j] ), type="l", col=cola[nb1], pch=20)
      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avv[peptide_all,j])*(-1), type="l", col=cola[nb1],pch=20)

      points(c(sh_avc$Start[peptide_index]+sh_avc$End[peptide_index])/2, c(sh_avc[peptide_index,j] ), type="p", col=cola[nb1], pch=20, lwd=2)
      points(c(sh_avc$Start[peptide_index]+sh_avc$End[peptide_index])/2, c(sh_avv[peptide_index,j])*(-1), type="p", col=cola[nb1],pch=20, lwd=2)
    }
    text(x=(min(thP$Start)+max(thP$End))/2,y=107.5, states[1], cex=0.7)
    text(x=(min(thP$Start)+max(thP$End))/2,y=-107.5, state, cex=0.7)
  }
  legend_tc_bottom(sh_avc, cola[2:length(cola)])

  }
