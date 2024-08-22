#' Returns a woods plot for comparisons of the timecourse samples
#'
#' Peptides are compared between the sets for the significance change between sets.
#' The y-axis is fraction deuteration: that is uptake divided by number of exchangeable protons
#' If there is significant change beteween sets peptides are plotted for all timepoints.
#' Significanty different timepoints for the peptides are colored.
#' Peptides ranges are plotted as a line at corresponding fraction deuteration values.
#'
#'
#' @param th output of output_tcourse() function. Raw data for uptake deuteration for time courses
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param states Protein states from the set. As default all states are chosen.
#' @param ylim y axis limit
#' @param method options either uptake or fraction
#' @param ... other variables
#' @param alpha critical interval, to have more restrictive use lower values, default=0.01
#' @return Woods plots with chosen statistically different peptides
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(file_nm)
#' woods_CI_plot_frac(th=a, pv_cutoff = 0.001, alpha = 0.01, replicates=3)
#' }
#' @export
#'
woods_CI_plot_frac<-function(th, replicates=3,
                        pv_cutoff=0.01, states, alpha=0.01, ylim=c(0,5),
                        method="fraction", ...){
  if(missing(states)) states=unique(th$Protein.State)

  nm<-colnames(ave_timepoint(th, replicates))
  nm1<-str_sub(nm[7:length(nm)], start=4, end=-10)
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  if (method=="fraction"){
  thP<-data.frame(th[,1:6],th[,7:dim(th)[2]]/nb_exch_deut(th))


  pl1f<-function(){
    xlim1=c(min(thP$Start), max(thP$End))
    plot(x=1, type = "n", xlim=xlim1, ylab="",
         xlab="", yaxt="n", ylim=ylim, ...)
    axis(1,at=pretty(xlim1, n=5), cex.axis=1, labels=FALSE,tcl=-0.2)
    axis(2, at=pretty(ylim, n=3), cex.axis=1, labels=pretty(ylim, n=3))
    axis(2, at=pretty(ylim, n=10), cex.axis=1, labels=FALSE,tcl=-0.2)
    exp_ddu<-expression('Fraction deuterated')
    mtext(c("Residue"),  c(SOUTH<-1),line=0.3, outer=TRUE, cex=0.8)
    mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE, cex=0.85)
  }} else {
    thP<-th


    pl1f<-function(){
      xlim1=c(min(thP$Start), max(thP$End))
      plot(x=1, type = "n", xlim=xlim1, ylab="",
           xlab="", yaxt="n", ylim=ylim, ...)
      axis(1, at=pretty(xlim1, n=5), cex.axis=1, labels=FALSE,tcl=-0.2)
      axis(2, at=pretty(ylim, n=3), cex.axis=1, labels=pretty(ylim, n=3))
      axis(2, at=pretty(ylim, n=10), cex.axis=1, labels=FALSE,tcl=-0.2)
      exp_ddu<-expression('Fraction deuterated')
      mtext(c("Residue"),  c(SOUTH<-1),line=0.3, outer=TRUE, cex=0.8)
      mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE, cex=0.85)
  }}


  par(mfcol=c(length(nm1),length(states)-1), mar = c(1.5, 1.5, 1.5, 1.5),
      oma=c(4,4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2,
      bg="white", font.lab=2, font.axis=2)





  for ( state in states[2:length(states)]) {

    control_df<- thP[thP$Protein.State==states[1],]
    variant_df<- thP[thP$Protein.State==state,]

    control_df_up<- th[th$Protein.State==states[1],]
    variant_df_up<- th[th$Protein.State==state,]

    pv1<-pv_timecourse(df_c = control_df_up, df_v=variant_df_up, replicates)
    lav.proc<-prep_timecourse_plot_ave(control_df, variant_df, replicates)
    lav.proc_up<-prep_timecourse_plot_ave(control_df_up, variant_df_up, replicates)

    sh_avc<-lav.proc[[1]]
    sh_avv<-lav.proc[[2]]
    sh_avc_up<-lav.proc_up[[1]]
    sh_avv_up<-lav.proc_up[[2]]

    CI_all<-prep_timecourse_plot_sd(control_df_up, variant_df_up, replicates=3,
                                    alpha=alpha)

    print(CI_all)

    cola<-(brewer.pal(n = 9, name = "Reds"))
    colg<-(brewer.pal(n = 9, name = "Blues"))

    peptide_all<-c()
    for (i in 7:dim(sh_avc)[2]) {
      peptide_all<-c(peptide_all, which(pv1[, i]<pv_cutoff &
                                          abs(sh_avc_up[, i]-sh_avv_up[, i]) > CI_all[i-6]))
    }
    peptide_all<-sort(unique(peptide_all))




    nb1=1
    for ( j in 7:dim(sh_avc)[2]){
      pl1f()

      peptide_index<-which(pv1[,j]<pv_cutoff & abs(sh_avc_up[,j]-sh_avv_up[,j]) > CI_all[j-6])
      nb1=nb1+1
      for ( i in peptide_all){
        points(c(sh_avc$Start[i], sh_avc$End[i]), c(sh_avc[i,j],sh_avc[i,j] ), type="l", col=cola[4])


        points(c(sh_avv$Start[i], sh_avv$End[i]), c(sh_avv[i,j],sh_avv[i,j]), type="l", col=colg[4])
      }

      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avc[peptide_all,j] ), type="p", col="grey25", pch=20, lwd=2)
      points(c(sh_avc$Start[peptide_index]+sh_avc$End[peptide_index])/2, c(sh_avc[peptide_index,j] ), type="p", col=cola[8], pch=20, lwd=2)


      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avv[peptide_all,j]), type="p", col="grey65",pch=20, lwd=2)

      points(c(sh_avc$Start[peptide_index]+sh_avc$End[peptide_index])/2, c(sh_avv[peptide_index,j]), type="p", col=colg[8],pch=20, lwd=2)
      text(x=(min(thP$Start)+max(thP$End))/2,y=ylim[2]-0.1, paste(states[1], "vs", state,";", nm1[nb1-1]), cex=0.7)}

  }
  legend_nm_bottom(c(states[1], "other state"), c(cola[8], colg[8]))
}






path="C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/mia/all_data.csv"

path="C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/heterooligomers/B1_heterooligomers_Allresults_culled.csv"


#path=system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

names_states<- nm_states(path)
th<-output_tc(path,replicates=3, states =names_states[c(1:3)], seq_match = F )
woods_CI_plot_frac(th, alpha=0.01, pv_cutoff = 0.01, ylim=c(0,1.1))
woods_CI_plot_frac(th, alpha=0.1, pv_cutoff = 0.01, method="uptake", ylim=c(0,16))


thP<-data.frame(th[,1:6],th[,7:dim(th)[2]]/nb_exch_deut(th))
h<-output_tp(path,replicates=3, states =names_states[c(1:3)], seq_match = F )


h1<-output_tp(path,replicates=3, states =names_states[c(3)], seq_match = F,
              times ="300.00s"  )
h2<-output_tp(path,replicates=3, states =names_states[c(2)], seq_match = F ,
              times ="300.00s" )



plot_peptide_sig_tp(h, pv_cutoff = 0.01)

plots_av_tp(h)
plots_diff_tp(h)

plot_heat_map_max_uptake_tp(h)


prep_timecourse_plot_sd(h1, h2, replicates=3,
                        alpha=0.01)

s1<-sd_timepoint(h1)
s2<-sd_timepoint(h2)

s<-sd_timepoint(h)




which(s1[,7]> 0.1)


s1[which(s1[,7]> 0.15),]
s2[which(s2[,7]> 0.1),]
