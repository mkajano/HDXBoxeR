#' Return woods plots for the timepoints
#'
#' All the peptides are plotted based on their uptake.
#'
#' @param input_data output from function output_tp(..., percent=TRUE)
#' @param cola colors, default NA
#' @param times Deuteration times, if missing all deuteration times used
#' @param replicates replicates
#' @param ylim y axis limits
#' @param ... other parameters
#' @return Woods plots for the timepoints
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm, percent=TRUE)
#' deuteration_woods_timepoints(a[1:12,])
#' }
#' @export
deuteration_woods_timepoints<-function(input_data,times, replicates=3,
                                       cola=NA, ylim=c(0,120), ...) {
  if(missing(times)) times=unique(input_data$Deut.Time)

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))

  a1<-ave_timepoint(input_data, replicates)
  s1<-sd_timepoint(input_data, replicates)

  if (is.na(cola[1])==FALSE){

  } else if (length(7:dim(a1)[2])<9){
    cola<-c(1,brewer.pal(9,"Set1"))
  } else{
    cola<-c(1,colorRampPalette(brewer.pal(9,"Set1"))(length(7:dim(a1)[2])))}


  par(mfrow=c(length(times), 1), mar = c(1.5, 1.5, 1.5, 1.5), oma=c(4,4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)


  for (time in times){
    indc<-c()
    indc=0
    pl_gen_ch2(a1, ddlab = 1,ylim=ylim, ...)

    a2<-c()
    a2<-a1[a1$Deut.Time==time,]
    s2<-c()
    s2<-s1[s1$Deut.Time==time,]

    for (j in 7:dim(a2)[2]){
      mtext(time,  c(North<-3), line=-1, outer=FALSE, cex=0.6)
      indc=indc+1
      for ( i in 1:dim(a2)[1]){

       suppressWarnings(arrows((a2$Start[i]+2 + a2$End[i])/2, a2[i, j]-s2[i,j],
               (a2$Start[i]+2 + a2$End[i])/2, a2[i, j]+s2[i,j], length=0.02,
               angle=90, code=3, col=cola[indc]))
        points(c(a2$Start[i]+2, a2$End[i]), c(a2[i, j], a2[i, j]), type = "l", col=cola[indc])
      }}}
  legend_states_PerD_bottom(input_data, cola[1:length(7:dim(a1)[2])])
}




#' Return woods plots for the timecourse
#'
#' All the peptides are plotted based on their uptake.
#'
#' @param input_data output from function output_tc(..., percent=TRUE)
#' @param states states, if missing all states used
#' @param replicates replicates
#' @param ylim y axis limits
#' @param ... other parameters
#' @return Woods plots for the timecourse
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(file_nm, percent=TRUE)
#' deuteration_woods_timecourse(a)
#' }
#' @export
deuteration_woods_timecourse<-function(input_data, states, replicates=3, ylim=c(0,120), ...) {
  if(missing(states)) states=unique(input_data$Protein.State)
  a1<-ave_timepoint(input_data, replicates)
  s1<-sd_timepoint(input_data, replicates)
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  cola<-colorRampPalette(brewer.pal(9,"YlGnBu")[4:9])(length(7:dim(a1)[2]))


  par(mfrow=c(length(states), 1), mar = c(1.5, 1.5, 1.5, 1.5), oma=c(4,4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)


  for (state in states){
    indc<-c()
    indc=0
    pl_gen_ch2(a1, ddlab = 1, ylim=ylim, ...)

    a2<-c()
    a2<-a1[a1$Protein.State==state,]
    s2<-c()
    s2<-s1[a1$Protein.State==state,]

    for (j in 7:dim(a2)[2]){
      mtext(state,  c(North<-3), line=-1, outer=FALSE, cex=0.6)
      indc=indc+1
      for ( i in 1:dim(a2)[1]){

        suppressWarnings(arrows((a2$Start[i]+2 + a2$End[i])/2, a2[i, j]-s2[i,j],
               (a2$Start[i]+2 + a2$End[i])/2, a2[i, j]+s2[i,j], length=0.02,
               angle=90, code=3, col=cola[indc]))
        points(c(a2$Start[i]+2, a2$End[i]), c(a2[i, j], a2[i, j]), type = "l", col=cola[indc])
      }}}
  legend_states_PerD_bottom(input_data, cola[1:length(7:dim(a1)[2])])
}



#' Returns a woods plot for comparisons of the timepoints samples
#'
#' Modification of butterfly plot. x axis residues.
#' y axis % deuteration for
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
#' @param ylim y axis limit
#' @param ... other variables
#' @param alpha critical interval, to have more restrictive use lower values, default=0.01
#' @return Woods plots with chosen statistically different peptides
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(file_nm)
#' b<-output_tc(file_nm, percent=TRUE)
#' woods_CI_plot(thP=b, th=a, pv_cutoff = 0.001, alpha = 0.01, replicates=3)
#' }
#' @export
woods_CI_plot<-function(thP, th, replicates=3,
                        pv_cutoff=0.01, states, alpha=0.01, ylim=c(0,120), ...){
  if(missing(states)) states=unique(th$Protein.State)

  nm<-colnames(ave_timepoint(th, replicates))
  nm1<-str_sub(nm[7:length(nm)], start=4, end=-10)
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  pl1f<-function(){
    plot(x=1, type = "n", xlim=c(min(thP$Start), max(thP$End)), ylab="",
         xlab="", yaxt="n", ylim=ylim, ...)
    axis(1, at=seq(0, 1000, by=10), cex.axis=1, labels=FALSE,tcl=-0.2)
    axis(2, at=seq(-1000, 1000, by=50), cex.axis=1, labels=c(rev(seq(50,1000, by=50)), seq(0,1000, by=50)))
    axis(2, at=seq(-1000, 1000, by=10), cex.axis=1, labels=FALSE,tcl=-0.2)
    exp_ddu<-expression('% Deuteration')
    mtext(c("Residue"),  c(SOUTH<-1),line=0.3, outer=TRUE, cex=0.8)
    mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE, cex=0.85)
  }


  par(mfcol=c(length(nm1),length(states)-1), mar = c(1.5, 1.5, 1.5, 1.5), oma=c(4,4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)


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
      text(x=(min(thP$Start)+max(thP$End))/2,y=115, paste(states[1], "vs", state,";", nm1[nb1-1]), cex=0.7)}

  }
  legend_nm_bottom(c(states[1], "other state"), c(cola[8], colg[8]))
}

