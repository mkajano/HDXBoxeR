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
#' @export
deuteration_woods_timepoints_frac<-function(input_data,times, replicates=3,
                                       cola=NA, ylim=c(0,1.1), ...) {
  if(missing(times)) times=unique(input_data$Deut.Time)

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))

  input_data2<-data.frame(input_data[,1:6],
                                               input_data[,7:dim(input_data)[2]]/nb_exch_deut(input_data))
  a1<-ave_timepoint(input_data2, replicates)

  s1<-sd_timepoint(input_data2, replicates)

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
    pl_gen_ch_frac(a1, ddlab = 1,ylim=ylim, ...)

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


#' Return woods plots for the timepoints
#'
#' All the peptides are plotted based on their uptake.
#'
#' @param list_input_data list of multiple outputs from function
#' output_tp(..., percent=TRUE)
#' @param cola colors, default NA
#' @param times Deuteration times, if missing all deuteration times used
#' @param replicates replicates
#' @param ylim y axis limits
#' @param ... other parameters
#' @return Woods plots for the timepoints
#' @export


deuteration_woods_timepoints_frac_list<-function(list_input_data,times, replicates=3,
                                            cola=NA, ylim=c(0,1.1), ...) {
  if(missing(times)) times=unique(input_data$Deut.Time)

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))

  ls_input_data<-list()
  input_data<-c()
  ls_a1<-list()
  ls_s1<-list()
  Start_seq<-c()
  End_seq<-c()
  nb_datasets=0
  nms_dt<-c()
  for (len_list in 1: length(list_input_data)){
    input_data<-c()
    input_data<-data.frame(list_input_data[[len_list]][,1:6],
                                   list_input_data[[len_list]][,7:dim(list_input_data[[len_list]])[2]]/nb_exch_deut(list_input_data[[len_list]]))
    ls_a1[[len_list]]<-ave_timepoint(input_data, replicates)

    ls_s1[[len_list]]<-sd_timepoint(input_data, replicates)
    Start_seq<-c(Start_seq, min(list_input_data[[len_list]]$Start))
    End_seq<-c(End_seq, max(list_input_data[[len_list]]$End))
  nb_datasets=nb_datasets+dim(ls_a1[[len_list]])[2]-6
  nms_dt<-c(nms_dt, unique(str_sub(colnames(input_data[7:dim(input_data)[2]]), start = 1, end = -11)))
  }


  if (is.na(cola[1])==FALSE){

  } else if (nb_datasets<9){
    cola<-c(1,brewer.pal(9,"Set1"))
  } else{
    cola<-c(1,colorRampPalette(brewer.pal(9,"Set1"))(nb_datasets))}


  par(mfrow=c(length(times), 1), mar = c(1.5, 1.5, 1.5, 1.5), oma=c(4,4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)



  for (time in times){
    indc<-c()
    indc=0
    pl_gen_ch_list_frac(min(Start_seq)[1], max(End_seq)[1],
                        ddlab = 3,ylim=ylim, ...)
    s1<-c()
    a1<-c()
for (list_nb in 1: length(ls_a1)){

  a1<-ls_a1[[list_nb]]

  s1<-ls_s1[[list_nb]]


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
      }}}}
legend_nm_bottom(nms_dt, cola[1:nb_datasets])
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
deuteration_woods_timecourse_frac<-function(input_data, states,
                                            replicates=3, ylim=c(0,1.2), ...) {
  if(missing(states)) states=unique(input_data$Protein.State)

  input_data2<-data.frame(input_data[,1:6],
                          input_data[,7:dim(input_data)[2]]/nb_exch_deut(input_data))
  a1<-ave_timepoint(input_data2, replicates)

  s1<-sd_timepoint(input_data2, replicates)
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))


  cola<-colorRampPalette(brewer.pal(9,"YlGnBu")[4:9])(length(7:dim(a1)[2]))


  par(mfrow=c(length(states), 1), mar = c(1.5, 1.5, 1.5, 1.5), oma=c(4,4,2,2), cex.axis=1,
      cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)


  for (state in states){
    indc<-c()
    indc=0
    pl_gen_ch_frac(a1, ddlab = 1, ylim=ylim, ...)

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

