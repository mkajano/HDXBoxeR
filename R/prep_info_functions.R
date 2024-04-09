######
#' Provides summary table for all data.sets.
#'
#' Returns data frame sumamrizing general information about the data sets.
#' Function returns: Protein states, timepoints, number of replicates,  # peptides, % coveregae, average peptide length and redundancy.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return Returns summary table.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- general_info(file_nm)
#' @export
general_info<- function(filepath){
  a<-read.csv(file=filepath,  header = FALSE, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = TRUE, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)]
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),]
    a<-a[-which(a$Deut.Time == c('FD')),]}

  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows

  ##loop below will go through Protein states, timepoints and Experiments to get replicates
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"

  ##creates temporary df, temp1, with Protein states going through all unique protein states

  summ<-c()
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    summ<-c(summ, state)
    st_l<-c()
    nbs=0
    tmps<-paste(as.vector(unique(temp1$Deut.Time)),collapse=" ")
    summ<-c(summ, tmps)

    summ<-c(summ, length(unique(temp1$Experiment))/length(as.vector(unique(temp1$Deut.Time))))
    for (time in unique(temp1$Deut.Time)){##

      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()

      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] ,
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3)
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs)
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp)
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)

    summ<-c(summ, dim(bp2)[1])

    cvr<-coverage_residue(bp2,start_col = 2, end_col = 3)[min(bp2$Start):max(bp2$End)]
    cvrLong<-coverage_residue(bp2,start_col = 2, end_col = 3)
    ppt_cov<-round((length(cvr)-length(which(cvr==0)))/length(cvr)*100, digits=2)
    av<-round(mean(bp2$End-bp2$Start),digits=2)

    redund<-c()
    for (nb1 in 1:dim(bp2)[1]){
      r1<-sum(cvrLong[bp2[nb1,2]:bp2[nb1,3]])/(bp2[nb1,3]-bp2[nb1,2]+1)
      redund<-c(redund, r1)
    }
    summ<-c(summ,  ppt_cov, av, round(mean(redund), digits=2))
    b=rbind(b, bp2)
  }

  sum1<-data.frame(matrix(summ,nrow=length(unique(a$Protein.State)), byrow = TRUE))
  names(sum1)<-c("Protein.State", "Timepoints", "# replicates", "# peptides", "peptide coverage %", "<Peptide Length>", "<Redundancy>")
  return(sum1)}




