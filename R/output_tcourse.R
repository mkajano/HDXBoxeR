
######
#' Prepares output for HDX-MS for the timecourses for deuteration uptake data.
#'
#' Returns a dataframe organized for additional analysis.
#' In reprocessed data frame columns have time course data.
#' Function takes raw data from HDXexaminer and organize data into columns with uptake deuteration data.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return Returns reprocessed data. Data.frame with reorganized data where in columns are deuterium uptake for different timepoints in the analysis.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tcourse(file_nm)
#' @export
output_tcourse<- function(filepath){

  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
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

  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
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
    b=rbind(b, bp2)

  } ## b has all information bound together again to have all information df
  b<-arrange(b, Protein.State,Start, End,  Charge)
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}

######
#' Prepares output for procent deuteration of HDX-MS for the timecourses.
#'
#' Returns a dataframe organized for additional analysis.
#' In reprocessed data frame columns have time course data for procent deuteration.
#' Function takes raw data from HDXexaminer and organize data into columns with % deuteration
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return Returns reprocessed data. Data.frame with reorganized data where in columns are deuterium uptake for different timepoints in the analysis.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tcourse_proc(file_nm)
#' @export
output_tcourse_proc<- function(filepath){

  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
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
  b<-c()

  ##creates temporary df, temp1, with Protein states going through all unique protein states

  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
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
    b=rbind(b, bp2)

  } ## b has all information bound together again to have all information df
  b<-arrange(b, Protein.State, Start, End,  Charge)
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}
