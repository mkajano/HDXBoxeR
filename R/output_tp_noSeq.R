####################
#' Prepares output for HDX-MS for the timepoints without matching sequences.
#'
#' Returns a dataframe organized for additional analysis.
#' Function does not match for the sequence. It can be used to match mutant sequences.
#' In columns are uptake data for loaded protein states.
#' Number of replicates must be the same for all protein states.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return data frame with reorganized data where in columns is uptake data for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp_noSeq(file_nm)
#' @export
output_tp_noSeq<- function(filepath){
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
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in unique(temp1$Protein.State)){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:5] ,
               paste(st,n_tmp[6],"_",nb,sep=""),n_tmp[7:8] ,
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3)
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs)
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp)
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))

  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]

  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )


  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}

####################
#' Prepares output for HDX-MS for the timepoints without matching sequences.
#'
#' Returns a dataframe organized for additional analysis.
#' Function does not match for the sequence. It can be used to match mutant sequences.
#' In columns are procent deuteration data for loaded protein states.
#' Number of replicates must be the same for all protein states.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return data frame with reorganized data where in columns is procent deuteration data for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp_noSeq_proc(file_nm)
#' @export
output_tp_noSeq_proc<- function(filepath){
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
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in unique(temp1$Protein.State)){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:5] ,
               paste(st,n_tmp[6],"_",nb,sep=""),n_tmp[7:8] ,
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3)
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs)
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp)
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))

  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]

  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )


  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}

####################
#' Prepares output for HDX-MS for the timepoints without matching sequences for selected states.
#'
#' Returns a dataframe organized for additional analysis.
#' Function does not match for the sequence. It can be used to match mutant sequences.
#' In columns are uptake data for selected protein states.
#' Number of replicates must be the same for all protein states.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param states list of states that will be chosen for analysis
#' @return data frame with reorganized data where in columns is uptake data for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' names_states<- nm_states(file_nm)
#' a<- output_tp_noSeq_states(file_nm, states=names_states)
#' @export
output_tp_noSeq_states<- function(filepath, states){

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
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in states){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:5] ,
               paste(st,n_tmp[6],"_",nb,sep=""),n_tmp[7:8] ,
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3)
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs)
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp)
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))

  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]

  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )


  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}

#' Prepares output for HDX-MS for the timepoints without matching sequences for selected states.
#'
#' Returns a dataframe organized for additional analysis.
#' Function does not match for the sequence. It can be used to match mutant sequences.
#' In columns are procent deuteration for selected protein states.
#' Number of replicates must be the same for all protein states.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param states list of states that will be chosen for analysis
#' @return data frame with reorganized data where in columns is procent deuteration for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' names_states<- nm_states(file_nm)
#' a<- output_tp_noSeq_states_proc(file_nm, states=names_states)
#' @export
output_tp_noSeq_states_proc<- function(filepath, states){

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
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in states){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:5] ,
               paste(st,n_tmp[6],"_",nb,sep=""),n_tmp[7:8] ,
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3)
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs)
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp)
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))

  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]

  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )


  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}
