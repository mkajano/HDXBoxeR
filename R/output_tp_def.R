####################
#' Prepares output for HDX-MS for the deuteration uptake or percent deuteration for the time points.
#'
#' Returns a data frame organized for additional analysis.
#' In columns are deuteration uptake or percent deuteration data for the given protein states.
#' Function allows for writing csv with data, matching sequences of peptide.
#' Protein.States, Deut.times, or number of replicates can be specified.
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param replicates number of replicates to be used in analysis. The function takes number of replicates up to specified number. If no argument provided number maximal common number of replicates it used.
#' @param states function allows to choose what states should be used for analysis. Default all states are used.
#' @param times lists the deuteration times to be used in analysis. Default all states used.
#' @param seq_match Flag allows to choose if the peptide sequences should be matched between states. seq_match=FALSE signifies no sequence matching, seq_match=T states that the sequences are matched between the sets.
#' @param csv Flag allowing saving the output as csv. With default csv="NA", data is not saved. If csv output is desided, provide output name.
#' @param percent Flag allowing to choose output as deteuration uptake (FALSE) or percent deuteration (TRUE). Default deuteration uptake.
#' @return data frame with reorganized data where in columns is the deuteration uptake for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(filepath=file_nm) ###all default parameters used
#'
#'
#' # all possible flags listed & percent deuteration output,
#' # with sequences matching for protein states.
#'
#' a<-output_tp(filepath=file_nm, replicates=3, states=c("bound", "Unbound"),
#' times=c("3.00s", "72000.00s"), seq_match=TRUE, csv="NA", percent=TRUE)
#' @export
output_tp<- function(filepath, replicates, states, times, seq_match=FALSE, csv="NA", percent=FALSE){
  if(missing(states)) { states=arguments_call1(filepath); message(cat("Protein.States used:", states, sep=" "))}
  if(missing(times)) times=arguments_call2(filepath, states); message(cat("Deut.times used:", times, sep=" "))
  if(missing(replicates)) replicates=arguments_call3(filepath, states, times); message(cat("Number of replicates used:", replicates, sep=" "))

  a<-arg_df(filepath)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  ##creates temporary df, temp1, with Protein states going through all unique protein states

  uptake_seq_matchF<-function(a, percent){
    b<-c()
    ##creates temporary df, temp1, with Protein states going through all unique protein states
    for (time in times){
      temp1<-a[which(a$Deut.Time ==time),]
      st_l<-c()
      nbs=0
      for (state in states){##
        temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
        nb=0
        nbs=nbs+1
        st<-gsub(c(' '),'',state)
        df_nm_st<-paste(st, "_", nbs,sep="")
        st_l<-c(st_l, df_nm_st)
        bs<-c()
        for (exp in unique(temp2$Experiment)[1:replicates]){
          nb=nb+1
          df_nm<-paste("b",nb,sep="")
          temp3<-temp2[which(temp2$Experiment == exp),]
          n_tmp<-names(temp3)
          nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] ,
                 paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
          colnames(temp3)<-nms
          assign(df_nm, temp3)
          bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
          df_List<-mget(bs)
          ##will merge all the replicates dataframes to bp dataframe
          bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                       'Charge')), df_List)}
        assign(df_nm_st, bp)
      }
      df_List2<-mget(st_l)
      bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                    'Charge')), df_List2)
      b=rbind(b, bp2)
    } ## b has all information bound together again to have all information df
    b<-arrange(b, Deut.Time, Start, End, Charge)
    ord=b$Deut.Time
    ord=as.numeric(str_sub(ord, end=-2))
    bp1<-data.frame(b, ord)
    bp1<-arrange(bp1, ord)
    b<-bp1[,-dim(bp1)[2]]


    if (percent==FALSE){
      tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
    } else if (percent==T){
      tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
    }

    return(tp)}
  uptake_seq_matchT<-function(a, percent){
    b<-c()
    ##creates temporary df, temp1, with Protein states going through all unique protein states
    for (time in times){
      temp1<-a[which(a$Deut.Time ==time),]
      st_l<-c()
      nbs=0
      for (state in states){##
        temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
        nb=0
        nbs=nbs+1
        st<-gsub(c(' '),'',state)
        df_nm_st<-paste(st, "_", nbs,sep="")
        st_l<-c(st_l, df_nm_st)
        bs<-c()
        for (exp in unique(temp2$Experiment)[1:replicates]){
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
          bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End',  'Search.RT',
                                                       'Charge')), df_List)}
        assign(df_nm_st, bp)
      }
      df_List2<-mget(st_l)
      bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End',  'Search.RT',
                                                    'Charge')), df_List2)
      b=rbind(b, bp2)
    } ## b has all information bound together again to have all information df


    seq_df<- data.frame(b[,grep("Sequence", colnames(b))])
    seq_df<-data.frame(seq_df[,grep("Sequence_1", colnames(seq_df))])

    seq_list<-c()
    for (k in 1:dim(seq_df)[1]){
      if (length(unique(c(seq_df[k,])))==1){
        seq_list<-c(seq_list, seq_df[k,1])
      } else {
        seq1<-c()
        for (m in 1:dim(seq_df)[2]){
          seq1<-c(seq1, seq_df[k,m])}
        seq_list<-c(seq_list,paste(seq1, sep = "' '", collapse = "."))
      }
    }


    cl_nb1<-grep("Sequence", colnames(b))
    b<-b[,-cl_nb1]
    b<-data.frame(seq_list, b)
    colnames(b)[1]<- "Sequence"

    b<-arrange(b, Deut.Time, Start, End, Charge)
    ord=b$Deut.Time
    ord=as.numeric(str_sub(ord, end=-2))
    bp1<-data.frame(b, ord)
    bp1<-arrange(bp1, ord)
    b<-bp1[,-dim(bp1)[2]]

    df_description <- data.frame(b[, 1:6])

    if (percent==FALSE){
      tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
    } else if (percent==T){
      tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
    }

    return(tp)

  }



  if (seq_match==FALSE){
    all1<-uptake_seq_matchF(a, percent)
  } else if (seq_match==T){
    all1<-uptake_seq_matchT(a, percent)
  } else {warning("incorrect seq_match argument provided, halted")
  }

  if (csv=="NA"){
    message("no csv file written")
  } else {
    write.csv(all1, file=csv)}

  for ( i in 7:dim(all1)[2]){
    if (is.character(all1[,i])==TRUE){
      all1[,i]<-as.numeric(all1[,i])
    }
  }

  return(all1)}
