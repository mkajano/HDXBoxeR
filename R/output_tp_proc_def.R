

####################
#' Prepares output for HDX-MS for the percent deuteration for the timepoints.
#'
#' Returns a data frame organized for additional analysis.
#' In columns are uptake data for loaded protein states.
#' In the columns are uptake data for all protein state.
#' Number of replicates, protein states and Deut.Times can be specified
#'
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param replicates number of replicates to be used in analysis. The function takes number of replicates up to specified number. If no argument provided number maximal common number of replicates it used.
#' @param states function allows to choose what states should be used for analysis. Default all states are used.
#' @param times lists the deuteration times to be used in analysis. Default all states used.
#' @param seq_match Flag allows to choose if the peptide sequences should be matched between states. seq_match=F signifies no sequence matching, seq_match=T states that the sequences are matched between the sets.
#' @param csv Flag allowing saving the output as csv. With default csv="NA", data is not saved. If csv output is desided, provide output name.
#' @return data frame with reorganized data where in columns is the percent deuteration for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' #a<- output_tp_proc_def(filepath=file_nm) ###all default parameters used
#' #a<-output_tp_proc_def(filepath=file_nm, replicates=3,
#' states=c("Bound", "Unbound"), times=c("3.00s", "72000.00s"), seq_match=T, csv="NA")
#' @export
output_tp_proc_def<-function(filepath, replicates, states, times, seq_match=F, csv="NA"){
  if(missing(states)) { states=arguments_call(filepath)[[1]]; print(c("Protein.States used:", states))}
  if(missing(times)) times=arguments_call(filepath)[[2]]; print(c("Deut.times used:", times))
  if(missing(replicates)) replicates=arguments_call(filepath)[[3]]; print(c("Number of replicates used:", replicates))



  a <- read.csv(file = filepath, header = F, skip = 1)
  nm <- read.csv(file = filepath, header = T, row.names = NULL,
                 nrows = 1)
  nm1 <- colnames(nm)
  colnames(a) <- c(nm1)

  a <- a[, c(1:6, 8, 9, 17)]


  undeut_theor_cent_av<-function(a){


    if (all(a$Deut.Time == "0s") == FALSE & length(unique(a$Deut.Time ==
                                                          "0s")) == 2) {
      undeut <- a[which(a$Deut.Time == "0s"), ]}


    b <- c()

    st_l <- c()
    nbs = 0
    for (state in unique(undeut$Protein.State)) {
      un1 <- undeut[which(undeut$Protein.State == state), ]
      nb = 0
      nbs = nbs + 1
      st <- gsub(c(" "), "", state)

      df_nm_st <- paste(st, "_", nbs, sep = "")

      st_l <- c(st_l, df_nm_st)
      bs <- c()
      for (exp in unique(un1$Experiment)[1:length(unique(un1$Experiment))]) {

        nb = nb + 1
        df_nm <- paste("b", nb, sep = "")
        un3 <- un1[which(un1$Experiment == exp),]
        n_tmp <- names(un3)

        nms <- c(paste(st, n_tmp[1], "_", nb, sep = ""),
                 n_tmp[2], paste(st, n_tmp[3], "_", nb,
                                 sep = ""), n_tmp[4:8], paste(st, "_",
                                                              n_tmp[9], "_UN_", nb, sep = ""))
        colnames(un3) <- nms

        assign(df_nm, un3)
        bs <- c(bs, df_nm)
        df_List <- mget(bs)
        bp <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                       "Start", "End", "Sequence",
                                                       "Search.RT", "Charge"),all=T), df_List)


      }
      assign(df_nm_st, bp)


    }
    df_List2 <- mget(st_l)
    bp2 <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                    "Start", "End", "Sequence", "Search.RT",
                                                    "Charge"), all=T), df_List2)
    b = rbind(b, bp2)

    b <- arrange(b, Deut.Time, Start, End, Charge)
    ord = b$Deut.Time
    ord = as.numeric(str_sub(ord, end = -2))
    bp1 <- data.frame(b, ord)
    bp1 <- arrange(bp1, ord, Start, End)
    b <- bp1[, -dim(bp1)[2]]

    un.centr <- b[, c(2:6,  grep("Exp.Cent", colnames(b)))]

    av.u<-c()
    un_nm_c<-c()
    for ( i in paste(gsub(c(" "), "", unique(undeut$Protein.State)), "_Exp.Cent", sep="")){
      cm.1<-un.centr[, grep(i, colnames(un.centr))]
      if (is.data.frame(cm.1)==TRUE){
        av.u<- cbind(av.u,  rowMeans(cm.1, na.rm = T))} else(
          av.u<-cbind(av.u, cm.1))

      un_nm_c <- c(un_nm_c,paste(i, "_Exp.Cent_UN", sep = ""))
    }

    colnames(av.u)<-un_nm_c
    av.u<-data.frame(b[, 1:6], av.u)
    return(av.u)
  }
  FDdeut_theor_cent_av<-function(a){


    if (all(a$Deut.Time == "FD") == FALSE ) {
      fd <- a[which(a$Deut.Time == "FD"), ]
    }


    b <- c()

    st_l <- c()
    nbs = 0
    for (state in unique(fd$Protein.State)) {

      un1 <- fd[which(fd$Protein.State == state), ]
      nb = 0
      nbs = nbs + 1
      st <- gsub(c(" "), "", state)

      df_nm_st <- paste(st, "_", nbs, sep = "")

      st_l <- c(st_l, df_nm_st)
      bs <- c()
      for (exp in unique(un1$Experiment)[1:length(unique(un1$Experiment))]) {

        nb = nb + 1
        df_nm <- paste("b", nb, sep = "")
        un3 <- un1[which(un1$Experiment == exp),]
        n_tmp <- names(un3)

        nms <- c(paste(st, n_tmp[1], "_", nb, sep = ""),
                 n_tmp[2], paste(st, n_tmp[3], "_", nb,
                                 sep = ""), n_tmp[4:8], paste(st, "_",
                                                              n_tmp[9], "_FD_", nb, sep = ""))
        colnames(un3) <- nms

        assign(df_nm, un3)
        bs <- c(bs, df_nm)
        df_List <- mget(bs)
        bp <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                       "Start", "End", "Sequence",
                                                       "Search.RT", "Charge"), all=T), df_List)


      }
      assign(df_nm_st, bp)


    }
    df_List2 <- mget(st_l)
    bp2 <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                    "Start", "End", "Sequence", "Search.RT",
                                                    "Charge"), all=T), df_List2)
    b = rbind(b, bp2)

    b <- arrange(b, Deut.Time, Start, End, Charge)

    bp1 <- arrange(b, Start, End)
    b <- bp1[, -dim(bp1)[2]]

    un.centr <- b[, c(2:6,  grep("Exp.Cent", colnames(b)))]

    av.u<-c()
    un_nm_c<-c()
    for ( i in paste(gsub(c(" "), "", unique(undeut$Protein.State)), "_Exp.Cent", sep="")){
      cm.1<-un.centr[, grep(i, colnames(un.centr))]
      if (is.data.frame(cm.1)==TRUE){
        av.u<- cbind(av.u,  rowMeans(cm.1, na.rm = T))} else(
          av.u<-cbind(av.u, cm.1))

      un_nm_c <- c(un_nm_c,paste(i, "_Exp.Cent_FD", sep = ""))
    }

    colnames(av.u)<-un_nm_c
    av.u<-data.frame(b[, 1:6], av.u)

    return(av.u)
  }

  undeut_theor_cent_avS<-function(a){


    if (all(a$Deut.Time == "0s") == FALSE & length(unique(a$Deut.Time ==
                                                          "0s")) == 2) {
      undeut <- a[which(a$Deut.Time == "0s"), ]
    }


    b <- c()

    st_l <- c()
    nbs = 0
    for (state in unique(undeut$Protein.State)) {
      un1 <- undeut[which(undeut$Protein.State == state), ]
      nb = 0
      nbs = nbs + 1
      st <- gsub(c(" "), "", state)

      df_nm_st <- paste(st, "_", nbs, sep = "")

      st_l <- c(st_l, df_nm_st)
      bs <- c()
      for (exp in unique(un1$Experiment)[1:length(unique(un1$Experiment))]) {

        nb = nb + 1
        df_nm <- paste("b", nb, sep = "")
        un3 <- un1[which(un1$Experiment == exp),
        ]
        n_tmp <- names(un3)

        nms <- c(paste(st, "_",n_tmp[1], "_", nb, sep = ""),
                 n_tmp[2], paste(st, n_tmp[3], "_", nb,  sep = ""),
                 n_tmp[4:5],
                 paste(st, "_", n_tmp[6], "_", nb, sep = ""),
                 n_tmp[7:8],
                 paste(st, "_", n_tmp[9], "_UN_", nb, sep = ""))
        colnames(un3) <- nms

        assign(df_nm, un3)
        bs <- c(bs, df_nm)
        df_List <- mget(bs)
        bp <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                       "Start", "End",
                                                       "Search.RT", "Charge"),all=T), df_List)


      }
      assign(df_nm_st, bp)


    }
    df_List2 <- mget(st_l)
    bp2 <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                    "Start", "End",  "Search.RT",
                                                    "Charge"), all=T), df_List2)
    b = rbind(b, bp2)

    b <- arrange(b, Deut.Time, Start, End, Charge)
    ord = b$Deut.Time
    ord = as.numeric(str_sub(ord, end = -2))
    bp1 <- data.frame(b, ord)
    bp1 <- arrange(bp1, ord, Start, End)
    b <- bp1[, -dim(bp1)[2]]

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


    un.centr <- b[, c(c(1, 3:6),  grep("Exp.Cent", colnames(b)))]

    av.u<-c()
    un_nm_c<-c()
    for ( i in paste(gsub(c(" "), "", unique(undeut$Protein.State)), "_Exp.Cent", sep="")){
      cm.1<-un.centr[, grep(i, colnames(un.centr))]
      if (is.data.frame(cm.1)==TRUE){
        av.u<- cbind(av.u,  rowMeans(cm.1, na.rm = T))} else(
          av.u<-cbind(av.u, cm.1))

      un_nm_c <- c(un_nm_c,paste(i, "_Exp.Cent_UN", sep = ""))
    }

    colnames(av.u)<-un_nm_c
    av.u<-data.frame(b[, 1:6], av.u)
    return(av.u)
  }
  FDdeut_theor_cent_avS<-function(a){


    if (all(a$Deut.Time == "FD") == FALSE ) {
      fd <- a[which(a$Deut.Time == "FD"), ]
    }


    b <- c()

    st_l <- c()
    nbs = 0
    for (state in unique(fd$Protein.State)) {

      un1 <- fd[which(fd$Protein.State == state), ]
      nb = 0
      nbs = nbs + 1
      st <- gsub(c(" "), "", state)

      df_nm_st <- paste(st, "_", nbs, sep = "")

      st_l <- c(st_l, df_nm_st)
      bs <- c()
      for (exp in unique(un1$Experiment)[1:length(unique(un1$Experiment))]) {

        nb = nb + 1
        df_nm <- paste("b", nb, sep = "")
        un3 <- un1[which(un1$Experiment == exp),
        ]
        n_tmp <- names(un3)

        nms <- c(paste(st, "_",n_tmp[1], "_", nb, sep = ""),
                 n_tmp[2], paste(st, n_tmp[3], "_", nb,  sep = ""),
                 n_tmp[4:5],
                 paste(st, "_", n_tmp[6], "_", nb, sep = ""),
                 n_tmp[7:8],
                 paste(st, "_", n_tmp[9], "_FD_", nb, sep = ""))


        colnames(un3) <- nms

        assign(df_nm, un3)
        bs <- c(bs, df_nm)
        df_List <- mget(bs)
        bp <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                       "Start", "End",
                                                       "Search.RT", "Charge"), all=T), df_List)


      }
      assign(df_nm_st, bp)


    }
    df_List2 <- mget(st_l)
    bp2 <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                    "Start", "End", "Search.RT",
                                                    "Charge"), all=T), df_List2)
    b = rbind(b, bp2)

    b <- arrange(b, Deut.Time, Start, End, Charge)
    b <- arrange(b, Start, End)

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



    un.centr <- b[, c(c(1, 3:6),  grep("Exp.Cent", colnames(b)))]

    av.u<-c()
    un_nm_c<-c()
    for ( i in paste(gsub(c(" "), "", unique(undeut$Protein.State)), "_Exp.Cent", sep="")){
      cm.1<-un.centr[, grep(i, colnames(un.centr))]
      if (is.data.frame(cm.1)==TRUE){
        av.u<- cbind(av.u,  rowMeans(cm.1, na.rm = T))} else(
          av.u<-cbind(av.u, cm.1))

      un_nm_c <- c(un_nm_c,paste(i, "_Exp.Cent_FD", sep = ""))
    }

    colnames(av.u)<-un_nm_c
    av.u<-data.frame(b[, 1:6], av.u)


    return(av.u)
  }




  #######

  theor_cent_av_SeqF<-function(a){
  prep_all_theor_cent_av<-function(a){
    if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
      a<-a[-which(a$Deut.Time == c('0.00s')),]}
    if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
      a<-a[-which(a$Deut.Time == c('FD')),]}

    a <- na.omit(a)
    ###define what to do with missing parameters


    rownames(a) <- 1:dim(a)[1]
    b <- c()
    for (time in times) {


      temp1 <- a[which(a$Deut.Time == time), ]
      print(head(temp1))
      st_l <- c()
      nbs = 0
      for (state in states) {
        print(state)

        temp2 <- temp1[which(temp1$Protein.State == state),]
        nb = 0
        nbs = nbs + 1
        st <- gsub(c(" "), "", state)
        df_nm_st <- paste(st, "_", nbs, sep = "")
        st_l <- c(st_l, df_nm_st)
        bs <- c()
        for (exp in unique(temp2$Experiment)[1:replicates]) {

          nb = nb + 1
          df_nm <- paste("b", nb, sep = "")
          temp3 <- temp2[which(temp2$Experiment == exp), ]
          n_tmp <- names(temp3)
          nms <- c(paste(st, "_",n_tmp[1], "_", nb, sep = ""),
                   n_tmp[2], paste(st, n_tmp[3], "_", nb,  sep = ""),
                   n_tmp[4:8], paste(st, "_", n_tmp[9], "_", nb, sep = ""))
          colnames(temp3) <- nms
          assign(df_nm, temp3)
          bs <- c(bs, df_nm)
          df_List <- mget(bs)
          bp <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                         "Start", "End", "Sequence",
                                                         "Search.RT", "Charge"), all=T), df_List)
        }
        assign(df_nm_st, bp)
      }
      df_List2 <- mget(st_l)
      bp2 <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                      "Start", "End", "Sequence", "Search.RT",
                                                      "Charge"), all=T), df_List2)
      b = rbind(b, bp2)
    }
    b <- arrange(b, Deut.Time, Start, End, Charge)
    ord = b$Deut.Time
    ord = as.numeric(str_sub(ord, end = -2))
    bp1 <- data.frame(b, ord)
    bp1 <- arrange(bp1, ord, Start, End)
    b <- bp1[, -dim(bp1)[2]]

    th.cent <- data.frame(b[, 1:6], b[, grep("Exp.Cent",
                                             colnames(b))])


    return(th.cent)
  }



  un.a<-undeut_theor_cent_av(a)[,-1]
  fd.a<-FDdeut_theor_cent_av(a)[,-1]
  all.tc<-prep_all_theor_cent_av(a)



  df_centr<-list(all.tc, un.a, fd.a)

  big_cm<-Reduce(function(x, y) merge(x, y, by = c(
    "Start", "End", "Sequence", "Search.RT",
    "Charge"), all=T), df_centr)


  all1<-big_cm[,1:6]

  for ( i in gsub(c(" "), "", states)){
    temp4<-c()
    temp4<-big_cm[, grep(i,colnames(big_cm))]

    all1<-
      data.frame(all1, (temp4[,1:replicates]-temp4[,replicates+1])/
                   (temp4[, replicates+2]-temp4[,replicates+1])*100)}



  nm_all1<-str_replace_all(colnames(all1), "Exp.Cent", "Deut..")
  colnames(all1)<-nm_all1
  all1 <- na.omit(all1)

  all1<-arrange(all1, Deut.Time, Start, End, Charge)
  ord = all1$Deut.Time
  ord = as.numeric(str_sub(ord, end = -2))
  all2 <- data.frame(all1, ord)
  all2 <- arrange(all2, ord)
  all1 <- all2[, -dim(all2)[2]]
  return(all1)}

  theor_cent_av_SeqT<-function(a){
    prep_all_theor_cent_avS<-function(a){
      if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
        a<-a[-which(a$Deut.Time == c('0.00s')),]}
      if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
        a<-a[-which(a$Deut.Time == c('FD')),]}

      a <- na.omit(a)
      ###define what to do with missing parameters


      rownames(a) <- 1:dim(a)[1]
      b <- c()
      for (time in times) {

        temp1 <- a[which(a$Deut.Time == time), ]
        st_l <- c()
        nbs = 0
        for (state in states) {

          temp2 <- temp1[which(temp1$Protein.State == state), ]
          nb = 0
          nbs = nbs + 1
          st <- gsub(c(" "), "", state)
          df_nm_st <- paste(st, "_", nbs, sep = "")
          st_l <- c(st_l, df_nm_st)
          bs <- c()
          for (exp in unique(temp2$Experiment)[1:replicates]) {

            nb = nb + 1
            df_nm <- paste("b", nb, sep = "")
            temp3 <- temp2[which(temp2$Experiment == exp), ]

            n_tmp <- names(temp3)
            nms <- c(paste(st, "_",n_tmp[1], "_", nb, sep = ""),
                     n_tmp[2], paste(st, n_tmp[3], "_", nb,  sep = ""),
                     n_tmp[4:5],
                     paste(st, "_", n_tmp[6], "_", nb, sep = ""),
                     n_tmp[7:8],
                     paste(st, "_", n_tmp[9], "_", nb, sep = ""))
            colnames(temp3) <- nms
            assign(df_nm, temp3)
            bs <- c(bs, df_nm)
            df_List <- mget(bs)
            bp <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                           "Start", "End",
                                                           "Search.RT", "Charge"), all=T), df_List)
          }
          assign(df_nm_st, bp)
        }
        df_List2 <- mget(st_l)
        bp2 <- Reduce(function(x, y) merge(x, y, by = c("Deut.Time",
                                                        "Start", "End",  "Search.RT",
                                                        "Charge"), all=T), df_List2)
        b = rbind(b, bp2)
      }

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

      b <- arrange(b, Deut.Time, Start, End, Charge)

      df_description <- data.frame(b[, 1:6])

      th.cent <- data.frame(df_description, b[, grep("Exp.Cent",colnames(b))])

      return(th.cent)
    }



    un.a<-undeut_theor_cent_avS(a)[,-c(2)]
    fd.a<-FDdeut_theor_cent_avS(a)[,-c(2)]
    all.tc<-prep_all_theor_cent_avS(a)


    df_centr<-list(all.tc, un.a, fd.a)

    big_cm<-Reduce(function(x, y) merge(x, y, by = c(
      "Start", "End", "Sequence", "Search.RT",
      "Charge")), df_centr)

    all1<-big_cm[,1:6]


    for ( i in gsub(c(" "), "", states)){
      temp4<-c()
      temp4<-big_cm[, grep(i,colnames(big_cm))]

      all1<-
        data.frame(all1, (temp4[,1:replicates]-temp4[,replicates+1])/
                     (temp4[, replicates+2]-temp4[,replicates+1])*100)}

    nm_all1<-str_replace_all(colnames(all1), "Exp.Cent", "Deut..")
    colnames(all1)<-nm_all1
    all1 <- na.omit(all1)


    all1<-arrange(all1, Deut.Time, Start, End, Charge)
    ord = all1$Deut.Time
    ord = as.numeric(str_sub(ord, end = -2))
    all2 <- data.frame(all1, ord)
    all2 <- arrange(all2, ord)
    all1 <- all2[, -dim(all2)[2]]

    return(all1)}


  if (seq_match==F){
    all1<-theor_cent_av_SeqF(a)
  } else if (seq_match==T){
    all1<-theor_cent_av_SeqT(a)
  } else {print("incorrect seq_match argument provided, not matching sequence")
    all1<-theor_cent_av_SeqF(a)
    }

all1<-data.frame(all1$Deut.Time, all1$Start, all1$End, all1$Sequence, all1$Search.RT, all1$Charge,
                 all1[,7:dim(all1)[2]])

colnames(all1)[1:6]<-c("Deut.Time", "Start", "End", "Sequence", "Search.RT", "Charge")

  if (csv=="NA"){
    print("no csv file written")
  } else {
    write.csv(all1, file=csv)}
  return(all1)
}

