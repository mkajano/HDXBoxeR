#' Prepares output with HDX-MS data for publications
#'
#' Format prepared based of example from:
#' Masson, G.R., Burke, J.E., Ahn, N.G. et al. Recommendations for performing, interpreting and reporting hydrogen deuterium exchange mass spectrometry (HDX-MS) experiments. Nat Methods 16, 595–602 (2019). https://doi.org/10.1038/s41592-019-0459-y
#' It generates csv file in format ready for publication of the data.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param output_name Name of output file. It has to be csv file
#' @param replicates number of replicates to be used in analysis. The function takes number of replicates up to specified number. If no argument provided number maximal common number of replicates it used.
#' @param states function allows to choose what states should be used for analysis. Default all states are used.
#' @param times lists the deuteration times to be used in analysis. Default all states used.
#' @param percent return either uptake or percent deuteration, default=FALSE, return uptake
#' @return Returns&saves data.frame in format that is accepted for the publications.
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom graphics abline axis box boxplot legend mtext par plot points polygon rect text arrows
#' @importFrom stats df na.omit qt sd t.test
#' @importFrom utils head read.csv tail type.convert write.csv write.table
#' @importFrom dplyr arrange %>%
#' @importFrom stringr str_sub
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr separate
#' @importFrom methods rbind2
#' @importFrom wrapr orderv
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' output_prep(filepath=file_nm, output_name=tempfile())
#' }
#' @export
output_prep<-function(filepath, output_name, states, replicates, times, percent=FALSE){
  if(missing(states)) { states=arguments_call1(filepath); message(cat("Protein.States used:", states))}
  if(missing(times)) times=arguments_call2(filepath, states); message(cat("Deut.times used:", times))
  if(missing(replicates)) replicates=arguments_call3(filepath, states, times); message(cat("Number of replicates used:", replicates))

  undeut<-arg_UN_FD(filepath)[[1]]
  FD<-arg_UN_FD(filepath)[[2]]
  a<-arg_UN_FD(filepath)[[3]]

  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  b<-c()
  for (state in states){##
    temp1<-a[which(a$Protein.State ==state ),] ##creates temporary df, temp1, with Protein states going through all unique protein states
    for (time in times){
      temp2<-temp1[which(temp1$Deut.Time ==time),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      bs<-c()
      for (exp in unique(temp2$Experiment)[1:replicates]){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1:2],paste(n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , paste(n_tmp[9],"_",nb,sep=""), paste(n_tmp[10],"_",nb,sep=""),
               paste(n_tmp[11],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3)
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs)
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      b=rbind(b, bp)}} ## b has all information bound together again to have all information df

  ###order columns
  #### names(b)[grep("Exp.Cent", colnames(b))] gives position of all variable with pattern and return position
  ord1<-c('Protein.State', 'Deut.Time', 'Start', 'End', 'Sequence', 'Search.RT', 'Charge',
          names(b)[grep("Experiment", colnames(b))],
          names(b)[grep("Exp.Cent", colnames(b))], names(b)[grep("X..Deut", colnames(b))], names(b)[grep("Deut.._", colnames(b))] )
  b<-b[,ord1]



  ###calculate means +sd and write to out dataframe. Later assign names, choose only important

  FD2<-data.frame(FD[,c(1,2,4:6,8,7,9)], rep(0,times=dim(FD)[1]), rep(0,times=dim(FD)[1]), rep(0,times=dim(FD)[1]))
  undeut2<-data.frame(undeut[,c(1,2,4:6,8,7,9)], rep(0,times=dim(undeut)[1]),  rep(0,times=dim(undeut)[1]), rep(0,times=dim(undeut)[1]))

  if (percent==FALSE){
    col_deut<-c(grep("X..Deut", colnames(b)))
    chars <- sapply(b[, col_deut], is.character)
    if (chars[1] == TRUE){
      b[ , col_deut[which(chars==TRUE)]] <- as.data.frame(apply(b[ , col_deut[which(chars==TRUE)]], 2, as.numeric))}

  ###prepare Full deuteration data.frame to be bound without data.frame.
    out<-data.frame(b[,1:7], rowMeans(b[,grep("Exp.Cent", colnames(b))]), apply(b[,grep("Exp.Cent", colnames(b))],1,sd),
                    round(rowMeans(b[,grep("X..Deut", colnames(b))]), digits = 4),
                    round(apply(b[,grep("X..Deut", colnames(b))],1,sd), digits=4))

  nm_final<-c(names(out)[1:5], "Retention_Time_[min]","Charge",
              "Mean.Peptide_Mass_[Da]", "st.dev_Peptide_Mass",
              "Mean.Deut.Uptake_[Da]", "st.dev_Deut_Uptake")

  colnames(FD2)<-nm_final
  colnames(undeut2)<-nm_final
  colnames(out)<-nm_final
  } else { ###percent deuteration
    col_deut<-c(grep("Deut.._", colnames(b)))
    chars <- sapply(b[, col_deut], is.character)
    if (chars[1] == TRUE){
      b[ , col_deut[which(chars==TRUE)]] <- as.data.frame(apply(b[ , col_deut[which(chars==TRUE)]], 2, as.numeric))}


    out<-data.frame(b[,1:7], rowMeans(b[,grep("Exp.Cent", colnames(b))]), apply(b[,grep("Exp.Cent", colnames(b))],1,sd),
                    round(rowMeans(b[,grep("Deut.._", colnames(b))]), digits = 4),
                    round(apply(b[,grep("Deut.._", colnames(b))],1,sd), digits=4))

    nm_final<-c(names(out)[1:5], "Retention_Time_[min]","Charge",
                "Mean.Peptide_Mass_[Da]", "st.dev_Peptide_Mass",
                "Mean.Percent_Deut_[Da]", "st.dev_Percent_Deut")

    colnames(FD2)<-nm_final
    colnames(undeut2)<-nm_final
    colnames(out)<-nm_final
  }



  out<-rbind(out, FD2)
  out<-rbind(out, undeut2)
  out<-out[,c(1:8,10:11)]


  out<-(arrange(out, Start, End, Protein.State, Charge))
  ##find a way to order correctly assign negative time to zero and
  ##very long time for FD and order based of these columns
  ##if bp does not work use out
  ord=out$Deut.Time
  ord=type.convert(ord, as.is = TRUE)
  ord[ord=="0s"]<- c("-5.0")
  ord[ord== "FD"]<- c("10000000000000.0s")
  ord=as.numeric(str_sub(ord, end=-2))

  bp<-data.frame(out, ord)
  bp<-arrange(bp, Start, End, Protein.State, Charge, ord)
  bp<-bp[,-dim(bp)[2]]
  ###write output
  write.csv(bp, output_name, row.names = FALSE)
  return()}




#' Returns initially processed data.frame from the export from the HDXExaminer
#'
#' Function used as internal function
#'
#' @param filepath input file location
#' @return Data.frame for further processing
#' @export
arg_UN_FD<-function(filepath){
  a <- read.csv(file = filepath, header = F, skip = 1)
  nm <- read.csv(file = filepath, header = T, row.names = NULL,
                 nrows = 1)
  nm1 <- colnames(nm)

  if(length(unique(nm1=="row.names"))==2){
    nm1<-nm1[-which(nm1=="row.names")]
    nm1<-c(nm1, "row.names")}

  colnames(a) <- c(nm1)

  dif_col<- setdiff(c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                      "Charge", "Search.RT", "X..Deut", "Deut.."), nm1)

  if (length(dif_col)!=0)
  {stop(paste("input file missing column(s): ", dif_col, sep=""))}


  undeut <- a[a$Deut.Time=="0s",c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                                  "Charge", "Search.RT","Exp.Cent", "X..Deut", "Deut..")]
  FD <- a[a$Deut.Time=="FD",c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                              "Charge", "Search.RT","Exp.Cent", "X..Deut", "Deut..")]


  a <- a[,c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
            "Charge", "Search.RT","Exp.Cent", "X..Deut", "Deut..")]
  if (all(a$Deut.Time == "0.00s") == FALSE & length(unique(a$Deut.Time ==
                                                           "0.00s")) == 2) {
    a <- a[-which(a$Deut.Time == c("0.00s")), ]  }
  if (all(a$Deut.Time == "FD") == FALSE & length(unique(a$Deut.Time ==
                                                        "FD")) == 2) {
    a <- a[-which(a$Deut.Time == c("FD")), ]}
  if (all(a$Deut.Time == "0s") == FALSE & length(unique(a$Deut.Time ==
                                                        "0s")) == 2) {
    a <- a[-which(a$Deut.Time == c("0s")), ] }
  a <- na.omit(a)


  return(list(undeut, FD, a))
}


