####################
#' Prepares output for HDX-MS Full deuteration data
#'
#' Returns a data frame for Full deuteration set
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return data frame with reorganized data where in columns is uptake data for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<-output_FD(file_nm)
#'
#' @export
  output_FD<- function(filepath){


    arg_df_FD<-function(filepath){
      a <- read.csv(file = filepath, header = FALSE, skip = 1)
      nm <- read.csv(file = filepath, header = TRUE, row.names = NULL,
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


      a <- a[,c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                "Charge", "Search.RT", "X..Deut", "Deut..")]

      FD<-a[which(a$Deut.Time == c("FD")), ]

      return(FD)
    }

FD<-arg_df_FD(filepath)



  a<-na.omit(FD)
  rownames(a)<-1:dim(a)[1] ##name rows
  #### Full deuteration to analysis.
  nbs=0
  fd<-c()

  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(FD$Protein.State)){
    st<-gsub(c(' '),'',state)
    df_nm_st<-paste(st, "_", nbs,sep="")
    nbs=nbs+1
    temp1<-FD[which(FD$Protein.State ==state),]
    st_l<-c()
    nb=0
    fs<-c()
    for (exp in unique(temp1$Experiment)){
      nb=nb+1
      df_nm<-paste("FD_",nb,sep="")
      temp2<-temp1[which(temp1$Experiment == exp),]
      n_tmp<-names(temp2)
      nms<-c(n_tmp[1],paste("t_FD",n_tmp[2],"_",nb,sep=""),paste("t_FD",n_tmp[3],"_",nb,sep=""),n_tmp[4:8] ,
             paste("t_FD", "_",n_tmp[9],"_",nb,sep=""), paste("t_FD", "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
      colnames(temp2)<-nms
      assign(df_nm, temp2)
      fs<-c(fs, df_nm) ##creates number of data.frames that equals to number of replicates
      df_List<-mget(fs)
      ##will merge all the replicates dataframes to bp dataframe
      fp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                   'Charge')), df_List)}
    assign(df_nm_st, fp)
    }

  nb_X..deut=grep("X..Deut", colnames(fp))
  uptake.fp<-data.frame(fp[,c(1, 4:8, nb_X..deut)])
  nb_deut=grep("Deut.._", colnames(fp))
  procent.fp<-data.frame(fp[,c(1, 4:8, nb_deut)])

  return(uptake.fp)}

####################
#' Prepares output for HDX-MS  Undeuterated sample data.
#'
#' Returns a data frame for Full deuteration set
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return data frame with reorganized data where in columns is uptake data for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_UD(file_nm)
#' @export
output_UD<- function(filepath){
  a<-read.csv(file=filepath,  header = FALSE, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = TRUE, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9, 21, 22)] ##choose only useful columns
  rownames(a)<-1:dim(a)[1] ##name rows
  a<-na.omit(a) ##remove missing values, removes non-deuterated state
  zero<-a[which(a$Deut.Time == '0.00s'),] ##this line can be used to take advantage of control sample
  nbs=0
  ud<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(zero$Protein.State)){
    st<-gsub(c(' '),'',state)
    df_nm_st<-paste(st, "_", nbs,sep="")
    nbs=nbs+1
    temp1<-zero[which(zero$Protein.State ==state),]
    st_l<-c()
    nb=0
    us<-c()
    for (exp in unique(temp1$Experiment)){
      nb=nb+1
      df_nm<-paste("zero_",nb,sep="")
      temp2<-temp1[which(temp1$Experiment == exp),]
      n_tmp<-names(temp2)
      nms<-c(n_tmp[1],paste("t_zero",n_tmp[2],"_",nb,sep=""),paste("t_zero",n_tmp[3],"_",nb,sep=""),n_tmp[4:8] ,
             paste("t_zero", "_",n_tmp[9],"_",nb,sep=""), paste("t_zero", "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
      colnames(temp2)<-nms
      assign(df_nm, temp2)
      us<-c(us, df_nm) ##creates number of data.frames that equals to number of replicates
      df_List<-mget(us)
      ##will merge all the replicates dataframes to bp dataframe
      ud<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                   'Charge')), df_List)}
    assign(df_nm_st, ud) }
  nb_X..deut=grep("X..Deut", colnames(ud))
  uptake.ud<-data.frame(ud[,c(1, 4:8, nb_X..deut)])
  nb_deut=grep("Deut.._", colnames(ud))
  procent.ud<-data.frame(ud[,c(1, 4:8, nb_deut)])
  return(uptake.ud)}

####################
#' Prepares output for HDX-MS Full deuteration data  for procent deuteration.
#'
#' Returns a data frame for Full deuteration set
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return data frame with reorganized data where in columns is procent deuteration for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_FD_proc(file_nm)
#' @export
output_FD_proc<- function(filepath){
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
  #### Full deuteration to analysis.
  nbs=0
  fd<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(FD$Protein.State)){
    temp1<-FD[which(FD$Protein.State ==state),]
    st<-gsub(c(' '),'',state)
    df_nm_st<-paste(st, "_", nbs,sep="")
    nbs=nbs+1
    st_l<-c()
    nb=0
    fs<-c()
    for (exp in unique(temp1$Experiment)){
      nb=nb+1
      df_nm<-paste("FD_",nb,sep="")
      temp2<-temp1[which(temp1$Experiment == exp),]
      n_tmp<-names(temp2)
      nms<-c(n_tmp[1],paste("t_FD",n_tmp[2],"_",nb,sep=""),paste("t_FD",n_tmp[3],"_",nb,sep=""),n_tmp[4:8] ,
             paste("t_FD", "_",n_tmp[9],"_",nb,sep=""), paste("t_FD", "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
      colnames(temp2)<-nms
      assign(df_nm, temp2)
      fs<-c(fs, df_nm) ##creates number of data.frames that equals to number of replicates
      df_List<-mget(fs)
      ##will merge all the replicates dataframes to bp dataframe
      fp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                   'Charge')), df_List)}
    assign(df_nm_st, fp) }

  nb_X..deut=grep("X..Deut", colnames(fp))
  uptake.fp<-data.frame(fp[,c(1, 4:8, nb_X..deut)])
  nb_deut=grep("Deut.._", colnames(fp))
  procent.fp<-data.frame(fp[,c(1, 4:8, nb_deut)])

  return(procent.fp)}

####################
#' Prepares output for HDX-MS Undeuterated data for procent deuteration.
#'
#' Returns a data frame for Undeuterated control set
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @return data frame with reorganized data where in columns is procent deuteration for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_UD_proc(file_nm)
#' @export
output_UD_proc<- function(filepath){
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
  zero<-a[which(a$Deut.Time == '0.00s'),]
  nbs=0
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ud<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(zero$Protein.State)){
    st<-gsub(c(' '),'',state)
    df_nm_st<-paste(st, "_", nbs,sep="")
    nbs=nbs+1
    temp1<-zero[which(zero$Protein.State ==state),]
    st_l<-c()
    nb=0
    us<-c()
    for (exp in unique(temp1$Experiment)){
      nb=nb+1
      df_nm<-paste("zero_",nb,sep="")
      temp2<-temp1[which(temp1$Experiment == exp),]
      n_tmp<-names(temp2)
      nms<-c(n_tmp[1],paste("t_zero",n_tmp[2],"_",nb,sep=""),paste("t_zero",n_tmp[3],"_",nb,sep=""),n_tmp[4:8] ,
             paste("t_zero", "_",n_tmp[9],"_",nb,sep=""), paste("t_zero", "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
      colnames(temp2)<-nms
      assign(df_nm, temp2)
      us<-c(us, df_nm) ##creates number of data.frames that equals to number of replicates
      df_List<-mget(us)
      ##will merge all the replicates dataframes to bp dataframe
      ud<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                   'Charge')), df_List)}
    assign(df_nm_st, ud) }
  nb_X..deut=grep("X..Deut", colnames(ud))
  uptake.ud<-data.frame(ud[,c(1, 4:8, nb_X..deut)])
  nb_deut=grep("Deut.._", colnames(ud))
  procent.ud<-data.frame(ud[,c(1, 4:8, nb_deut)])

  return(procent.ud)}
